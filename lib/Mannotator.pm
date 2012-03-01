###############################################################################
#
#    Copyright (C) 2011 2012 Michael Imelfort Florent Angly Connor Skennerton
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
package Mannotator;

use strict;
use warnings;
use threads;
use base 'Exporter';
use File::Path qw( remove_tree make_path ) ;
use File::Basename;
use File::Spec::Functions qw( :ALL ); # catfile
use Bio::SeqIO;
use Bio::Seq;
use Bio::SearchIO;
use Bio::Tools::GFF;
use Class::Struct;


our @EXPORT = qw(gff2fasta generateAnnotations loadU2A splitGffs splitFasta combineGffs blastUnknowns splitBlastResults 
                cleanTmps annotate createFlatFile insertFeature debugFeature nextInList recombineGff3);
our $VERSION = '0.1';

sub prodigal_seqname {
    my ($comment) = @_;
    # Prodigal's comment line may contain the sequence's true name. Try to get it!
    #   # Sequence Data: seqnum=1;seqlen=16231;seqhdr="4 len=16231"
    #   Prodigal_Seq_1  Prodigal_v2.50  CDS     3       1592    158.2   -       0       ID=1_1;partial=10;start_type=ATG;rbs_m
    # Should be
    #   4  Prodigal_v2.50  CDS     3       1592    158.2   -       0       ID=1_1;partial=10;start_type=ATG;rbs_m
    my $seqname = undef;
    if ($comment =~ m/Sequence Data:.*seqhdr="(\S+) /i ) {
        $seqname = $1;
    }
    return $seqname;
}

sub gff2fasta {
    # this is a simple way of getting the orfs from a single
    # gff file. Use it when there is no need to do any fancy 
    # spliting or recombining
    my ($gff, $fasta, $tmp_fasta_ref, $tmp_gff) = @_;
    open my $gff_fh, "<", $gff or die "Error: Could not read file $gff\n$!\n";
    open my $tmp_gff_fh, '>', $tmp_gff or die "ERROR: Could not read file $tmp_gff\n$!\n";
    print $tmp_gff_fh "##gff-version 3\n";
    my $gff_def = <$gff_fh>;
    if ($gff_def !~ m/^##gff-version\s+3/i) {
        die "Error: File $gff does not seem to be a valid GFF file\nDid not see header line: ##gff-version 3\n";
    }
    my $true_fasta_header = undef;
    my %gff_hash;
    while(<$gff_fh>){
          # ##FASTA section indicates end of ##gff section
            if ($_ =~ m/^##FASTA/i) {
               last;
            }

            # comment lines split sequences or include prodigal information
            if ($_ =~ m/^#(.*)$/) {
                my $prodigal_seqname = prodigal_seqname($1);
                $true_fasta_header = $prodigal_seqname if defined $prodigal_seqname;
                next;
            }
      my ($seqid, $field_1, $field_2, $start, $end,
          $field_5, $field_6, $field_7, $attrs) = split;
      if (defined $true_fasta_header) {
          $seqid = $true_fasta_header;
      }
      print $tmp_gff_fh "$seqid\t$field_1\t$field_2\t$start\t$end\t$field_5\t$field_6\t$field_7\t$attrs\n";
      push @{$gff_hash{$seqid}}, [$start, $end, $attrs];
    }
    close $gff_fh;
    open my $out_fasta_fh, '>', ${$tmp_fasta_ref} or die "ERROR: could not read file $fasta\n$1\n";

    my $seqio = Bio::SeqIO-> new( -file => $fasta, -format => 'fasta' );
     while(my $sobj = $seqio->next_seq){
         my $seqid = $sobj->id;

        next if not defined $gff_hash{$seqid};
        my $seq = $sobj->seq;

        for(@{$gff_hash{$seqid}}){
            my ($start, $end, $attrs) = @$_;

            print $out_fasta_fh ">$seqid"."_".$start."_".$end."\n";

            print $out_fasta_fh substr($seq, $start, $end-$start+1), "\n";
        }
    }
}
sub splitGffs {
    #----
    # split the gffs into multiple files, one for each sequence
    # and place them in separate directories
    #
    my ($gff_files, $tmp_folders_ref, $tmp_prefix_ref) = @_;
    my @gffs = split /,/, $gff_files;
    foreach my $gff (@gffs)
    {
        print "Parsing GFF file $gff\n";
        my $current_fh;
        my $current_fasta_header = '__DOOF';
        my $true_fasta_header = undef;
        open my $gff_fh, "<", $gff or die "Error: Could not read file $gff\n$!\n";

        my $gff_def = <$gff_fh>;
        if ($gff_def !~ m/^##gff-version\s+3/i) {
           die "Error: File $gff does not seem to be a valid GFF file\nDid not see header line: ##gff-version 3\n";
        }

        while(<$gff_fh>)
        {
            # ##FASTA section indicates end of ##gff section
            if ($_ =~ m/^##FASTA/i) {
               last;
            }

            # comment lines split sequences or include prodigal information
            if ($_ =~ m/^#(.*)$/) {
                my $prodigal_seqname = prodigal_seqname($1);
                $true_fasta_header = $prodigal_seqname if defined $prodigal_seqname;
                next;
            }
            my @bits = split /\t/, $_;
            my $fasta_header = $bits[0];

            if (defined $true_fasta_header)
            {
                $fasta_header = $true_fasta_header;
                s/^(.+?)\t/$true_fasta_header\t/;
            }

            if($fasta_header ne $current_fasta_header)
            {
                if($current_fasta_header ne '__DOOF')
                {
                    # this is not the first run.
                    # so we have a file to close...
                    close $current_fh;
                }
               
                $current_fasta_header = $fasta_header;
                  
                # make a new temporary directory
                my $tmp_folder = ${$tmp_prefix_ref}.'_'.$current_fasta_header;
                make_path( $tmp_folder ) if not -e $tmp_folder;
                $tmp_folders_ref->{$tmp_folder} = 1;
                    
                # write in output file
                my $out_file = catfile( $tmp_folder, basename($gff) );
                if (not -e $out_file) {
                    # Start a new file and print the GFFheader information
                    open $current_fh, '>', $out_file or die "Error: Could not write file $out_file\n$!\n";
                    print $current_fh $gff_def;

                } else {
                    # Add to an exiting file (if all the features of a sequence are not contiguous in a GFF file)
                    open $current_fh, '>>', $out_file or die "Error: Could not write file $out_file\n$!\n";
                }

            }
            print $current_fh $_;
        }
        if("__DOOF" ne $current_fasta_header)
        {
            # we have a file to close...
            close $current_fh;
        }
    }
}


sub splitFasta {
    #-----
    # Simple script to split a fasta sequence and put into multiple folders
    # If the sequence had some annotations, then a folder with a name based
    # on the sequence fasta header should have been created.
    #
    my ($contig_file, $tmp_prefix_ref) = @_;
    my $seqio_object = Bio::SeqIO->new(-file => $contig_file, -format => "fasta");
    while (my $seq = $seqio_object->next_seq)
    {
        # Skip sequence if it has no annotation (and hence no folder)
        my $tmp_folder = ${$tmp_prefix_ref}.'_'.$seq->display_id;
        next if not -e $tmp_folder;
        # Write sequence in folder
        my $seq_fn = catfile( $tmp_folder, 'sequence.fa' );
        open my $ffh, '>', $seq_fn or die "Error: Could not write file $seq_fn\n$!\n";
        print $ffh ">" . $seq->display_id() . "\n";
        print $ffh $seq->seq(). "\n";
        close $ffh;
    }
}


sub combineGffs {
    #-----
    # Wrapper to combine multiple orf calls into one gff
    # by calling an external script
    #
    my ($gff_files, $tmp_folders_ref, $min_len, $tmp_fasta_ref) = @_;
    my @gffs = split /,/, $gff_files;
    foreach my $current_folder (keys %{$tmp_folders_ref})
    {

        # string of gff files
        my $gff_str = "";
        foreach my $gff (@gffs)
        {
            my $gff_path = catfile( $current_folder, basename($gff) );
            $gff_str .= $gff_path.",";

            # If that contig had no annotations, the file does not exist. Just create an empty one.
            touch($gff_path) if (not -e $gff_path);
        }
        $gff_str =~ s/,$//; # take off the last comma
        
        # run the external script!
        my $sequence_file = catfile( $current_folder, "sequence.fa" );
        my $unknowns_file = catfile( $current_folder, "unknowns.fa" );
        my $combined_file = catfile( $current_folder, "combined.gff3" );
        my $cmd = "combineGffOrfs -c $sequence_file -g $gff_str -o $combined_file -a $unknowns_file -m $min_len";
        run($cmd);

        # move the unknowns onto the pile
        concat( $unknowns_file, $$tmp_fasta_ref );
    }
}


sub blastWorker {
    my ($in_fasta, $out_blast, $db, $new_blast, $blast_program) = @_;
    if( $new_blast) {
        run("$blast_program -query $in_fasta -db $db -out $out_blast -outfmt 6");
    } else {
        run("blastall -p $blast_program -i $in_fasta -d $db -o $out_blast -m 8");
    }
}


sub blastUnknowns {
    #-----
    # BLAST unknowns against the Uniref or Nr protein database
    #

    # Skip BLAST if a BLAST result file was provided
    my ($blast_sims, $tmp_fasta_ref, $tmp_fasta_prefix_ref, $threads, $blast_program, $db_ref, $new_blast) = @_;
    return 1 if $blast_sims;

    # BLAST them
    my $num_seq = count_fasta_sequences($tmp_fasta_ref);
    print "$num_seq sequence(s) to BLAST in ${$tmp_fasta_ref}\n";
    my $blast_file = "${$tmp_fasta_prefix_ref}.$blast_program";
    if ($threads > 1)
    {
        my $num_seq_per_file = int ($num_seq / $threads);
        my $seqio_global = Bio::SeqIO->new(-file => ${$tmp_fasta_ref}, -format => 'fasta');
        print "Splitting ${$tmp_fasta_ref} into $threads parts, $num_seq_per_file sequence(s) per file\n";

        # Open output files
        my @out_fhs;
        my @out_files;
        my @tmp_blast_files;
        for my $i (1 .. $threads)
        {
            my $tmp_fasta = ${$tmp_fasta_prefix_ref}.'_'.$i.'.fna';
            push @out_files, $tmp_fasta;
            my $seqio_out = Bio::SeqIO->new(
                -file   => '>'.$tmp_fasta,
                -format => 'fasta',
                -flush  => 0, # go as fast as we can!
            );
            push @out_fhs, $seqio_out;
            my $tmp_blast = ${$tmp_fasta_prefix_ref}.'_'.$i.'.'.$blast_program;
            push @tmp_blast_files, $tmp_blast;
        }

        # Distribute sequences equally in output files
        my $i = 1;
        while (my $seq = $seqio_global->next_seq())
        {
            my $out_fh = $out_fhs[$i-1];
            $out_fh->write_seq($seq);
            $i = ($i == $threads) ? 1 : $i + 1;
        }

        # Close output files
        for my $out_fh (@out_fhs)
        {
            $out_fh->close;
        }

        # BLAST the sequences
        for my $i (1 .. scalar @out_files) 
        {
            print "Spawning BLAST worker $i (out of $threads)\n";
            my $q = $out_files[$i-1];
            my $o = $tmp_blast_files[$i-1];
            threads->new(\&blastWorker, $q, $o, ${$db_ref}, $new_blast, $blast_program );
        }
            
        $_->join() for threads->list();

        # Put all blast results in a unique file
        for my $tmp_blast_file (@tmp_blast_files)
        {
            concat($tmp_blast_file, $blast_file);
        }

    }
    else
    {
        blastWorker(${$tmp_fasta_ref}, $blast_file, ${$db_ref}, $new_blast, $blast_program);
    }
}


sub splitBlastResults {
    # now split BLAST Results across multiple folders...
    my ($blast_sims, $tmp_fasta_prefix_ref, $blast_program, $tmp_prefix_ref, $tmp_folders_hash_ref) = @_;
    my $current_dir_name = "__DOOF";
    my $current_file_handle;
    my $blast_results;
    
    # Determine input file
    if ($blast_sims)
    {
        open $blast_results, "<", $blast_sims or die "Error: Could not read file $blast_results\n$!\n";
    }
    else
    {
        my $blast_file = "${$tmp_fasta_prefix_ref}.$blast_program";
        open $blast_results, "<", $blast_file or die "Error: Could not read file $blast_file\n$!\n";
    }

    # Let the splitting occur
    my $unknown_filename = "unknowns.$blast_program";
    while(<$blast_results>)
    {
        # split the line, we need to know where to put this guy
        my @bits = split /\t/, $_;
        my @bits_bits = split /_/, $bits[0];
        
        # all but the last 2 underscores are the name!
        my $dir_name = join "_", @bits_bits[0..$#bits_bits-2];
        
        # time to write to a new file!
        if($dir_name ne $current_dir_name)
        {
            if("__DOOF" ne $current_dir_name)
            {
                # we have a file to close
                close $current_file_handle;
            }

            my $unknown_filepath = catfile( ${$tmp_prefix_ref}.'_'.$dir_name, $unknown_filename );
            open $current_file_handle, ">", $unknown_filepath or die "Error: Could not write file $unknown_filepath\n$!\n";
            $current_dir_name = $dir_name;
        }
        print $current_file_handle $_;
    }
    if("__DOOF" ne $current_dir_name)
    {
        # we have a file to close
        close $current_file_handle;
    }


    # If there are not BLAST hits, just create an empty file
    for my $tmp_folder (keys %{$tmp_folders_hash_ref}) {
        my $unknown_filepath = catfile( $tmp_folder, $unknown_filename );
        if (not -e $unknown_filepath) {
            touch($unknown_filepath);
        }
    }

}


sub touch {
    my ($file) = @_;
    open my $fh, '>', $file or die "Error: Could not write file $file\n$!\n";
    close $fh;
    return 1;
}


sub annotate {
    my ($blast_program, $seq_embed, $u2a_file, $u2a_hash_ref, $tmp_folder_ref,
        $contig_file, $evalue, $max_results_cutoff, $ann_hash_ref, $output_file_name ) = @_;

    # load the databases
    &loadU2A($u2a_file, $u2a_hash_ref);
    
    my $counter = 0;
    my $big_counter = 0;
    # first do the call out
    foreach my $current_folder (keys %{$tmp_folder_ref})
    {
        
        # get the blast results

        my $pattern = catfile( $current_folder, "*.$blast_program" );
        my @blast_files = glob $pattern;

        # Even if there were no BLAST hits, we created an empty file
        if (scalar @blast_files == 0)
        {
            die "Error: No BLAST file found in $current_folder\n";
        }

        # parse the blast results and report the predictions
        &generateAnnotations(\@blast_files, $u2a_hash_ref, $evalue, $max_results_cutoff, $ann_hash_ref);

        # insert these new values into the GFF3 file
        my $combined_file  = catfile( $current_folder, "combined.gff3"  );
        my $annotated_file = catfile( $current_folder, "annotated.gff3" );
        &recombineGff3($combined_file, $annotated_file, $ann_hash_ref);

        $counter++;
        if($counter == 10)
        {
            $big_counter += $counter;
            print "Parsed: $counter ... \n";
            $counter = 0;
        }
    }
    
    # then recombine it all!
    print "Recombining results\n";
    open my $out_fh, ">", $output_file_name or die "Error: Could not write file $output_file_name\n$!\n";
    print $out_fh "##gff-version 3\n";
    foreach my $current_folder (keys %{$tmp_folder_ref})
    {
        my $annotated_file = catfile( $current_folder, "annotated.gff3" );
        open my $this_gff3, '<', $annotated_file or die "Error: Could not read file $annotated_file\n$!\n";
        while(<$this_gff3>)
        {
            next if $_ =~ /^#/;
            print $out_fh $_;
        }
        close $this_gff3;
    }

    # finally, embed the FASTA sequences in the GFF file
    if ($seq_embed) {
        print $out_fh "##FASTA\n";
        my $in  = Bio::SeqIO->new( -file => $contig_file  , -format => 'fasta' ); 
        my $out = Bio::SeqIO->new( -fh   => $out_fh, -format => 'fasta' );
        while (my $seq = $in->next_seq)
        {
            $out->write_seq($seq);
        }
        $in->close;
        $out->close;
    }

    close $out_fh;
}


sub cleanTmps {
    #-----
    # clean up all the tmp files we added
    #
    my ($blast_program, $keep_bx, $tmp_folders_ref, $tmp_fasta_ref, $option_sims, $threads, $tmp_fasta_prefix) = @_;
    foreach my $current_folder (keys %{$tmp_folders_ref})
    {
        remove_tree( $current_folder );
    }

    unlink ${$tmp_fasta_ref} or die "Error: Could not delete file ${$tmp_fasta_ref}\n$!\n";
   
    if ( (not $option_sims) && ($threads > 1) )
    {
        for (my $i = 1; $i <= $threads; $i++)
        {
            my $fasta_chunk = $tmp_fasta_prefix.'_'.$i.'.fna';
            unlink $fasta_chunk or die "Error: Could not delete file $fasta_chunk\n$!\n";
            my $other = $tmp_fasta_prefix.'_'.$i.'.'.$blast_program;
            unlink $other or die "Error: Could not delete file $other\n$!\n";
        }
    }
    
    if($keep_bx == 0)
    {
       my $file = "$tmp_fasta_prefix.$blast_program";
       unlink $file or die "Error: Could not delete file $file\n$!\n";
    }
}

sub loadU2A() {
    #-----
    # load the annotation association file.
    # 
    # File must look like this:
    # 
    # Q17FA7^inNOG11764  aag:AaeL_AAEL003451 
    # Q0K0A5^COG0583  Transcriptional regulator   reh:H16_B1787   
    # D3FA93^NO_COG   cwo:Cwoe_0729   
    # Q1GAF8^COG1028  Dehydrogenases with different specificities (related to short-chain alcohol dehydrogenases) ldb:Ldb0903 path:ldb00061   ko:K00059   
    
    my ($Uniprot2ANN_file, $U2A_hash_ref) = @_;
    
    
    open my $U2A_fh, "<", $Uniprot2ANN_file or die "Error: Could not read U2A file $Uniprot2ANN_file\n$!\n";
    while (<$U2A_fh>)
    {
        chomp $_;
        my @data = split(/\^/, $_);
        $U2A_hash_ref->{$data[0]} = $data[1];
    }
    close $U2A_fh;
    print scalar(keys(%{$U2A_hash_ref}))." records loaded\n";
}


sub generateAnnotations() {
    #-----
    # parse the blast results and report the predictions
    #
    my ($blast_files_ref, $U2A_hash_ref, $evalue, $max_results_cutoff, $ann_hash_ref) = @_;

    foreach my $blast_file (@{$blast_files_ref}) 
    {
        # Load the result into memory
        my $in = new Bio::SearchIO(
            -format => 'blastTable',
            -file => $blast_file
        );

        # get all the accessions with a significance less than $global_evalue_cutoff
        while( my $result = $in->next_result ) 
        {
            my $tmp_ann = "__DOOF";
            my $num_hits_done = 0;
            my %ogs_done = ();
            while( my $hit = $result->next_hit ) 
            {
                # Assume that the BLAST results are sorted by increasing E-value
                # i.e. that the first hit we encounter is the best one.
                last if $num_hits_done > $max_results_cutoff;
                my $hit_name = $hit->name;
                $hit_name =~ s/UniRef90_([^ ]*).*/$1/;
                next if $hit_name eq 'Unknown';
                next if $hit->significance > $evalue;
                next if not exists $U2A_hash_ref->{$hit_name};
                $tmp_ann = $U2A_hash_ref->{$hit_name};
                $num_hits_done++;
            }
            $ann_hash_ref->{$result->query_name} = gff_escape_chars($tmp_ann);
        }
    }
}

sub gff_escape_chars {
    # GFF3 basics:
    #    = separates tag from values: tag1=value
    #    ; separates tag and value pairs: tag1=value1;tag2=value2
    #    , separates multiple values that a tag can have: tag2=value2,value3
    # As a consequence =;, need to be escaped in hexadecimal if their occur in
    # a tag or value:
    #    ; (semicolon) - %3B
    #    = (equals) - %3D
    #    , (comma) - %2C
    # Here, we deal with the ANN mapping file, which uses unescaped commas in
    # tag values, e.g.:
    #    A8MM54^;Ontology_term=KEGG_ENZYME:ftsX, cell division transport permease
    # To get valid GFF3 files, we need to escape commas.
    # Getting rid of the commas is easier and more fun!
    my $ann_string = shift;
    $ann_string =~ s/,/ /g;
    return $ann_string;
}
sub gff_collapse_tags {
    # Collapse multiple identical tags into a single one, e.g. transform:
    #    tag2=value2;tag2=value3
    # into:
    #    tag2=value2,value3
    my $string = shift;
    my @pairs = split ';', $string;
    my %tag_hash;
    $string = '';
    for my $pair (@pairs) {
       my ($tag, $value) = split '=', $pair;
       if ( (defined $value) and ($value !~ m/^\s*$/) ) {
          if (exists $tag_hash{$tag}) {
             $tag_hash{$tag} .= ',';
          }
          $tag_hash{$tag} .= $value;
       }
    }
    while ( my ($tag, $values) = each %tag_hash ) {
       $string .= ';' if $string;
       $string .= "$tag=$values";
    }
    return $string;
}


sub recombineGff3() {
    #-----
    # blend the new annotations into the old
    #
    my ($gff3_in, $gff3_out, $ann_hash_ref ) = @_;

    # specify input via -fh or -file
    my $gffio_in = Bio::Tools::GFF->new(-file => $gff3_in, -gff_version => 3);
    open my $ann_fh, '>', $gff3_out or die "Error: Could not write file $gff3_out\n$!\n";
    print $ann_fh "##gff-version 3\n";
    
    # loop over the input stream
    while(my $feat = $gffio_in->next_feature()) 
    {
        my $gff_string = $feat->gff_string($gffio_in);
        my @gff_bits = split /\t/, $gff_string;
        my $feat_key = $gff_bits[0];#."_".$feat->start."_".$feat->end;
        my $annotation = $ann_hash_ref->{$feat_key};
        if( (defined $annotation) && ($annotation ne '__DOOF') )
        {
            # this is a bit dodge... earlier, I put __DOOF in place of any null annotations
            $gff_bits[8]  =~ s/: hypothetical protein//;
            $gff_bits[8] .= ';' if $annotation !~ m/^;/;
            $gff_bits[8] .= $annotation;
            $gff_bits[8]  = gff_collapse_tags($gff_bits[8]);
            $gff_string   = join "\t", @gff_bits;
        }
        print $ann_fh "$gff_string\n";
    }
    # clean up
    $gffio_in->close();
    close $ann_fh;
}


sub createFlatFile {
    print "generating genbank files for contigs...";
    my ($contig_file_name, $output_file_ref) = @_;
    my $cmd = run("gff2genbank $contig_file_name ${$output_file_ref}");
    print "done\n";
}


sub run {
    # Run a command, check that it completed successfully and return its output
    my ($cmd) = @_;
    my $results = `$cmd`;
    # In theory, a return status of -1 is an error, but in practice, other values
    # are also errors. A return status of 0 seems to be ok though.
    die "Error: Command '$cmd' failed with return status $?\n$!\n" if ($? != 0);
    return $results;
}


sub concat {
    # Concatenate or append file content (a single file name or an
    # arrayref of filenames $in_files) into another file (scalar
    # $out_file). If mode is '>', create/overwrite the output file,
    # but if mode is '>>' (default), create/append. If del_input is 0,
    # delete the input files after concatenation.
    my ($in_files, $out_file, $mode, $del_input) = @_;
    if (not defined $mode) {
        $mode = '>>';
    }
    if (not defined $del_input) {
        $del_input = 0;
    } 
    if ($mode !~ /^>{1,2}$/) {
        die "Error: Invalid cat mode $mode\n";
    }
    if (not ref $in_files) {
        # Put single file into an array
        $in_files = [ $in_files ];
    }
    open my $ofh, $mode, $out_file or die "Error: Could not write file $out_file\n$!\n";
    for my $in_file (@$in_files) {
        open my $ifh, '<', $in_file or die "Error: Could not open file $in_file\n$!\n";
        while (my $line = <$ifh>) {
            print $ofh $line;
        }
        close $ifh;
        if ($del_input == 1) {
            unlink $in_file or die "Error: Could not remove file $in_file\n$!\n";
        }
    }
    close $ofh;
    return 1;
}


sub count_fasta_sequences {
    my ($fasta_file) = @_;
    my $count = 0;
    open my $fh, '<', $fasta_file or die "Error: Could not read file $fasta_file\n$!\n";
    while (my $line = <$fh>) {
        $count++ if $line =~ m/^>/;
    }
    close $fh;
    return $count;
}


##################################################################
# Following Code originally from combineGffOrfs.pl
##################################################################
struct GffObject =>
{
    GO_prevGO_ref => '$',   # next / prev in list
    GO_nextGO_ref => '$',       
    GO_gffRef =>'$',        # reference to the gff object
    GO_start => '$',        # we use these A LOT, so just save them twice...
    GO_end => '$',
    GO_next_start => '$',   # how much "free" space lies beyond this end...? Where does the next orf start?
};

sub debugFeature {
    # display why a feature was kept or not
    my ($msg, $details) = @_;
    print "Debug: ".$msg.".\n";
    print "Reason: $details\n" if defined $details;
    return 1;
}


sub feat2str {
    my $feat = shift;
    return $feat->source_tag." ".$feat->primary_tag." ".$feat->start."-".$feat->end;
}


sub insertFeature {
    #-----
    # insert a feature, provided that it does not displace any existing feature...
    # start_ref: the list of features
    # feat_ref: the feature that we try to put in the list
    my ($start_ref, $feat_ref, $gff_list_ref, $debug, $shared_olap_cutoff, $end_ref) = @_;

    my $feat = ${$feat_ref};
    my $feat_start  = $feat->start;
    my $feat_end    = $feat->end;
    my $feat_length = $feat->length; # start <= end , always
    my $first_loop_done = 0;
    my $list_handle_ref = $start_ref;
    do
    {
        my $current_node = ${$list_handle_ref};
        my $current_start = $current_node->GO_start;
        my $current_end = $current_node->GO_end;
        my $current_next_start = $current_node->GO_next_start;

        if($gff_list_ref == $list_handle_ref)
        {
            $current_start = 0;
            $current_end = 0;
        }
        
        # sometimes we start a little too far in...
        if($first_loop_done == 0)
        {
            if($current_start >= $feat_end)
            {
                my $prev_ref = $current_node->GO_prevGO_ref;
                my $prev_node = ${$prev_ref};
                $current_next_start = $current_start;
                $current_node = $prev_node;
                $list_handle_ref = \$current_node;
                if($gff_list_ref == $prev_ref)
                {
                    # the prev node is the list head...
                    $current_end = 0;
                    $current_start = 0;
                }
                else
                {
                    # things are relatively normal...
                    $current_start = $current_node->GO_start;
                    $current_end = $current_node->GO_end;
                }
            }
            $first_loop_done = 1;
        }
        
        # if we haven't inserted it by the time we have moved past it,
        # there is no need to keep searching the list
        if($current_start > $feat_start)
        {
            debugFeature("Not keeping feature ".feat2str($$feat_ref), 'no need to keep searching the list') if $debug;
            return $start_ref;
        }
        
        # so, $current_start is less than or equal to
        # the start of the feature being added...
        if($current_end >= $feat_end)
        {
            # the feature to be added is completely contained
            # in an existing feature. no go!
            debugFeature("Not keeping feature ".feat2str($$feat_ref), 'entirely contained in another feature') if $debug;
            return $start_ref;
        }

        # we would like to know if the new feature overlaps with the 
        # end of the current feature
        my $should_add = 0;
        if($feat_end > $current_next_start)
        {
            if($feat_start < $current_end) 
            { 
                # new feature overlaps with only the next feature
                #
                #  CURRENT             NEXT
                # XXXXXXXXX          XXXXXXX
                #        NNNNNNNNNNNNN
                if( ( ( ($feat_end - $current_next_start) + ($current_end - $feat_start) ) / $feat_length ) < $shared_olap_cutoff ) { $should_add = 1; }
                
            }
            else
            {
                # new feature overlaps with two existing features
                #
                #  CURRENT             NEXT
                # XXXXXXXXX          XXXXXXX
                #            NNNNNNNNNN
                if( ( ($feat_end - $current_next_start) / $feat_length ) < $shared_olap_cutoff) { $should_add = 1; }
            }
        }
        else
        {
            if($feat_start < $current_end) 
            { 
                # new feature overlaps with only the current feature
                #
                #  CURRENT             NEXT
                # XXXXXXXXX          XXXXXXX
                #        NNNNNNNNNN
                if( ( ($current_end - $feat_start) / $feat_length ) < $shared_olap_cutoff) { $should_add = 1; }
            }
            else
            {
                # new feature floats in the ether...
                #
                #  CURRENT             NEXT
                # XXXXXXXXX           XXXXXX
                #           NNNNNNN
                $should_add = 1;
            }
        }
        if($should_add == 1)
        {
            debugFeature("Keeping feature ".feat2str($$feat_ref)) if $debug;


            my $new_obj = GffObject->new();
            $new_obj->GO_gffRef($feat_ref);
            $new_obj->GO_start($feat_start);
            $new_obj->GO_end($feat_end);
            insertInList($list_handle_ref, \$new_obj);
            # return a reference to the inserted object to speed up the next insertion...
            return \$new_obj;
        }
    } while(nextInList(\$list_handle_ref, $end_ref ) == 1);

    # nothing!
    return $gff_list_ref;
}


sub insertInList {
    #-----
    # Given a node in the LL and a ref to a node to add
    # Insert the new ref AFTER the given LL ref
    #
    my ($ll_ref, $new_ref) = @_;

    my $ll_obj = ${$ll_ref};
    my $new_obj = ${$new_ref};
    
    my $tmp_nxt_ref = $ll_obj->GO_nextGO_ref;
    my $tmp_nxt_obj = ${$tmp_nxt_ref};
    
    $ll_obj->GO_nextGO_ref($new_ref);
    $ll_obj->GO_next_start($new_obj->GO_start);
    $tmp_nxt_obj->GO_prevGO_ref($new_ref);
    
    $new_obj->GO_nextGO_ref($tmp_nxt_ref);
    $new_obj->GO_prevGO_ref($ll_ref);
    $new_obj->GO_next_start($tmp_nxt_obj->GO_start);
}


sub nextInList {
    #-----
    # given a ref to a linked list node, update to reflect the next node in the list
    # return 1 or 0 accordingly
    #
    my ($ll_ref_ref, $end_ref) = @_;
    my $ll_ref = ${$ll_ref_ref};
    my $ll_obj = ${$ll_ref};
    if($ll_obj->GO_nextGO_ref == $end_ref) { return 0; }
    # not at the end of the list!
    $$ll_ref_ref = $ll_obj->GO_nextGO_ref;
    return 1;
}


1;

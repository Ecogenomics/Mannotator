#! /usr/bin/env perl

###############################################################################
#
#    mannotator.pl
#    
#    Annotate a genome!
#
#    Copyright (C) 2011 Michael Imelfort
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

#pragmas
use strict;
use warnings;

#perl modules
use File::Path qw( remove_tree make_path ) ;
use File::Basename;
use File::Spec::Functions qw( :ALL ); # catfile
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SearchIO;
use threads;
use Data::Dumper;
#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
printAtStart();
my $options = checkParams();

######################################################################
# CODE HERE
######################################################################
# globals
my %global_U2A_hash = ();
my %global_annotations_hash = ();

my %global_tmp_folders = ();
my $global_tmp_fasta = "mannotator_unknowns_".time.".fasta";

my $global_output_file = "mannotatored.gff3";
if(exists $options->{'out'}) { $global_output_file = $options->{'out'}; }

# evalue cutoff
my $global_evalue_cutoff = 0.000000001;
if(exists $options->{'evalue'}) { $global_evalue_cutoff = $options->{'evalue'}; }

# results cut off
my $global_results_max_cut_off = 1; # use the top hit only

#blast program to run
my $blast_program = "blastx";
if (exists $options->{'blast_prg'}) { $blast_program = $options->{'blast_prg'}; }

my $threads = 1;
if (exists $options->{'threads'}) {$threads = $options->{'threads'}; }

#
# Step 1. Split down the GFF3 and FASTA files into multiple folders
#
unless ($options->{'one'})
{
    print "Step 1: Splitting GFF and FASTA files into temporary directories\n";
    splitGffs();
    splitFasta();
}

#
# Step 2. For each folder, combine the gff files and produce a list of sequences to be blasted!
#
unless($options->{'two'})
{
    print "Step 2: Combining Gffs\n";
    combineGffs();
}

#
# Step 3. Blast the unknowns against UniRef or Nr proteic database
#
unless ($options->{'three'})
{
    print "Step 3: Blasting unknown ORFs against protein database... This could take some time. Perhaps you need a coffee?\n";
    blastUnknowns();
}

#
# Step 4. Splitting the BLAST results
#
unless ($options->{'four'})
{
    print "Step 4: Splitting BLAST results\n";
    splitBlastResults();
}

#
# Step 5. Annotate the GFFs using BLAST results and re-combine!
#
unless ($options->{'five'})
{
    print "Step 5: Annotating using the BLAST results\n";
    annotate();
}

#
# Finally, remove all the temp directories (if we need to)
#
if(exists $options->{'blastx'})
{
    print "Removing all but the blastx file\n";
    cleanTmps(1);
}
elsif(!exists $options->{'keep'})
{
    print "Cleaning up tmp folders\n";
    cleanTmps(0);
}

if(exists $options->{'flatfile'})
{
    createFlatFile();
}


######################################################################
# CUSTOM SUBS
######################################################################

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


sub splitGffs {
    #----
    # split the gffs into multiple files, one for each sequence
    # and place them in separate directories
    #
    my @gffs = split /,/, $options->{'gffs'};
    foreach my $gff (@gffs)
    {
        print "Parsing GFF file $gff\n";
        my $current_fh;
        my $current_fasta_header = "__DOOF";
        my $true_fasta_header = undef;
        open my $gff_fh, "<", $gff or die "Error: Could not read file $gff\n$!\n";
        my $gff_def = <$gff_fh>;
        while(<$gff_fh>)
        {

            # ##FASTA section indicates end of ##gff section
            if ($_ =~ m/^##FASTA/i) {
               last;
            }

            # comment lines split sequences
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
                if("__DOOF" ne $current_fasta_header)
                {
                    # this is not the first run.
                    # so we have a file to close...
                    close $current_fh;
                }
               
                $current_fasta_header = $fasta_header;
                  
                # make a new file and directory
                make_path( $current_fasta_header );
                $global_tmp_folders{$current_fasta_header} = 1;
                    
                # make a new file handle
                my $out_file = catfile( $current_fasta_header, basename($gff) );
                open $current_fh, ">", $out_file or die "Error: Could not write file $out_file\n$!\n";
                    
                # print back in the header information
                print $current_fh $gff_def;
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
    # Folders should be made by now and should be named according to the 
    # fasta_headers...
    #
    my $seqio_object = Bio::SeqIO->new(-file => $options->{'contigs'}, -format => "fasta");
    while (my $seq = $seqio_object->next_seq)
    {
        my $seq_fn = catfile( $seq->display_id(), "sequence.fa" );
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
    my @gffs = split /,/, $options->{'gffs'};
    foreach my $current_folder (keys %global_tmp_folders)
    {

        my $gff_str = "";
        foreach my $gff (@gffs)
        {
            $gff_str .= catfile( $current_folder, basename($gff) ) . ",";
        }

        # take off the last comma
        $gff_str =~ s/,$//;
        
        # run the script!
        my $sequence_file = catfile( $current_folder, "sequence.fa" );
        my $unknowns_file = catfile( $current_folder, "unknowns.fa" );
        my $combined_file = catfile( $current_folder, "combined.gff3" );
        run("combineGff3.pl -c $sequence_file -g $gff_str -o $combined_file -a $unknowns_file");

        # move the unknowns onto the pile
        concat( $unknowns_file, $global_tmp_fasta );
    }
}


sub worker {
    my ($in_fasta, $out_blast) = @_;
    run("blastall -p $blast_program -i $in_fasta -d $options->{'protdb'} -o $out_blast -m 8");
}


sub blastUnknowns {
    #-----
    # BLAST unknowns against the Uniref or Nr protein database
    #

    # Skip BLAST if a BLAST result file was provided
    return 1 if $options->{'sims'};

    # BLAST them
    my $num_seq = count_fasta_sequences($global_tmp_fasta);
    print "$num_seq sequences to BLAST in $global_tmp_fasta\n";
    if ($threads > 1)
    {
        my $num_seq_per_file = int ($num_seq / $threads);
        my $seqio_global = Bio::SeqIO->new(-file => $global_tmp_fasta, -format => 'fasta');
        print "Splitting $global_tmp_fasta into $threads parts, $num_seq_per_file sequence(s) per file\n";

        # Open output files
        my @out_fhs;
        for my $i (1 .. $threads)
        {
            my $seqio_out = Bio::SeqIO->new(
                -file   => '>'.$global_tmp_fasta.'_'.$i,
                -format => 'fasta',
                -flush  => 0, # go as fast as we can!
            );
            push @out_fhs, $seqio_out;
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
        my @blast_files;
        for my $i (1 .. $threads)
        {
            print "Spawning BLAST worker $i (out of $threads)\n";
            my $q = $global_tmp_fasta."_".$i;
            my $o = "$global_tmp_fasta.$i.$blast_program";
            threads->new(\&worker, $q, $o);
            push @blast_files, $o;
        }
            
        $_->join() for threads->list();

        # Put all blast results in a unique file
        for my $blast_file (@blast_files)
        {
            concat( $blast_file, "$global_tmp_fasta.$blast_program" );
        }

    }
    else
    {
        worker($global_tmp_fasta, "$global_tmp_fasta.$blast_program");
    }
}


sub splitBlastResults {
    # now split BLAST Results across multiple folders...
    my $current_dir_name = "__DOOF";
    my $current_file_handle;
    my $blast_results;
    
    # Determine input file
    if ($options->{'sims'})
    {
        open $blast_results, "<", $options->{'sims'} or die "Error: Could not read file $blast_results\n$!\n";
    }
    else
    {
        my $blast_file = "$global_tmp_fasta.$blast_program";
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
            my $unknown_filepath = catfile( $dir_name, $unknown_filename );
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
    for my $tmp_folder (keys %global_tmp_folders) {
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
    #-----
    # call out to blast2ann.pl
    #
    my ($seq_embed) = shift;

    # load the databases
    &loadU2A($options->{'i2a'});
    
    my $counter = 0;
    my $big_counter = 0;
    # first do the call out
    foreach my $current_folder (keys %global_tmp_folders)
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
        &generateAnnotations(@blast_files);

        # insert these new values into the GFF3 file
        my $combined_file  = catfile( $current_folder, "combined.gff3"  );
        my $annotated_file = catfile( $current_folder, "annotated.gff3" );
        &recombineGff3($combined_file, $annotated_file);

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
    open my $out_fh, ">", $global_output_file or die "Error: Could not write file $global_output_file\n$!\n";
    print $out_fh "##gff-version  3\n";
    foreach my $current_folder (keys %global_tmp_folders)
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
    if (exists $options->{'seq_embed'}) {
        print $out_fh "##FASTA\n";
        my $file = $options->{'contigs'};
        my $in  = Bio::SeqIO->new( -file => $file  , -format => 'fasta' ); 
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
    my ($keep_bx) = @_;
    foreach my $current_folder (keys %global_tmp_folders)
    {
        remove_tree( $current_folder );
    }

    unlink $global_tmp_fasta or die "Error: Could not delete file $global_tmp_fasta\n$!\n";
   
    if ( (not $options->{'sims'}) && ($threads > 1) )
    {
        for (my $i = 1; $i <= $threads; $i++)
        {
            my $global_fasta_chunk = $global_tmp_fasta."_".$i;
            unlink $global_fasta_chunk or die "Error: Could not delete file $global_fasta_chunk\n$!\n";
            my $other = "$global_tmp_fasta.$i.$blast_program";
            unlink $other or die "Error: Could not delete file $other\n$!\n";
        }
    }
    
    if($keep_bx == 0)
    {
       my $file = "$global_tmp_fasta.$blast_program";
       unlink $file or die "Error: Could not delete file $file\n$!\n";
    }
}


# stolen from blast2ann.pl
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
    
    my ($Uniprot2ANN_file) = @_;
    
    print "Loading ID to annotations mapping file... ";
    
    open my $U2A_fh, "<", $Uniprot2ANN_file or die "Error: Could not read U2A file $Uniprot2ANN_file\n$!\n";
    while (<$U2A_fh>)
    {
        chomp $_;
        my @data = split(/\^/, $_);
        $global_U2A_hash{$data[0]} = $data[1];
    }
    close $U2A_fh;
    print scalar(keys(%global_U2A_hash))." records loaded\n";
}


sub generateAnnotations() {
    #-----
    # parse the blast results and report the predictions
    #
    my (@blast_files) = @_;
    
    foreach my $blast_file (@blast_files) 
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
                last if $num_hits_done > $global_results_max_cut_off;
                my $hit_name = $hit->name;
                $hit_name =~ s/UniRef90_([^ ]*).*/$1/;
                next if $hit_name eq 'Unknown';
                next if $hit->significance > $global_evalue_cutoff;
                next if not exists $global_U2A_hash{$hit_name};
                $tmp_ann = $global_U2A_hash{$hit_name};
                $num_hits_done++;
            }
            $global_annotations_hash{$result->query_name} = gff_escape_chars($tmp_ann);
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
    my ($gff3_in, $gff3_out) = @_;

    # specify input via -fh or -file
    my $gffio_in = Bio::Tools::GFF->new(-file => $gff3_in, -gff_version => 3);
    open my $ann_fh, ">", $gff3_out or die "Error: Could not write file $gff3_out\n$!\n";
    print $ann_fh "##gff-version 3\n";
    
    # loop over the input stream
    while(my $feat = $gffio_in->next_feature()) 
    {
        my $gff_string = $feat->gff_string($gffio_in);
        my @gff_bits = split /\t/, $gff_string;
        my $feat_key = $gff_bits[0]."_".$feat->start."_".$feat->end;
        my $annotation = $global_annotations_hash{$feat_key};
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
    my $cmd = run("gff2genbank.pl $options->{'c'} $global_output_file");
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


######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = (
         "help|h+",
         "gffs|g:s",
         "keep|k+",
         "contigs|c:s",
         "i2a|i:s",
         "protdb|p:s",
         "out|o:s",
         "blastx|x:+",
         "sims|s:s",
         "evalue|e:s",
         "blast_prg|p:s",
         "flatfile|f:+",
         "threads|t:s",
         "one|1+",
         "two|2+",
         "three|3+",
         "four|4+",
         "five|5+",
         "seq_embed|d",
    );
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};
    
    if(!exists $options{'gffs'})
    {
        print "ERROR: You need to input some gff3 files to continue!\n";
        exec("pod2usage $0");
    }

    if(!exists $options{'contigs'})
    {
        print "ERROR: You need to give me the location of the contigs to continue!\n";
        exec("pod2usage $0");
    }

    if(!exists $options{'protdb'})
    {
        print "ERROR: You need to give me the location of a proteic BLAST database (UniRef or Nr) to continue!\n";
        print "Perhaps, try this one (on EURY): /Volumes/Biodata/BLAST_databases/UniProt/UniRef90/uniref90_micro.fasta\n";
        exec("pod2usage $0");
    }
    
    if(!exists $options{'i2a'})
    {
        print "ERROR: You need to give me the location of the ID to annotations mapping file!\n";
        print "Perhaps, try this one (on EURY): /Volumes/Biodata/ANN_mappings/ANN_mappings.txt\n";
        exec("pod2usage $0");
    }

    return \%options;
}

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2011 Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    mannotator.pl

=head1 COPYRIGHT

   copyright (C) 2011 Michael Imelfort

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

    
    STEPS:

    BEFORE YOU RUN THE MANOTATOR!
    
    1a. Send entire genome out to RAST server.
    1b. Run entire genome though prodigal
    
    You now have 2 gff3 files...

    NOW, UNLEASH THE MANOTATOR!


=head1 SYNOPSIS

    mannotator.pl -gffs|g GFF_FILE1[,GFF_FILE2...] -contigs|c FILE -blast_prg|p BLAST_TYPE -protdb|p LOCATION -i2a|i FILE

       -gffs      -g FILE[,FILE]   List of gff3 files in order of trust!
       -contigs   -c FILE          FASTA file of contigs to be annotated...
       -protdb    -p LOCATION      Location of the UniRef or Nr BLAST database
       -i2a       -i FILE          ID to annotations mapping :file

      [-seq_embed -d           ]   Embed sequences into GFF files (useful to view the annotation in Artemis)     
      [-threads   -t INTEGER   ]   Number of blast jobs to run [default: 1]
      [-flatfile  -f           ]   Optionally create multiple genbank files for your contigs [default: do not create]
      [-blast_prg -b BLAST TYPE]   The type of blast to run [default: blastx]
      [-evalue    -e DECIMAL   ]   E-value cut off for blastx [default: 0.000000001]
      [-out       -o FILE      ]   Filename of the final gff3 file [default: mannotatored.gff3]
      [-keep      -k           ]   Keep all the tmp directories
      [-blastx    -x           ]   Keep only the blastx file (overrides -k option)
      [-sims      -s FILE      ]   Use BLAST similarities given in this file
      [-one                    ]   Skip step 1 of Mannotator
      [-two                    ]   Skip step 2 of Mannotator
      [-three                  ]   Skip step 3 of Mannotator
      [-four                   ]   Skip step 4 of Mannotator
      [-five                   ]   Skip step 5 of Mannotator
      [-help      -h           ]   Displays basic usage information
         
=cut


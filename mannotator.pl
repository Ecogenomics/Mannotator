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
use File::Basename;
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

#results cut off
my $global_results_max_cut_off = 1;

#blast program to run
my $blast_program = "blastx";
if (exists $options->{'blast_prg'}) { $blast_program = $options->{'blast_prg'}; }

my $threads = 1;
if (exists $options->{'threads'}) {$threads = $options->{'threads'}; }
#
# Step 1a. Split down the Gff3 files into multiple folders
#
unless ($options->{'one'})
{
print "Step 1a: Splitting gffs and making tmp directories\n";
splitGffs();

#
# Step 1b. Split down the fasta file and also put into the same folders
#
print "Step 1b: Splitting fasta\n";
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
unless ($options->{'four'})
{
print "Step 4: Splitting Blast results\n";
splitBlastResults();
}
#
# Step 4. Annotate the GFFs using blast results and re-combine!
#
unless ($options->{'five'})
{
print "Step 5: Annotating using the blast results\n";
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
sub splitGffs {
    #----
    # split the gffs into multiple files, one for each sequence
    # and place them in separate directories
    #
    my @gffs = split /,/, $options->{'gffs'};
    foreach my $gff (@gffs)
    {
        print "Parsing: $gff \n";
        my $current_fh;
        my $current_fasta_header = "__DOOF";
        open my $gff_fh, "<", $gff or die "Error: Could not read file $gff\n$!\n";
        my $gff_def = <$gff_fh>;
        while(<$gff_fh>)
        {
            # comment lines split sequences
            next if($_ =~ /^#/);
            my @bits = split /\t/, $_;
            if($bits[0] ne $current_fasta_header)
            {
                if("__DOOF" ne $current_fasta_header)
                {
                    # this is not the first run.
                    # so we have a file to close...
                    close $current_fh;
                }
                
                $current_fasta_header = $bits[0];
                  
                # make a new directory.
                run("mkdir -p $bits[0]");
                $global_tmp_folders{$bits[0]} = 1;
                    
                # make a new file handle
                my $out_file = $bits[0].'/'.basename($gff);
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
    # Simple script to split a fasta sequence ansd put into multiple folders
    # Folders should be made by now and should be named according to the 
    # fasta_headers...
    #
    my $seqio_object = Bio::SeqIO->new(-file => $options->{'contigs'}, -format => "fasta");
    while (my $seq = $seqio_object->next_seq)
    {
        my $seq_fn = $seq->display_id(). "/sequence.fa";
        if(open my $ffh, ">", $seq_fn)
        {
            print $ffh ">" . $seq->display_id() . "\n";;
            print $ffh $seq->seq(). "\n";
            close $ffh;
        }
    }
}

sub combineGffs {
    #-----
    # Wrapper to combine multiple orf calls nto one gff
    # calls external script
    #
    my @gffs = split /,/, $options->{'gffs'};
    foreach my $current_folder (keys %global_tmp_folders)
    {
        my $gff_str = "";
        foreach my $gff (@gffs)
        {
            $gff_str .= "$current_folder/$gff,";
        }
        
        # take off the last comma
        $gff_str =~ s/,$//;
        
        # run the script!
        run("combineGff3.pl -c $current_folder/sequence.fa -g $gff_str -o $current_folder/combined.gff3 -a $current_folder/unknowns.fa");
        
        # move the unknowns onto the pile
        run("cat $current_folder/unknowns.fa >> $global_tmp_fasta");

    }
}

sub worker {
	my ($chunk, $n) = @_;
	run("blastall -p $blast_program -i $chunk -d $options->{'protdb'} -o $global_tmp_fasta.$n.$blast_program -m 8");

}


sub blastUnknowns {
    #-----
    # Blast unknowns against the Uniref or Nr protein database
    #
    # first blast them
    my $num_seq = run("grep -c '>' $global_tmp_fasta");

    print "total sequences to blast: $num_seq\n";
    if ($threads > 1)
    {
    	my $num_seq_per_file = int ($num_seq / $threads);
    	my $seqio_global = Bio::SeqIO->new(-file => $global_tmp_fasta, -format => 'fasta');
    	print "splitting $global_tmp_fasta into $threads parts, $num_seq_per_file sequences per file\n";
    	my $i = 1;
    	my $j = 0;
    	while ($i < $threads)
    	{
    		# open a file to hold a chunk
                my $out_file = $global_tmp_fasta."_".$i;
		open (CH, ">", $out_file) or die "Error: Could not write file $out_file\n$!\n";
		while(my $fasta = $seqio_global->next_seq())
		{
              last if $j > $num_seq_per_file;
              print CH ">".$fasta->primary_id."\n".$fasta->seq()."\n";
     		$j++;
		}
            close CH;
            $i++;
    	}
        my $out_file = $global_tmp_fasta."_".$i;
    	open (CH, ">",$out_file) or die "Error: Could not write file $out_file\n$!\n";
    	while(my $fasta = $seqio_global->next_seq())
	{
        print CH ">".$fasta->primary_id."\n".$fasta->seq()."\n";
	}
   	close CH;

   	for (my $x = 1; $x <= $threads; $x++) 
   	{
     	print "spawning thread $i\n";
     	my $q = $global_tmp_fasta."_".$i;
     	threads->new(\&worker, $q, $i);
   	}
            
        $_->join() for threads->list();
   		
		
		for (my $i = 1; $i <= $threads; $i++)
		{
                        run("cat $global_tmp_fasta.$i.$blast_program >> $global_tmp_fasta.$blast_program");
		}
		
    }
    else
    {
    	run("blastall -p $blast_program -i $global_tmp_fasta -d $options->{'protdb'} -o $global_tmp_fasta.$blast_program -m 8");
    }
}


sub splitBlastResults {

    # now split them across multiple folders...
    my $current_dir_name = "__DOOF";
    my $current_file_handle;
    my $blast_results;
    
    if ($options->{'sims'})
    {
    open $blast_results, "<", $options->{'sims'} or die "Error: Could not read file $blast_results\n$!\n";
    }
    else
    {
        my $blast_file = "$global_tmp_fasta.$blast_program";
    	open $blast_results, "<", $blast_file or die "Error: Could not read file $blast_file\n$!\n";
    }
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
            my $unknown_file = "$dir_name/unknowns.$blast_program";
            open $current_file_handle, ">", $unknown_file or die "Error: Could not write file $unknown_file\n$!\n";
            $current_dir_name = $dir_name;
        }
        print $current_file_handle $_;
    }
    if("__DOOF" ne $current_dir_name)
    {
        # we have a file to close
        close $current_file_handle;
    }
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
        
        # get the blast results:
        my @blast_files = ();
        if (-d $current_folder) {
            opendir(INDIR, $current_folder)
            or die "Error: Could not read directory $current_folder\n$!\n";
            @blast_files = grep /\.$blast_program$/, readdir (INDIR);
            closedir (INDIR);
        }
        next if($#blast_files == -1);
		if (scalar @blast_files == 0)
		{
			warn "Warning: hmm... this is strange...\nthere does not seem to be any blast output files in $current_folder\n";
		}
        # parse the blast results and report the GO predictions
        &generateAnnotations($current_folder, @blast_files);
		
        # insert these new values into the GFF3 file
        &recombineGff3("$current_folder/combined.gff3", "$current_folder/annotated.gff3");

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
        if(open my $this_gff3, "<", "$current_folder/annotated.gff3")
        {
            while(<$this_gff3>)
            {
                next if $_ =~ /^#/;
                print $out_fh $_;
            }
            close $this_gff3;
        }
    }

    # finally, embed the FASTA sequences in the GFF file
    if (exists $options->{'seq_embed'}) {
        print $out_fh "##FASTA\n";
        my $file = $options->{'contigs'};
        my $in  = Bio::SeqIO->new( -file => $file  , -format => 'fasta' ); 
        my $out = Bio::SeqIO->new( -fh   => $out_fh, -format => 'fasta' );
        while (my $seq = $in->next_seq) {
            $out->write_seq($seq);
        }
        $in->close;
        $out->close;
    }

    close $out_fh;

}

sub cleanTmps {
    #-----
    # clean up all the tmp files we maded
    #
    my ($keep_bx) = @_;
    foreach my $current_folder (keys %global_tmp_folders)
    {
        run("rm -rf $current_folder");
    }
   
    run("rm $global_tmp_fasta");
    
    if ($threads > 1)
    {
    	for (my $i = 1; $i <= $threads; $i++)
    	{
    	my $global_fasta_chunk = $global_tmp_fasta."_".$i;
        run("rm $global_fasta_chunk $global_tmp_fasta.$i.$blast_program");
    	}
    }
    
    if(0 == $keep_bx)
    {
       run("rm $global_tmp_fasta.$blast_program");
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
    
    print "Loading ID to annotations mapping file...";
    
    open my $U2A_fh, "<", $Uniprot2ANN_file or die "Error: Could not read U2A file $Uniprot2ANN_file\n$!\n";
    while (<$U2A_fh>) {
        chomp $_;
        my @data = split(/\^/, $_);
        $global_U2A_hash{$data[0]} = $data[1];
    }
    close $U2A_fh;
    print "done.  Loaded ".keys(%global_U2A_hash)." ID to annotations records\n";
}

sub generateAnnotations() {
    #-----
    # parse the blast results and report the predictions
    #
    my ($current_folder, @blast_files) = @_;
    
    foreach my $blast_file (@blast_files) 
    {
        # Load the result into memory
        my $in = new Bio::SearchIO(-format => 'blastTable',
            -file => "$current_folder/$blast_file")
            or die "Error: Could not read file $blast_file\n$!\n";
        
        my $query_name = undef;
        
        # get all the accessions with a significance less than $global_evalue_cutoff
        while( my $result = $in->next_result ) 
        {
            my $tmp_ann = "__DOOF";
            my $num_hits_done = 0;
            my %ogs_done = ();
            while( my $hit = $result->next_hit ) 
            {
                last if($num_hits_done > $global_results_max_cut_off);
                my $hit_name = $hit->name;
                $hit_name =~ s/UniRef90_([^ ]*).*/$1/;
                if("Unknown" ne $hit_name)
                {
                    if ($hit->significance <= $global_evalue_cutoff) 
                    {
                        if(exists $global_U2A_hash{$hit_name})
                        {
                            $tmp_ann = $global_U2A_hash{$hit_name};
                            $num_hits_done++;
                        }
                    }
                }
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
        if(exists $global_annotations_hash{$feat_key})
        {
            # this is a bit dodge...
            # earlier, I put __DOOF in place of any null annotations
            if("__DOOF" ne $global_annotations_hash{$feat_key})
            {
                $gff_bits[8] =~ s/: hypothetical protein//;
                $gff_bits[8] = $gff_bits[8].$global_annotations_hash{$feat_key};
                $gff_bits[8] = gff_collapse_tags($gff_bits[8]);
                print $ann_fh (join "\t", @gff_bits)."\n";
            }
            else
            {
                print $ann_fh "$gff_string\n";
            }
        }
        else
        {
            print $ann_fh "$gff_string\n";
        }
    }
    # clean up
    $gffio_in->close();
    close $ann_fh;
}


sub createFlatFile
{
	print "generating genbank files for contigs...";
	my $cmd = run("gff2genbank.pl $options->{'c'} $global_output_file");
	print "done\n";
}


sub run {
   # Run a command, check that it completed successfully and return its output
   my ($cmd) = @_;
   my $results = `$cmd`;
   die "Error: Command '$cmd' failed\n$!\n" if ($? == -1);
   return $results;
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
         "seq_embed|s",
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

#    
#    STEPS:
#
#    BEFORE YOU RUN THE MANOTATOR!
#    
#    1a. Send entire genome out to RAST server.
#    1b. Run entire genome though prodigal
#    
#    You now have 2 gff3 files...
#
#    NOW, UNLEASH THE MANOTATOR!
#

=head1 SYNOPSIS

    mannotator.pl -gffs|g GFF_FILE1[,GFF_FILE2...] -contigs|c FILE -blast_prg|p BLAST_TYPE -protdb|p LOCATION -i2a|i FILE

      -gffs -g FILE[,FILE]         List of gff3 files in order of trust!
      -contigs -c FILE             Contigs to be annotated...
      -protdb -p LOCATION          Location of the UniRef or Nr BLAST database
      -i2a -i FILE                 ID to annotations mapping :file

      [-seq_embed -s]              Embed sequences into GFF files (useful to view the annotation in Artemis)     
      [-threads -t]                Number of blast jobs to run [default: 1]
      [-flatfile -f]               Optionally create multiple genbank files for your contigs [default: do not create]
      [-blast_prg -b BLAST TYPE]   The type of blast to run [default: blastx]
      [-evalue -e DECIMAL]         E-value cut off for blastx [default: 0.000000001]
      [-out -o FILE]               Filename of the final gff3 file [default: mannotatored.gff3]
      [-keep -k]                   Keep all the tmp directories
      [-blastx -x]                 Keep only the blastx file (overrides -k option)
      [-sims -s FILE]              Use BLAST similarities given in this file
      [-one]                       Skip step 1 of Mannotator
      [-two]                       Skip step 2 of Mannotator
      [-three]                     Skip step 3 of Mannotator
      [-four]                      Skip step 4 of Mannotator
      [-five]                      Skip step 5 of Mannotator
      [-help -h]                   Displays basic usage information
         
=cut


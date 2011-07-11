#!/usr/bin/perl
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
print "Splitting gffs and making tmp directories\n";
splitGffs();

#
# Step 1b. Split down the fasta file and also put into the same folders
#
print "Splitting fasta\n";
splitFasta();

#
# Step 2. For each folder, combine the gff files and produce a list of sequences to be blasted!
#
print "Combining Gffs\n";
combineGffs();

#
# Step 3. Blast the unknowns against UniRef
#
print "Blasting unknown ORFs against UniRef... This could take some time. Perhaps you need a coffee?\n";
blastUnknowns();

#
# Step 4. Annotate the GFFs using blast results and re-combine!
#
print "Annotating using the blast results\n";
annotate();

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
        open my $gff_fh, "<", $gff or die $!;
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
                my $cmd = "mkdir -p $bits[0]";
                `$cmd`;
                $global_tmp_folders{$bits[0]} = 1;
                    
                # make a new file handle
                open $current_fh, ">", $bits[0]."/$gff" or die $!;
                    
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
        my $cmd = "combineGff3.pl -c $current_folder/sequence.fa -g $gff_str -o $current_folder/combined.gff3 -a $current_folder/unknowns.fa";
        `$cmd`;
        
        # move the unknowns onto the pile
        `cat $current_folder/unknowns.fa >> $global_tmp_fasta`;
    }
}
sub worker {
	my ($chunk, $n) = @_;
	my $cmd = "blastall -p $blast_program -i $chunk -d $options->{'uniref'} -o $global_tmp_fasta.$n.$blast_program -m 8";
    `$cmd`;
}


sub blastUnknowns {
    #-----
    # Blast unknowns against the UniRef90 DB
    #
    # first blast em'
    my $num_seq = `grep -c ">" $global_tmp_fasta`;
    print "total sequences to blast: $num_seq\n";
    if ($threads > 1)
    {
    	my $num_seq_per_file = int ($num_seq / $threads);
    	my $seqio_global = Bio::SeqIO->new(-file => $global_tmp_fasta, -format => 'fasta');
    	
    	for (my $i = 1; $i <= $threads; $i++)
    	{
    		# open a file to hold a chunk
			open (CH, ">","$global_tmp_fasta_".$i) or die $!;
            my $j = 0;
			while(my $fasta = $seqio_global->next_seq())
			{
                last if $j > $num_seq_per_file;
                print CH ">".$fasta->primary_id."\n".$fasta->seq()."\n";	
				$j++;
			}
            close CH;
    	}
    	
   		for (my $i = 1; $i <= $threads; $i++) 
   		{
     		print "spawning thread $i\n";
     		my $q = "$global_tmp_fasta_$i";
     		threads->new(\&worker, $q, $i);
   		}
            
        $_->join() for threads->list();
   		
		
		for (my $i = 1; $i <= $threads; $i++)
		{
			`cat global_tmp_fasta.$i.$blast_program >> $global_tmp_fasta.$blast_program`
		}
		
    }
    else
    {
    	my $cmd = "blastall -p $blast_program -i $global_tmp_fasta -d $options->{'uniref'} -o $global_tmp_fasta.$blast_program -m 8";
    	`$cmd`;
    }
    
    # now split them across multiple folders...
    my $current_dir_name = "__DOOF";
    my $current_file_handle;
    open my $blast_results, "<", "$global_tmp_fasta.$blast_program" or die $!;
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
            open $current_file_handle, ">", "$dir_name/unknowns.$blast_program" or die $!;
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
    
    # load the databases
    &loadU2A($options->{'u2a'});
    
    my $counter = 0;
    my $big_counter = 0;
    # first do the call out
    foreach my $current_folder (keys %global_tmp_folders)
    {
        # get the blast results:
        my @blast_files = ();
        if (-d $current_folder) {
            opendir(INDIR, $current_folder)
            or die "Failed to read from directory $current_folder:$!\n";
            @blast_files = grep /\.[t]blast[pnx]$/, readdir (INDIR);
            closedir (INDIR);
        }

        next if($#blast_files == -1);

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
    open my $out_fh, ">", $global_output_file or die $!;
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
    close $out_fh;
}

sub cleanTmps {
    #-----
    # clean up all the tmp files we maded
    #
    my ($keep_bx) = @_;
    foreach my $current_folder (keys %global_tmp_folders)
    {
        my $cmd = "rm -rf $current_folder";
        `$cmd`;
    }
    
    `rm $global_tmp_fasta`;
    
    if(0 == $keep_bx)
    {
       `rm $global_tmp_fasta.$blast_program`;
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
    
    print "Loading Uniref => Annotations file....";
    
    open my $U2A_fh, "<", $Uniprot2ANN_file or die "Cannot read from U2A file: $Uniprot2ANN_file:$!\n";
    while (<$U2A_fh>) {
        chomp $_;
        my @data = split(/\^/, $_);
        $global_U2A_hash{$data[0]} = $data[1];
    }
    close $U2A_fh;
    print "done.  Loaded ".keys(%global_U2A_hash)." Uniprot 2 annotation refs\n";
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
            or die "Failed to read from $blast_file: $!\n";
        
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
            $global_annotations_hash{$result->query_name} = $tmp_ann;
        }
    }    
}

sub recombineGff3() {
    #-----
    # blend the new annotations into the old
    #
    my ($gff3_in, $gff3_out) = @_;
    
    # specify input via -fh or -file
    my $gffio_in = Bio::Tools::GFF->new(-file => $gff3_in, -gff_version => 3);
    open my $ann_fh, ">", $gff3_out or die $!;
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
	my $cmd = "gff2genbank.pl $options->{'c'} $global_output_file";
	`$cmd`;
	print "done\n";
}


######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "gffs|g:s", "keep|k+", "contigs|c:s", "u2a|a:s", "uniref|u:s", "out|o:s", "blastx|x:+", "evalue|e:s","blast_prg|p:s", "flatfile|f:+", "threads|t:s");
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

    if(!exists $options{'uniref'})
    {
        print "ERROR: You need to give me the location of the UniRef blast database to continue!\n";
        print "Perhaps, try this one (on EURY): /Volumes/Biodata/BLAST_databases/UniProt/UniRef90/uniref90_micro.fasta\n";
        exec("pod2usage $0");
    }
    
    if(!exists $options{'u2a'})
    {
        print "ERROR: You need to give me the location of the Uniref to Annotations file!\n";
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

    mannotator.pl -gffs|g GFF_FILE1[,GFF_FILE2...] -contigs|c FILE -blast_prg|p BLAST TYPE -uniref|u LOCATION -u2a|a FILE

      -gffs -g FILE[,FILE]         List of gff3 files in order of trust!
      -contigs -c FILE             Contigs to be annotated...
      -uniref -u LOCATION          Location of UniRef blast database
      -u2a -a FILE                 UniRef to annotations file
      [-threads -t]                Number of blast jobs to run [default: 1]
      [-flatfile -f]               Optionally create multiple genbank files for your contigs [default: do not create]
      [-blast_prg -p BLAST TYPE]   The type of blast to run [default: blastx]
      [-evalue -e DECIMAL]         E-value cut off for blastx [default: 0.000000001]
      [-out -o FILE]               Filename of the final gff3 file [default: mannotatored.gff3]
      [-keep -k]                   Keep all the tmp directories
      [-blastx -x]                 Keep only the blastx file (overrides -k option)
      [-help -h]                   Displays basic usage information
         
=cut


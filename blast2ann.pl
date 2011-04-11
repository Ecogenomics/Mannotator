#!/usr/bin/perl -w
###############################################################################
#
#    blast2ann.pl
#    
#    Convert blast results to GO annotations
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

#core Perl modules
use Getopt::Long;
use Bio::SearchIO;
use Bio::Tools::GFF;

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

# default working directory
my $global_working_dir = ".";
if(exists $options->{'working_dir'}) { $global_working_dir = $options->{'working_dir'}; }

# default output file
my $global_output_filename = "annotated.gff3";
if(exists $options->{'out'}) { $global_output_filename = $options->{'out'}; }

# evalue cutoff
my $global_evalue_cutoff = 0.000000001;
if(exists $options->{'evalue'}) { $global_evalue_cutoff = $options->{'evalue'}; }

#results cut off
my $global_results_max_cut_off = 1;

# get the blast results:
my @blast_files = ();
if (-d $global_working_dir) {
    opendir(INDIR, $global_working_dir)
    or die "Failed to read from directory $global_working_dir:$!\n";
    @blast_files = grep /\.blast[nx]$/, readdir (INDIR);
    closedir (INDIR);
}
else {
    push (@blast_files, $global_working_dir);
}

if($#blast_files == -1)
{
    print "ERROR: no blast files found in current working directory: [$global_working_dir]\n";
    exit;
}
print "parsing ".@blast_files." blast results\n";

# load the databases
&loadU2A($options->{'u2a'});

# parse the blast results and report the GO predictions
&generateAnnotations(@blast_files);

# insert these new values into the GFF3 file
&recombineGff3($options->{'gff3'}, $global_output_filename);

######################################################################
# CUSTOM SUBS
######################################################################

sub loadU2A() {
    #-----
    # load the eggNOG association file.
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
    print "done.  Loaded ".keys(%global_U2A_hash)." Uniprot 2 eggNOG refs\n";
}

sub generateAnnotations() {
    #-----
    # parse the blast results and report the predictions
    #
    my @blast_files = @_;
    
    foreach my $blast_file (@blast_files) 
    {
        # Load the result into memory
        print "Now parsing: $blast_file\n";
        my $in = new Bio::SearchIO(-format => 'blast',
                                -file => "$global_working_dir/$blast_file")
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
    print "Parsed all blast files!\n";
}

sub recombineGff3() {
    #-----
    # blend the new annotations into the old
    #
    my ($gff3_in, $gff3_out) = @_;
    print "Merging new annotations with: $gff3_in. Printing to $gff3_out\n";
    
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
    
    print "Done!\n";
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "evalue|e:i", "working_dir|w:s", "out|o:s", "u2a|u:s", "gff3|g:s", );
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

    if(!exists $options{'u2a'})
    {
        print "ERROR: You forgot to supply the location of the UniRef 2 annotations file\n";
        print "Perhaps, try this one (on EURY): /Volumes/Biodata/ANN_mappings/ANN_mappings.txt\n";
        exec("pod2usage $0");
    }
    
    if(!exists $options{'gff3'})
    {
        print "ERROR: You need to supply a gff3 file\n";
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

    blast2ann.pl

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

   Insert detailed description here

=head1 SYNOPSIS

        blast2ann.pl -u2a|u ANN_file -gff3|g gff3_file

        -u2a -u FILE                    UniRef to annotations file
        -gff3 -g FILE                   Gff3 file to add annotations to
        [-working_dir -w DIRECTORY]     location of the blast output files [default: current directory]
        [-out -o FILE]                  file to write results to [default: annotated.gff3]
        [-evalue -e NUMBER]             evalue cut off for blast results [default: 0.000000001]
        [-help -h]                      Displays basic usage information
         
=cut


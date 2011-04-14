#!/usr/bin/perl
###############################################################################
#
#    idMerge.pl
#    
#    Merge UniProt, KEGG, COG-eggNOG ids and info into one big kludge
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

#CPAN modules

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
my $global_outfile_name = "ANN_mappings.txt";
if(exists $options->{'o'})
{
    $global_outfile_name = $options->{'o'};
}

my %global_seenUP_hash = ();
my %global_U2K_hash = ();
my %global_KEGG_p_hash = ();
my %global_KEGG_k_hash = ();
my %global_KEGG_e_hash = ();
my %global_UPN_hash = ();
my %global_COG2Txt_hash = ();

# open all the files!
open my $COG_fh, "<", $options->{'Uc'} or die $!;
open my $KEGGU_fh, "<", $options->{'Ku'} or die $!;

open my $KEGGP_fh, "<", $options->{'Kp'} or die $!;
open my $KEGGK_fh, "<", $options->{'Kk'} or die $!;
open my $KEGGE_fh, "<", $options->{'Ke'} or die $!;

open my $N2A_fh, "<", $options->{'Cd'} or die $!;

open my $OUT_fh, ">", $global_outfile_name or die $!;

# parse parse away!

# first load all the UniProt ID's vs KEGG
print "Loading UniProt Vs KEGG...";
while(<$KEGGU_fh>)
{
    #
    # hsa:1000        up:P19022       original
    #
    chomp $_;
    my @data = split /\t/, $_;
    $data[1] =~ s/up://;
    
    {
        if(!exists $global_U2K_hash{$data[1]})
        {
            $global_U2K_hash{$data[1]} = $data[0];
        }
        
        if(!exists $global_seenUP_hash{$data[1]})
        {
            $global_seenUP_hash{$data[1]} = 0;
        }
    }
}
close $KEGGU_fh;
print "done\n";

# make a hash (or two) of KEGG entry ID's versus pathway etc...
print "Loading KEGG pathways...";
while(<$KEGGP_fh>)
{
    chomp $_;
    my @data = split /\t/, $_;
    if(!exists $global_KEGG_p_hash{$data[0]})
    {
        $global_KEGG_p_hash{$data[0]} = $data[1];
    }
}
print "done\n";
print "Loading KEGG ontology...";
while(<$KEGGK_fh>)
{
    chomp $_;
    my @data = split /\t/, $_;
    if(!exists $global_KEGG_k_hash{$data[0]})
    {
        $global_KEGG_k_hash{$data[0]} = $data[1];
    }
}
print "done\n";
print "Loading KEGG enzyme...";
while(<$KEGGE_fh>)
{
    chomp $_;
    my @data = split /\t/, $_;
    if(!exists $global_KEGG_k_hash{$data[0]})
    {
        $global_KEGG_e_hash{$data[0]} = $data[1];
    }
}
print "done\n";

close $KEGGP_fh;
close $KEGGK_fh;
close $KEGGE_fh;

# load all the COG text
print "Loading COG/NOG descriptors...";
while (<$N2A_fh>) {
    chomp $_;
    my @data = split(/\t/, $_);
    if(!exists $global_COG2Txt_hash{$data[3]})
    {
        if("Annotation not available" ne $data[4])
        {
            $global_COG2Txt_hash{$data[3]} = $data[4];
        }
    }
}
print "done\n";
close $N2A_fh;

# Link UPIDs to COG IDs
print "Loading UniProt Vs COG/NOG...";
while(<$COG_fh>)
{
    chomp $_;
    my @data = split /\t/, $_;
    if(!exists $global_UPN_hash{$data[0]})
    {
        $global_UPN_hash{$data[0]} = $data[1];
    }
    if(!exists $global_seenUP_hash{$data[0]})
    {
        $global_seenUP_hash{$data[0]} = 1;
    }    
}
print "done\n";
close $COG_fh;

# print the results
print "Vomiting on the file system...\n";
foreach my $UPID (keys %global_seenUP_hash)
{
    my $found_one = 0;
    
    # UPID
    my $out_string = "$UPID^";
    
    # cog
    if(exists $global_UPN_hash{$UPID})
    {
        $found_one = 1;
        $out_string .= ";Ontology_term=COG_ID:".$global_UPN_hash{$UPID};
        if(exists $global_COG2Txt_hash{$global_UPN_hash{$UPID}})
        {
            $out_string .= ";Ontology_term=COG_DESC:".$global_COG2Txt_hash{$global_UPN_hash{$UPID}};
        }
    }
    
    # KEGG entry ID
    if(exists $global_U2K_hash{$UPID})
    {
        $found_one = 1;
        $out_string .= ";Ontology_term=KEGG_ID:".$global_U2K_hash{$UPID};
        
        # KEGG pathway ID
        if(exists $global_KEGG_p_hash{$global_U2K_hash{$UPID}})
        {
            $out_string .= ";Ontology_term=KEGG_PATHWAY:".$global_KEGG_p_hash{$global_U2K_hash{$UPID}};
        }
        
        # KEGG ko ID
        if(exists $global_KEGG_k_hash{$global_U2K_hash{$UPID}})
        {
            my $ko_ID = $global_KEGG_k_hash{$global_U2K_hash{$UPID}};
            $out_string .= ";Ontology_term=KEGG_ONTOLOGY:$ko_ID";
            
            # KEGG enzyme
            if(exists $global_KEGG_e_hash{$ko_ID})
            {
                $out_string .= ";Ontology_term=KEGG_ENZYME:".$global_KEGG_e_hash{$ko_ID};
            }
        }
    }
    
    if(0 != $found_one)
    {
        chomp $out_string;
        print $OUT_fh $out_string."\n";
    }
}
print "done\n";
close $OUT_fh;

######################################################################
# CUSTOM SUBS
######################################################################


######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "Ku:s", "Uc:s", "Kp:s", "Kk:s", "Ke:s", "Cd:s", "o:s", );
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
    
    if(!exists $options{'Uc'} || !exists $options{'Kp'} || !exists $options{'Ke'} || !exists $options{'Kk'} || !exists $options{'Cd'} || !exists $options{'Ku'})
    {
        print "ERROR: Check your input parameters!\n";
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

    idMerge.pl

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

    Merge UniProt, KEGG, COG-eggNOG ids and info into one big kludge
    UniProt ID is THE key. This script will keep all  UniProt IDs which
    have at least one of the following included as a reference...
    
    COG ID
    KEGG ID

=head1 SYNOPSIS

    idMerge.pl -Ku KEGG2Uniprot_file -Uc UniProt2Cog_file -Kp KEGG_pathways_file -Kk KEGG_ko_file -Ke KEGG2Enzyme_file -cD COG_text [-o OUTFILE] [-help|h]

        -Uc   UniProt2Cog_file      UniProt 2 COG File [http://eggnog.embl.de/cgi_bin/show_download_page.pl]
        -Ku   KEGG2Uniprot file     KEGG genes to Uniprot accessions [ftp://ftp.genome.jp/pub/kegg/linkdb/genes/]
        -Ke   KEGG2Enzyme file      KEGG genes to Uniprot accessions [munged from ftp://ftp.genome.jp/pub/kegg/brite/ko/]
        -Kp   KEGG_pathways_file    KEGG entry ID to pathway ID [ftp://ftp.genome.jp/pub/kegg/linkdb/genes/]
        -Kk   KEGG_ko_file          KEGG entry ID to KEGG Ontology ID [ftp://ftp.genome.jp/pub/kegg/linkdb/genes/]
        -Cd   COG_text              Links descriptions to COG/NOG IDs [http://eggnog.embl.de/cgi_bin/show_download_page.pl]
        
        
        [-o   OUTFILE]              File to write mappings to. [Default ANN_mappings.txt]
        [-help -h]                  Displays basic usage information
         
=cut


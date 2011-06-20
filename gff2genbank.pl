#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;


if(! $ARGV[0]){  die "\nUSAGE: ./gff2genbank <FASTA> <GFF>\n";}

my $subject = $ARGV[0];
my $fasta = Bio::SeqIO->new('-file' => $subject,
                             	 '-format' => 'fasta');

my %features;


my $gffio = Bio::Tools::GFF->new(-file => $ARGV[1], -gff_version => 3);
    my $feature;
    # loop over the input stream
    while($feature = $gffio->next_feature()) 
    {
        push(@{$features{ $feature->seq_id }}, $feature);
    }
$gffio->close();

while(my $seqobj = $fasta->next_seq())
{
	if (exists $features{$seqobj->primary_id})
	{
		my $seqout = Bio::SeqIO->new(-file   => ">".$seqobj->primary_id().".gb",
                                 -format => 'genbank');
    	$seqobj->add_SeqFeature(@{$features{$seqobj->primary_id}});

		$seqout->write_seq($seqobj);
	}
	else
	{
		warn "there are no features for $seqobj->primary_id\n";
	}
}
exit;


# 
# my %genes;
# while( my $f = $gff->next_feature ) {
#    my ($group) = $feature->get_tag_values('Group'); # substitute group 
# #with whatever you have in the group field
#   push @{$gene{$group}}, $feature;
# }
# # get a Bio::Seq object called $seq somehow, either by reading in a 
# #fasta sequence file, etc...
# while( my ($gene,$features) = each %genes ) 
# {
#   my $location = Bio::Location::Split->new();
#   for my $f ( @$features ) {
#     $location->add_sub_Location($f->location);
#   }
#   my $genef = Bio::SeqFeature::Generic->new(-location =>$location, 
# -primary_tag => 'CDS');
#   $seq->add_SeqFeature($genef);
# }
# my $seqio = Bio::SeqIO->new(-format => 'genbank');
# $seqio->write_seq($seq);

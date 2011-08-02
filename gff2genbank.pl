#! /usr/bin/env perl


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
        # if ($feature->primary_tag()  eq 'CDS' )
#         {
#         	my $gene_feature = Bio::SeqFeature::Generic->new( -start => $feature->start(),
#         													  -end => $feature->end(),
#         													  -strand => $feature->strand(),
#         													  -primary_tag => 'gene')
#         	$feature->add_tag_feature('gene','')
#         }
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
		warn "there are no features for ",$seqobj->primary_id, "\n";
	}
}
exit;


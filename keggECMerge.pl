#! /usr/bin/env perl


use warnings;
use strict;

my $fh;

open $fh, "<", $ARGV[0] or die $!;

my %hash;
while(<$fh>)
{
	chomp $_;
	if ($_ =~ /\>(\w+)\<\/a\>(.+)\[EC.*(\d+\.(\d+|-)\.(\d+|-)\.(\d+|-))\]$/)
	{
		my ($gene_symbol, $gene_name) = split (/;/, $2);
		$hash{$1} = $3;
		print "$1\t$3";
		if (defined $gene_name)
		{
			$hash{$1} .= "\t$gene_name";
		}
	}
}

while(my ($key, $value) = each(%hash))
{
	print "$key\t$value\n";
}
close $fh;
exit;

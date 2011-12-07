#! /usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;

if($#ARGV != 2) { print "Usage: $0 uniref_fasta ANN_mappings.txt out_file\n"; exit; }

my %UR_exists_hash = ();

open my $AM_fh, "<", $ARGV[1] or die $!;
open my $out_fh, ">", $ARGV[2] or die $!;

# load all "seen" Uniref IDs
print "loading all 'seen' Uniref IDs\n";
my $count = 0;
my $tmp_count = 0;
while(<$AM_fh>)
{
    chomp $_;
    my @fields = split /\^/, $_;
    $UR_exists_hash{$fields[0]} = 1;
    $tmp_count++;
    if($tmp_count == 100000)
    {
        $count+=100000;
        print "Parsed $count\n";
        $tmp_count = 0;
    }
}

print "done!\nLoading Uniref strings and printing!\n";
$count = 0;
$tmp_count = 0;
my $dun_count = 0;
my $seqio_object = Bio::SeqIO->new(-file => $ARGV[0], -format => "fasta");
while (my $seq = $seqio_object->next_seq)
{
    my $id = $seq->display_id();
    $id =~ s/UniRef90_//;
    if(exists $UR_exists_hash{$id})
    {
        print $out_fh ">UniRef90_" . $id . "\n";
        print $out_fh $seq->seq(). "\n";
        $dun_count++;
    }
    $tmp_count++;
    if($tmp_count == 100000)
    {
        $count+=100000;
        print "Parsed $count : retained : $dun_count\n";
        $tmp_count = 0;
    }
}
print "done!\n";

close $seqio_object;
close $AM_fh;
close $out_fh;

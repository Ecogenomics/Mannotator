#! /usr/bin/env perl

=head1 NAME

combineGffAnnotations.pl

=head1 COPYRIGHT

Copyright (C) 2011 Florent Angly

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

Combine two Mannotator GFF files, one that have been assigned a function, the
other that has been assigned a taxonomy, into a single GFF file. This is done
by merge the taxonomic and functional annotations of ORFs that are common
between the two files.

=head1 REQUIRED ARGUMENTS

=over

=item -f <function_gff>

GFF file that contains functional annotations

=for Euclid:
     function_gff.type: readable

=item -t <taxonomy_gff>

GFF file that contains taxonomic annotations

=for Euclid:
     taxonomy_gff.type: readable

=back

=head1 OPTIONAL ARGUMENTS

=over

=item -o <output_gff>

Location where to write the combined GFF. Default: output_gff.default

=for Euclid:
     output_gff.type:    writable
     output_gff.default: 'combinotatored.gff'

=cut


use strict;
use warnings;
use Getopt::Euclid;
use Bio::SeqIO;
use Bio::Tools::GFF;

combineGffAnnotations($ARGV{-f}, $ARGV{-t}, $ARGV{-o});

exit;



sub combineGffAnnotations {
   my ($function_gff, $taxonomy_gff, $output_gff) = @_;

   #1# Deal with annotations

   # Record functional annotations
   my %functions;
   my $function_io = Bio::Tools::GFF->new(-file => $function_gff , -gff_version => 3);
   while ( my $ffeat = $function_io->next_feature() ) {
      $functions{feat_id($ffeat)} = $ffeat;
   }

   # Read taxonomic annotations, merge them with functional and write them 
   my $taxonomy_io = Bio::Tools::GFF->new(-file => $taxonomy_gff , -gff_version => 3);
   my $output_io   = Bio::Tools::GFF->new(-file => ">$output_gff", -gff_version => 3);
   while ( my $tfeat = $taxonomy_io->next_feature() ) {
      my $tfeat_id = feat_id($tfeat);
      my $ofeat;
      if (exists $functions{$tfeat_id}) {
         # Merge functional and taxonomic annotation
         my $ffeat = $functions{$tfeat_id};
         $ofeat = merge_feats($ffeat, $tfeat);
         # Delete functional annotation
         delete $functions{$tfeat_id};
      } else {
         # Just keep the taxonomic annotation
         $ofeat = $tfeat;
      }
      # Write annotation
      $output_io->write_feature($ofeat);
   }

   # Now write remaining functional annotations, if any
   while (my (undef, $ffeat) = each %functions) {
      $output_io->write_feature($ffeat);
   }
   $output_io->close;
   %functions = ();


   #2# Now write all sequences, if any

   my @fseqs = $function_io->get_seqs();
   $function_io->close;
   my @tseqs = $function_io->get_seqs();
   $taxonomy_io->close;

   if ( (scalar @fseqs > 0) || (scalar @tseqs > 0) ) {
   
      # Open filehandle
      open my $output_fh, '>>', $output_gff or die "Error: Could not write file $output_gff\n$!\n";
      print $output_fh "##FASTA\n";
      $output_io = Bio::SeqIO->new( -format => 'fasta', -fh => $output_fh );

      # Record functional sequences
      for my $fseq ( @fseqs ) {
         $functions{$fseq->id} = $fseq;
      }

      # Write taxonomic sequences
      for my $tseq ( @tseqs ) {
         $output_io->write_seq($tseq);
         delete $functions{$tseq->id};
      }
 
      # Write remaining functional sequences, if any
      while ( my (undef, $fseq) = each %functions ) {
         $output_io->write_seq($fseq);
      }
      %functions = ();

      $output_io->close;
   }

   return 1;
}


sub feat_id {
   # Give an ID to a feature
   my ($feat) = @_;
   my $id = $feat->start.'-'.$feat->end.'/'.$feat->strand;
   return $id;
}


sub merge_feats {
   # Merge a functional feature and taxonomic feature into an output feature
   my ($ffeat, $tfeat) = @_;
   my $ofeat = $ffeat;
   my @tags_to_import = ('Note', 'Dbxref');
   for my $tag_name (@tags_to_import) {
      if ($tfeat->has_tag($tag_name)) {
         $ffeat->add_tag_value( $tag_name, $tfeat->get_tag_values($tag_name) );
      }
   }
   return $ofeat;
}


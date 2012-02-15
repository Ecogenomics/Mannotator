#! /usr/bin/env perl

=head1 NAME

idMergeNr.pl

=head1 COPYRIGHT

Copyright (C) 2011 2012 Florent Angly Connor Skennerton

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

Prepare a subset of the NCBI nr proteic database (non-redundant) for use 
in Mannotator by merging GI number, sequence description and taxonomy into
a single file.  Preparing these files takes a while (1-2 days)!

CAVEAT: The nr database defline (FASTA header) for each sequence is the
concatenation of the headers of multiple identical protein sequences. Here,
only the information of the first sequence is used.

=head1 REQUIRED ARGUMENTS

=over

=item -Fa <fasta>

FASTA file that contains the Nr database of proteic sequences (ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz)

=for Euclid:
     fasta.type: readable

=item -No <nodes>

NCBI taxonomy nodes file (nodes.dmp) download it from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

=for Euclid:
     nodes.type: readable

=item -Na <names>

NCBI taxonomy names files (names.dmp) downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

=for Euclid:
     names.type: readable

=item -Gt <gitaxid>

NCBI GI (genbank ID) to taxid (taxon ID) file downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz

=back

=head1 OPTIONAL ARGUMENTS

=over

=item -o <outfile>

File to write mappings to. Default: outfile.default

=for Euclid:
     outfile.type:    writable
     outfile.default: 'Nr_mappings.txt'

=item -Fo <outfasta>

Output fasta file containing sequences of the subset

=for Euclid:
    outfasta.type: writable
    outfasta.default: 'Nr_subset.fa'

=cut


use strict;
use warnings;
use Cwd;
use File::Spec::Functions;
use Getopt::Euclid;
use Bio::DB::Taxonomy;
use URI::Escape;
use DB_File;
idMergeNr($ARGV{-Fa}, $ARGV{-No}, $ARGV{-Na}, $ARGV{-Gt}, $ARGV{-o}, $ARGV{-Fo});
exit;


sub idMergeNr {
   my ($fasta, $nodes, $names, $gitaxid, $outmappings, $fasta_subset) = @_;

   # Reading and indexing taxonomy files
   print "Reading taxonomy files...\n";
   my $tax_db = create_taxonomy_db($nodes, $names);
   print "done\n";

   # Read GI to taxid file
   print "Reading gi2taxid file...\n";
   my ($gitaxid_db, $gitaxid_db_file) = create_gitaxid_db($gitaxid);
   print "done\n";

   # Reading sequences
   print "Reading sequences file and writing mappings...\n";
   process_sequences($fasta, $outmappings, $gitaxid_db, $tax_db, $fasta_subset);
   print "done\n";

   print "Mapping file is ready at $outmappings\n\n";

   print "You can now delete the following database files unless you plan on running $0 again:\n";
   print "  id2names, names2id, nodes, parents, gi2taxid\n\n";
 
   return 1;
}


sub create_taxonomy_db {
  my ($nodes, $names) = @_;
  my $db = Bio::DB::Taxonomy->new(
     -source    => 'flatfile',
     -nodesfile => $nodes    ,
     -namesfile => $names    ,
     -directory => getcwd    , # directory for index file storage (default: /tmp)
     #-force     => 1         , # force the creation of the index file
  );
  return $db;
}


sub create_gitaxid_db {
   my ($gitaxid) = @_;
   my $name = 'gi2taxid';
   my ($gitaxid_db, $gitaxid_db_file, $gitaxid_db_exists) = tie_hash( getcwd, $name);
   if (not $gitaxid_db_exists) {
      # Populate tied hash with gi2taxid entries
      my $rec_no = 0;
      open my $in, '<', $gitaxid or die "Error: Could not read file '$gitaxid'\n$!\n";
      while ( my $line = <$in> ) {
         $rec_no++;
         print "   $rec_no\n" if ($rec_no % 1E5 == 0);
         chomp $line;
         my ($gi, $taxid) = split /\t/, $line;
         $$gitaxid_db{$gi} = $taxid;
      }
      close $in;
   }
   return $gitaxid_db, $gitaxid_db_file;
}


sub delete_gitaxid_db {
   my ($hashref, $file) = @_;
   untie_hash($hashref, $file);
   return 1;
}


sub tie_hash {
   # Create a single-level hash tied to a Berkeley database file
   # Input:  directory where to put database file (optional)
   #         prefix for database file name (optional)
   # Output: full path of database file
   my ($dir, $file) = @_;
   my $hashref = {};
   if (not defined $dir) {
      $dir = getcwd;
   }
   if (not defined $file) {
      $file = File::Temp->new(TEMPLATE => 'temp_XXXXXXXXX', SUFFIX => '.db')->filename;
   }
   $file = catfile($dir, $file);
   my $exists = 0;
   if ( -e $file ) {
      # Re-use a previously created database
      $exists = 1;
      my $err_msg = "Error: Could not read tied hash file '$file'\n"; 
      tie %$hashref, 'DB_File', $file, O_RDONLY, 0644, $DB_BTREE or die $err_msg."$!\n";
   } else { 
      # Create a BerkeleyDB Btree database. Btrees are fast to build and fast to
      # search. Also the keys are ordered!
      my $err_msg = "Error: Could not write tied hash file '$file'\n";
      tie %$hashref, 'DB_File', $file, O_CREAT|O_RDWR, 0644, $DB_BTREE or die $err_msg."$!\n";
   }
   return $hashref, $file, $exists;
}


sub untie_hash {
  # Remove a tied hash and its Berkeley database file
  # Input:  ref of variable to untie and file location
  # Output: 1 for success
  my ($var, $file) = @_;
  my $ref = ref $var;
  if (not defined $ref) {
    die "Error: untie_hash() expects a reference to a variable but got none.\n";
  }
  if ($ref eq 'HASH') {
    untie %$var;
  } else {
    die "Error: Cannot untie variable of type '$ref'\n";
  }
  unlink $file or die "Error: Could not remove temporary file '$file'\n$!\n";
  undef %$var;
  return 1;
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


sub process_sequences {
   my ($fasta, $outmappings, $gitaxid_db, $tax_db, $outFasta) = @_;
   open my $in, '<', $fasta or die "Error: Could not read file '$fasta'\n$1\n";
   open my $outmap, '>', $outmappings or die "Error: Could not write file '$outmappings'\n$1\n";
   open my $outfastafh, '>', $outFasta or die "Error: Could not write output fasta file '$outFasta'\n$1\n";
   my $esc = ',;='; # characters to escape in GFF tags or values
   my $rec_no = 0;
   my @aux = undef;
    my ($name, $seq, $qual);
    while (($name, $seq, $qual) = readfq($in, \@aux)) {
      # In nr, multiple identical protein sequences are merged and their defline
      # is concatenated. Example of nr sequence defline:
      # >gi|66816243|ref|XP_642131.1| hypothetical protein DDB_G0277827 [Dictyostelium discoideum AX4]^Agi|1705556|sp|P54670.1|CAF1_DICDI RecName: Full=Calfumirin-1; Short=CAF-1^Agi|793761|dbj|BAA06266.1| calfumirin-1 [Dictyostelium discoideum]
      $rec_no++;
      print "   $rec_no\n" if ($rec_no % 1E5 == 0);
      my ($id, $gi, $desc) = ( $name =~ m/(gi\|(\d+)\S+)\s+(.*?)(\[|\001)?/ );
      my $gff_blurb = "Name=".uri_escape($desc, $esc).";Dbxref=NCBI_gi:$gi;";
      my $taxid = $$gitaxid_db{$gi};
      if (defined $taxid) { # some sequences do not have a taxid
         my $taxon = $tax_db->get_taxon( -taxonid => $taxid );
         if ($taxon) { # some GIs have a taxon ID of zero, i.e. unknown source organism
            my $taxo_str = $taxon->node_name;
            while ( $taxon = $taxon->ancestor ) {
               $taxo_str = $taxon->node_name . '/' . $taxo_str;
            }
            $gff_blurb .= "Dbxref=taxon:$taxid;Note=".uri_escape($taxo_str, $esc).";";
         }
      }
      print $outmap "$id^".gff_collapse_tags($gff_blurb)."\n";
   }
   close $in;
   close $outmap;
   return 1;
}

sub readfq {
	my ($fh, $aux) = @_;
	@$aux = [undef, 0] if (!defined(@$aux));
	return if ($aux->[1]);
	if (!defined($aux->[0])) {
		while (<$fh>) {
			chomp;
			if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
				$aux->[0] = $_;
				last;
			}
		}
		if (!defined($aux->[0])) {
			$aux->[1] = 1;
			return;
		}
	}
	my $name = /^.(\S+)/? $1 : '';
	my $seq = '';
	my $c;
	$aux->[0] = undef;
	while (<$fh>) {
		chomp;
		$c = substr($_, 0, 1);
		last if ($c eq '>' || $c eq '@' || $c eq '+');
		$seq .= $_;
	}
	$aux->[0] = $_;
	$aux->[1] = 1 if (!defined($aux->[0]));
	return ($name, $seq) if ($c ne '+');
	my $qual = '';
	while (<$fh>) {
		chomp;
		$qual .= $_;
		if (length($qual) >= length($seq)) {
			$aux->[0] = undef;
			return ($name, $seq, $qual);
		}
	}
	$aux->[1] = 1;
	return ($name, $seq);
}


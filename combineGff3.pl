#! /usr/bin/env perl

###############################################################################
#
#    combineGff3.pl
#    
#    parse a list of gff3 files in order, 
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

#Perl modules
use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Class::Struct;

#locally-written modules

# implement a linked list...
struct GffObject =>
{
    GO_prevGO_ref => '$',   # next / prev in list
    GO_nextGO_ref => '$',       
    GO_gffRef =>'$',        # reference to the gff object
    GO_start => '$',        # we use these A LOT, so just save them twice...
    GO_end => '$',
    GO_next_start => '$',   # how much "free" space lies beyond this end...? Where does the next orf start?
};

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
my %global_gff_used_list = ();
my $global_any_inserted = 0;

# turn debugging on (1) or off (0)
my $debug = 0;

# algorithm based parameters...
# reject all orfs less than this amount!
my $global_min_orf_cutoff = 100;
# reject all orfs which overlap with already accepted orfs by
# at least this much
my $global_shared_olap_cutoff = 0.1;

# output gff3 file
my $default_feat_file = 'parsed.gff3';
if(exists $options->{'out'}) { $default_feat_file = $options->{'out'}; }
open my $feat_fh, ">", $default_feat_file or die "Error: Could not write file $default_feat_file\n$!\n";

# output todo annotation file
my $default_seq_file = 'todo.fa';
if(exists $options->{'ann'}) { $default_seq_file = $options->{'ann'}; }
open my $seq_fh, ">", $default_seq_file or die "Error: Could not write file $default_seq_file\n$!\n";

# first parse the fasta file of contigs to get headers and sequence lengths...
# really there should only be one guy here...
my $seqio_object = Bio::SeqIO->new(-file => $options->{'contigs'}, -format => 'fasta');
my $global_seq = $seqio_object->next_seq;
my $global_seq_length = $global_seq->length;
$seqio_object->close;

# initialise the linked list
my $global_gff_list = GffObject->new();
$global_gff_list->GO_prevGO_ref(\$global_gff_list);
$global_gff_list->GO_nextGO_ref(\$global_gff_list);
$global_gff_list->GO_next_start($global_seq_length);
$global_gff_list->GO_start($global_seq_length);
$global_gff_list->GO_end(0);

#nd value used so we don't loop the loop!
my $global_end_ref = \$global_gff_list;

# get all the gff file names!
my @gff_fns = split /,/, $options->{'gffs'};
foreach my $gff3 (@gff_fns)
{
    if (not -e $gff3) {
        die "Error: Could not read GFF file $gff3\n$!\n";
    }

    # specify input via -fh or -file
    my $gffio = Bio::Tools::GFF->new(-file => $gff3, -gff_version => 3);
    my $features_used = 0;

    my $ins_ref = \$global_gff_list;
    # loop over the input stream
    while(my $feat = $gffio->next_feature()) {
        
        # filter out orf if it is too short
        if ($feat->length < $global_min_orf_cutoff) {
            debugFeature('Not keeping feature '.feat2str($feat), "less than $global_min_orf_cutoff bp") if ($debug);
            next;
        }

        # then insert
        $ins_ref = insertFeature($ins_ref, \$feat);
    }
    
    # clean up
    $gffio->close();
    
    # record keeping
    $global_gff_used_list{$gff3} = $features_used;
}



# print to the output file
print $feat_fh "##gff-version 3\n";
my $gffio = Bio::Tools::GFF->new(-gff_version => 3);
my $list_handle_ref = \$global_gff_list;
while(nextInList(\$list_handle_ref) == 1)
{
    my $current_node = ${$list_handle_ref};
    my $current_feat_ref = $current_node->GO_gffRef;
    my $current_feature = ${$current_feat_ref};
    my $gff_string = $current_feature->gff_string($gffio);
    my @gff_bits = split /\t/, $gff_string;

    # Write sequence except if it is a RAST annotated gene
    unless( ($gff_bits[1] eq 'FIG') && ($gff_bits[8] !~ /hypothetical/) ) {
        my $seq_entry = '>'.$gff_bits[0].'_'.$current_node->GO_start.'_'.$current_node->GO_end."\n".
                        $global_seq->subseq($current_node->GO_start, $current_node->GO_end)."\n";
        print $seq_fh $seq_entry;
    }
    
    # check to see if we'll need to do some annotation afterwards...
    print $feat_fh "$gff_string\n";
}


######################################################################
# CUSTOM SUBS
######################################################################

sub debugFeature {
    # display why a feature was kept or not
    my ($msg, $details) = @_;
    print "Debug: ".$msg.".\n";
    print "Reason: $details\n" if defined $details;
    return 1;
}


sub feat2str {
    my $feat = shift;
    return $feat->source_tag." ".$feat->primary_tag." ".$feat->start."-".$feat->end;
}


sub insertFeature {
    #-----
    # insert a feature, provided that it does not displace any existing feature...
    # start_ref: the list of features
    # feat_ref: the feature that we try to put in the list
    my ($start_ref, $feat_ref) = @_;

    ####
    #use Data::Dumper;
    #print "FEATURE TO INSERT: ".Dumper($feat_ref);
    #print "LIST: ".Dumper($start_ref);
    #print "\n";
    #die "Temp exit\n";
    ####

    my $feat = ${$feat_ref};
    my $feat_start  = $feat->start;
    my $feat_end    = $feat->end;
    my $feat_length = $feat->length; # start <= end , always
    my $first_loop_done = 0;
    my $list_handle_ref = $start_ref;
    do
    {
        my $current_node = ${$list_handle_ref};
        my $current_start = $current_node->GO_start;
        my $current_end = $current_node->GO_end;
        my $current_next_start = $current_node->GO_next_start;

        if(\$global_gff_list == $list_handle_ref)
        {
            $current_start = 0;
            $current_end = 0;
        }
        
        # sometimes we start a little too far in...
        if($first_loop_done == 0)
        {
            if($current_start >= $feat_end)
            {
                my $prev_ref = $current_node->GO_prevGO_ref;
                my $prev_node = ${$prev_ref};
                $current_next_start = $current_start;
                $current_node = $prev_node;
                $list_handle_ref = \$current_node;
                if(\$global_gff_list == $prev_ref)
                {
                    # the prev node is the list head...
                    $current_end = 0;
                    $current_start = 0;
                }
                else
                {
                    # things are relatively normal...
                    $current_start = $current_node->GO_start;
                    $current_end = $current_node->GO_end;
                }
            }
            $first_loop_done = 1;
        }
        
        # if we haven't inserted it by the time we have moved past it,
        # there is no need to keep searching the list
        if($current_start > $feat_start)
        {
            debugFeature("Not keeping feature ".feat2str($$feat_ref), 'no need to keep searching the list') if $debug;
            return $start_ref;
        }
        
        # so, $current_start is less than or equal to
        # the start of the feature being added...
        if($current_end >= $feat_end)
        {
            # the feature to be added is completely contained
            # in an existing feature. no go!
            debugFeature("Not keeping feature ".feat2str($$feat_ref), 'entirely contained in another feature') if $debug;
            return $start_ref;
        }

        # we would like to know if the new feature overlaps with the 
        # end of the current feature
        my $should_add = 0;
        if($feat_end > $current_next_start)
        {
            if($feat_start < $current_end) 
            { 
                # new feature overlaps with only the next feature
                #
                #  CURRENT             NEXT
                # XXXXXXXXX          XXXXXXX
                #        NNNNNNNNNNNNN
                if( ( ( ($feat_end - $current_next_start) + ($current_end - $feat_start) ) / $feat_length ) < $global_shared_olap_cutoff ) { $should_add = 1; }
                
            }
            else
            {
                # new feature overlaps with two existing features
                #
                #  CURRENT             NEXT
                # XXXXXXXXX          XXXXXXX
                #            NNNNNNNNNN
                if( ( ($feat_end - $current_next_start) / $feat_length ) < $global_shared_olap_cutoff) { $should_add = 1; }
            }
        }
        else
        {
            if($feat_start < $current_end) 
            { 
                # new feature overlaps with only the current feature
                #
                #  CURRENT             NEXT
                # XXXXXXXXX          XXXXXXX
                #        NNNNNNNNNN
                if( ( ($current_end - $feat_start) / $feat_length ) < $global_shared_olap_cutoff) { $should_add = 1; }
            }
            else
            {
                # new feature floats in the ether...
                #
                #  CURRENT             NEXT
                # XXXXXXXXX           XXXXXX
                #           NNNNNNN
                $should_add = 1;
            }
        }
        if($should_add == 1)
        {
            debugFeature("Keeping feature ".feat2str($$feat_ref)) if $debug;


            my $new_obj = GffObject->new();
            $new_obj->GO_gffRef($feat_ref);
            $new_obj->GO_start($feat_start);
            $new_obj->GO_end($feat_end);
            insertInList($list_handle_ref, \$new_obj);
            # return a reference to the inserted object to speed up the next insertion...
            return \$new_obj;
        }
    } while(nextInList(\$list_handle_ref) == 1);

    # nothing!
    return \$global_gff_list;
}


sub insertInList {
    #-----
    # Given a node in the LL and a ref to a node to add
    # Insert the new ref AFTER the given LL ref
    #
    my ($ll_ref, $new_ref) = @_;

    my $ll_obj = ${$ll_ref};
    my $new_obj = ${$new_ref};
    
    my $tmp_nxt_ref = $ll_obj->GO_nextGO_ref;
    my $tmp_nxt_obj = ${$tmp_nxt_ref};
    
    $ll_obj->GO_nextGO_ref($new_ref);
    $ll_obj->GO_next_start($new_obj->GO_start);
    $tmp_nxt_obj->GO_prevGO_ref($new_ref);
    
    $new_obj->GO_nextGO_ref($tmp_nxt_ref);
    $new_obj->GO_prevGO_ref($ll_ref);
    $new_obj->GO_next_start($tmp_nxt_obj->GO_start);
}


sub nextInList {
    #-----
    # given a ref to a linked list node, update to reflect the next node in the list
    # return 1 or 0 accordingly
    #
    my ($ll_ref_ref) = @_;
    my $ll_ref = ${$ll_ref_ref};
    my $ll_obj = ${$ll_ref};
    if($ll_obj->GO_nextGO_ref == $global_end_ref) { return 0; }
    # not at the end of the list!
    $$ll_ref_ref = $ll_obj->GO_nextGO_ref;
    return 1;
}


######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "gffs|g:s", "contigs|c:s", "out|o:s", "ann|a:s", );
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

    # check that the user gave us some GFF3 files to parse
    if(!exists($options{'gffs'}))
    {
        print "ERROR: You need to supply at least one gff file to continue!";
        exec("pod2usage $0");
    }

    # check that the user gave us a contig file
    if(!exists($options{'contigs'}))
    {
        print "ERROR: You need to supply the contigs file!";
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

    combineGff3.pl

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

    combineGff3.pl -gffs|g gff_file1[,gff_file2[, ... ]] -contigs|c contigs_file

     -gffs    -g FILE    Gff3 files to parse (comma separated list, in order of "trust"
     -contigs -c FILE    Contigs which were annotated to produce the gff3 files. (Headers must match!)
    [-out     -o FILE]   Output file to write to [default: parsed.gff3]
    [-ann     -a FILE]   File to place fasta type sequences for further annotation [default: todo.fa]
    [-help    -h     ]   Displays basic usage information
      
         
=cut


#!/usr/bin/perl
###############################################################################
#
#    glim2gff3.pl
#
#    Convert glimmer output to a gff3 feature file
#
#    Copyright (C) 2010 Michael Imelfort
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


open(my $GLIM_FILE, "<", $options->{'in'}) or die "Could not open ".$options->{'in'}." $!\n";
open(my $GFF3_FILE, ">", $options->{'out'}) or die "Could not open ".$options->{'out'}." $!\n";
print $GFF3_FILE "##gff-version 3\n";

my $con_header = "";
while(<$GLIM_FILE>)
{
    chomp $_;
    if($_ =~ /^>/)
    {
        # header line
        $_ =~ s/>//;
        $con_header = $_;
    }
    else
    {
        #data line
        $_ =~ s/(\d) /\1\t/g;
        $_ =~ s/ //g;
        my @fields = split /\t/, $_;
        
        next if(int($fields[1]) > int($fields[2]));
        
        my $strand = '-';
        if($fields[3] =~ /^\+/) { $strand = '+'};
        my $phase = abs(int($fields[3]))  - 1;
        $phase = 0;
        
        print $GFF3_FILE "$con_header\t.\torf\t$fields[1]\t$fields[2]\t$fields[4]\t$strand\t$phase\tname=$fields[0]\n";
    }
}

close $GLIM_FILE;
close $GFF3_FILE;

sub checkParams {
    my @standard_options = ( "help+", "in:s", "out:s" );
    my %options;

    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # user must set input and output files
    # 
    exec("pod2usage $0") if (!exists $options{'in'} || !exists $options{'out'});

    return \%options;
}


sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2010 Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    glim2gff3.pl

=head1 COPYRIGHT

   copyright (C) 2010 Michael Imelfort

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

    glim2gff3.pl

=head1 SYNOPSIS

    glim2gff3.pl -in GLIM_FILE -out GFF3_FILE [-help]

      -in GLIM_FILE      Output from glimmer3
      -out GFF3_FILE     gff3 file to create
      [-help]            Displays basic usage information
         
=cut

#!/usr/bin/perl
###############################################################################
#
#    extractECnumbers.pl
#    
#    parse mannortator output for genes containing EC numbers
#
#    Copyright (C) 2011 Connor Skennerton
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

use warnings;
use strict;
use Pod::Usage;

pod2usage if (!$ARGV[0]);
open(IN, "<", $ARGV[0]) or die $!;

while(<IN>)
{
    if ($_ =~ /KEGG_ENZYME:(\d+\.\d+\.\d+\.\d+)/ )
    {
        print "$1\n";
        next;
        
    }
    elsif ($_ =~ /\[EC.*(\d+\.\d+\.\d+\.\d+)\]$/)
    {
        print "$1\n";
        next;
    }
    elsif ($_ =~ /\(EC.*(\d+\.\d+\.\d+\.\d+)\)/)
    {
        print "$1\n";
        next;
    }
}

exit;

__DATA__

=head1 NAME

    extractECnumbers.pl

=head1 COPYRIGHT

   copyright (C) 2011 Connor Skennerton

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

	Get your mannotator output file, probably called "mannotated.gff3".
	Runs a regular expression on the data to extract all of the EC numbers.
	Prints to stdout a list of one EC number per line.  Can then use the output 
	in the script "markKeggPathway.rb

=head1 SYNOPSIS

	extractECnumbers.pl mannotated.gff3|any_other.gff3 
         
=cut
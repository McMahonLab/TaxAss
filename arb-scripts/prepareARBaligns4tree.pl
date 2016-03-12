#!/opt/local/bin/perl

#


use strict;
use Benchmark qw(:all);
use Data::Dumper;


my $USAGE=<<USAGE;

 * prepareARBaligns4tree.pl *
 
   Replaces the U's in a sequence with T's and "."s with "-"s, without touching the defline.
   Useful for ARB alignments that need to go into tree-building programs or to use in 16S 
   taxonomy trainingsets. 
 
   usage:

   cat myseqs.fa | prepareARBaligns4tree.pl > mycleanseqs.fa

    - for getting this help -

    replaceUs.pl -h
    replaceUs.pl --help


USAGE

my $informat   = shift;

if ($informat eq "-h")     { die $USAGE; }
if ($informat eq "--help") { die $USAGE; }

while (my $line = <STDIN>) {
	
	

	if ($line !~ /^(\>)/) {
	
	$line =~ s/U/T/g;
	$line =~ s/\./-/g;
	}

        print "$line";
}

exit; 


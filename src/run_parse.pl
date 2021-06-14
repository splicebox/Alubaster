#!/usr/bin/perl
use strict;

($#ARGV==6) or die "Usage: $0 set genome genomehdrs alufasta matesfasta infile outfile \n";
my $set = shift;
my $GENOME = shift;
my $GENOMEHDRS = shift;
my $ALUFA = shift;
my $MATESFA = shift;
my $infile = shift;
my $outfile = shift; 

my $scripts = $ENV{SCRIPTS} ;
 $| = 1;

# reset outfile
open(F, ">$outfile") or die "could not open outfile $outfile.\n";
close(F);

open(F, "<$infile") or die "could not open screened file $infile.\n";
while (<F>) {
   chomp; `echo \"$_\" > tmp.test`;
   `${scripts}/parse_screened.pl $set $GENOME $GENOMEHDRS $ALUFA $MATESFA < tmp.test > tmp.out; cat tmp.out >> $outfile`;
}
close(F);


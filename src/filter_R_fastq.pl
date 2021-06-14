#!/usr/bin/perl
use strict;

($#ARGV==0) or die "Usage: $0 mates\n";
my $matesFile = shift;

# read in a set of >readid and print out only those that are not simultaneosuly listed in both fastq files
# (i.e., anchor does not contain repeat); to simplify, simpy eliminate from the Yes list those reads that appear in the mates file.

my %RMate;
open(F, "<$matesFile") or die "Could not open $matesFile.\n"; # xxx.nonconcordant.mates.fastq";
while (<F>) {
   next if !/length/;
   /^@(\S+) / or die "died. $_";
   $RMate{$1} = 1;
}
close(F);

# read in list of SIGNAL mates
while (<>) {
   /^(\S+).(\d):\S$/ or die "died. $_";
   my ($readid,$idx) = ($1,$2);
   defined($RMate{$1 . "/" . (3-$idx)}) or print $_;
}

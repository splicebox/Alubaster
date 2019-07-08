#!/usr/bin/perl
use strict;

my %Seen;
while (<>) { 
  chomp;
  # ENST00000291107.5;chr17:1091665-1091850	2	SRR1284895.48542266/2+,SRR1284895.180677625/1+
  /^\S+;(\S+:\d+\-\d+)\s/ or die "died. $_";
  next if defined($Seen{$1});
  print $_, "\n";
  $Seen{$1} = 1;
}

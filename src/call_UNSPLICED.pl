#!/usr/bin/perl
use strict;

($#ARGV==1) or die "Usage: $0 cov pctid\n";

my $COV = shift;
my $PCTID = shift;

my %Seen; # 1 if good match found ('Yes')

my $last_key;
# read in a sorted UNSPLICED stats file
while (<>) {
  next if /^#/;


  # >ERR188081.33635998/1:+ >ENSG000001:ENST00000613788.1;chr1:160684887-160685189 75 803 m1fwd + 1 38 21 57 38 92 1-38 AU 0-37 37-0 Not valid: wrong region, correct ori
  # >ERR188081.28372058/2:- >ENSG000001:ENST00000613788.1;chr1:160684887-160685189 75 803 m1rev - 1 28 593 620 28 96 1-28 UA 0-28 0-28 Not valid: wrong region, correct ori
  /^(\S+) (\S+) (\d+) \d+ (\S+) (\S) \d+ \d+ \d+ \d+ (\d+) (\d+) (\S+) (\S+) \S+ \S+$/ or die "died. $_";
  my ($readid,$anchorid,$readlen,$m1ori,$ori,$cov,$pctid,$exons,$type) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
  $readid =~/^\S+(\S)$/ or die "died readid ($readid).\n";
  my $readori = $1;

  next if (defined($Seen{"$readid $anchorid"}));
  if (defined($last_key) && ($last_key ne "$readid $anchorid")) {
     if (!defined($Seen{$last_key})) {
        print "$last_key No\n";
     }
  }
  $last_key = "$readid $anchorid";

  # NOTE: m1ori and type are redundant; m1ori refers to the type of region, not to the original read orientation!!!
  # i.e., m1rev == UA; m1fwd = AU
  # eliminate the wrong set of regions (facing in the wrong direction)
  next if (($readori eq "+") && !($m1ori eq "m1rev"));
  next if (($readori eq "-") && !($m1ori eq "m1fwd"));

  # Good matches: $readori = '+'; m1rev; UA; ori =='+';
  #               $readori = '-'; m1fwd; AU; ori =='-';
  if ((($cov/(1.0*$readlen))>=$COV) && ($pctid>=$PCTID) && 
      ((($type eq "AU") && ($ori eq "-")) || (($type eq "UA") && ($ori eq "+"))) &&
      !($exons=~/,/)) {
      print "$readid $anchorid Yes\n";
      $Seen{"$readid $anchorid"} = 1;
  }
# $last_key = "$readid $anchorid";
}
if (defined($last_key) && !defined($Seen{$last_key})) {
  print "$last_key No\n";
}


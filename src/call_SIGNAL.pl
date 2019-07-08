#!/usr/bin/perl
use strict;

($#ARGV==2) or die "Usage: $0 cov pctid rmin\n";

my $COV = shift;    # 80
my $PCTID = shift;  # 85
my $RMIN = shift;   # 10
my $XPCTID = 93;

my $ACOVMIN = 0.7;
my $APCTIDMIN = 92.0;

# NOTE: m1ori and type are redundant; m1ori refers to the type of region, not to the original read orientation!!! Exception: ARA

my %isGood;
my %isExon;

my $last_key;
# read in a sorted SIGNAL stats file
while (<>) {
  next if /^#/;

  # edef ddef elen glen ori efrom eto dfrom dto cov pctid exons type readsR-A-C reads-by-type
  # For region, compare readori and m1{fwd,rev} - should be opposite; for ori, compare readori with matchori, shoudl be the same
  # >ERR188081.26338794/1:+ >ENSG000001:ENST00000367042.4;chr1:207793519-207795513 75 2371 m1rev - 13 75 2312 2371 63 87 13-75 CAr 60-0-0 0-0-60 Bad match: correct region, wrong ori
  # >ERR188081.35572929/2:+ >ENSG000001:ENST00000367042.4;chr1:207793519-207795513 75 2371 m1rev + 1 39 1195 1232 39 89 1-39 CAr 0-38-0 0-38-0 Good match: correct region, correct ori
  # >ERR188081.38553975/2:- >ENSG000001:ENST00000367042.4;chr1:207793519-207795513 75 2371 m1rev - 1 43 1264 1304 43 95 1-43 CAr 0-41-0 0-41-0 Bad match: wrong region, correct ori
  # >ERR188081.16021143/2:- >ENSG000001:ENST00000367042.4;chr1:207793519-207795513 75 2370 m1rev + 1 15 1910 1926 15 88 1-15 CAr 0-17-0 0-17-0 Bad match: wrong region, wrong ori
  # >ERR188081.24688841/2:+ >ENSG000001:ENST00000367042.4;chr1:207793519-207795513 75 2370 m1rev - 23 75 2317 2370 53 92 23-75 CAr 54-0-0 0-0-54 Bad match: correct region, wrong ori

  (/^(\S+) (\S+) (\d+) \d+ (\S+) (\S) \d+ \d+ \d+ \d+ (\d+) (\d+) (\S+) (\S+) (\d+)\-(\d+)\-(\d+) \S+$/) or die "died. $_";
  my ($readid,$anchorid,$readlen,$m1ori,$ori,$cov,$pctid,$exons,$type,$Rnum,$Anum,$Cnum) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12);
  $readid =~ /^\S+(\S)$/ or die "readid ($readid).\n";
  my $readori = $1;

  if (defined($last_key) && ($last_key ne "$readid $anchorid")) {
     if (defined($isExon{$last_key})) {
        print "$last_key No\n";
     } elsif (defined($isGood{$last_key})) {
        print "$last_key Yes\n";
     } else {
        print "$last_key No\n";
     }
  }

  $last_key = "$readid $anchorid";
  next if (!/A[Rr]C/ && !/C[Rr]A/);

  # NOTE: m1ori and type are redundant; m1ori refers to the type of region, not to the original read orientation!!!
  # i.e., m1fwd == ARC,RAC,ARA; m1rev = CRA,CAR,(ARA)
  # eliminate the wrong set of regions (facing in the wrong direction)
  next if (($readori eq "+") && !($type=~/A[rR]A/) && !($m1ori eq "m1rev"));
  next if (($readori eq "-") && !($type=~/A[rR]A/) && !($m1ori eq "m1fwd"));
  if ((($Anum>=$ACOVMIN*$readlen) || ($Cnum>=$ACOVMIN*$readlen)) && ($pctid>=$APCTIDMIN)) { $isExon{"$readid $anchorid"} = 1; }

  # Good matches: $readori = '+'; m1rev; ARC,ARA,RAC; ori =='+';
  #               $readori = '-'; m1fwd; AG; ori =='-';
  if ((($cov/(1.0*$readlen))>=$COV) && ($pctid>=$PCTID) && ($Rnum>=$RMIN) && 
       ((($type=~/C[rR]A/) && ($ori eq "+")) || (($type=~/A[rR]C/) && ($ori eq "-")) ||    # next exon, Alu exonization
        (($type=~/CA[rR]/) && ($ori eq "+")) || (($type=~/CA[rR]/) && ($ori eq "-")) ||
         ($type=~/A[rR]A/))) {
      if ($Rnum>=$RMIN) { $isGood{"$readid $anchorid"} = 1; }
      elsif ($pctid>=$XPCTID) { $isExon{"$readid $anchorid"} = 1; }
  }
# $last_key = "$readid $anchorid";
}
if (defined($last_key)) {
  if (defined($isGood{$last_key}) && !defined($isExon{$last_key})) {
     print "$last_key Yes\n";
   } else {
     print "$last_key No\n";
   }
}

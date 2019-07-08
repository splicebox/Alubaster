#!/usr/bin/perl
use strict;

($#ARGV==1) or die "Usage: $0 cov pctid [-g]\n";

my $COV = shift;
my $PCTID = shift;

my %Seen;

# NOTE: m1ori and type are redundant; m1ori refers to the type of region, not to the original read orientation!!!

my $last_key;
# read in a sorted REGION stats file
while (<>) {
  next if /^#/;


  # >ERR188081.22436574/2:+ >ENSG0000001:ENST00000629041.1;chr1:45581219-45581278 75 5059 m1rev + 1 75 4933 5007 75 93 1-75 GA 68-7 68-7  Valid: correct region, correct ori
  # >ERR188081.907544/2:+ >ENSG0000001:ENST00000433834.4;chr1:19344358-19344434 75 930 m1fwd + 24 43 128 147 20 85 24-43 AG 20-0 0-20 Not valid: wrong region, correct ori
  /^(\S+) (\S+) (\d+) \d+ (\S+) (\S) \d+ \d+ \d+ \d+ (\d+) (\d+) (\S+) (\S+) \S+ \S+$/ or die "died. $_";
  my ($readid,$anchorid,$readlen,$m1ori,$ori,$cov,$pctid,$exons,$type) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
  $readid =~/^\S+:(\S)$/ or die "died (readid). $readid";
  my $readori = $1;
  next if (defined($Seen{"$readid $anchorid"}));
  if (defined($last_key) && ($last_key ne "$readid $anchorid")) {
     if (!defined($Seen{$last_key})) {
        print "$last_key No .\n";
     }
  }

  $last_key = "$readid $anchorid";
  # NOTE: m1ori and type are redundant; m1ori refers to the type of region, not to the original read orientation!!!
  # i.e., m1rev == GA; m1fwd = AG
  # eliminate the wrong set of regions (facing in the wrong direction)
  next if (($readori eq "+") && !($m1ori eq "m1rev"));
  next if (($readori eq "-") && !($m1ori eq "m1fwd"));
  
  # Good matches: $readori = '+'; m1rev; GA; ori =='+';
  #               $readori = '-'; m1fwd; AG; ori =='-';
  if ((($cov/(1.0*$readlen))>=$COV) && ($pctid>=$PCTID) &&
      ((($type eq "AG") && ($ori eq "-")) || (($type eq "GA") && ($ori eq "+")))) {
      print "$readid $anchorid Yes";
      print (($exons=~/,/) ? " Spliced\n" : " Unspliced\n");
      $Seen{"$readid $anchorid"} = 1;
  }
# $last_key = "$readid $anchorid";
}
if (defined($last_key) && !defined($Seen{$last_key})) {
  print "$last_key No .\n";
}


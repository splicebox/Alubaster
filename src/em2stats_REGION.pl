#!/usr/bin/perl
use strict;

#sim4begin
#0[75-0-0] 1[0-5059] <70-0-93-forward-unknown>
#edef=>ERR188081.22436574/2:- HWI-962:71:D0PEYACXX:3:1308:18505:84334/2
#ddef=>"ENST00000629041.1";chr1:45581219-45581278 m1rev:(m1-5000):45576219 (GA)
#1-75 (4933-5007) <70-0-93>
#caggGccccaccaccatgtccggctaCtttttgaatttttagtagagatgggAtttCaccatgttggctaggaGg
#caggCccccaccaccatgtccggctaAtttttgaatttttagtagagatgggTtttTaccatgttggctaggaTg
#sim4end

print "#edef ddef elen glen m1ori ori efrom eto dfrom dto cov pctid exons type readsRA reads-by-type\n";

# read in a sim4db file
while (<>) {
  /^sim4begin/ or die "died (sim4begin). $_";
  $_ = <>; chomp; 
  /^\d+\[(\d+)\-\d+\-\d+\] \d+\[(\d+)\-(\d+)\] <\d+\-\d+\-(\d+)\-(\S+)\-\S+>/ or die "died (header). $_";
  my ($elen,$GF,$GT,$pctid,$ori) = ($1,$2,$3,$4,$5);
  my $glen = $GT-$GF;
  $ori = ($ori eq "forward") ? "+" : "-";

  $_ = <>; chomp;
  /^edef=(\S+)\s/ or die "died (edef). $_";
  my $tid = $1;
  $_ = <>; chomp;
  my $dline = $_;
  # ddef=>"ENST00000629041.1";chr1:45581219-45581278 m1rev:(m1-5000):45576219 (GA)
# /^ddef=(\S+) (m1\w\w\w):\(\S+\):(\d+) \((\S\S)\)$/ or die "died (ddef). $_";     HERE, last parenthesis
  /^ddef=(\S+) (m1\w\w\w):\(\S+\):(\d+) \((\S\S)/ or die "died (ddef). $_";
  my ($ganame,$m1ori,$pos,$type) = ($1,$2,$3,$4);

  my $cov = 0;
  my ($efrom,$eto,$gfrom,$gto);
  my @exons = ();
  my @gexons = ();
  while (<>) {
    last if /^sim4end/; 
    next if /^[\-actgnACGTN]/;

    /^(\d+)\-(\d+) \((\d+)\-(\d+)\)/ or die "died (cDNA). $_";
    my ($f,$t,$F,$T) = ($1,$2,$3,$4);
    if (!defined($efrom)) {
       $efrom = $f; $gfrom = $GF+$F;
    }
    $eto = $t; $gto = $GF+$T;
    if ($ori eq "-") { ($f,$t) = ($elen-$t+1,$elen-$f+1); }
    push (@exons,"$f-$t");
    push (@gexons, "$F-$T");

    $cov += ($t-$f+1);
  }
  if ($ori eq "-") {
     my $tmp = $elen-$efrom+1;
     $efrom = $elen-$eto+1;
     $eto = $tmp;
     @exons = reverse(@exons);
  }

  # parse these from the dline
  # ddef=>"ENST00000629041.1";chr1:45581219-45581278 m1rev:(m1-5000):45576219 (GA)
  # ddef=>ENSG00000280279.1:ENST00000623180.1;chr17:65830-65887 m1fwd:(m1+71366+100):71466 (AG)
# $dline =~ /\S+;\S+:(\d+)\-(\d+) m1\S+:\S+:(\d+) / or die "died (dline). $dline";  HERE, last parenthesis
  $dline =~ /\S+;\S+:(\d+)\-(\d+) m1\S+:\S+:(\d+) / or die "died (dline). $dline";
  my ($L1,$L2);
  if ($type eq "AG") {
     # L1 is for A
     $L1 = $2-$1+1; $L2 = $3-$2;
  } elsif ($type eq "GA") {
     # L2 is for A
     $L1 = $1-$3; $L2 = $2-$1+1;
  } else {
     die "Unrecognized type. $type";
  }

  # Intervals: [1 .. L1]; [L1+1 .. L1+L2];
  my ($reads1,$reads2) = (0,0);
  foreach my $x (@gexons) {
     $x =~ /^(\d+)\-(\d+)$/ or die "died (exon). $x";
     my ($b,$e) = ($1,$2);
     $reads1 += &overlap($b,$e,1,$L1); 
     $reads2 += &overlap($b,$e,$L1+1,$L1+$L2); 
  }

  # Types: 
  my $rstring = "";  #GA
  if ($type =~ /^GA$/) { $rstring = "$reads1-$reads2"; }
  elsif ($type =~ /^AG$/) { $rstring = "$reads2-$reads1"; }
  else {
    die "unrecognized type. $type\n";
  }

  # edef ddef elen glen m1ori ori efrom eto dfrom dto cov pctid exons type readsRA reads-by-type
  print "$tid $ganame $elen $glen $m1ori $ori $efrom $eto $gfrom $gto $cov $pctid ", join(',', @exons), " $type $rstring $reads1-$reads2\n";
}


sub overlap {
   my ($a,$b,$u,$v) = @_;
   my $t = &min($b,$v) - &max($a,$u) +1;
   return ($t>=0) ? $t : 0;
}

sub min { return ($_[0]<=$_[1]) ? $_[0] : $_[1]; }

sub max { return ($_[0]>=$_[1]) ? $_[0] : $_[1]; }

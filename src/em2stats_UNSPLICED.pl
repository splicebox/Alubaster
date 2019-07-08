#!/usr/bin/perl
use strict;

#sim4begin
#0[75-30-0] 0[0-2495] <15-0-88-forward-unknown>
#edef=>ERR188081.16021143/2:+ HWI-962:71:D0PEYACXX:3:1302:10933:7241/2
#ddef=>"ENST00000367042.4";chr1:207793519-207795513 m1fwd:(m1+500):207796013 (AU)
#1-15 (1817-1833) <15-0-88>
#tt-tc-ggaagtttgtt
#ttCtcTggaagtttgtt
#sim4end


print "#edef ddef elen glen m1ori ori efrom eto dfrom dto cov pctid exons type readsU-A reads-by-type\n";

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
  # ddef=>"ENST00000367042.4";chr1:207793519-207795513 m1fwd:(m1+500):207796013 (AU)
# /^ddef=(\S+) (m1\S+):\S+:\d+ \((\S\S)\)$/ or die "died (ddef). $_";  HERE, last parenthesis
  /^ddef=(\S+) (m1\S+):\S+:\d+ \((\S\S)/ or die "died (ddef). $_";
  my ($ganame,$m1ori,$type) = ($1,$2,$3);

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
  # ddef=>"ENST00000629041.1";chr1:45581219-45581278 m1fwd:(m1+500):45581778 (AU)
  # ddef=>ENSG00000280279.1:ENST00000623180.1;chr17:65830-65887 m1fwd:(m1+500):66387 (AU)

  my ($L1,$L2);
  $dline =~ / m1\S+:\S+:(\d+) / or die "died. (dline 1). $dline";
  my $tmpto = $1;
  $dline =~ /\S+;\S+:(\d+)\-(\d+) / or die "died (dline 2). $dline";
  if ($type eq "AU") {
     $L1 = $2-$1+1; $L2 = $tmpto-$2;
  } elsif ($type eq "UA") {
     $L1 = $1-$tmpto; $L2 = $2-$1+1;
  } else {
     die "unrecognized type. $type";
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
  my $rstring = "";  #UA
  if ($type =~ /^AU$/) { $rstring = "$reads2-$reads1"; }
  elsif ($type =~ /^UA$/) { $rstring = "$reads1-$reads2"; }
  else {
    die "unrecognized type. $type\n";
  }

  # edef ddef elen glen m1ori ori efrom eto dfrom dto cov pctid exons type readsU-A reads-by-type
  print "$tid $ganame $elen $glen $m1ori $ori $efrom $eto $gfrom $gto $cov $pctid ", join(',', @exons), " $type $rstring $reads1-$reads2\n";
}


sub overlap {
   my ($a,$b,$u,$v) = @_;
   my $t = &min($b,$v) - &max($a,$u) +1;
   return ($t>=0) ? $t : 0;
}

sub min { return ($_[0]<=$_[1]) ? $_[0] : $_[1]; }

sub max { return ($_[0]>=$_[1]) ? $_[0] : $_[1]; }

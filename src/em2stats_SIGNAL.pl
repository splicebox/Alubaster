#!/usr/bin/perl
use strict;

#sim4begin
#0[75-0-0] 1[0-381] <57-0-83-complement-unknown>
#edef=>ERR188081.22436574/2:- HWI-962:71:D0PEYACXX:3:1308:18505:84334/2
#ddef=>"ENST00000629041.1";chr1:45581219-45581278 m1fwd:I:45581283-45581321 (ARC)
#3-69 (146-213) <57-0-83>
#tcctAgcCaacaTggtgaaaTcccAtctctactaaaaatTcaaaaaG-tagccggAcAtggtggTggg
#tcctGgcTaacaCggtgaaaCcccGtctctactaaaaatAcaaaaaATtagccggGcGtggtggCggg
#sim4end

#sim4begin
#0[75-0-0] 1[0-381] <57-0-83-complement-unknown>
#edef=>ERR188081.22436574/2:- HWI-962:71:D0PEYACXX:3:1308:18505:84334/2
#ddef=>"ENST00000629041.1";chr1:45581219-45581278 m1fwd:I:45581283-45581321 (ARC)
#3-69 (146-213) <57-0-83>
#tcctAgcCaacaTggtgaaaTcccAtctctactaaaaatTcaaaaaG-tagccggAcAtggtggTggg
#tcctGgcTaacaCggtgaaaCcccGtctctactaaaaatAcaaaaaATtagccggGcGtggtggCggg
#sim4end

print "#edef ddef elen glen m1ori ori efrom eto dfrom dto cov pctid exons type readsR-A-C reads-by-type\n";

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
  # ddef=>"ENST00000629041.1";chr1:45581219-45581278 m1fwd:I:45581283-45581321 (ARC)
# /^ddef=(\S+) (m1\S\S\S):\S+ \((\S+)\)$/ or die "died (ddef). $_";  HERE, last parenthesis
  /^ddef=(\S+) (m1\S\S\S):\S+ \(([^\)]+)/ or die "died (ddef). $_";
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
  # ddef=>"ENST00000629041.1";chr1:45581219-45581278 m1fwd:I:45581283-45581321 (ARC)
  # NEW: ddef=>ENSG00000280279.1:ENST00000623180.1;chr17:65830-65887 m1fwd:I:71366-71556 (ARC)
  my ($L1,$L2,$L3);
  $dline =~ />\S+;\S+:(\d+)\-(\S+) m1\w\w\w:\S:(\d+)\-(\d+) / or die "died (dline). $dline";
  my ($from,$to,$F,$T) = ($1,$2,$3,$4);
  if ($type =~ /A[Rr]C/) { # I
     $L1 = $to-$from+1; $L3 = $T-$F+1; $L2 = $glen-$L1-$L3;
  } elsif ($type =~ /C[rR]A/) { # I
     $L1 = $T-$F+1; $L3 = $to-$from+1; $L2 = $glen-$L1-$L3;
  } elsif ($type =~ /CA[Rr]/)  { # X
     $L1 = $T-$F+1; $L2 = $to-$from+1; $L3 = $glen-$L1-$L2; 
  } elsif ($type =~ /[Rr]AC/) { # X
     $L2 = $to-$from+1; $L3 = $T-$F+1; $L1 = $glen-$L1-$L2;
  } elsif ($type =~ /A[rR]A/) { # X
     $L1 = $L3 = $to-$from+1; $L2 = $glen-$L1-$L3;
  } else {
     die "unrecognized type. $type\n";
  }

  # Intervals: [1 .. L1]; [L1+1 .. L1+L2]; [L1+L2+1 .. L1+L2+L3]
  my ($reads1,$reads2,$reads3) = (0,0,0);
  foreach my $x (@gexons) {
     $x =~ /^(\d+)\-(\d+)$/ or die "died (exon). $x";
     my ($b,$e) = ($1,$2);
     $reads1 += &overlap($b,$e,1,$L1); 
     $reads2 += &overlap($b,$e,$L1+1,$L1+$L2); 
     $reads3 += &overlap($b,$e,$L1+$L2+1,$L1+$L2+$L3);
  }

  # Types: 
  my $rstring = "";  #RAC
  if ($type =~ /^C[Rr]A$/) { $rstring = "$reads2-$reads3-$reads1"; }
  elsif ($type =~ /^A[Rr]C$/) { $rstring = "$reads2-$reads1-$reads3"; }
  elsif ($type =~ /^CA[Rr]$/) { $rstring = "$reads3-$reads2-$reads1"; }
  elsif ($type =~ /^[Rr]AC$/) { $rstring = "$reads1-$reads2-$reads3"; }
  elsif ($type =~ /^A[rR]A$/) { $rstring = "$reads2-$reads1-$reads3"; }

  else {
    die "unrecognized type. $type\n";
  }

  # edef ddef elen glen m1ori ori efrom eto dfrom dto cov pctid exons type readsR-A-C reads-by-type
  print "$tid $ganame $elen $glen $m1ori $ori $efrom $eto $gfrom $gto $cov $pctid ", join(',', @exons), " $type $rstring $reads1-$reads2-$reads3\n";
}


sub overlap {
   my ($a,$b,$u,$v) = @_;
   my $t = &min($b,$v) - &max($a,$u) +1;
   return ($t>=0) ? $t : 0;
}

sub min { return ($_[0]<=$_[1]) ? $_[0] : $_[1]; }

sub max { return ($_[0]>=$_[1]) ? $_[0] : $_[1]; }

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
  # Update these for unmapped_realign to allow ids that have only 1 field
  my $tid;
  if (/^edef=(\S+)\s/) { $tid = $1; }
  elsif (/^edef=(\S+)$/) { $tid = $1; }
  else { die "died (edef). $_"; }
  $_ = <>; chomp;
  my $dline = $_;
  # ddef=>"ENST00000629041.1";chr1:45581219-45581278 m1fwd:I:45581283-45581321 (ARC)
# /^ddef=(\S+) (m1\S\S\S):\S+ \((\S+)\)$/ or die "died (ddef). $_";  
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
  # ddef="ENST00000629041.1";chr1:45581219-45581278 m1fwd:I:45581283-45581321 (ARC)
  # NEW: ddef=ENSG00000280279.1:ENST00000623180.1;chr17:65830-65887 m1fwd:I:71366-71556 (ARC)
  # Note: L1, L2 and L3 are the lengths of the three segments of geomic sequence, in that order, from left-to-right (irrespective of their label)
  my ($L1,$L2,$L3);
  my @LFrom; my @LTo; my @GFrom; my @GTo; my @GgaID;
  $dline =~ /\S+;(\S+):(\d+)\-(\S+) m1\w\w\w:\S:(\d+)\-(\d+) / or die "died (dline). $dline";
  my ($chrom,$from,$to,$F,$T) = ($1,$2,$3,$4,$5);
  if ($type =~ /A[Rr]C/) { # I
     $L1 = $to-$from+1; $L3 = $T-$F+1; $L2 = $glen-$L1-$L3;
     $GFrom[0] = $from; $GTo[0] = $to; $GgaID[0] = $chrom;
     $GFrom[1] = 1; $GTo[1] = $L2; $GgaID[1] = "Alu";
     $GFrom[2] = $F; $GTo[2] = $T; $GgaID[2] = $chrom;
  } elsif ($type =~ /C[rR]A/) { # I
     $L1 = $T-$F+1; $L3 = $to-$from+1; $L2 = $glen-$L1-$L3;
     $GFrom[0] = $F; $GTo[0] = $T; $GgaID[0] = $chrom;
     $GFrom[1] = 1; $GTo[1] = $L2; $GgaID[1] = "Alu";
     $GFrom[2] = $from; $GTo[2] = $to; $GgaID[2] = $chrom;
  } elsif ($type =~ /CA[Rr]/)  { # X
     $L1 = $T-$F+1; $L2 = $to-$from+1; $L3 = $glen-$L1-$L2; 
     $GFrom[0] = $F; $GTo[0] = $T; $GgaID[0] = $chrom;
     $GFrom[1] = $from; $GTo[1] = $to; $GgaID[1] = $chrom;
     $GFrom[2] = 1; $GTo[2] = $L3; $GgaID[2] = "Alu";
  } elsif ($type =~ /[Rr]AC/) { # X
     $L2 = $to-$from+1; $L3 = $T-$F+1; $L1 = $glen-$L1-$L2;
     $GFrom[0] = 1; $GTo[0] = $L1; $GgaID[0] = "Alu";
     $GFrom[1] = $from; $GTo[1] = $to; $GgaID[1] = $chrom;
     $GFrom[2] = $F; $GTo[2] = $T; $GgaID[2] = $chrom;
  } elsif ($type =~ /A[rR]A/) { # X
     $L1 = $L3 = $to-$from+1; $L2 = $glen-$L1-$L3;
     $GFrom[0] = $from; $GTo[0] = $to; $GgaID[0] = $chrom;
     $GFrom[1] = 1; $GTo[1] = $L2; $GgaID[1] = "Alu";
     $GFrom[2] = $from; $GTo[2] = $to; $GgaID[2] = $chrom;
  } else {
     die "unrecognized type. $type\n";
  }
  @LFrom = (1,$L1+1,$L1+$L2+1);
  @LTo = ($L1,$L1+$L2,$L1+$L2+$L3);

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

  # now calculate the breakdown of the read sequence by chromosome versus Alu portions
  my $readComposition = "";
  foreach my $x (@gexons) {
    $x =~ /^(\d+)\-(\d+)$/ or die "died. $_"; 
    my ($xf,$xt) = ($1,$2);
    # avoid a false start when the first component is empty (1-0)
    my $i0 = &find_idx($xf,$L1,$L2,$L3);
    my $i1 = &find_idx($xt,$L1,$L2,$L3);


    $readComposition .= $GgaID[$i0] . ":" . ($xf-$LFrom[$i0]+$GFrom[$i0]);
    foreach (my $i=$i0; $i<$i1; $i++) {
       $readComposition .=  "-" . $GTo[$i] . ";" . $GgaID[$i+1] . ":" . $GFrom[$i+1]; 
    }
    $readComposition .= "-" . ($xt-$LFrom[$i1]+$GFrom[$i1]) . ";";
  }

  # edef ddef elen glen m1ori ori efrom eto dfrom dto cov pctid exons type readsR-A-C reads-by-type
  print "$tid $ganame $elen $glen $m1ori $ori $efrom $eto $gfrom $gto $cov $pctid ", join(',', @exons), " $type $rstring $reads1-$reads2-$reads3 $readComposition\n";
}

sub find_idx { # val L1 L2 L3
   my ($val,$l1,$l2,$l3) = @_;

   if (1<=$val && $val<=$l1) { return 0; }
   if ($l1+1<=$val && $val<=$l1+$l2) { return 1; }
   if ($l1+1<=$val && $val<=$l1+$l2+$l3) { return 2; }

   die "Illegal coordinate in find_idx. $val ($l1,$l2,$l3)\n";
}

sub overlap {
   my ($a,$b,$u,$v) = @_;
   my $t = &min($b,$v) - &max($a,$u) +1;
   return ($t>=0) ? $t : 0;
}

sub min { return ($_[0]<=$_[1]) ? $_[0] : $_[1]; }

sub max { return ($_[0]>=$_[1]) ? $_[0] : $_[1]; }



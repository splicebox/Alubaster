#!/usr/bin/perl
use strict;

#sim4begin
#120[510-0-0] 4[677-5187] <510-0-100-forward-unknown>
#edef=>gi|18401728|ref|NM_103637.1| Arabidopsis thaliana uncharacterized protein mRNA, complete cds
#ddef=>scaffold_4
#1-510 (2001-2510) <510-0-100>
#sim4end


# read in a polished-good file, filtered by %id
while (<>) {
  /^sim4begin/ or die "begin .$_";
  $_ = <>; chomp; 
  /^\d+\[(\d+)\-\d+\-\d+\] \S+ <\d+\-\d+\-(\d+)\-(\S+)\-/ or die "died (header). $_";
  my ($len,$pctid,$relori) = ($1,$2,$3);
  $relori = ($relori eq "complement") ? "-" : "+";
  $_ = <>; chomp;
  if ( ($_ !~ m/^edef=(\S+)\s/) and ($_ !~ m/^edef=(\S+)$/) ){die "died (edef). $_";}
 
  my $tid = $1;
  $_ = <>; chomp;
  my $ganame;
  if (/^ddef=(\S+)$/) { $ganame = $1; }
  elsif (/^ddef=(\S+)\s/) { $ganame = $1; } 
  else { die "died (ddef). $_"; }
  my $cov = 0;
  my ($beg,$end,$gabeg,$gaend);
  while (<>) {
    last if /^sim4end/; 
    next if /^[\-actgnACGTN]/;
    /^(\d+)\-(\d+) \((\d+)\-(\d+)\) / or die "died (cDNA). $_";
    my ($a,$b,$c,$d) = ($1,$2,$3,$4);
    $cov += ($b-$a+1);
    if (!defined($beg)) { $beg = $a; $gabeg = $c; }
    $end = $b;
    $gaend = $d;
  }
  if ($relori eq "-") { my $tmp = $len-$beg+1; $beg = $len-$end+1; $end = $tmp; }
  print "$tid $ganame $beg $end $gabeg $gaend $cov $len $relori ";
  printf "%1.2f %d\n", 100*$cov/(1.0*$len), $pctid;
}

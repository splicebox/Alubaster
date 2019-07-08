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
  /^edef=>(\S+)\s/ or die "died (edef). $_";
  my $tid = $1;
  $_ = <>; chomp;
  /^ddef=>(\S+)$/ or die "died (ddef). $_";
  my $ganame = $1;
  my $cov = 0;
  while (<>) {
    last if /^sim4end/; 
    next if /^[\-actgnACGTN]/;
    /^(\d+)\-(\d+) / or die "died (cDNA). $_";
    $cov =+ ($2-$1+1);
  }
  print "$tid $ganame $cov $len $relori ";
  printf "%1.2f %1.2n\n", 100*$cov/(1.0*$len), $pctid;
}

    


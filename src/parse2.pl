#!/usr/bin/perl
use strict;

my ($aln,$good);
while (<>) {
    $good = 1;
   /^sim4begin/ or die "died (sim4begin).$_";
   $aln = $_;

   $_ = <>;
   /^\d+\[/ or die "died (header). $_";
   $aln .= $_;

   $_ = <>;
   /^edef/ or die "died (edef). $_";
   $aln .= $_;

   $_ = <>; 
   /^ddef/ or die "died (ddef). $_";
   $aln .= $_;
   ((/m1fwd/ || /m1rev/) && / \(\S\S\S.?$/) or $good = 0;
#  ((/m1fwd/ || /m1rev/) && / \(\S\S.?$/) or $good = 0;
   
   $_ = <>;
   while (!/^sim4end/) {
      $aln .= $_;
      $_ = <>;
   }
   $aln .= $_;

   if ($good) { print $aln; }
   else { print STDERR $aln; }
}
   

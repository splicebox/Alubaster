#!/usr/bin/perl
use strict;


# read one alignment record at a time (sim4begin,sim4end) and parse/correct

my $aln;

while (<>) {

   next if !/^sim4begin/;

   if (/^sim4begin/) {
      $aln = "";
      while (!/^sim4end/) { $aln .= $_; $_ = <>; }
      $aln .= $_;
      print &parse_aln($aln);
   }
}


sub parse_aln {# aln
   my $aln = shift;
   my $alnout = "";

   my $failCheck;

   my @lines = split '\n', $aln;

   while (@lines) {

      my $l = shift @lines;
      $l =~ /^sim4begin/ or die "Something wrong.$_\n";
      $alnout = "sim4begin\n";

      # exits when it finds the header line or at the end of the list
      $l = shift @lines;
      while (defined($l) && !($l=~/^\d+\[\S+\] \d+\[\S+\] <\S+>/)) { $l = shift @lines; }
      if ($l) { $alnout .= "$l\n"; }

      # exits when it finds the edef line or at the end of the list
      $l = shift @lines;
      while (defined($l) && !($l=~/^edef/)) { $l = shift @lines; }
      if ($l) { $alnout .= "$l\n"; }

      # exits when it finds the ddef line or at the end of the list
      $l = shift @lines;
      while (defined($l) && !($l=~/^ddef/)) { $l = shift @lines; }
      ((($l=~/m1fwd/) || ($l=~/m1rev/)) && ($l=~/ \(\S\S\S.?$/)) or $failCheck = 1;
      if ($l) {$alnout .= "$l\n"; }

      # exits when it finds an alignment coords line or at the end of the list
      # skips aligbment )text) lines, not useful anyway
      $l = shift @lines;
      while (defined($l) && !($l=~/^sim4end/)) {
        if ($l=~/^\d+\-\d+ \(\d+\-\d+\)/) { $alnout .= "$l\n"; }
        $l = shift @lines;
      }
      
      # skip this altogether, not useful for stats
      if ($l=~/^sim4end/ && !defined($failCheck)) {
         $alnout .= "sim4end\n";
      } else {
         $alnout = "";
      }

      $l = shift @lines;
   }
 
   return $alnout;
}

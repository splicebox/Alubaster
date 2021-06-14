#!/usr/bin/perl
use strict;

my $Usage = "Usage: $0 < em > BED\n";
my $tmp = shift;
if ($tmp eq "-h") { print $Usage; exit(0); }

#sim4begin
#3[100-0-0] 0[0-5059] <94-0-94-complement-reverse>
#edef=SRR1284895.171640632/2:- HWI-ST1022:122:D1689ACXX:4:2308:8456:72068 length=100
#ddef=ENST00000623180.1;chr17:65829-65887 m1fwd:(m1+5000):70887 (AG)
#1-25 (36-59) <24-0-96> <-
#26-100 (2120-2194) <70-0-93>
#Ntttcctgcaactggatggtcccac
#-tttcctgcaactggatggtcccac
#gtttTgctgtcacccctctGgaatccacgtataTgaaaattccaaatAttagttgggcatGgtggcaagcacctg
#gtttCgctgtcacccctctAgaatccacgtataCgaaaattccaaatGttagttgggcatAgtggcaagcacctg
#sim4end

while (<>) {
  chomp;

  /^sim4begin/ or die "Error: Expected sim4begin, got $_.\n";

  $_ = <>; chomp;
  #0[75-0-0] 0[0-1496] <51-0-100-complement-unknown>
  /^(\d+)\[(\d+)\-\d+\-\d+] (\d+)\[(\d+)\-(\d+)\] <(\d+)\-\d+\-(\d+)\-(\S+)\-(\S+)>$/ or die "Error: Expected sim4db header, got $_.\n";
  my ($cdnaid,$cdnalen,$gaid,$gafrom,$gato,$nmatches,$pctid,$matchori,$strand) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
  $gafrom=0; #no need for the new version of sim4db
  if ($matchori eq "forward") {
     $matchori = "+";
  }  elsif ($matchori eq "complement") {
     $matchori = "-";
  } else {
     die "Error: Unrecognized matchori $matchori.\n";
  }

  $_ = <>; chomp;
  my $cdnaname;
  /^edef=/ or die "Error: Expected edef, got $_.\n";
  if (/edef=(\S+)$/) { $cdnaname = $1; } elsif (/edef=(\S+)\s/) { $cdnaname = $1; }

  $_ = <>; chomp;
  my $ganame;
  /^ddef=/ or die "Error: Expected ddef, got $_.\n";
  if (/ddef=(\S+)$/) { $ganame = $1; } elsif (/ddef=(\S+)\s/) { $ganame = $1; }


  $_ = <>; chomp;
  my @cdnaBeg = (); my @cdnaEnd = (); my @gaBeg = (); my @gaEnd = (); my @Xpctid = (); my @Xnmatches = ();
  my $i = 0;
  while (/^\d/) {
     /^(\d+)\-(\d+) \((\d+)\-(\d+)\) <(\d+)\-\d+\-(\d+)>/ or die "Error: Expected match coordinate line, got $_.\n";
     $cdnaBeg[$i] = $1; $cdnaEnd[$i] = $2;
     $gaBeg[$i] = $gafrom+$3; $gaEnd[$i] = $gafrom+$4;
     $Xnmatches[$i] = $5; $Xpctid[$i] = $6; 
     $i++;

     $_ = <>; chomp;
  }
  while (!/sim4end/) {
     $_ = <>;
  }
  
  my $nx = scalar(@Xpctid);
  my ($From,$To) = ($gaBeg[0],$gaEnd[$nx-1]);
  $From--;   # First base in chromosome in BED format is 0

  print "$ganame\t$From\t$To\t$cdnaname\t$pctid\t$matchori\t$From\t$To\t255,0,0\t$nx\t";
  print ($gaEnd[0]-$gaBeg[0]+1);
  for ($i=1; $i<$nx; $i++) {
     print "," . ($gaEnd[$i]-$gaBeg[$i]+1);
  }
  print "\t0";
  for ($i=1; $i<$nx; $i++) {
     print "," . ($gaBeg[$i]-$gaBeg[0]);
  }
  print "\n";


}
exit(0);

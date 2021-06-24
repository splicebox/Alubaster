#!/usr/bin/perl
use strict;

my $Usage = "$0 txpt2gene signalnr unsplicednr regionnr\n";
($#ARGV==3) or die $Usage;
my $ANNOT = shift;
my $SIGNALfile = shift;
my $UNSPLICEDfile = shift;
my $REGIONfile = shift;

my %isREGION;
my %isUNSPLICED;

my %Txpt2Gene;
open(T, "<$ANNOT") or die "txpt2gene.\n";
while (<T>) {
  # "ENSG00000000003.13"; "ENST00000373020.7"; "TSPAN6"; chr1  1623413265 1623416666
  /^\"(\S+)\"; \"(\S+)\"; \"(\S+)\";/ or die "died. $_";
  $Txpt2Gene{"$1:$2"} = $3;
}
close(T);

open(R, "<$REGIONfile") or die "REGION: $REGIONfile\n";
while (<R>) {
   chomp;

   # ERR188081.25054577 ENSG00000000457.12:ENST00000367771.9;chr1:169823606-169827746 No .
   next if / No/;
   /^(\S+ \S+) Yes/ or die "died. $_";
   $isREGION{$1} = 1;       
}
close(R);

open(U, "<$UNSPLICEDfile") or die "UNSPLICED: $UNSPLICEDfile\n";
while (<U>) {
   chomp;

   # ERR188081.25054577 ENSG00000000457.12:ENST00000367771.9;chr1:169823606-169827746 No
   next if / No/;
   /^(\S+ \S+) Yes/ or die "died. $_";     
   $isUNSPLICED{$1} = 1;
}
close(U);

open(S, "<$SIGNALfile") or die "SIGNAL: $SIGNALfile\n";
while (<S>) {
   chomp;
   # ERR188081.10000020 ENSG00000000457.12:ENST00000367771.9;chr1:169823606-169827746 No
   next if /No/;
   /^(\S+ \S+) / or die "died. $_";
   my $key = $1;
   /^\S+ (\S+);/ or die "died. $_"; 
   next if (defined($isREGION{$key}) || defined($isUNSPLICED{$key}));
   defined($Txpt2Gene{$1}) or die "undef txpt2gene for $1.\n";
   print $_, " ", $Txpt2Gene{$1}, "\n";
}
close(S);

#!/usr/bin/perl
use strict;

#Usage: fastq2fasta.pl [-i idx] < fastqfile

my $idx;

if ($#ARGV==1) {
   my $tmparg = shift;
   ($tmparg eq "-i") or die "Usage: fastq2fasta.pl [-i idx] < fastqfile\n";
   $tmparg = shift;
   ($tmparg==1 || $tmparg==2) or die "Usage: fastq2fasta.pl [-i [1,2]] < fastqfile\n";
   $idx = $tmparg;
} elsif ($#ARGV==-1) {
   ;
} else {
   die "Usage: fastq2fasta.pl [-i idx] < fastqfile\n";
}

while (<>) {
   /^@(.+)$/ or die "died. $_";
   print (defined($idx) ? ">$1/$idx\n" : ">$1\n");
   $_ = <>; print $_;
   $_ = <>; $_ = <>;
}

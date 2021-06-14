#!/usr/bin/perl
use strict;

while (<>) {
   chomp;
   #(2) NT5C2 ENSG00000076685.17 chr10 103174913 103181267 [2 15 0 0] [,ERR127302.20865209/2,ERR127302.23911463/2,ERR127302.458172/2,ERR127302.11171813/2,ERR127302.27249525/2,ERR127302.5490442/2] [103174857-103174913,103181267-103181327]
    /^\S+ (\S+) \S+ (\S+) \S+ \S+ \[\d+ \d+ \d+ \d+\] \[\S+\] \[(\S+)\]$/ or die "died. $_";
    my ($gene,$chrom,$rest) = ($1,$2,$3);

    my @elems = split ',', $rest;
    foreach my $x (@elems) {
       $x =~ /^(\d+)\-(\d+)$/ or die "died. $_";
       print "$chrom\t$1\t$2\t$gene\n";
    }
}

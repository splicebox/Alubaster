#!/usr/bin/perl
use strict;

($#ARGV==0) or die "Usage: $0 setaluovlp_prefix < selection\n";

my $SETPATH = shift;

my %Seen;
#$SET/$SET.selection_Alu_overlaps.bed
open(F, "<$SETPATH" .".selection_Alu_overlaps.bed") or die "could not open $SETPATH.selection_Alu_overlaps.bed";
while (<F>) {
  # chr10   102901245       102901319       AS3MT   chr10   102901177       102901502       h38_mk_AluY     0       -       74 
  /^(\S+)\t(\d+)\t(\d+)\t/ or die "died. $_"; 
  $Seen{"$1:$2-$3"} = 1;
}
close(F);


# read in a selection file
while (<>) {
   chomp;
   #(2) NT5C2 ENSG00000076685.17 chr10 103174913 103181267 [2 15 0 0] [,ERR127302.20865209/2,ERR127302.23911463/2,ERR127302.458172/2,ERR127302.11171813/2,ERR127302.27249525/2,ERR127302.5490442/2] [103174857-103174913,103181267-103181327]
    /^(\S+ \S+ \S+) (\S+) (\S+ \S+ \[\d+ \d+ \d+ \d+\] \[\S+\]) \[(\S+)\]$/ or die "died. $_";
    my ($left,$chrom,$middle,$rest) = ($1,$2,$3,$4);
    print "$left $chrom $middle ";

    my @elems = split ',', $rest;
    print "[";
    while (@elems) {
       my $x = shift @elems;
       $x =~ /^(\d+)\-(\d+)$/ or die "died. $_";
       if (defined($Seen{"$chrom:$1-$2"})) { print "R $x"; } else { print $x; }
       !scalar(@elems) or print ",";
    }
    print "]\n";
}

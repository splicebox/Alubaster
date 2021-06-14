#!/usr/bin/perl
use strict;

($#ARGV==0) or die "Usage: $0 txpt2gene\n";
my $txpt2geneFile = shift;

my %geneId2geneName;
open(F, "<$txpt2geneFile") or die "could not open $txpt2geneFile";  #gencode.v22.Txpt2Gene
while (<F>) {
  # "ENSG00000000003.13"; "ENST00000373020.7"; "TSPAN6";
  /^\"(\S+)\"; \S+ \"(\S+)\";$/ or die "died (1) add_gene_name at $_";
  $geneId2geneName{$1} = $2; 
}
close(F);

# read in the matches file and add the gene information
while (<>) {
   if (/^::/) { print $_; }
   else {
     #>ERR188058.23232999/1:- >ENSG00000090686.14:ENST00000526044.4;chr1:21729704-21729829 75 588 m1fwd - 9 75 91 176 67 85 9-41,42-75 ARC 30-33-0 33-30-0 chr1:21729794-21729826;Alu:21-50;50-47
     chomp;
     /^(\S+) (\S+):(\S+;\S+ .*)$/ or die "died (2) add_gene_name at $_";
     print "$1 " . $geneId2geneName{$2} . ":$2:$3\n";
   }
}


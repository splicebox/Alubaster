#!/usr/bin/perl
use strict;

($#ARGV==0) or die "Usage: $0 annotfile < genes.bed > groups\n";
my $annotFile = shift;

my %Txpt2Gene;

open(F, "<$annotFile") or die "died. $_";
while (<F>) {
   next if /^#/;
   next if !/transcript_id/;

   m/transcript_id \"(\S+)\";/ or die "died (1). $_";
   my $tid = $1;
   m/gene_id \"(\S+)\";/ or die "died (2). $_";
   my $geneid = $1;
   $Txpt2Gene{$tid} = $geneid;
}
close(F);

# read in a genes.bed file and group readids by anchor exons; "NA12814.nonconcordant.genes.last-exon.bed.\n";

my %Seen;

while (<>) {
  # chr10   100150269       100150370       HWI-D00372:257:C5KMHANXX:2:2209:12402:85389/2   50      -       100150269       100150370       0,0,0   1       101,    0,      chr10   100150093       100152352       ENST00000421367.5       0       -       101
  # chr17	10634024	10634124	SRR1284895.65279512/2	50	+	10634024	10634124	0,0,0	1	100,	0,	chr17	10634016	10634182	ENST00000583535.4	0	-	100

  /^\S+\t\d+\t\d+\t(\S+)\t\S+\t(\S)\t\d+\t\d+\t\S+\t\S+\t\S+\t\S+\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t/ or die "died (1). $_";
  my ($readid,$ori,$chrom,$from,$to,$tid) = ($1,$2,$3,$4,$5,$6);
  $from++;    # Adjust for -1 in bed first coordinate
  my $key = $Txpt2Gene{$tid} . ":$tid;$chrom:$from-$to";
  #"ENSG00001234.5:ENST00000466430.4;chr1:89295:91629 ERR188081.10684278/1
  if (defined($Seen{$key})) {
     $Seen{$key} .= "\t$readid$ori";
  } else {
     $Seen{$key} = "$readid$ori";
  }
}

#HWI-D00372:257:C5KMHANXX:2:2209:12402:85389/2
foreach my $k (keys %Seen) {
   $k =~ /^\S+;(\S+):(\d+)\-(\d+)$/ or die "died. (2). $k\n";
   my ($chrom,$f,$t) = ($1,$2,$3);
   my @elems = split '\t', $Seen{$k};
   print "$chrom\t$f\t$t\t";
   print "$k\t", scalar(@elems), "\t", join(',',@elems), "\n";
}


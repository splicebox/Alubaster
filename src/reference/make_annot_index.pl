#!/usr/bin/perl
use strict;

## Read in a GENCODE gene annotations file and extract gene name, transcript name and transcript ID
#
#Example input line:
#chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; level 2; tag "basic"; transcript_support_level "1"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
 
while (<>) {
  next if /^#/;
  next if !/\ttranscript\t/;
  /^(\S+)\t\S+\ttranscript\t(\d+)\t(\d+)\t\S+\t\S+\t\S+\t/ or die "died. $_";
  my ($chrom,$from,$to) = ($1,$2,$3);
  /gene_id (\"\S+\";)/ or die "no geneid? $_";
  my $geneid = $1;
  /transcript_id (\"\S+\";)/ or die "no transcriptid? $_";
  my $tid = $1;
  my $genename = (/gene_name (\"\S+\";)/) ? $1 : "-";

  #"ENSG00000000003.13"; "ENST00000373020.7"; "TSPAN6";
  print "$geneid $tid $genename $chrom $from $to\n";
}

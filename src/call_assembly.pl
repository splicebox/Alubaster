#!/usr/bin/perl
use strict;

#export OMP_THREAD_LIMIT=16
#export OMP_NUM_THREADS=10
# OUTPUT FILES: 
#   prefix.contigs.fa
#   prefix.scr


my $Usage = "Usage: $0 selectionfile reads mates -gi geneindex -ge genomeindex -gen genome -p prefix [k] \n";

($#ARGV>=11) or die "5", $Usage;;
my $selectionFile = shift;
my $READSFA = shift;
my $MATESFA = shift;
my $GINDEX = shift;
($GINDEX eq "-gi") or die "1", $Usage; $GINDEX = shift;
my $GENDEX = shift;
($GENDEX eq "-ge") or die "2", $Usage; $GENDEX = shift;
my $GENOME = shift;
($GENOME eq "-gen") or die "3", $Usage; $GENOME = shift;
my $prefix = shift;
($prefix eq "-p") or die "4", $Usage; $prefix = shift;
my $K = shift;
(defined($K)) or $K = 23;

my $SIM4DB = $ENV{SIM4DB};
(defined($SIM4DB)) or die "can't find sim4db";
my $SCRIPTS = $ENV{SCRIPTS};
(defined($SCRIPTS)) or die "can't find scripts directory";
my $OASES = $ENV{OASES};
(defined($OASES)) or die "can't find OASES !";
my $VELVETH = $ENV{VELVETH};
(defined($VELVETH)) or die "can't find VELVETH !";
my $VELVETG = $ENV{VELVETG};
(defined($VELVETG)) or die "can't find VELVETG !";


my %GeneFrom; my %GeneTo; my %GeneChrom;
#### Read and parse the gene coordinate index
#
open(F, "<$GINDEX") or die "Could not open gene indexfile $GINDEX.\n";
while (<F>)  {
  #"ENSG00000243485.3"; "ENST00000469289.1"; "RP11-34P13.3"; chr1 30267 31109
  /^\"(\S+)\"; \"(\S+)\"; \"(\S+)\"; (\S+) (\d+) (\d+)$/ or die "died. $_";
  my ($geneid,$tid,$genename,$chrom,$from,$to) = ($1,$2,$3,$4,$5,$6);
  if (!defined($GeneFrom{$geneid}) || ($GeneFrom{$geneid}>$from)) {
     if (!defined($GeneTo{$geneid}) || ($GeneTo{$geneid}<$to)) {
        $GeneTo{$geneid} = $to;
        $GeneTo{$genename} = $to;
     }
     $GeneFrom{$geneid} = $from;
     $GeneChrom{$geneid} = $chrom;
     $GeneFrom{$genename} = $from;
     $GeneChrom{$genename} = $chrom;
  }
}
close(F);

my %GAid;
my $count = 0;
### Read and convert the chromsome id
#
open(F, "<$GENDEX") or die "Could not open gene indexfile $GENDEX.\n";
while (<F>)  {
   chomp;
   if (/^>(\S+)$/) { $GAid{$1} = $count; }
   elsif (/^>(\S+)\s/) { $GAid{$1} = $count; }
   $count++;
} 
close(F);

### Collect all contigs and all sim4db alignments in BIG genome-wide files; even if we might need to split them later
# 
my $nContigs = 0;

open(C, ">$prefix.contigs.fa") or die "could not open big contigs file for writing.\n";
open(R, ">$prefix.scr") or die "could not open big scrips file for writing.\n";

### Parse the selection file and process each line (to be changed to process one gene at a time)
#
#(2) ZNF175 ENSG00000105497.6 chr19 51591925 51592058 [5 10 0 0] [,ERR204931.3192493/1,ERR204931.21777482/2,ERR204931.13385252/1,ERR204931.16811222/2,ERR204931.14571901/1,ERR204931.9473943/2,ERR204931.21128203/1,ERR204931.9449828/2,ERR204931.3839247/1,ERR204931.14163241/1;,ERR204931.13690777/2,ERR204931.16245050/1,ERR204931.20786637/2]
#
open(F, "<$selectionFile") or die "Could not open selection file $selectionFile.\n";
while (<F>) {
  /^\S+ (\S+) (\S+) (\S+) (\S+) (\S+) \S+ \S+ \S+ \S+ \[(\S+)\]$/ or die "died. $_";
  my ($genename,$geneid,$chrom,$iFrom,$iTo,$reads) = ($1,$2,$3,$4,$5,$6);
  $reads =~ s/;,/,/; $reads =~ s/^,//; $reads =~ tr/,/\n/;
  my @read_list = split '\n', $reads;

  if ($iFrom =~ /^\d+$/) { $iFrom -= 200; } else { $iFrom = $GeneFrom{$geneid}-200; }
  if ($iTo =~ /^\d+$/) { $iTo += 200; } else { $iTo = $GeneTo{$geneid}+200; }

  my $gaid;
  if (defined($GAid{$chrom})) {
     $gaid = $GAid{$chrom};
  } else {
     die "Error: Unknown gaid for $chrom\n";
  }
 
  my $st = 0;

`cat <<"EOF" > $prefix.tmpreads.ids
$reads
EOF
`;
  if (-s "$prefix.tmpreads.ids") {
    `${SIM4DB}/leaff -F $READSFA -q $prefix.tmpreads.ids > $prefix.tmpreads_1.fa`;
  } else {
     $st = 1;
     print STDERR "WARNING: Read_1 ID collection failed for $genename. Skipping assembly.\n";
  } 

  # modify the read list to show the mate ids
  foreach my $x (@read_list) {
       $x =~ /^(\S+)\/(\d)$/ or die "died at x: $x\n";
       $x = "$1/" . (3-$2);
  }
  $reads = join "\n", @read_list;
`cat <<"EOF" > $prefix.tmpreads.ids
$reads
EOF
`;
  if (-s "$prefix.tmpreads.ids") {
  # $st = system("echo \"$reads\" > $prefix.tmpreads.ids");
  # if ($st) { die "st is $st\n"; }
    `${SIM4DB}/leaff -F $MATESFA -q $prefix.tmpreads.ids > $prefix.tmpreads_2.fa`;
  } else {
     $st = 1;
     print STDERR "WARNING: Read_2 ID collection failed for $genename. Skipping assembly.\n";
  }

  if (!$st) {
    `$VELVETH ${prefix}_oases_$K $K -fmtAuto -separate -shortPaired $prefix.tmpreads_1.fa $prefix.tmpreads_2.fa`;
    `$VELVETG ${prefix}_oases_$K -read_trkg yes`;
#print STDERR "Run OASES now\n";
    `$OASES ${prefix}_oases_$K -unused_reads yes -alignments yes`;
#print STDERR "Done OASES \n";
    open(G, "< ${prefix}_oases_$K/contigs.fa") or die "could not open contigs file ${prefix}_oases_$K/contigs.fa.\n";
    while (<G>) {
        if (/^>([\S \t]+)$/) {
           print C ">$genename:$1\n";
           print R "-f -e $nContigs -D $gaid $iFrom $iTo\n";
           print R "-r -e $nContigs -D $gaid $iFrom $iTo\n";
           $nContigs++;
        } else {
           print C $_;
        }
    }
    close(G);
  }
 
# `rm -r ${prefix}_oases_$K`;
}
close(F);
close(C);
close(R);

### Now call sim4db on the big contig fasta file
#print STDERR "simdb4 -- output = $prefix.contigs.sim4db\n";

print STDERR "comanda: ${SIM4DB}/sim4db -cdna $prefix.contigs.fa -gen $GENOME -script $prefix.scr -aligns -o $prefix.contigs.sim4db \n";
`${SIM4DB}/sim4db -cdna $prefix.contigs.fa -gen $GENOME -script $prefix.scr -aligns -o $prefix.contigs.sim4db`;

#print STDERR "em2BED -- output = $prefix.contigs.bed\n";
print STDERR "comanda: ${SCRIPTS}/em2BED.pl < $prefix.contigs.sim4db > $prefix.contigs.bed \n";
`${SCRIPTS}/em2BED.pl < $prefix.contigs.sim4db > $prefix.contigs.bed`;

#print STDERR "DONE call_assembly\n";

exit(0);

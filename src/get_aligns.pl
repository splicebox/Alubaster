#!/usr/bin/perl
use strict;

# Note: At least one of -k or -d must be specified, possibly both
($#ARGV>=2) or die "Usage: $0 -k kraken_prefix out_prefix [-u] < samfile\n";
my ($krakenPrefix,$outPrefix,$uniq);
$krakenPrefix = shift;
if ($krakenPrefix eq "-k") {
   $krakenPrefix = shift;
} else {
    die "Usage: $0 -k kraken_prefix out_prefix [-u] < samfile\n";
}
$outPrefix = shift;
$uniq = shift;
!defined($uniq) or ($uniq eq "-u") or die "Usage: $0 -k kraken_prefix out_prefix [-u] < samfile\n";

my $L =  10000;

# Procedure:
# read in a SAM file sorted by read id, one block (readid=pairid:idx) at a time
# store it, and determine whether it is concordant or not, and whether the mate is repeat or not
#      attributes: per align: concordant_align && (Alu(m1) ||  Alu(m2)), add alignment to CLines
#                  per align: discordant_align && ((m1 && Alu(m2)) || (m2 && Alu(m1)) && isPrimary (optional), add alignment to NCLines
#                  per block: print CLines; if (!hasConcordant) { print $NCLines; }

my %Alu;

# read in the two (paired) lists of kraken files
if (defined($krakenPrefix)) {
   my $tmpkFile = $krakenPrefix . "_1.kraken_out";
   open(F, "<$tmpkFile") or die "could not open kraken file $tmpkFile.\n";
   while (<F>) {
      next if !/^C\s/;
      chomp;
      /^C\s(\S+)\s/ or die "died (kraken). $_";
      $Alu{$1} = 1;
   }
   close(F);

   my $tmpkFile = $krakenPrefix . "_2.kraken_out";
   open(F, "<$tmpkFile") or die "could not open kraken file $tmpkFile.\n";
   while (<F>) {
      next if !/^C\s/;

      chomp;
      /^C\s(\S+)\s/ or die "died (kraken). $_";
      if (defined($Alu{$1})) { $Alu{$1} = 3; } else { $Alu{$1} = 2; }
   }
   close(F);
}

open(C, ">$outPrefix.concordant.sam") or die "could not open concordant for writing.\n";
open(NC, ">$outPrefix.nonconcordant.sam") or die "could not open nonconcordant for writing.\n";

# read in a SAM file, sorted by read pair id
my $prev_readid;
my $cLines = "";
my $ncLines = "";
my $hasConcordant = 0;
my $hasUniq = 1;

my %Seen;
while (<>) {
   if (/^@/) { print C $_; print NC $_; next; }
   
   /^(\S+)\t(\d+)\t/ or die "died (readid). $_";
   my ($pairid,$flag) = ($1,$2);
   my $idx = (($flag & 0x40) ? 1 : (($flag & 0x80) ? 2 : 0));
   my $readid = "$pairid:$idx";

   if (!defined($prev_readid) || ($readid ne $prev_readid)) {
      if (defined($prev_readid) && ($readid lt $prev_readid)) {
         print STDERR "WARNING: File does not appear to be sorted: $readid<$prev_readid\n";
      } 
      if (defined($prev_readid)) {
         print C $cLines;
         if (!$hasConcordant && (!defined($uniq) || ($hasUniq==1))) {
            print NC $ncLines;
         }
      }
      $prev_readid = $readid;
      $hasConcordant = 0;
      %Seen = ();
      $hasUniq = 1;
      $cLines = ""; $ncLines = "";
   }
   
   # determine whether concordant or discordant for the current line

   # ERR188081.3807988  417     chr1    26530   0       75M     =       197302  170847  CCCGTTCAAGAATGGGACTGAATACACCTGATGAGTGGTTTACTTTCTGTCTGCAAACATCTACTGATCATCTGT     CCCFFFFFHHHHHJJJJJJIJJJJIGJJJJJJFGICGIBGGIIJJJJJIJHIIIIIIIJJJJJJIJJJJIJJJJJ     AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YT:Z:UU NH:i:10 CC:Z:=  CP:i:26530      HI:i:0
   /^\S+\t\d+\t(\S+)\t(\d+)\t\d+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/ or die "died. $_";
   my ($chrom,$pos,$cigar,$next_chrom,$next_pos,$dist) = ($1,$2,$3,$4,$5,$6);

   my $ori1 = ($flag & 0x10) ? 1 : 0;
   my $ori2 = ($flag & 0x20) ? 1 : 0;

   my $isPrimary = (!($flag & 0x100) && !($flag & 0x900));

   if (defined($Seen{$readid}) && ($Seen{$readid} ne "$chrom:$pos:$cigar")) {
      $Seen{$readid} .= ";$chrom:$pos:$cigar";
      $hasUniq = 0;
   } elsif (!defined($Seen{$readid})) { # otherwise stays as is, equal to the current chrom:pos:cigar
      $Seen{$readid} = "$chrom:$pos:$cigar";
   }

   if (defined($Alu{$pairid}) &&                         # condition: one of the reads in the pair has Alu match
       (($chrom ne "*") && ($next_chrom eq "=")) &&      # condition: same chromosome
       ((!$ori1 && $ori2) || (!$ori2 && $ori1)) &&       # condition: opposite orientations
#      ((!$ori1 && $dist>=0) || ($ori1 && $dist<=0)) &&  # condition: correct order (see ppp.pl and ggg: m1fwd, m2fwd (!ori1): most dist>0; m1rev, m2rev (ori1): most dist are <0))
       ((!$ori1 && ($next_pos>=$pos-10)) ||              # condition: correct order, m1(m2)fwd
        ($ori1 && ($next_pos<=$pos+10))) &&              # condition: correct order, m1(m2)rev
       (abs($dist) <= $L)) {                             # insert size
       $cLines .= $_;
       $hasConcordant = 1;
   } else {
       my $isAlu = ($Alu{$pairid}==3) || ($idx==1 && $Alu{$pairid}==2) || ($idx==2 && $Alu{$pairid}==1);
       if ($isAlu && $isPrimary) { $ncLines .= $_; }   # change to isALu && isPrimary
   }
}
if (defined($prev_readid)) {
   print C $cLines;
   if (!$hasConcordant && (!defined($uniq) || ($hasUniq==1))) {
      print NC $ncLines;
   }
}
close(C); close(NC);

exit(0);

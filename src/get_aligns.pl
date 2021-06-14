#!/usr/bin/perl
use strict;

my $MINALULEN = 20;

# Note: At least one of -k or -d must be specified, possibly both
($#ARGV>=2) or die "Usage(1): $0 [-k kraken_prefix] [-d Alusim4dbprefix] [-o out_prefix] [-u] < samfile\n";
my ($krakenPrefix,$alusim4dbPrefix,$outPrefix,$uniq);
my $argtmp;
while (@ARGV) {
   $argtmp = shift;
   if ($argtmp eq "-k") {
       $krakenPrefix = shift;
   } elsif ($argtmp eq "-d") {
       $alusim4dbPrefix = shift;
   } elsif ($argtmp eq "-o") {
       $outPrefix = shift;
   } elsif ($argtmp eq "-u") {
       $uniq = 1; 
   } else {
       die "Usage(2): $0 [-k kraken_prefix] [-d alusim4dbprefix] [-o out_prefix] [-u] < samfile\n";
   }
}

#!defined($uniq) or die "Usage(3): $0 [-k kraken_prefix] [-d alusim4dbprefix] [-o out_prefix] [-u] < samfile\n";
defined($krakenPrefix) or defined($alusim4dbPrefix) or die "Did not specify a classification file. \nUsage(4): $0 [-k kraken_prefix] [-d alusim4dbprefix] [-o out_prefix] [-u] < samfile\n";

my $L =  10000;

# Procedure:
# read in a SAM file sorted by read id, one block (readid=pairid:idx) at a time
# store it, and determine whether it is concordant or not, and whether the mate is repeat or not
#      attributes: per align: concordant_align && (Alu(m1) ||  Alu(m2)), add alignment to CLines
#                  per align: discordant_align && ((m1 && Alu(m2)) || (m2 && Alu(m1)) && isPrimary (optional), add alignment to NCLines
#                  per block: print CLines; if (!hasConcordant) { print $NCLines; }

### Allow Alu information from two sources (either or both); sim4db ALU alignments and Kraken classification
#
#
#### Kraken information:
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

my %Alu4;

# read in the two (paired) lists of sim4db ALU alignment files; Example: SRR1284895sim_1.Alu.sim4db.stats
# SRR1284895.5018 ALU 30 100 + 30.00 - OLD; see below
if (defined($alusim4dbPrefix)) {
   my $tmp4File = $alusim4dbPrefix . "_1.Alu.sim4db.coords.stats";
   open(F, "<$tmp4File") or die "could not open sim4db Alu alignments file $tmp4File.\n";
   while (<F>) {
      chomp;
      #tid ganame beg end gabeg gaend cov len relori pctcov pctid
      #ERR127302.129 ALU 17 72 1 56 56 72 + 77.78 89
      /^(\S+)\sALUY?\s\d+\s\d+\s\d+\s\d+\s(\d+)\s\d+\s\S\s\S+\s\d+$/ or die "died (alusim4db). $_";
      if ($2>=$MINALULEN) { $Alu4{$1} = 1; }
   }
   close(F);

   $tmp4File = $alusim4dbPrefix . "_2.Alu.sim4db.coords.stats";
   open(F, "<$tmp4File") or die "could not open sim4db Alu alignments file $tmp4File.\n";
   while (<F>) {
      chomp;
      /^(\S+)\sALUY?\s\d+\s\d+\s\d+\s\d+\s(\d+)\s\d+\s\S\s\S+\s\d+$/ or die "died (alusim4db). $_";
      if ($2>=$MINALULEN) {
         if (defined($Alu{$1})) { $Alu{$1} = 3; } else { $Alu{$1} = 2; }
      }
   }
   close(F);
}

### Now reconcile all information about ALU match in the %Alu data structure; when both criteria used, then both must show match
if (defined($krakenPrefix) && defined($alusim4dbPrefix)) {
   foreach my $k (keys %Alu) {
      next if ($Alu{$k} == $Alu4{$k});
      if (!defined($Alu4{$k})) { undef $Alu{$k}; next; }   # strip the label if not confirmed by sim4db
      # strip the label if incompatible
      if ($Alu{$k}+$Alu4{$k} == 3) { undef $Alu{$k}; next; }
      if ($Alu{$k}>$Alu4{$k}) { $Alu{$k} = $Alu4{$k}; }  # Basically, here Alu{$K} is 3, and restrict to what sim4db says
      # the only remaining cases are 1,3 and 1,2 which stay as they are; no correction
   }
} elsif (defined($alusim4dbPrefix)) {
   foreach my $k (keys %Alu4) {
      $Alu{$k} == $Alu4{$k};
   }
} # otherwise leave Alu as is, as the only criterion (kraken-based)

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

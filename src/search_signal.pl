#!/usr/bin/perl
use strict;

($#ARGV==4) or die "Usage: $0 prefix groupfile matesfastafile annotfile genomefile\n";
my $SET     = shift; # ./chr17.test
my $GROUPS  = shift; # "/scratch0/igm3/florea/repeat/Spliced/NA12814.nonconcordant.R.genes.last-exon.nr.sorted.groups";
my $MATES   = shift; # "/scratch0/igm3/florea/repeat/Spliced/NA12814.nonconcordant.R.mates.fa";
my $ANNOT   = shift; # "/home/florea/gencode.v22.annotation.gtf";
my $GENOME  = shift; # "/scratch0/igm3/florea/repeat/hg38c.fa";
my $HDRS    = $GENOME . ".hdrs"; # "/scratch0/igm3/florea/repeat/hg38c.fa.hdrs";

# debug hg38c.Sim
#$SET =  "ZFX.left";
#$GROUPS = "/ccb/salz4-1/florea/repeat/Pipeline/SRR1284895sim/ZFX.left.groups";
#$MATES = "/ccb/salz4-1/florea/repeat/Pipeline/SRR1284895sim/SRR1284895sim.nonconcordant.mates.fa";
#$ANNOT = "/ccb/salz4-1/florea/repeat/Data/hg38c.Sim/hg38sim/transcriptome/gencode.v22.annotation.sim.gtf";
#$GENOME = "/ccb/salz4-1/florea/repeat/Data/hg38c.Sim/hg38sim/hg38sim.fa";
#$HDRS = $GENOME . ".hdrs"; 



# $SET = "test.groups1";

# DEBUG ONLY
#$GROUPS = "./groups"; # "./chr17.test/chr17.test.nonconcordant.genes.nr.sorted.groups";
#$MATES  = "./chr17.test/chr17.test.nonconcordant.mates.fa";
#$ANNOT  = "chr17.test/anno"; # "/home/florea/gencode.v22.annotation.gtf";
#$GENOME = "/scratch0/igm3/florea/repeat/hg38c.fa";
#$HDRS    = $GENOME . ".hdrs";

my $ALU     = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGAGGATTGCTTGAGC" .
              "CCAGGAGTTCGAGACCAGCCTGGGCAACATAGCGAGACCCCGTCTCTACAAAAAATACAAAAATTAGCCG" .
              "GGCGTGGTGGCGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGGATCGCTTGAGCCCAGG" .
              "AGTTCGAGGCTGCAGTGAGCTATGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACCCTGTC" .
              "TCA";
my $revALU  = "TGAGACAGGGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCGCGATCATAGCTCACTGCAGCCTCGA" .
              "ACTCCTGGGCTCAAGCGATCCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCGCGCCACCAC" .
              "GCCCGGCTAATTTTTGTATTTTTTGTAGAGACGGGGTCTCGCTATGTTGCCCAGGCTGGTCTCGAACTCC" .
              "TGGGCTCAAGCAATCCTCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCG" .
              "GCC";
my $ALUY    = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTC" .
              "AGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGG" .
              "GCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGA" .
              "GGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCT" .
              "CA";
my $revALUY = "TGAGACGGAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCGCGATCTCGGCTCACTGCAAGCTCCG" .
              "CCTCCCGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCAC" .
              "GCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTC" .
              "CTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGG" .
              "CC";

my $ScaffSeq;
my $revScaffSeq;
my $crtGAid;

my %GaId;
my $count = 0;
open(F, "<$HDRS") or die "$HDRS";
while (<F>) {
  chomp;
  if (/^>(\S+)$/) { $GaId{$1} = $count; }
  elsif (/^>(\S+)\s/) { $GaId{$1} = $count; }
  else { die "died. $_"; }
  $count++;
}
close(F);


#######   CALCULATE NEXT_EXON AND PREV_EXON, FOR 'REGION'

my %TxptOri;          # Now indexed by gene_id:txpt_id
my %Seen;             # pairs of consecutive exons that have already been recorded for Next, Prev
my %ExonFarNextPos;   # farthest neighbor to the right
my %ExonFarPrevPos;   # farthest neighbor to the left
my %ExonNextPos;      # list of neighboring next exons (on the right), separated by ';'
my %ExonPrevPos;      # list of neighboring previous exons (on the left), separated by ';'

my @annot = ();
my ($prev_gtid,$prev_from,$isFlipped);
open(F, "<$ANNOT") or die "$ANNOT";
while (<F>) {
  next if (/^#/ || !/\texon\t/);
   /gene_id \"(\S+)\";/ or die "died. $_";
   my $geneid = $1;
   /transcript_id \"(\S+)\";/ or die "died. $_";
   my $tid = $1;
   /^\S+\t\S+\t\S+\t(\d+)\t\d+\t\S+\t(\S)+/ or die "died. $_";
   my ($from,$ori) = ($1,$2);
   $TxptOri{"$geneid:$tid"} = ($ori eq "+") ? 1 : -1; 
   if (defined($prev_gtid) && ($prev_gtid ne "$geneid:$tid")) {
      if (defined($isFlipped)) {
         @annot = reverse(@annot);
      }
      &process_Next_Prev(scalar(@annot),@annot);
      undef $isFlipped;
      @annot = ();
   } elsif (defined($prev_gtid) && ($prev_from>$from)) {
      $isFlipped = 1;
   }
   push @annot, $_;
   $prev_gtid = "$geneid:$tid";
   $prev_from = $from;
}
if (defined($prev_gtid)) {
   if (defined($isFlipped)) {
         @annot = reverse(@annot);
   }
   &process_Next_Prev(scalar(@annot),@annot);
}
close(F);

sub process_Next_Prev # nannot annotlist
{
   my $n = shift;
   my @annot = @_;

   my ($chrom,$from,$to,$prev_from,$prev_to);

   # first line:
   my $l = shift @annot;
   $l =~ /^(\S+)\t\S+\t\S+\t(\d+)\t(\d+)\t\S+\t\S+/ or die "died. $_";
   my ($chrom,$prev_from,$prev_to) = ($1,$2,$3);
 
   foreach my $l (@annot) {
      $l =~ /^(\S+)\t\S+\t\S+\t(\d+)\t(\d+)\t\S+\t\S+/ or die "died. $_";
      ($chrom,$from,$to) = ($1,$2,$3);    
      my $x = "$chrom:$prev_from-$prev_to";
      if (!defined($ExonFarNextPos{$x}) || ($ExonFarNextPos{$x}<$from)) { $ExonFarNextPos{$x} = $from; }
      if (!defined($ExonNextPos{$x})) {
         $ExonNextPos{$x} = "$from-$to";
      } elsif (!defined($Seen{"$prev_from-$prev_to;$from-$to"})) {
         $ExonNextPos{$x} .= ";$from-$to";
      }

      $x = "$chrom:$from-$to";
      if (!defined($ExonFarPrevPos{$x}) || ($ExonFarPrevPos{$x}>$prev_to)) { $ExonFarPrevPos{$x} = $prev_to; }
      if (!defined($ExonPrevPos{$x})) {
         $ExonPrevPos{$x} = "$prev_from-$prev_to";
      } elsif (!defined($Seen{"$prev_from-$prev_to;$from-$to"})) {
         $ExonPrevPos{$x} .= ";$prev_from-$prev_to";
      }
      $Seen{"$prev_from-$prev_to;$from-$to"} = 1;
      $prev_from = $from; $prev_to = $to;
   }
}

 
#my ($prev_tid,$prev_from,$prev_to);
#open(F, "<$ANNOT") or die "$ANNOT";
#while (<F>) {
#   next if (/^#/ || !/\texon\t/);
#   /transcript_id (\"\S+\");/ or die "died. $_";
#   my $tid = $1;
#   /^(\S+)\t\S+\t\S+\t(\d+)\t(\d+)\t\S+\t(\S)+/ or die "died. $_";
#   my ($chrom,$from,$to,$ori) = ($1,$2,$3,$4);
#   $TxptOri{$tid} = ($ori eq "+") ? 1 : -1;
#   if (defined($prev_tid) && ($prev_tid eq $tid)) {
#      my $x = "$chrom:$prev_from-$prev_to";
#      if (($TxptOri{$tid}>0) && ($from>$prev_from)) {
#         if (!defined($ExonFarNextPos{$x}) || ($ExonFarNextPos{$x}<$from)) { $ExonFarNextPos{$x} = $from; }
#         if (!defined($ExonNextPos{$x})) {
#            $ExonNextPos{$x} = "$from-$to";
#         } elsif (!defined($Seen{"$prev_from-$prev_to;$from-$to"})) {
#            $ExonNextPos{$x} .= ";$from-$to";
#         }
#         $x = "$chrom:$from-$to";
#         if (!defined($ExonFarPrevPos{$x}) || ($ExonFarPrevPos{$x}>$prev_to)) { $ExonFarPrevPos{$x} = $prev_to; }
#         if (!defined($ExonPrevPos{$x})) {
#            $ExonPrevPos{$x} = "$prev_from-$prev_to";
#         } elsif (!defined($Seen{"$prev_from-$prev_to;$from-$to"})) {
#            $ExonPrevPos{$x} .= ";$prev_from-$prev_to";
#         }
#         $Seen{"$prev_from-$prev_to;$from-$to"} = 1;
#      } else {
#         ($TxptOri{$tid}<0) or die "died (ori). $tid\n";
#         if (!defined($ExonFarPrevPos{$x}) || ($ExonFarPrevPos{$x}>$to)) { $ExonFarPrevPos{$x} = $to; }
#         if (!defined($ExonPrevPos{$x})) {
#            $ExonPrevPos{$x} = "$from-$to";
#         } elsif (!defined($Seen{"$from-$to;$prev_from-$prev_to"})) {
#            $ExonPrevPos{$x} .= ";$from-$to";
#         }
#         $x = "$chrom:$from-$to";
#         if (!defined($ExonFarNextPos{$x}) || ($ExonFarNextPos{$x}<$prev_from)) { $ExonFarNextPos{$x} = $prev_from; }
#         if (!defined($ExonNextPos{$x})) {
#            $ExonNextPos{$x} = "$prev_from-$prev_to";
#         } elsif (!defined($Seen{"$from-$to;$prev_from-$prev_to"})) {
#            $ExonNextPos{$x} .= ";$prev_from-$prev_to";
#         }
#         $Seen{"$from-$to;$prev_from-$prev_to"} = 1;
#      }
#   }
#   $prev_tid = $tid; $prev_from = $from; $prev_to = $to;
#}  
#close(F);
    
my %ReadMatchOri;

##### DETERMINE THE ORI OF THE MATCHED MATE;
# Now done from the groups file
# This will be input in the temporary fasta file for reference during the selection (for consistency with alignment)

#open(F, "<$OVERLAP") or die "died. $OVERLAP";
#while (<F>) {
#   next if /^#/;
#   # chr1	90937	91012	ERR188081.10684278/1	1	+	90937	91012	0,0,0	1	75,	0,	chr1	89295	91629	"ENST00000466430.4";	75
#   /^\S+\t\d+\t\d+\t(\S+)\t\S+\t(\S)\t\d+\t\d+\t\S+\t\d+\t\S+\t\S+\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t/ or die "died. $_";
#   my ($readid,$ori,$achr,$afrom,$ato,$atid) = ($1,$2,$3,$4,$5,$6);
#   # key is readid:anchor_exon
#   $ReadMatchOri{"$readid:$achr:$afrom-$ato"} = $ori;   # + or -
##  print "<<<$readid:$achr:$afrom-$ato : " . $ReadMatchOri{"$readid:$achr:$afrom-$ato"} . ">>>\n";
#}
#close(F);
 

##### GENERATE THE GENOMIC INTERVAL  FILES FOR SEARCHES
# read in the output from bedtools, Grouped by Anchor exon
# foreach anchor exon, extract the list of mate reads (leaff)
# 1. sim4db-search the list of reads against (anchor exon, ALU) - check orientation, order;  filter matches that cover >90% of the read
# 2. sim4db-search the list of reads against (anchor exon, +500 bp) - check orientation, order; filter matches that cover >90% at >9X% identity
# 3. sim4db-search the list of exons against (genomic region: anchor exon .. next exon)


open(F, "<$GROUPS") or die "$GROUPS";
while (<F>) {
   print STDERR $_;

   # ENSG00000197381.14:ENST00000348831.7;chr21:45184621-45184789    5       SRR1284895.127662969/1-,SRR1284895.124641841/1-,SRR1284895.12560659/2+,SRR1284895.165006817/2+,SRR1284895.159902753/2-
   /^(\S+);(\S+):(\d+)\-(\d+)\t\d+\t(\S+)$/ or die "died. $_";
   my ($txptid,$chrom,$from,$to,$reads) = ($1,$2,$3,$4,$5);

   my $tori = $TxptOri{$txptid};
#  print "!!!$txptid $tori!!!\n";
   my $gaid = $GaId{$chrom};

   if (!defined($crtGAid) || ($gaid!=$crtGAid)) {
      $ScaffSeq = `leaff -F $GENOME -s $gaid | tail -1`;
      $ScaffSeq =~ tr/acgtn/ACGTN/;
      $revScaffSeq = reverse $ScaffSeq;
      $revScaffSeq =~ tr/acgtn/ACGTN/;
      $crtGAid = $gaid;
   }

   # extract read sequences
   # insert the opposite of matched mate's ori -- this is the valid ori for an alignment

   my @readids = split ',', $reads;
   open(S, ">$SET.tmpreads.ids") or die "$SET.tmpreads.ids";
   foreach my $r (@readids) { 
      $r =~ /^(\S+).([12])([+-])$/ or die "died. $_";
      my ($readid,$idx,$ori) = ($1,$2,$3);
      $ReadMatchOri{"$readid/$idx:$chrom:$from-$to"} = $ori;   # FIXED
      $r = "$1/" . (3-$2);
      print S $r . "\n";
   }
   close(S);
   open(S, ">$SET.tmpreads.fa") or die "$SET.tmpreads.fa";
   open(M, "leaff -Fdc $MATES -q $SET.tmpreads.ids |") or die "could not extract the reads.";
   while (<M>) {
      chomp;
      /^>(\S+).([12]) ([\S \t]+)$/ or die "died. $_";
      my ($readid,$idx,$rest) = ($1,$2,$3);
      my $otheridx = 3-$idx;
      defined($ReadMatchOri{"$readid/$otheridx:$chrom:$from-$to"}) or die "undef ori for $readid/$otheridx:$chrom:$from-$to";
      my $sign = ($ReadMatchOri{"$readid/$otheridx:$chrom:$from-$to"} eq "-") ? "+" : "-";
      print S ">$readid/$idx:$sign $rest\n";
      $_ = <M>; print S $_;
   }
   close(M);
   close(S);

   ######   MAIN - GENERATE THE REGIONS 
   ###### 


   ######  'SIGNAL'
   # generate signal sequences
   # order and ori:
   # m1 - anchor read; m2 - its Alu-labeled mate
   # m1 fwd = sequence match as is; m1 rev (flag & 0x10 == 1) = sequence matched in reverse complement
   # if gene on + strand: if m1 fwd, then anchor . reverse(ALU); else reverse(ALU) . anchor
   # if gene on - strand: if m1 fwd, then anchor . ALU; else ALU . anchor
   open(S, ">$SET.signal.fa") or die "$SET.signal.fa";

   my $pos = $from-1;
   my ($I_signal_seq_m1fwd,$I_signalY_seq_m1fwd,$I_signal_seq_m1rev,$I_signalY_seq_m1rev);
   my ($X_signal_seq_m1fwd,$X_signalY_seq_m1fwd,$X_signal_seq_m1rev,$X_signalY_seq_m1rev,$XA_signal_seq_m1fwd,$XA_signalY_seq_m1fwd,$XA_signal_seq_m1rev,$XA_signalY_seq_m1rev);

   # truly, should only be included when m1fwd
   my @elems = split ';', $ExonNextPos{"$chrom:$from-$to"};

   if (scalar(@elems)) {   # There is at least one 'next' exon

      foreach my $x (@elems) {
         #####  Intronic insertion: A(rR)C; m1fwd => m2 should be rev
         my ($F,$T) = ($x=~/^(\d+)\-(\d+)$/);
#        print "<<<$txptid $tori>>>\n";
         $I_signal_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:I:$F-$T (A" . (($tori>0) ? "r" : "R") . "C)\n";
         $I_signal_seq_m1fwd .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALU : $ALU) . substr($ScaffSeq,$F-1,$T-$F+1);
  
         $I_signalY_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:I:$F-$T (A" . (($tori>0) ? "r" : "R") . "C)\n";
         $I_signalY_seq_m1fwd .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALUY : $ALUY) . substr($ScaffSeq,$F-1,$T-$F+1);

         ##### Exonic insertion: (rR)AC or A(rR)A; m1fwd => m2 should be rev 
         $X_signal_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:X:$F-$T (" . (($tori>0) ? "r" : "R") . "AC)\n"; 
         $X_signal_seq_m1fwd .= (($tori>0) ? $revALU : $ALU) . substr($ScaffSeq,$pos,$to-$pos) . substr($ScaffSeq,$F-1,$T-$F+1);

         $X_signalY_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:X:$F-$T (" . (($tori>0) ? "r" : "R") . "AC)\n";
         $X_signalY_seq_m1fwd .= (($tori>0) ? $revALUY : $ALUY) . substr($ScaffSeq,$pos,$to-$pos) . substr($ScaffSeq,$F-1,$T-$F+1);
   
         $XA_signal_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:X:$from-$to (A" . (($tori>0) ? "r" : "R") . "A)\n"; 
         $XA_signal_seq_m1fwd .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALU : $ALU) . substr($ScaffSeq,$pos,$to-$pos);

         $XA_signalY_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:X:$from-$to (A" . (($tori>0) ? "r" : "R") . "A)\n";
         $XA_signalY_seq_m1fwd .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALUY : $ALUY) . substr($ScaffSeq,$pos,$to-$pos);


         print S $I_signal_seq_m1fwd, "\n";
         print S $I_signalY_seq_m1fwd, "\n";
         print S $X_signal_seq_m1fwd, "\n";
         print S $X_signalY_seq_m1fwd, "\n";
         print S $XA_signal_seq_m1fwd, "\n";
         print S $XA_signalY_seq_m1fwd, "\n";
      }
   } else { # No next exon; this is the last exon; ARC->AR
      #####  Intronic insertion: A(rR)C; m1fwd => m2 should be rev
      # print "<<<$txptid $tori>>>\n";
      $I_signal_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:I:1-0 (A" . (($tori>0) ? "r" : "R") . "C)\n";
      $I_signal_seq_m1fwd .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALU : $ALU);
      
      $I_signalY_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:I:1-0 (A" . (($tori>0) ? "r" : "R") . "C)\n";
      $I_signalY_seq_m1fwd .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALUY : $ALUY);
      
      ##### Exonic insertion: (rR)AC or A(rR)A; m1fwd => m2 should be rev 
      $X_signal_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:X:1-0 (" . (($tori>0) ? "r" : "R") . "AC)\n";
      $X_signal_seq_m1fwd .= (($tori>0) ? $revALU : $ALU) . substr($ScaffSeq,$pos,$to-$pos);
      
      $X_signalY_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:X:1-0 (" . (($tori>0) ? "r" : "R") . "AC)\n";
      $X_signalY_seq_m1fwd .= (($tori>0) ? $revALUY : $ALUY) . substr($ScaffSeq,$pos,$to-$pos);
      
      $XA_signal_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:X:$from-$to (A" . (($tori>0) ? "r" : "R") . "A)\n";
      $XA_signal_seq_m1fwd .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALU : $ALU) . substr($ScaffSeq,$pos,$to-$pos);
      
      $XA_signalY_seq_m1fwd  = ">$txptid;$chrom:$from-$to" . " m1fwd:X:$from-$to (A" . (($tori>0) ? "r" : "R") . "A)\n";
      $XA_signalY_seq_m1fwd .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALUY : $ALUY) . substr($ScaffSeq,$pos,$to-$pos);
      
      
      print S $I_signal_seq_m1fwd, "\n";
      print S $I_signalY_seq_m1fwd, "\n";
      print S $X_signal_seq_m1fwd, "\n";
      print S $X_signalY_seq_m1fwd, "\n";
      print S $XA_signal_seq_m1fwd, "\n";
      print S $XA_signalY_seq_m1fwd, "\n";
   }

   # truly, should only be included when m1rev
   @elems = split ';', $ExonPrevPos{"$chrom:$from-$to"};

   if (scalar(@elems)) {  # there is at least one previous exon
      foreach my $x (@elems) {
         my ($F,$T) = ($x=~/^(\d+)\-(\d+)$/);

         ##### Intronic insertion: C(rR)A
         $I_signal_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:I:$F-$T (C" . (($tori>0) ? "r":"R") . "A)\n";
         $I_signal_seq_m1rev .= substr($ScaffSeq,$F-1,$T-$F+1) . (($tori>0) ? $revALU : $ALU) . substr($ScaffSeq,$pos,$to-$pos);

         $I_signalY_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:I:$F-$T (C" . (($tori>0) ? "r":"R") . "A)\n";
         $I_signalY_seq_m1rev .= substr($ScaffSeq,$F-1,$T-$F+1) . (($tori>0) ? $revALUY : $ALUY) . substr($ScaffSeq,$pos,$to-$pos);

         ##### Exonic insertion: CA(rR) or A(rR)A (the latter may have been calculated above, but we add it here so we can show the m1rev)

         $X_signal_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:X:$F-$T (CA" . (($tori>0) ? "r":"R") . ")\n";
         $X_signal_seq_m1rev .= substr($ScaffSeq,$F-1,$T-$F+1) . substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALU : $ALU);

         $X_signalY_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:X:$F-$T (CA" . (($tori>0) ? "r":"R") . ")\n";
         $X_signalY_seq_m1rev .= substr($ScaffSeq,$F-1,$T-$F+1) . substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALUY : $ALUY);

         $XA_signal_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:X:$from-$to (A" . (($tori>0) ? "r" : "R") . "A)\n";
         $XA_signal_seq_m1rev .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALU : $ALU) . substr($ScaffSeq,$pos,$to-$pos);

         $XA_signalY_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:X:$from-$to (A" . (($tori>0) ? "r" : "R") . "A)\n";
         $XA_signalY_seq_m1rev .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALUY : $ALUY) . substr($ScaffSeq,$pos,$to-$pos);

         print S $I_signal_seq_m1rev, "\n";
         print S $I_signalY_seq_m1rev, "\n";
         print S $X_signal_seq_m1rev, "\n";
         print S $X_signalY_seq_m1rev, "\n";
         print S $XA_signal_seq_m1rev, "\n";
         print S $XA_signalY_seq_m1rev, "\n";
      }
   } else {  # There is no previous exon; this is the first exon; CRA -> RA, CAR->AR
      ##### Intronic insertion: C(rR)A
      $I_signal_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:I:1-0 (C" . (($tori>0) ? "r":"R") . "A)\n";
      $I_signal_seq_m1rev .= (($tori>0) ? $revALU : $ALU) . substr($ScaffSeq,$pos,$to-$pos);
      
      $I_signalY_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:I:1-0 (C" . (($tori>0) ? "r":"R") . "A)\n";
      $I_signalY_seq_m1rev .= (($tori>0) ? $revALUY : $ALUY) . substr($ScaffSeq,$pos,$to-$pos);
      
      ##### Exonic insertion: CA(rR) or A(rR)A (the latter may have been calculated above, but we add here too so we can have m1rev)
      $X_signal_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:X:1-0 (CA" . (($tori>0) ? "r":"R") . ")\n";
      $X_signal_seq_m1rev .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALU : $ALU);
      
      $X_signalY_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:X:1-0 (CA" . (($tori>0) ? "r":"R") . ")\n";
      $X_signalY_seq_m1rev .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALUY : $ALUY);

      $XA_signal_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:X:$from-$to (A" . (($tori>0) ? "r" : "R") . "A)\n";
      $XA_signal_seq_m1rev .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALU : $ALU) . substr($ScaffSeq,$pos,$to-$pos);

      $XA_signalY_seq_m1rev  = ">$txptid;$chrom:$from-$to" . " m1rev:X:$from-$to (A" . (($tori>0) ? "r" : "R") . "A)\n";
      $XA_signalY_seq_m1rev .= substr($ScaffSeq,$pos,$to-$pos) . (($tori>0) ? $revALUY : $ALUY) . substr($ScaffSeq,$pos,$to-$pos);

      
      print S $I_signal_seq_m1rev, "\n";
      print S $I_signalY_seq_m1rev, "\n";
      print S $X_signal_seq_m1rev, "\n";
      print S $X_signalY_seq_m1rev, "\n";
      print S $XA_signal_seq_m1rev, "\n";
      print S $XA_signalY_seq_m1rev, "\n";
   }
   close(S);

   my $from1 = $from-1;

   ##### 'UNSPLICED'
   #####
   $pos = $to+500;
   my $unspliced_seq_m1fwd = ">$txptid;$chrom:$from-$to m1fwd:(m1+500):$pos (AU)\n";
   $unspliced_seq_m1fwd .= substr($ScaffSeq,$from1,$pos-$from+1);
   $pos = $from1-500;
   my $unspliced_seq_m1rev = ">$txptid;$chrom:$from-$to m1rev:(m1-500):$pos (UA)\n";
   $unspliced_seq_m1rev .= substr($ScaffSeq,$pos,$to-$pos);

   open(S, ">$SET.unspliced.fa") or die "$SET.unspliced.fa";
   print S $unspliced_seq_m1fwd, "\n";
   print S $unspliced_seq_m1rev, "\n";
   close(S);
   
   $from1 = $from-1;
   # generate region_seq sequences from the annotation, next coord + 100

   my ($region_seq_m1fwd,$region_seq_m1rev);
   if (defined($ExonFarNextPos{"$chrom:$from-$to"})) {
      $region_seq_m1fwd = ">$txptid;$chrom:$from-$to" . " m1fwd:(m1+" . $ExonFarNextPos{"$chrom:$from-$to"} . "+100):" . ($ExonFarNextPos{"$chrom:$from-$to"} +100) . " (AG)\n";
      $pos = $ExonFarNextPos{"$chrom:$from-$to"} + 100;
   } else {
      $region_seq_m1fwd = ">$txptid;$chrom:$from-$to" . " m1fwd:(m1+5000):" . ($to+5000) . " (AG)\n";
      $pos = $to + 5000;
   }
   $region_seq_m1fwd .= substr($ScaffSeq,$from1,$pos-$from1);

   if (defined($ExonFarPrevPos{"$chrom:$from-$to"})) {
      $region_seq_m1rev = ">$txptid;$chrom:$from-$to" . " m1rev:(m1-" . $ExonFarPrevPos{"$chrom:$from-$to"} . "-100):" . ($ExonFarPrevPos{"$chrom:$from-$to"}-100) . " (GA)\n";
      $pos = $ExonFarPrevPos{"$chrom:$from-$to"} - 100 -1;
   } else {
      $region_seq_m1rev = ">$txptid;$chrom:$from-$to" . " m1rev:(m1-5000):" . ($from-5000) . " (GA)\n";
      $pos = $from-5000;
   }
   $region_seq_m1rev .= substr($ScaffSeq,$pos,$to-$pos);

   open(S, ">$SET.region.fa") or die "$SET.region.fa";
   print S $region_seq_m1fwd, "\n";
   print S $region_seq_m1rev, "\n";
   close(S);
   

   # run the searches
   # option: -minid 80 -minlength 22
   `sim4db -gen $SET.signal.fa -cdna $SET.tmpreads.fa -aligns -o - > $SET.tmpreads.signal.sim4db`;
   `sim4db -gen $SET.unspliced.fa -cdna $SET.tmpreads.fa -aligns -o - > $SET.tmpreads.unspliced.sim4db`;
   `sim4db -gen $SET.region.fa -cdna $SET.tmpreads.fa -aligns -o - > $SET.tmpreads.region.sim4db`;
   `cat $SET.tmpreads.signal.sim4db >> $SET.SIGNAL.sim4db`;
   `cat $SET.tmpreads.unspliced.sim4db >> $SET.UNSPLICED.sim4db`;
   `cat $SET.tmpreads.region.sim4db >> $SET.REGION.sim4db`;
#  `rm $SET.signal.fa $SET.region.fa $SET.unspliced.fa $SET.signal.faidx $SET.region.faidx $SET.unspliced.faidx $SET.tmpreads.fa $SET.tmpreads.faidx`;

   # analyze the searches
}

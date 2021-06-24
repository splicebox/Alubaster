#!/usr/bin/perl
use strict;

my $Usage =
"Usage: $0 [options] annotindex overlaps selection SIGNAL REGION UNSPLICED
    [options]
    -enum | -b mins minsn mins2sn mins2sru mins2ru (default: -enum)
               mutually exclusive; list current values (-enum) vs apply filter (-b)
    -debug     print signal, signanone, regionunspliced information in BED files
               (default: off)\n
    -p prefix  used with -debug; use this as prefix for debug file (default:AA)\n
    [required]
    annotindex txpt2gene file
    overlaps   overlaps file
    selection  selection file
    SIGNAL     SIGNAL nr file
    REGION     REGION nr file
    UNSPLICED  UNSPLICED nr file
";
($#ARGV>2) or die $Usage;

my $isEnum = 1;
my $isDebug = 0;
my $prefix = "AA";

my ($MIN_SIGNALS,$MIN_SIGNALNONE,$MIN_S2SN,$MIN_S2SRU,$MIN_S2RU);
my $prefix;

my $argtmp = shift;
while ($argtmp=~/^\-/) {
   if ($argtmp eq "-enum") {
      $isEnum = 1;
      ($MIN_SIGNALS,$MIN_SIGNALNONE,$MIN_S2SN,$MIN_S2SRU,$MIN_S2RU) = (0,0,0,0,0);
   } elsif ($argtmp eq "-b") {
      $isEnum = 0;
      $MIN_SIGNALS = shift; 
      $MIN_SIGNALNONE = shift;
      $MIN_S2SN = shift;
      $MIN_S2SRU = shift;
      $MIN_S2RU = shift;
   } elsif ($argtmp eq "-p") {
      $prefix = shift;
   } elsif ($argtmp eq "-debug") {
      $isDebug = 1;
   }
   $argtmp = shift; 
}
my $Txpt2GeneFile = $argtmp;
my $OverlapsFile = shift;
my $SelectionFile = shift;
my $SignalFile = shift;
my $RegionFile = shift;
my $UnsplicedFile = shift;


# sort by chrom,genename: groups, bed, Selection, SIGNAL.nr, UNSPLICED.nr, REGION.nr

if ($isDebug) {
   open(DBG, ">$prefix.debug") or die "$prefix.debug.\n";
}


my %isSignal;
my %isCandidate;
my %isRegionUnspliced;
my %GeneId2Name;

### Record 'passing' signals: Y,N,N
open(F, "<$SelectionFile") or die "Could not open selection file $SelectionFile.\n";
while (<F>) {
   chomp;

   #SRR1284895.19496361/2:- ENSG00000000003.13:ENST00000612152.3;chrX:100620211-100623088 Yes TSPAN6
   /^(\S+)(\d):\S (\S+):\S+;\S+ Yes (\S+)$/ or die "died. $_";
   my ($readid,$idx,$geneid,$genename)= ($1,$2,$3,$4);
   $readid = $readid . (3-$idx); 
   $isSignal{"$geneid:$readid"} = 1;
#  print "Signal: <<<", "$geneid:$readid", ">>>\n";
}
close(F);

# Record candidates: Y,*,*
open(F, "<$SignalFile") or die "Could not open SIGNAL file $SignalFile.\n";
while (<F>) {
   chomp;

   # SRR1284895.123549302/2:+ ENSG00000000003.13:ENST00000612152.3;chrX:100620211-100623088 No 
   next if !/Yes/;
   /^(\S+)(\d):\S (\S+):\S+;\S+/ or die "died. $_";
   my ($readid,$idx,$geneid) = ($1,$2,$3);
   $readid = $readid . (3-$idx);
   $isCandidate{"$geneid:$readid"} = 1;
}
close(F);
# End Record candidates

### Record "region-unspliced": *,Y,* or *,*,Y
open(F, "<$RegionFile") or die "Could not open REGION file $RegionFile.\n";
while (<F>) {
   chomp;

   # SRR1284895.161017365/2:- ENSG00000000003.13:ENST00000612152.3;chrX:100620211-100623088 Yes Unspliced
   next if !/ Yes/;
   /^(\S+)(\d):\S (\S+):\S+;\S+/ or die "died. $_";
   my ($readid,$idx,$geneid) = ($1,$2,$3);
   $readid = $readid . (3-$idx);
   $isRegionUnspliced{"$geneid:$readid"} = 1;
}
close(F);

open(F, "<$UnsplicedFile") or die "Could not open UNSPLICED file $UnsplicedFile.\n";
while (<F>) {
   chomp;

   # SRR1284895.161017365/2:- ENSG00000000003.13:ENST00000612152.3;chrX:100620211-100623088 Yes
   next if !/ Yes/;
   /^(\S+)(\d):\S (\S+):\S+;\S+/ or die "died. $_";
   my ($readid,$idx,$geneid) = ($1,$2,$3);
   $readid = $readid . (3-$idx);
   $isRegionUnspliced{"$geneid:$readid"} = 1;
}
close(F);
# End Record "region-unspliced"

my %Txpt2Gene;
open(F, "<$Txpt2GeneFile") or die "Could not open Txpt2Gene file $Txpt2GeneFile";
while (<F>) {
   # ENSG00000000003.13"; "ENST00000373020.7"; "TSPAN6"; chr1 61324347 61330001
   chomp;
   /^\"(\S+)\"; \"(\S+)\"; \"(\S+)\";/ or die "died. $_";
   $Txpt2Gene{$2} = $1;
   $GeneId2Name{$1} = $3;
}
close(F);

my %Seen;
my $lastgene;
my $lastchrom;
my $Lines;
my @signals_plus;
my @signalnone_plus;
my @regionunspliced_plus;
my @signalregionunspliced_plus;
my @signals_minus;
my @signalnone_minus;
my @regionunspliced_minus;
my @signalregionunspliced_minus;

#$prefix = "AA"; #"SRR1284895sim";
open(S, ">$prefix.signals.bed") or die "signals.bed";
open(SN, ">$prefix.signalnone.bed") or die "signalnone.bed";
open(RU, ">$prefix.regionunspliced.bed") or die "regionunspliced.bed";
open(SRU, ">$prefix.signalregionunspliced.bed") or die "signalregionunspliced.bed";

######### Now read in the overlaps bed file sorted by gene, and break into genes
open(F, "<$OverlapsFile") or die "Could not open OverlapsFile $OverlapsFile\n";
while (<F>) {
    chomp;
    # chr10   100146454       100147759       SRR1284895.134062683/2  50      -       100146454       100147759       0,0,0   2       73,27,  0,1278, chr10   100146448       100146527       ENST00000421367.5       0       -       73
    # chr10   100142146       100142246       SRR1284895.22012/1      50      +       100142146       100142246       0,0,0   1       100,    0,      chr10   100141682       100143940       ENST00000421367.5       0       -       100

   /^(\S+)\t(\d+)\t(\d+)\t(\S+)\t\S+\t(\S)\t\d+\t\d+\t\S+\t(\d+)\t(\S+)\t(\S+)\t\S+\t\d+\t\d+\t(\S+)\t/ or die "died. $_";
   my ($chrom,$from,$to,$readid,$ori,$nexons,$xlen,$xoffset,$tid) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
   if ($Txpt2Gene{$tid} ne $lastgene) {
      if (defined($lastgene)) {
          my $tmp_clusters_signal_plus = &cluster_counts(scalar(@signals_plus), @signals_plus);
          my $tmp_clusters_signalnone_plus = &cluster_counts(scalar(@signalnone_plus),@signalnone_plus);
          my $tmp_clusters_regionunspliced_plus = &cluster_counts(scalar(@regionunspliced_plus),@regionunspliced_plus);
          my $tmp_clusters_signalregionunspliced_plus = &cluster_counts(scalar(@signalregionunspliced_plus),@signalregionunspliced_plus);
          my $tmp_clusters_signal_minus = &cluster_counts(scalar(@signals_minus), @signals_minus);
          my $tmp_clusters_signalnone_minus = &cluster_counts(scalar(@signalnone_minus),@signalnone_minus);
          my $tmp_clusters_regionunspliced_minus = &cluster_counts(scalar(@regionunspliced_minus),@regionunspliced_minus);
          my $tmp_clusters_signalregionunspliced_minus = &cluster_counts(scalar(@signalregionunspliced_minus),@signalregionunspliced_minus);

          if ($isDebug) {
             &print_debug("$lastgene", "S", "+", $chrom, $tmp_clusters_signal_plus);
             &print_debug("$lastgene", "SN", "+", $chrom, $tmp_clusters_signalnone_plus);
             &print_debug("$lastgene", "RU", "+", $chrom, $tmp_clusters_regionunspliced_plus);
             &print_debug("$lastgene", "SRU", "+", $chrom, $tmp_clusters_signalregionunspliced_plus);
             &print_debug("$lastgene", "S", "-", $chrom, $tmp_clusters_signal_minus);
             &print_debug("$lastgene", "SN", "-", $chrom, $tmp_clusters_signalnone_minus);
             &print_debug("$lastgene", "RU", "-", $chrom, $tmp_clusters_regionunspliced_minus);
             &print_debug("$lastgene", "SRU", "+", $chrom, $tmp_clusters_signalregionunspliced_minus);
          }

          &find_insertion_loc_filter($lastgene,$lastchrom,$tmp_clusters_signal_plus,$tmp_clusters_signalnone_plus,$tmp_clusters_regionunspliced_plus,$tmp_clusters_signalregionunspliced_plus,$tmp_clusters_signal_minus,$tmp_clusters_signalnone_minus,$tmp_clusters_regionunspliced_minus,$tmp_clusters_signalregionunspliced_minus,$Lines);
          # print "$lastgene $lastchrom $from $to\n";
      }
      $lastgene = $Txpt2Gene{$tid};
      $lastchrom = $chrom;
      $Lines = "";
      @signals_plus = (); @signalnone_plus = (); @regionunspliced_plus = (); @signalregionunspliced_plus = ();
      @signals_minus = (); @signalnone_minus = (); @regionunspliced_minus = (); @signalregionunspliced_minus = ();
   }
   next if (defined($Seen{"$lastgene:$readid"}));
   # Note: the above would remove duplicate entries that come from the same interval (read) in set1 
   # overlapping multiple entries in set 2 (exons); however, the code below still produced multiple
   # lines for the same alignment - namely, one for each 'block' in the alignment. Changed that
   # by printing only when i==0 (except for signal, which is only printed one, for either the 1st or the last exon..

   # add read to the clusters, if not Seen; mark them as Seen;
   # if SIgnal read, only record the last interval in the signals array, and both intervals in teh signalnone array; all others, record all intervals
   $Lines .= $_ . "\n";

   my @xlens = split ',', $xlen;
   my @xoffsets =  split ',', $xoffset;
   if ($ori eq "+") {
      foreach (my $i=0; $i<$nexons; $i++) {
      
         if (defined($isRegionUnspliced{"$lastgene:$readid"})) {
            push @regionunspliced_plus, ($from+$xoffsets[$i]) . ":" . ($from+$xoffsets[$i]+$xlens[$i]-1);
            if ($isDebug) { $i or print RU $_, "\t$lastgene:", $GeneId2Name{$lastgene}, "\t+\n"; }

            if (defined($isCandidate{"$lastgene:$readid"})) {
               push @signalregionunspliced_plus, ($from+$xoffsets[$i]) . ":" . ($from+$xoffsets[$i]+$xlens[$i]-1);
               if ($isDebug) { $i or print SRU $_, "\t$lastgene:", $GeneId2Name{$lastgene}, \"t+\n"; }
            }
         } elsif (defined($isSignal{"$lastgene:$readid"})) {
            if ((($ori eq "+") && ($i==$nexons-1)) || (($ori eq "-") && ($i==0))) {
               push @signals_plus, ($from+$xoffsets[$i]) . ":" . ($from+$xoffsets[$i]+$xlens[$i]-1);
               if ($isDebug) { print S $_, "\t$lastgene:", $GeneId2Name{$lastgene}, "\t+\n"; }
            }
            push @signalnone_plus, ($from+$xoffsets[$i]) . ":" . ($from+$xoffsets[$i]+$xlens[$i]-1);
            if ($isDebug) { $i or print SN $_, "\t$lastgene:", $GeneId2Name{$lastgene}, "\t+\n"; }
         } else { # SignalNone
            push @signalnone_plus, ($from+$xoffsets[$i]) . ":" . ($from+$xoffsets[$i]+$xlens[$i]-1);
            if ($isDebug) { $i or print SN $_, "\t$lastgene:", $GeneId2Name{$lastgene}, "\t+\n"; }
         }
      }
   } elsif ($ori eq "-") {
      foreach (my $i=0; $i<$nexons; $i++) {
      
         if (defined($isRegionUnspliced{"$lastgene:$readid"})) {
            push @regionunspliced_minus, ($from+$xoffsets[$i]) . ":" . ($from+$xoffsets[$i]+$xlens[$i]-1);
            if ($isDebug) { $i or print RU $_, "\t$lastgene:", $GeneId2Name{$lastgene}, "\t-\n"; }

            if (defined($isCandidate{"$lastgene:$readid"})) {
               push @signalregionunspliced_minus, ($from+$xoffsets[$i]) . ":" . ($from+$xoffsets[$i]+$xlens[$i]-1);
               if ($isDebug) { $i or print SRU $_, "\t$lastgene:", $GeneId2Name{$lastgene}, "\t-\n"; }
            }
         } elsif (defined($isSignal{"$lastgene:$readid"})) {
            if ((($ori eq "+") && ($i==$nexons-1)) || (($ori eq "-") && ($i==0))) {
               push @signals_minus, ($from+$xoffsets[$i]) . ":" . ($from+$xoffsets[$i]+$xlens[$i]-1);
               if ($isDebug) { print S $_, "\t$lastgene:", $GeneId2Name{$lastgene}, "\t-\n"; }
            }
            push @signalnone_minus, ($from+$xoffsets[$i]) . ":" . ($from+$xoffsets[$i]+$xlens[$i]-1);
            if ($isDebug) { $i or print SN $_, "\t$lastgene:", $GeneId2Name{$lastgene}, "\t-\n"; }
         } else { # SignalNone
            push @signalnone_minus, ($from+$xoffsets[$i]) . ":" . ($from+$xoffsets[$i]+$xlens[$i]-1);
            if ($isDebug) { $i or print SN $_, "\t$lastgene:", $GeneId2Name{$lastgene}, "\t-\n"; }
         }
      }
   } else {
      die "unrecognized ori. $ori\n";
   }
   $Seen{"$lastgene:$readid"} = 1;
}
if (defined($lastgene)) {
   my $tmp_clusters_signal_plus = &cluster_counts(scalar(@signals_plus),@signals_plus);
   my $tmp_clusters_signalnone_plus = &cluster_counts(scalar(@signalnone_plus),@signalnone_plus);
   my $tmp_clusters_regionunspliced_plus = &cluster_counts(scalar(@regionunspliced_plus),@regionunspliced_plus);
   my $tmp_clusters_signalregionunspliced_plus = &cluster_counts(scalar(@signalregionunspliced_plus),@signalregionunspliced_plus);
   my $tmp_clusters_signal_minus = &cluster_counts(scalar(@signals_minus),@signals_minus);
   my $tmp_clusters_signalnone_minus = &cluster_counts(scalar(@signalnone_minus),@signalnone_minus);
   my $tmp_clusters_regionunspliced_minus = &cluster_counts(scalar(@regionunspliced_minus),@regionunspliced_minus);
   my $tmp_clusters_signalregionunspliced_minus = &cluster_counts(scalar(@signalregionunspliced_minus),@signalregionunspliced_minus);

   if ($isDebug) {
      &print_debug("$lastgene", "S", "+", $lastchrom, $tmp_clusters_signal_plus);
      &print_debug("$lastgene", "SN", "+", $lastchrom, $tmp_clusters_signalnone_plus);
      &print_debug("$lastgene", "RU", "+", $lastchrom, $tmp_clusters_regionunspliced_plus);
      &print_debug("$lastgene", "SRU", "+", $lastchrom, $tmp_clusters_signalregionunspliced_plus);
      &print_debug("$lastgene", "S", "-", $lastchrom, $tmp_clusters_signal_minus);
      &print_debug("$lastgene", "SN", "-", $lastchrom, $tmp_clusters_signalnone_minus);
      &print_debug("$lastgene", "RU", "-", $lastchrom, $tmp_clusters_regionunspliced_minus);
      &print_debug("$lastgene", "SRU", "-", $lastchrom, $tmp_clusters_signalregionunspliced_minus);
   }

   my ($from,$to) = &find_insertion_loc_filter($lastgene,$lastchrom,$tmp_clusters_signal_plus,$tmp_clusters_signalnone_plus,$tmp_clusters_regionunspliced_plus,$tmp_clusters_signalregionunspliced_plus,$tmp_clusters_signal_minus,$tmp_clusters_signalnone_minus,$tmp_clusters_regionunspliced_minus,$tmp_clusters_signalregionunspliced_minus,$Lines);
}
close(F); 

close(S); close(SN); close(RU);

##### Subroutines
my $W = 0;
sub cluster_counts {
   my $n = shift @_; $n or return "";
   my @I = sort byPair @_;
   
   my ($b,$e);
   my $count = 0;
   my $clusters = "";
    
   foreach my $x (@I) {
      my ($from,$to) = split ':', $x;
      if (!defined($b)) {
         $b = $from; $e = $to; $count = 0;
         $clusters = "";
      } elsif ($from>$e+$W) {
         if ($clusters eq "") {
            $clusters = "$b:$e:$count";
         } else {
            $clusters .= ",$b:$e:$count";
         }
         $count = 0; $b = $from; $e = $to;
      } else {
         $e = &max($to,$e);
      }
      $count++;
   }
   if ($clusters eq "") {
      $clusters = "$b:$e:$count";
   } else {
      $clusters .= ",$b:$e:$count";
   }

   return $clusters;
}

sub byPair {
   my @temp_a = split ':', $a;
   my @temp_b = split ':', $b;

   if ($temp_a[0]!=$temp_b[0]) { $temp_a[0] <=> $temp_b[0]; }
   else { $temp_a[1] <=> $temp_b[1]; }
}

sub max { return (($_[0]>=$_[1]) ? $_[0] : $_[1]); }



sub find_insertion_loc_filter { # geneid chrom $tmp_clusters_signal+,$tmp_clusters_signalnone+,$tmp_clusters_regionunspliced+,$tmp_clusters_signalregionunspliced+,$tmp_clusters_signal-,$tmp_clusters_signalnone-,$tmp_clusters_regionunspliced-,$tmp_clusters_signalregionunspliced-,Lines
    my ($From,$To);

    my $geneid = shift;
    my $chrom = shift;
    my $tmp_clusters_signal_plus = shift;
    my $tmp_clusters_signalnone_plus = shift;
    my $tmp_clusters_regionunspliced_plus = shift;
    my $tmp_clusters_signalregionunspliced_plus = shift;
    my $tmp_clusters_signal_minus = shift;
    my $tmp_clusters_signalnone_minus = shift;
    my $tmp_clusters_regionunspliced_minus = shift;
    my $tmp_clusters_signalregionunspliced_minus = shift;
    my $Lines = shift;

    my @signalFplus = (); my @signalTplus = (); my @signalCplus = ();
    my @signalFminus = (); my @signalTminus = (); my @signalCminus = ();
    my @signalnoneFplus = (); my @signalnoneTplus = (); my @signalnoneCplus = ();
    my @signalnoneFminus = (); my @signalnoneTminus = (); my @signalnoneCminus = ();
    my @regionunsplicedFplus = (); my @regionunsplicedTplus = (); my @regionunsplicedCplus = ();
    my @regionunsplicedFminus = (); my @regionunsplicedTminus = (); my @regionunsplicedCminus = ();
    my @signalregionunsplicedFplus = (); my @signalregionunsplicedTplus = (); my @signalregionunsplicedCplus = ();
    my @signalregionunsplicedFminus = (); my @signalregionunsplicedTminus = (); my @signalregionunsplicedCminus = ();

    my @signalIDplus = (); my @signalIDminus = ();
    my @signalnoneIDplus = (); my @signalnoneIDminus = ();
    my @regionunsplicedIDplus = (); my @regionunsplicedIDminus = ();

    my @elems = split ',', $tmp_clusters_signal_plus;
    foreach my $x (@elems) {
       $x =~ /^(\d+):(\d+):(\d+)$/ or die "died. $x ";
       push @signalFplus, $1; push @signalTplus, $2; push @signalCplus, $3;
    }
    @elems = split ',', $tmp_clusters_signal_minus;
    foreach my $x (@elems) {
       $x =~ /^(\d+):(\d+):(\d+)$/ or die "died. $x ";
       push @signalFminus, $1; push @signalTminus, $2; push @signalCminus, $3;
    }

    @elems = split ',', $tmp_clusters_signalnone_plus;
    foreach my $x (@elems) {
       $x =~ /^(\d+):(\d+):(\d+)$/ or die "died. $x ";
       push @signalnoneFplus, $1; push @signalnoneTplus, $2; push @signalnoneCplus, $3;
    }
    @elems = split ',', $tmp_clusters_signalnone_minus;
    foreach my $x (@elems) {
       $x =~ /^(\d+):(\d+):(\d+)$/ or die "died. $x ";
       push @signalnoneFminus, $1; push @signalnoneTminus, $2; push @signalnoneCminus, $3;
    }

    @elems = split ',', $tmp_clusters_regionunspliced_plus;
    foreach my $x (@elems) {
       $x =~ /^(\d+):(\d+):(\d+)$/ or die "died. $x ";
       push @regionunsplicedFplus, $1; push @regionunsplicedTplus, $2; push @regionunsplicedCplus, $3;
    }
    @elems = split ',', $tmp_clusters_regionunspliced_minus;
    foreach my $x (@elems) {
       $x =~ /^(\d+):(\d+):(\d+)$/ or die "died. $x ";
       push @regionunsplicedFminus, $1; push @regionunsplicedTminus, $2; push @regionunsplicedCminus, $3;
    }

    @elems = split ',', $tmp_clusters_signalregionunspliced_plus;
    foreach my $x (@elems) {
       $x =~ /^(\d+):(\d+):(\d+)$/ or die "died. $x ";
       push @signalregionunsplicedFplus, $1; push @signalregionunsplicedTplus, $2; push @signalregionunsplicedCplus, $3;
    }
    @elems = split ',', $tmp_clusters_signalregionunspliced_minus;
    foreach my $x (@elems) {
       $x =~ /^(\d+):(\d+):(\d+)$/ or die "died. $x ";
       push @signalregionunsplicedFminus, $1; push @signalregionunsplicedTminus, $2; push @signalregionunsplicedCminus, $3;
    }

    # now assign readids to intervals
    my @elems = split '\n', $Lines;
#   print "<<<$Lines>>>\n";
    foreach my $l (@elems) {
       $l =~ /^(\S+)\t(\d+)\t(\d+)\t(\S+)\t\S+\t(\S)\t\d+\t\d+\t\S+\t(\d+)\t(\S+)\t(\S+)\t\S+\t\d+\t\d+\t(\S+)\t/ or die "died. $l";
       my ($chrom,$from,$to,$readid,$ori,$nexons,$xlen,$xoffset,$tid) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);

       my @xlens = split ',', $xlen;
       my @xoffsets =  split ',', $xoffset;
       if ($ori eq "+") {
         # foreach exon, find (an) interval
         for (my $j=0; $j<$nexons; $j++) {

           # find the interval
           if (defined($isRegionUnspliced{"$geneid:$readid"})) {
              for (my $i=0; $i<scalar(@regionunsplicedFplus); $i++) {
                 if (($regionunsplicedFplus[$i]<=$from+$xoffsets[$j]) && ($from+$xoffsets[$j]<=$regionunsplicedTplus[$i])) {
                    $regionunsplicedIDplus[$i] .= ",$readid";
                    last;
                 }
              }
           } elsif (defined($isSignal{"$geneid:$readid"})) {
              for (my $i=0; $i<scalar(@signalFplus); $i++) { 
                 if (($signalFplus[$i]<=$from+$xoffsets[$j]) && ($from+$xoffsets[$j]<=$signalTplus[$i])) {
                    $signalIDplus[$i] .= ",$readid";
                    last;
                 }
              }
              for (my $i=0; $i<scalar(@signalnoneFplus); $i++) {         
                 if (($signalnoneFplus[$i]<=$from+$xoffsets[$j]) && ($from+$xoffsets[$j]<=$signalnoneTplus[$i])) {
                    $signalnoneIDplus[$i] .= ",$readid";
                    last;
                 }
              }
           } else { 
              for (my $i=0; $i<scalar(@signalnoneFplus); $i++) {
                 if (($signalnoneFplus[$i]<=$from+$xoffsets[$j]) && ($from+$xoffsets[$j]<=$signalnoneTplus[$i])) {
                    $signalnoneIDplus[$i] .= ",$readid";
                    last;
                 }
              }
           }
         }
       # minus strand
       } elsif ($ori eq "-") {
         for (my $j=0; $j<$nexons; $j++) {

           #find the interval
           if (defined($isRegionUnspliced{"$geneid:$readid"})) {
             for (my $i=0; $i<scalar(@regionunsplicedFminus); $i++) {
               if (($regionunsplicedFminus[$i]<=$from+$xoffsets[$j]) && ($from+$xoffsets[$j]<=$regionunsplicedTminus[$i])) {
                  $regionunsplicedIDminus[$i] .= ",$readid";
                  last;
               }
             }
           } elsif (defined($isSignal{"$geneid:$readid"})) {
             for (my $i=0; $i<scalar(@signalFminus); $i++) {
               if (($signalFminus[$i]<=$from+$xoffsets[$j]) && ($from+$xoffsets[$j]<=$signalTminus[$i])) {
                  $signalIDminus[$i] .= ",$readid";
                  last;
               }
             }
             for (my $i=0; $i<scalar(@signalnoneFminus); $i++) {
               if (($signalnoneFminus[$i]<=$from+$xoffsets[$j]) && ($from+$xoffsets[$j]<=$signalnoneTminus[$i])) {
                  $signalnoneIDminus[$i] .= ",$readid";
                  last;
               }
             }
           } else {
               for (my $i=0; $i<scalar(@signalnoneFminus); $i++) {
                  if (($signalnoneFminus[$i]<=$from+$xoffsets[$j]) && ($from+$xoffsets[$j]<=$signalnoneTminus[$i])) {
                     $signalnoneIDminus[$i] .= ",$readid";
                     last;
                  }
               }
           }
         } # for j
       } # ori eq "-:
    } # top

    @elems = ();

    #### start with + signals
    foreach (my $i=0; $i<scalar(@signalFplus); $i++) {
       my ($j,$pi);

       my ($snum,$snnum,$runum,$srunum) = (0,0,0,0); 

       $snum = $signalCplus[$i];
       # find out if it overlaps with a signalnone - it should!
       foreach ($j=0; $j<scalar(@signalnoneFplus); $j++) {
          if (&overlap($signalFplus[$i],$signalTplus[$i],$signalnoneFplus[$j],$signalnoneTplus[$j])>=10) { $snnum += $signalnoneCplus[$j]; }
       }
       # find out if it overlaps with a regionunspliced 
       foreach ($j=0; $j<scalar(@regionunsplicedFplus); $j++) {
          if (&overlap($signalFplus[$i],$signalTplus[$i],$regionunsplicedFplus[$j],$regionunsplicedTplus[$j])>=10) { $runum += $regionunsplicedCplus[$j]; }
          # should I add overlap with region minus?  LLL
       }
       # find out if it overlaps with a signalregionunspliced 
       foreach ($j=0; $j<scalar(@signalregionunsplicedFplus); $j++) {
          if (&overlap($signalFplus[$i],$signalTplus[$i],$signalregionunsplicedFplus[$j],$signalregionunsplicedTplus[$j])>=10) { $srunum += $signalregionunsplicedCplus[$j]; }
       }

       # search for nearest opposite sign signal_none cluster 
       foreach ($pi=0; ($pi<scalar(@signalnoneFminus)) && ($signalnoneFminus[$pi]<$signalTplus[$i]-10); $pi++) { ; }
       my $isPaired = ($pi==scalar(@signalnoneFminus)) ? 0 : 1;
       
       # if paired, then add the counts
       if ($isPaired) {
          # find out if it overlaps with a signal
          foreach ($j=0; $j<scalar(@signalFplus); $j++) {
             if (&overlap($signalnoneFminus[$pi],$signalnoneTminus[$pi],$signalFminus[$j],$signalTminus[$j])>=10) { $snum += $signalCminus[$j]; }
          }
          # find out if it overlaps with a signal none - it should!
          foreach ($j=0; $j<scalar(@signalnoneFplus); $j++) {
             if (&overlap($signalnoneFminus[$pi],$signalnoneTminus[$pi],$signalnoneFminus[$j],$signalnoneTminus[$j])>=10) { $snnum += $signalnoneCminus[$j]; }
          }
          # find out if it overlaps with a regionunspliced
          foreach ($j=0; $j<scalar(@regionunsplicedFplus); $j++) {
             if (&overlap($signalnoneFminus[$pi],$signalnoneTminus[$pi],$regionunsplicedFminus[$j],$regionunsplicedTminus[$j])>=10) { $runum += $regionunsplicedCminus[$j]; }
          }
          # find out if it overlaps with a signalregionunspliced
          foreach ($j=0; $j<scalar(@signalregionunsplicedFplus); $j++) {
             if (&overlap($signalnoneFminus[$pi],$signalnoneTminus[$pi],$signalregionunsplicedFminus[$j],$signalregionunsplicedTminus[$j])>=10) { $srunum += $signalregionunsplicedCminus[$j]; }
          }
       }

       my $st = &filter($geneid,$snum,$snnum,$runum,$srunum);

       if ($st) {
          my $j;
          foreach ($j=0; $j<scalar(@signalnoneFplus); $j++) {
            last if (($signalnoneFplus[$j]<=$signalFplus[$i]) && ($signalTplus[$i]<=$signalnoneTplus[$j]));
          }
          ($j<scalar(@signalnoneFplus)) or die "Unexpected (1-2): $chrom $j>=", scalar(@signalnoneFplus), "<<<", $signalFplus[$i],":",$signalTplus[$i],">>>!!!",
                                               join(',', @signalnoneFplus), ";", join(',', @signalnoneTplus),"!!!\n"; 


          if (!$isPaired) { # report end of plus cluster, +inf
             print "(1) ", $GeneId2Name{$geneid}, " $geneid $chrom ", $signalTplus[$i], " ", "+inf [$snum $snnum $runum $srunum] ";
             # print "[", $signalIDplus[$i], "]\n"; 
             print "[", $signalnoneIDplus[$j], "] ";
             print "[", $signalFplus[$i], "-", $signalTplus[$i], "]\n";   # interval
          } else {  # found
             # what if too far? LLL - might need to look at next exon
             if ($signalnoneFminus[$pi]-$signalTplus[$i]<25000) {
                 print "(2) ", $GeneId2Name{$geneid}, " $geneid $chrom ", $signalTplus[$i], " ", $signalnoneFminus[$pi], " [$snum $snnum $runum $srunum] ";
                 #print "[", $signalIDplus[$i], ";",$signalnoneIDminus[$pi], "]\n";
                 print "[", $signalnoneIDplus[$j], ";",$signalnoneIDminus[$pi], "] ";
                 print "[", $signalFplus[$i], "-", $signalTplus[$i], ",", $signalnoneFminus[$pi], "-", $signalnoneTminus[$pi], "]\n";  # interval
             }
          }
       }
   }

   #### Now the - strand signals - might report duplicates of above, but ok, willsolve later LLL
   foreach (my $i=0; $i<scalar(@signalFminus); $i++) {
      my ($j,$pi);
        
      my ($snum,$snnum,$runum,$srunum) = (0,0,0,0);
       
      $snum = $signalCminus[$i];
       
      # find out if it overlaps with a signalnone - it should!
      foreach ($j=0; $j<scalar(@signalnoneFminus); $j++) {
         if (&overlap($signalFminus[$i],$signalTminus[$i],$signalnoneFminus[$j],$signalnoneTminus[$j])>=10) { $snnum += $signalnoneCminus[$j]; }
      }
      # find out if it overlaps with a regionunspliced 
      foreach ($j=0; $j<scalar(@regionunsplicedFminus); $j++) {
         if (&overlap($signalFminus[$i],$signalTminus[$i],$regionunsplicedFminus[$j],$regionunsplicedTminus[$j])>=10) { $runum += $regionunsplicedCminus[$j]; }
         # should I add overlap with region plus?  LLL
      }
      # find out if it overlaps with a regionunspliced
      foreach ($j=0; $j<scalar(@signalregionunsplicedFminus); $j++) {
         if (&overlap($signalFminus[$i],$signalTminus[$i],$signalregionunsplicedFminus[$j],$signalregionunsplicedTminus[$j])>=10) { $srunum += $signalregionunsplicedCminus[$j]; }
         # should I add overlap with region plus?  LLL
      }
      
       
      # search for nearest opposite sign signal_none cluster 
      foreach ($pi=scalar(@signalnoneFplus)-1; ($pi>=0) && ($signalnoneTplus[$pi]>$signalFminus[$i]+10); $pi--) { ; }
      my $isPaired = ($pi<0) ? 0 : 1;

      # if paired, then add the counts
      if ($isPaired) {
         # find out if it overlaps with a signal
         foreach ($j=0; $j<scalar(@signalFplus); $j++) {
            if (&overlap($signalnoneFplus[$pi],$signalnoneTplus[$pi],$signalFplus[$j],$signalTplus[$j])>=10) { $snum += $signalCplus[$j]; }
         }
         # find out if it overlaps with a signal none - it should!
         foreach ($j=0; $j<scalar(@signalnoneFplus); $j++) {
            if (&overlap($signalnoneFplus[$pi],$signalnoneTplus[$pi],$signalnoneFplus[$j],$signalnoneTplus[$j])>=10) { $snnum += $signalnoneCplus[$j]; }
         }
         # find out if it overlaps with a regiounspliced
         foreach ($j=0; $j<scalar(@regionunsplicedFplus); $j++) {
            if (&overlap($signalnoneFplus[$pi],$signalnoneTplus[$pi],$regionunsplicedFplus[$j],$regionunsplicedTplus[$j])>=10) { $snnum += $regionunsplicedCplus[$j]; }
         }
         # find out if it overlaps with a signalregiounspliced
         foreach ($j=0; $j<scalar(@signalregionunsplicedFplus); $j++) {
            if (&overlap($signalnoneFplus[$pi],$signalnoneTplus[$pi],$signalregionunsplicedFplus[$j],$signalregionunsplicedTplus[$j])>=10) { $snnum += $signalregionunsplicedCplus[$j]; }
         }
      }

      my $st = &filter($geneid,$snum,$snnum,$runum,$srunum);

      if ($st) {
          my $j;
          foreach ($j=0; $j<scalar(@signalnoneFminus); $j++) {
            last if (($signalnoneFminus[$j]<=$signalFminus[$i]) && ($signalTminus[$i]<=$signalnoneTminus[$j]));
          }  
          ($j<scalar(@signalnoneFminus)) or die "Unexpected (3-4): $chrom $j>=", scalar(@signalnoneFminus), "<<<",$signalFminus[$i],":",$signalTminus[$i],">>>!!!",
                                               join(':', @signalnoneFminus), ",", join(':', @signalnoneTminus),"!!!\n";

          if (!$isPaired) { # report -inf, start of minus cluster
             print "(3) ", $GeneId2Name{$geneid}, " $geneid $chrom -inf ", $signalFminus[$i], " [$snum $snnum $runum $srunum] ";
             #print "[", $signalIDminus[$i], "]\n";
             print "[", $signalnoneIDminus[$j], "] ";
             print "[", $signalFminus[$i], "-", $signalTminus[$i], "]\n";   # interval
          } else {  # pair cluster found
             # what if too far? LLL - might need to look at next exon
             if ($signalFminus[$i]-$signalnoneTplus[$pi]<25000) {
                print "(4) ", $GeneId2Name{$geneid}, " $geneid $chrom ", $signalnoneTplus[$pi], " ", $signalFminus[$i], " [$snum $snnum $runum $srunum] ";
                #print "[", $signalnoneIDplus[$pi], ";", $signalIDminus[$i], "]\n";
                print "[", $signalnoneIDplus[$pi], ";", $signalnoneIDminus[$j], "] ";
                print "[", $signalnoneFplus[$pi], "-", $signalnoneTplus[$pi], ",", $signalFminus[$i], "-", $signalTminus[$i],"]\n";  # interval
             } 
          }
      }
   }
}

sub filter { # snum snnum runum srunum
   my ($geneid,$snum,$snnum,$runum,$srunum) = @_;

   if ($isEnum) {
      print $GeneId2Name{$geneid}, " $snum $snnum $runum $srunum\n";
      return 0;
   }

   if (($snum>=$MIN_SIGNALS) && ($snnum>=$MIN_SIGNALNONE) && (!$snnum || ($snum>=$snnum*$MIN_S2SN)) && (!$srunum || ($snum>=$srunum*$MIN_S2SRU)) && (!$runum || ($snum>=$runum*$MIN_S2RU))) {
      return 1;
   }
   return 0;
}

sub overlap {
    my ($x,$y,$a,$b) = @_;

    return &min($y,$b)-&max($x,$a);
}

sub print_debug # $lastgene class ori chrom tmp_clusters_signal_plus
{
    my ($g,$c,$o,$chrom,$tmp) = @_;

    my @elems = split ',', $tmp;
    foreach my $x (@elems) {
       $x =~ /^(\d+):(\d+):(\d+)$/ or die "died. $x";
       print DBG  "$g\t", $GeneId2Name{$g}, "\t$c\t$o\t$chrom\t$1\t$2\t$3\n";
    }
}

sub min { ($_[0]<=$_[1]) ? $_[0] : $_[1]; }
sub max { ($_[0]>=$_[1]) ? $_[0] : $_[1]; }

close(DBG);

exit(0);

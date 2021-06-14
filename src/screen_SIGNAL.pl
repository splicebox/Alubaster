#!/usr/bin/perl
use strict;

#perl screen_SIGNAL.pl -s screen.unique < ERR188058.SIGNAL.screen.stats > ERR188058.SIGNAL.unique.screen.stats
($#ARGV==1 || $#ARGV==-1) or die "Usage: $0 [-s screenlist] < screenstats\n";

my $screenFile;
my %Screen;
if ($#ARGV==1) {
   my $tmparg = shift; ($tmparg eq "-s") or die "Usage: $0 [-s screenlist] < screenstats\n";
   $screenFile = shift;

   # note: the file lists the mate ids and desired orientation
   open(F, "<$screenFile") or die "died (1). $screenFile.\n";
   while (<F>) {
      #ERR188058.2207719/2:+ 
      /(\S+)$/ or die "died (2). $_";
      $Screen{$1} = 1;
   }
   close(F);
}
 
# Must ensure that the file is sorted by readid! It should be, from search_signal.pl.
#
my $MINLEN = 10;
my $MINOVL = 0.95;

my ($lastreadid,$lastreadlen);
my ($Lines,$LongestA,$LidxA,$LongestR,$LidxR,$count);

# read in the input stats file annotated with information about location
while (<>) {
   #>ERR188058.1350362/2:+ >ENSG00000237491.7:ENST00000585768.4;chr1:801607-801876 75 822 m1rev + 1 75 51 125 75 98 1-75 ArA 0-75-0 75-0-0 chr1:801657-801731;
   #>ERR188058.23964983/2:- >ENSG00000228794.7:ENST00000445118.5;chr1:827608-827775 75 553 m1fwd - 1 75 158 513 75 94 1-65,66-75 ArC 3-10-62 10-3-62 chr1:827765-827774;Alu:281-283;chr1:829003-829064;
   
   next if /^#/;

   /^(\S+) \S+ (\d+) \d+ (\S+) (\S) \d+ \d+ \d+ \d+ \d+ \d+ \S+ (\S+) (\d+)\-(\d+)\-(\d+) \S+ (\S+)$/ or die "died (3). $_"; 
   my ($readid,$readlen,$m1ori,$ori,$code,$R,$A,$C,$annot) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);

   if ($lastreadid ne $readid) {
      if (defined($lastreadid)) {
         # report
         if (!$count) {
            print "::Skipping $lastreadid (no_valid_aligns)\n";
         } elsif ($LongestA>=$MINOVL*$lastreadlen) {
            print "::Skipping $lastreadid (all_exon $LongestA)\n";
         } elsif ($LongestR>=$MINOVL*$lastreadlen) {
            print "::Skipping $lastreadid (all_repeat $LongestR)\n";
         } elsif ($LongestA<$MINLEN || $LongestR<$MINLEN) {
           print "::Skipping $lastreadid (not_enough_exon_repeat $LongestA $LongestR)\n";
         } elsif ($LidxA==0) {
           print "::Skipping $lastreadid (no_split_align $LongestA $LongestR)\n";
         } else {
           my @elems = split '\n', $Lines;
           print $elems[$LidxA], "$LongestA-$LongestR\n";
#          print "<<<", scalar(@elems), ">>>", $elems[$LidxA], "\n";
         }
      }

      # start a new one 
      $Lines = ""; $LidxA = $LidxR = -1; $LongestA = $LongestR = -1;
      $count = 0;
   }
      
   if (($m1ori eq "m1fwd") && ($ori eq "+")) { goto next_line; }
   if (($m1ori eq "m1rev") && ($ori eq "-")) { goto next_line; }

   if (defined($screenFile) && !defined($Screen{$readid})) { goto next_line; };

   $readid =~ /^\S+.(\d):(\S)$/ or die "died (4). $readid <<<$_>>>\n";
   my ($idx,$readori) = ($1,$2);
   if ($readori ne $ori) { goto next_line; };

   if ($LongestA<$A) { $LongestA = $A; }
   if ($LongestA<$C) { $LongestA = $C; }
   if ($LongestR<$R) { $LongestR = $R; } 

   if ((($A<$MINLEN) && ($C<$MINLEN)) || ($R<$MINLEN)) { goto next_line; }

   if ($LongestA==$A) { $LidxA = $count; }
   if ($LongestA==$C) { $LidxA = $count; }
   if ($LongestR==$R) { $LidxR = $count; }

   $Lines .= $_;

   $count++;

   next_line:
   $lastreadid = $readid;
   $lastreadlen = $readlen;
}
if (defined($lastreadid)) {
   if (!$count) {
      print "::Skipping $lastreadid (no_valid_aligns)\n";
   } elsif ($LongestA>=$MINOVL*$lastreadlen) {
      print "::Skipping $lastreadid (all_exon $LongestA)\n";
   } elsif ($LongestR>=$MINOVL*$lastreadlen) {
      print "::Skipping $lastreadid (all_repeat $LongestR)\n";
   } elsif ($LongestA<$MINLEN || $LongestR<$MINLEN) {
      print "::Skipping $lastreadid (not_enough_exon_repeat $LongestA $LongestR)\n";
   } elsif ($LidxA==0) {
      print "::Skipping $lastreadid (no_split_align $LongestA $LongestR)\n";
   } else {
      my @elems = split '\n', $Lines;
      print $elems[$LidxA], "\n";
#     print "<<<", scalar(@elems), ">>>", $elems[$LidxA], "\n";
  }
}

#!/usr/bin/perl
use strict;

# meant to be called from the working directory
#
($#ARGV==2) or die "Usage: SET krakenprefix unmappedpath\n";
my $SET = shift;
my $krakenPrefix = shift;
my $unmappedFile = shift;

my %Alu;

my $tmpFile;

# in the end, if we use the intersection between sim4db and kraken then
# there is no need/benefit to check kraken results (more inclusive), simply go
# with sim4db; but for the time bein gkeep both and see how many kraken
# reads remain unchanged (and therefore unmappable)
#
### read in the two (paired) lists of kraken files
if (defined($krakenPrefix)) {
   $tmpFile = $krakenPrefix . "_1.kraken_out";
   open(F, "<$tmpFile") or die "could not open kraken file $tmpFile.\n";
   while (<F>) {
     next if !/^C\s/; chomp;
     /^C\s(\S+)\s/ or die "died (kraken). $_";
     $Alu{$1} = 1;
   } close(F);

   $tmpFile = $krakenPrefix . "_2.kraken_out";
   open(F, "<$tmpFile") or die "could not open kraken file $tmpFile.\n";
   while (<F>) {
     next if !/^C\s/;
     chomp;
     /^C\s(\S+)\s/ or die "died (kraken). $_";
     if (defined($Alu{$1})) { $Alu{$1} = 3; } else { $Alu{$1} = 2; }
   }
   close(F);
}
#print STDERR "Finished reading Kraken Alu information...\n";


## Note: Must add position information in the ALU  sim4db stats file first (columns 3 and 4)
##
#Assume the following format (coordinates are in the original read sequence):
#tid ganame beg end gabeg gaend cov len relori pctcov pctid
#ERR127302.129 ALU 17 72 1 56 56 72 + 77.78 89
### read Alu match position from the sim4db ALu file; calculate the clear range
my %StartRange;
my %EndRange;

my $B = 10; # buffer; note: coordinates are in base 1 
foreach (my $i=1; $i<=2; $i++) {
  my $tmpFile = $SET . "_" . "$i.Alu.sim4db.coords.stats";
  open(F, "<$tmpFile") or die "could not open Alu sim4db file $tmpFile\n";
  while (<F>) {
    chomp;
    #ERR127302.129 ALU 17 72 1 56 56 72 + 77.78 89
    /^(\S+) \S+ (\d+) (\d+) \d+ \d+ \d+ (\d+) \S \S+ \d+$/ or die "died. $_";  
    my ($readid,$f,$t,$readlen) = ($1,$2,$3,$4);
    $readid .= ":$i";
    if ($f<$B) { $f = 1; }
    if ($t>$readlen-$B) { $t = $readlen; }
  
    if (!defined($StartRange{$readid})) {
       $StartRange{$readid} = 1; $EndRange{$readid} = $readlen;
    }

    next if ($StartRange{$readid}>$EndRange{$readid}); # this will get chopped entirely
    next if (($f!=1) && ($t!=$readlen));

    if (($f==1) && ($StartRange{$readid}<=$t)) {
       # chop at start
       $StartRange{$readid} = $t+1;
    }
    if (($t==$readlen) && ($EndRange{$readid}>=$f)) {
       # chop at end
       $EndRange{$readid} = $f-1;
    }
  }
  close(F);
}
# print STDERR "Finished reading Alu sim4db information and clear-range calculation...\n";

### Sweep 1: read in an unmapped file and first record the read ids
#   (could do this in one step, remembering the full read information,
#   but not sure how large the hash will be in memory).
#   
#   Further restrict these to the reads with kraken matches.
#   SRR1284895.22223173        69      *       0       0       *       *       0       0       CTAGTTTGTAAATGTGAACCCCCGGTGCCCTTCTTGGAGGAGGGGGCCAGGAGCAGCTTCAGGTAGGAAGGCACTGTGGCAGCCTCTTCAGGGATGANNN    ;;<44:AD?CFA>F,AEGIIA18@1?CFB*0?GD?BF0?/0-68BE'5;?@5=52;=?:>:A@5(:<::1((8?3>@A@@BBB8?9(4>:4:?#######

my %Unmapped;

open(F, "-|", "samtools view $unmappedFile") or die "could not open pipe from samtools view"; 
#open(F, "<UNMAPPED") or die "UNMAPPED";
while (<F>) {
   chomp;

   /^(\S+)\t(\d+)\t/ or die "died. $_";
   my ($readid,$flag) = ($1,$2);
   next if !defined($Alu{$readid});
   my $idx = ($flag & 0x40) ? 1 : 2;

   $Unmapped{$readid} += $idx; # unmapped reads are only listed once
}
close(F);

#print STDERR "Finished reading unmapped file (1st time). Start writing...\n";

# now read in the read information, including sequence
my %Reads1;
my %Reads2;


my %OrigReads1;
my %OrigReads2;


my $MINL = 25;
open(F, "-|", "samtools view $unmappedFile") or die "could not open pipe from samtools view"; 
#open(F, "<UNMAPPED") or die "UNMAPPED";
while (<F>) {
   chomp;
    
   # SRR1284895.67553584       69      *       0       0       *       *       0       0       ATAATATGCACAGCAGCCGTGAGGCACACAATGCCTTCAATATCATTAGAAGCCACACATGCCTTTAATTCATCCAAAAATGACCTGACCNNNNNNNNNN    @CBDFFFFHHHHHIJIIJJIJJGJIHJIIIIIGJJJJIJIBFGDHGIIJIJIIJJIJJIIIJ=CEHFEEE=AHBE@DEDDEEEDDDDCDCDDDDDDEDDB    ZT:A:N

   # SRR1284895.2      69      *       0       0       *       *       0       0       NTAAGATAGAGGAGACACCTGCTAGGTGTAAGGAGAAGANNNNTAGGTCTACGGAGGCTCCAGGNTGGGAGTAGTTCCCTGCTAAGGGAGGGTAGACNNN    #1:DDAEFHHFHHIIIIIIIIIEHHI<BFFGIIHGICHB####00?DBGHIIGA@GIIDCHG>D#(,=BA?;?@C@ACDCCDDCCDDAAB@@2>?CCC<C
   my ($readid,$flag,$seq,$quals);
   if (/^(\S+)\t(\d+)\t\S\t\S\t\S\t\S\t\S\t\S\t\S\t(\S+)\t(\S+)\t/) {
      ($readid,$flag,$seq,$quals) = ($1,$2,$3,$4);
   } elsif (/^(\S+)\t(\d+)\t\S\t\S\t\S\t\S\t\S\t\S\t\S\t(\S+)\t(\S+)$/) {
      ($readid,$flag,$seq,$quals) = ($1,$2,$3,$4);
   } else {
      die "died. $_";
   }
   next if ($Unmapped{$readid}!=3); 

   my $idx = ($flag & 0x40) ? 1 : 2;

   # Before any trimming, record th eoriginal reads
   if ($idx==1) { $OrigReads1{$readid} = "\@" . "$readid\n$seq\n+\n$quals\n"; }
   else { $OrigReads2{$readid} = "\@" . "$readid\n$seq\n+\n$quals\n"; }
   
   # now edit the sequence to remove the Alu matching portions
   if (defined($StartRange{"$readid:$idx"}) && (($StartRange{"$readid:$idx"}!=1) || ($EndRange{"$readid:$idx"}!=length($seq)))) {
      # need to modify the sequence and quals string
      if ($EndRange{"$readid:$idx"} - $StartRange{"$readid:$idx"} < $MINL) {
         $seq = "NNNNNNNNNN"; $quals = "IIIIIIIIII";
      } else {
         $seq = substr($seq,$StartRange{"$readid:$idx"}-1,$EndRange{"$readid:$idx"}-$StartRange{"$readid:$idx"}+1);
         $quals = substr($quals,$StartRange{"$readid:$idx"}-1,$EndRange{"$readid:$idx"}-$StartRange{"$readid:$idx"}+1);
      }
   } # otherwise no change to seq, quals

   if ($idx==1) { $Reads1{$readid} = "\@" . "$readid\n$seq\n+\n$quals\n"; }
   else { $Reads2{$readid} = "\@" . "$readid\n$seq\n+\n$quals\n"; }

}
close(F);

$tmpFile = $SET . "_unmapped_1.NEW.fastq";
open(R1, ">$tmpFile") or die "could not open $tmpFile for writing.\n";
$tmpFile = $SET . "_unmapped_2.NEW.fastq";
open(R2, ">$tmpFile") or die "could not open $tmpFile for writing.\n";

$tmpFile = $SET . "_unmapped_1.orig.fastq";
open(OR1, ">$tmpFile") or die "could not open $tmpFile for writing.\n";
$tmpFile = $SET . "_unmapped_2.orig.fastq";
open(OR2, ">$tmpFile") or die "could not open $tmpFile for writing.\n";


foreach my $k (keys %Reads1) {
   next if (($EndRange{"$k:1"}-$StartRange{"$k:1"}+1<$MINL) || ($EndRange{"$k:1"}-$StartRange{"$k:1"}+1<$MINL));
   print R1 $Reads1{$k};
   print R2 $Reads2{$k};
   print OR1 $OrigReads1{$k};
   print OR2 $OrigReads2{$k};
}
close(R1);
close(R2);

close(OR1);
close(OR2);

print STDERR " get_unmapped_readpairs --  done!\n";

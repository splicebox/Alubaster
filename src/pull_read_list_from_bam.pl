#!/usr/bin/perl
use strict;

my $Usage = "Usage: $0 [options] fastq_1 fastq_2 [--gzip] < samfile\n" .
            "           --gzip   input fastq files are compressed\n" .
            "           -o       output nonconcordant mates fastq file (required)\n"  .
            "           -d       output nonconcordant reads fastq file\n"; 
                   
my ($fqFile1,$fqFile2);
my ($outmateFile,$outncFile,$gzipFlag);

my $tmparg = shift;
while ($tmparg) {
  if ($tmparg eq "-o") {
     $outmateFile = shift;
  } elsif ($tmparg eq "-d") {
     $outncFile = shift;
  } elsif ($tmparg eq "--gzip") {
     $gzipFlag = 1;
  } else {
     $fqFile1 = $tmparg;
     $fqFile2 = shift;
  }
  $tmparg = shift;
}

(defined($outmateFile) && defined($fqFile1) && defined($fqFile2)) or die $Usage;

open(OM, ">$outmateFile") or die "could not open output mate file $outmateFile for writing.\n";
!defined($outncFile) or open(NC, ">$outncFile") or die "could not open output nonconcordant reads file $outncFile for writing.\n";


my %Seen;

# read in a nonconcordant sam file

while (<>) {
   /^(\S+)\t(\d+)\t/ or die "died (1). $_";
   my ($readid,$flag) = ($1,$2);
   my $idx = ($flag & 0x40) ? 1 : 2;
   $Seen{"$readid:$idx"} = 1;
}

my $fq;

# ../ERR188081_1.fastq : idx==1
my $idx = 1;
if (defined($gzipFlag)) {
   open (F1, "zcat $fqFile1 |") or die "$fqFile1.\n"; }
else {
   open (F1,"<$fqFile1") or die "$fqFile1.\n";
}
while (<F1>) {
   my ($readid,$rest) = ($1,$2);
   if (/^@(\S+)( .*)/) {
      ($readid,$rest) = ($1,$2);
   } elsif (/^@(\S+)$/) {
      ($readid,$rest) = ($1,$2);
   } else {
      die "died (2). $_";
   }

   $fq = "\@$1/" . $idx . "$rest\n";
   $_  = <F1>; $fq.= $_;
   $_  = <F1>; $fq.= "+\n";
   $_  = <F1>; $fq.= $_;
   my  $otheridx = (3-$idx);
   if (defined($Seen{"$readid:$otheridx"})) { print OM $fq; }
   if (defined($outncFile) && defined($Seen{"$readid:$idx"})) {
      print NC $fq;
   } 
}
close(F1);

$idx = 2;
if (defined($gzipFlag)) {
   open (F2, "zcat $fqFile2 |") or die "$fqFile2.\n"; }
else {
   open (F2,"<$fqFile2") or die "$fqFile2.\n";
}
while (<F2>) {
   my ($readid,$rest) = ($1,$2);
   if (/^@(\S+)( .*)/) { 
      ($readid,$rest) = ($1,$2);
   } elsif (/^@(\S+)$/) {
      ($readid,$rest) = ($1,$2);
   } else {
      die "died (3). $_";
   }

   $fq= "\@$1/" . $idx . "$rest\n";
   $_ = <F2>; $fq.= $_;
   $_ = <F2>; $fq.= "+\n";
   $_ = <F2>; $fq.= $_;
   my $otheridx = (3-$idx);
   if (defined($Seen{"$readid:$otheridx"})) { print OM $fq; }
   if (defined($outncFile) && defined($Seen{"$readid:$idx"})) {
      print NC $fq;  
   }
}
close(F2);

close(OM);
!defined($outncFile) or close(NC);

exit(0);

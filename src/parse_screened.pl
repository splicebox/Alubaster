#!/usr/bin/perl
use strict;

my $SIM4DB = $ENV{SIM4DB};
my $SCRIPTS = $ENV{SCRIPTS};


my $Usage = "Usage: $0 set genome genomehdrs alufasta matesfasta\n";
($#ARGV==4) or die $Usage;

my $SET = shift;
my $GENOME = shift;
my $GENOMEHDRS = shift;
my $ALU = shift;
my $MATESFA = shift;

my %GAid;
my $count = 0;
open (F, "<$GENOMEHDRS") or die "Could not open $GENOMEHDRS";
while  (<F>) {
   chomp;
   # >chr1
   if (/^>(\S+)$/) { $GAid{$1} = $count; } # print "(1) $1\n"; }
   elsif (/^>(\S+)\s/) { $GAid{$1} = $count; } # print "(2) $1\n";}
   elsif (/^(\S+)$/) { $GAid{$1} = $count; } # print "(3) $1\n"; }
   elsif (/^(\S+)\s/) { $GAid{$1} = $count; } # print "(4) $1\n"; }

   $count++;
}
close(F);

#>ERR188058.23232999/1:- >USP48:ENSG00000090686.14:ENST00000526044.4;chr1:21729704-21729829 75 588 m1fwd - 9 75 91 176 67 85 9-41,42-75 ARC 30-33-0 33-30-0 chr1:21729794-21729826;Alu:21-50;50-47
#>ERR188058.21978955/2:- >ZMYM6:ENSG00000163867.15:ENST00000460607.1;chr1:35017250-35019602 75 4988 m1fwd - 1 75 984 2464 75 92 1-47,48-75 ARA 47-28-0 28-47-0 chr1:35018233-35018260;Alu:65-111;33-47
#>ERR188058.25722012/2:- >NFYC:ENSG00000066136.18:ENST00000372652.4;chr1:40709357-40709616 75 655 m1fwd - 1 75 88 578 75 89 1-42,43-75 rAC 33-42-0 33-42-0 Alu:88-120;chr1:40709498-40709539;42-44
#>ERR188058.18404951/1:+ >RBMXL1:ENSG00000213516.8:ENST00000321792.5;chr1:88979456-88984066 75 4893 m1rev + 1 75 248 2542 75 86 1-29,30-40,41-75 CRA 31-49-0 0-31-49 Alu:248-278;chr1:88981176-88981187;chr1:88981679-88981715;50-31
#>ERR188058.10126975/1:- >DAP3:ENSG00000132676.14:ENST00000461479.1;chr1:155688839-155689337 75 781 m1fwd - 1 75 319 723 75 82 1-55,56-75 ArC 53-18-0 18-53-0 chr1:155689157-155689174;Alu:172-224;55-62

while (<>) {
   /^(\S+) (\S+) \d+ \d+ \S+ (\S) \d+ \d+ \d+ \d+ \d+ \d+ \S+ (\S+) \S+ \S+ (\S+)$/ or die "died (1). $_";  
   my ($readid,$info,$relori,$code,$tempaln) = ($1,$2,$3,$4,$5);

   $readid =~ /^(\S+):\S$/ or die "died (2). $readid";
   `echo $1 > tmp.id`;
   `${SIM4DB}/leaff -F $MATESFA -q tmp.id > tmp.est`;
   my $est = `tail -1 tmp.est`; chomp $est;

   my ($left,$alu,$right);

   my @elems = split ';', $tempaln;
   while (@elems) {
     my $x = shift @elems;
     if ($x=~/^chr/) {
       if (!defined($alu)) { # left
          if (defined($left)) { $left .= ";$x"; } else { $left = $x; }
       } else { # right
          if (defined($right)) { $right .= ";$x"; } else { $right = $x; }
       }
     } else {
       $alu = 1;
     }
   }
#  if ($tempaln =~ /^(\S+)Alu:\d+\-\d+;(\S+;)\d+\-\d+$/) { 
#     ($left,$right) = ($1,$2); 
#  } elsif ($tempaln =~ /^Alu:\d+\-\d+;(\S+;)\d+\-\d+$/) {
#     $right = $1;
#  } elsif ($tempaln =~ /^(\S+)Alu:\d+\-\d+;\d+\-\d+$/) {  
#     $left = $1;
#  } else {
#     die "died (3). $tempaln\n";
#  }
   
   print "$readid $est $info $code";
   my $scr;

   # left side
   if (defined($left)) {
     my @ints = split ';', $left;
     my ($chr,$b,$e);
     foreach my $x (@ints) {
        $x =~ /^(\S+):(\d+)\-(\d+)$/ or die "died (4). $x";
        $chr = $1;
        defined($b) or $b = $2;
        $e = $3;

        # generate a script; genomic coords are 0-based 
        # -f -e 61728 -D 0 2370482 2375224
        # -r -e 61730 -D 0 6723331 6757701
        my ($B,$E) = ($b-50-1,$e+50);
        $scr = ($relori eq "+") ? "-f" : "-r";
        $scr .= " -e 0";
        $scr .= " -D " . $GAid{$chr};
        $scr .= " $B $E";
     
        `echo \"$scr\" > tmp.scr`;
        my $aln = `${SIM4DB}/sim4db -cdna tmp.est -gen $GENOME -scr tmp.scr -o - `; #| perl ${SCRIPTS}/EM2cov_coords.pl;
        ##### Next process this into something printable 
        my $str = &process_aln($chr,$aln); print " $str";

        #print "$readid <<<< $aln\n";
     }
   }
   #Alu
   $scr = ($relori eq "+") ? (($code=~/R/) ? "-f" : "-r") : (($code=~/R/) ? "-r" : "-f"); 
   $scr .= " -e 0 -D 0 0 312";
   `echo \"$scr\" > tmp.scr`;
   my $aln = `${SIM4DB}/sim4db -cdna tmp.est -gen $ALU -scr tmp.scr -o -`; #| perl ${SCRIPTS}/EM2cov_coords.pl;
   ##### Next process this into something printable 
   my $str = &process_aln("ALU",$aln); print " $str";

   #print "$readid :::: $aln\n";
   
   # right side
   if (defined($right)) {
     #print "<<<$right>>>\n";
     my @ints = split ';', $right;
     my ($chr,$b,$e);
     foreach my $x (@ints) {
        $x =~ /^(\S+):(\d+)\-(\d+)$/ or die "died (5). $x";
        $chr = $1;
        defined($b) or $b = $2;
        $e = $3;
     
        # generate a script; genomic coords are 0-based 
        my ($B,$E) = ($b-50-1,$e+50);
        $scr = ($relori eq "+") ? "-f" : "-r";
        $scr .= " -e 0";
        $scr .= " -D " . $GAid{$chr};
        $scr .= " $B $E";

        `echo \"$scr\" > tmp.scr`;
        my $aln = `${SIM4DB}/sim4db -cdna tmp.est -gen $GENOME -scr tmp.scr -o - `; #| perl ${SCRIPTS}/EM2cov_coords.pl;
        ##### Next process this into something printable
        my $str = &process_aln($chr,$aln); print " $str";
        #print "$readid $aln>>>>\n";
     }
   }
   print "\n";
}

sub process_aln { # aln
  my ($ganame, $l) = @_;

  my @lines = split '\n', $l;

  my $str;
  my @genX = ();
  my @cdnaX = ();
  my ($len,$gaoffset,$pctid,$relori);

  while (@lines) {
    $_ = shift @lines; chomp;
    /^sim4begin/ or die "begin (6).$_";

    # header line
    $_ = shift @lines; chomp;
    /^\d+\[(\d+)\-\d+\-\d+\] \d+\[(\d+)\-\d+\] <\d+\-\d+\-(\d+)\-(\S+)\-/ or die "died (header). $_";
    ($len,$gaoffset,$pctid,$relori) = ($1,$2,$3,$4);
    $relori = ($relori eq "complement") ? "-" : "+";

    $gaoffset = 0;   ### HERE: Because the new versions of sim4db print the absolute values, not relative to offset

    # edef line - skip, we already know readid
    $_ = shift @lines; chomp;
    while (/^\s*$/) { $_ = shift @lines; chomp; }
    /^edef/ or die "died (edef). $_";

    #ddef line - skip, ewe already know chromosome/ALU
    $_ = shift @lines; chomp;
    while (/^\s*$/) { $_ = shift @lines; chomp; }
    /^ddef/ or die "died (ddef). $_";

    # start coordinate lines
    while (@lines) {
      $_ = shift @lines; chomp;
      last if /^sim4end/;
      next if /^[\-actgnACGTN]/;
      /^(\d+)\-(\d+) \((\d+)\-(\d+)\) / or die "died (cDNA). $_";
      
      my ($a,$b,$c,$d) = ($1,$2,$3,$4);
      $str = "";
      if ($relori eq "-") { $str = ($len-$b+1) . "-" . ($len-$a+1); push @cdnaX, $str; }
      else { push @cdnaX, "$a-$b"; }
      $str = ($gaoffset+$c) . "-" . ($gaoffset+$d);
      push @genX, $str;
    }
  }
  $str = join ',', @cdnaX;
  $str .= " $ganame:";
  $str .= join ',', @genX;
  $str .= " $relori $pctid";
}

exit(0);

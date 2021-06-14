#!/usr/bin/perl
use strict;

# sort by geneid, then by start (end coordinate)

## REGION:
##edef ddef elen glen m1ori ori efrom eto dfrom dto cov pctid exons type readsRA reads-by-type
#>SRR1284895.44160104/1:- >ENSG00000280279.1:ENST00000623180.1;chr17:65830-65887 100 5637 m1fwd - 1 92 5525 5618 92 90 1-92 AG 94-0 0-94
## UNSPLICED:
##edef ddef elen glen m1ori ori efrom eto dfrom dto cov pctid exons type readsU-A reads-by-type
#>SRR1284895.44160104/1:- >ENSG00000280279.1:ENST00000623180.1;chr17:65830-65887 100 558 m1fwd - 76 100 41 64 25 84 76-100 AU 6-18 18-6
## SIGNAL:
##edef ddef elen glen m1ori ori efrom eto dfrom dto cov pctid exons type readsR-A-C reads-by-type
#>SRR1284895.94379692/1:- >ENSG00000280279.1:ENST00000623180.1;chr17:65830-65887 100 532 m1fwd + 79 100 139 158 22 73 79-100 ARC 20-0-0 0-20-0

my %Lines;
my %GeneID;
my %ReadID;
my %From;
my %To;


while (<>) {
    if (/^#/) { print $_; next; }

    # >SRR1284895.80755432/1:+ >ENSG00000280071.2:ENST00000620015.3;chr21:5116341-5117231 100 1331 m1fwd + 1 89 1035 1121 89 79 1-89 ARC 87-0-0 0-87-0
    # >SRR1284895.44160104/1:- >ENSG00000280279.1:ENST00000623180.1;chr17:65830-65887 100 5637 m1fwd - 1 92 5525 5618 92 90 1-92 AG 94-0 0-94
    #/^>(\S+.[12]):\S >(\S+):\S+;\S+:(\d+)\-(\d+) / or die "died. $_";
    /^(\S+.[12]):\S (\S+):\S+;\S+:(\d+)\-(\d+) / or die "died. $_";
    my ($readid,$geneid,$from,$to) = ($1,$2,$3,$4);
    my $key = "$readid:$geneid:$from:$to";
    $Lines{$key} .= $_;  
    $GeneID{$key} = $geneid;
    $ReadID{$key} = $readid;
    $From{$key} = $from;
    $To{$key} = $to;
}


foreach my $k (sort mysort_readid keys %Lines) {
    print $Lines{$k};
}

sub mysort #
{
    if ($GeneID{$a} ne $GeneID{$b}) { return ($GeneID{$a} cmp $GeneID{$b}); }
    if ($From{$a}!=$From{$b}) { return ($From{$a} <=> $From{$b}); }
    if ($To{$a}!=$To{$b}) { return ($To{$a} <=> $To{$b}); }
    
    return ($ReadID{$a} cmp $ReadID{$b});
}

# alternative procedure - sort by geneid, then by readid
sub mysort_readid #
{
    if ($GeneID{$a} ne $GeneID{$b}) { return ($GeneID{$a} cmp $GeneID{$b}); }
    if ($ReadID{$a} ne $ReadID{$b}) { return ($ReadID{$a} cmp $ReadID{$b}); }

    if ($From{$a}!=$From{$b}) { return ($From{$a} <=> $From{$b}); }
    if ($To{$a}!=$To{$b}) { return ($To{$a} <=> $To{$b}); } 
}

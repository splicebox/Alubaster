#!/usr/bin/perl
use strict;

while (<>) {
   /^@(.+)$/ or die "died. $_";
   print ">$1\n";
   $_ = <>; print $_;
   $_ = <>; $_ = <>;
}

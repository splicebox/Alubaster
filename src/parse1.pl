#!/usr/bin/perl
use strict;

while (<>) {
   next if /^\d?\d?\d?$/;
   next if /^.$/;
   next if /^..\)/;
   next if /^>/;
   print $_;
}

while (<>) {
   next if /^$/;
   
   if (/sim4begin/) {
       print $_;
       $_ = <>; chomp;
       while (/^$/ || /^.$/) { $_ = <>; chomp; }
       /^\d+\[/ or die "expected header line. $_";
       print $_, "\n";

       $_ = <>; chomp; while (!/edef/) { $_ = <>; chomp; } 
       /^edef/ or die "expected edef line. $_";
       print $_, "\n";

       $_ = <>; chomp; while (!/ddef/) { $_ = <>; chomp; }
       /^ddef/ or die "expected ddef line. $_";
       print $_, "\n";

       $_ = <>; chomp; while (!/^\d+\-/) { $_ = <>; chomp; }
       /^\d+\-/ or die "expected alignment line. $_";
       print $_, "\n";

       $_ = <>; while (!/sim4end/) { if (!/^$/) { print $_; } $_ = <>; }
       print $_;
   }
}

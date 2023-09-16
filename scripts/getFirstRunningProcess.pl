#!/usr/bin/perl -w

use strict;
use Getopt::Std qw(getopt);

my ($topCmd) = "top -n 1 -b";
my ($start) = 0;
my (%count, $i);
getopt('c',\%count);
$count{c} = 1 if (! exists($count{c}));
open TOPCMD, "$topCmd |" or die "ERROR: Cannot fork: $!\n";
while (<TOPCMD>) {
    chomp;
    if ($_ =~ /^\s+PID\s+USER/) {
	$start = 1;
    }elsif ($start) {
	$i++;
	print "$_\n";
	last if ($i == $count{c});
    }
}
close TOPCMD;

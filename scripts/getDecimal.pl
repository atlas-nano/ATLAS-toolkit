#!/usr/bin/perl -w
use strict;
die "usage: $0 text\n" if (! @ARGV);

my ($intext) = $ARGV[0];
if ($intext =~ /(\d+\.?\d*)/) {
    print "$1\n";
} else {
    print "\n";
}

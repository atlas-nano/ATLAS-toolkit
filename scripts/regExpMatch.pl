#!/usr/bin/perl -w

use strict;
die "0\n" if ($#ARGV < 1);

my ($str, $match) = ($ARGV[0], $ARGV[1]);

if ($str =~ /$match/) {
    print "1\n";
} else {
    print "0\n";
}

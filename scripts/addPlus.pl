#!/usr/bin/perl -w
use strict;
die "usage: $0 number\n" if (! @ARGV);
my $num = $ARGV[0];
die "ERROR: Expected decimal for number, got \"$num\"!\n" if ($num !~ /^(\-?\d+\.?\d*)/);
if ($1 > 0) {
    print "+${1}\n";
} else {
    print "$1\n";
}

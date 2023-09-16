#!/usr/bin/perl

use Cwd 'abs_path';

die "usage: $0 relativepath\n" if (! @ARGV);
my $root = $ARGV[0];

$root = abs_path($root);

print "$root\n";

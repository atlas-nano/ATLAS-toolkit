#!/usr/bin/perl -w
use strict;

die "usage: $0 bgffile [outfile]\n" if (! @ARGV);

my ($infile, $outfile) = @ARGV;
my ($i, @outarray, $outline);

if (! $outfile) {
    $outfile = $infile;
}

die "Cannot locate $infile: $!\n" if (! -e $infile);

open INFILE, $infile or die "Cannot open $infile: $!\n";
while (<INFILE>) {
    chomp;
    $outline = $_;
    if ($_ =~ /(ATOM\s+\d+\s+\w+\s+)(\w+)\s\s(\s+\d+.+)$/) {
	$outline = $1 . $2 . " A" . $3;
    }
    push @outarray, $outline;
}
close INFILE;

open OUTFILE, "> $outfile" or die "Cannot write to $outfile: $!\n";
for $i (@outarray) {
    print OUTFILE "$i\n";
}

#!/usr/bin/perl -w

if (!@ARGV) {
    die "usage: add.pl datafile\n";
}

-e $ARGV[0] or die "Cannot open $ARGV[0]:$!\n";

open INFILE, $ARGV[0];
$instr = "";
$my_sum = 0.000000;
while (<INFILE>) {
    $instr = $_;
    chomp($instr);
    if ($instr =~ /^\s+\d+\s+(\-?\d+\.\d+)/) {
	$my_sum += $1;
    }
}

close INFILE;

print "Total: $my_sum\n";

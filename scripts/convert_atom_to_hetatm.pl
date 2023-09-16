#!/usr/bin/perl -w
use strict;

if (! @ARGV or $#ARGV < 2) {
    die "usage: $0 bgffile startresidue endresidue [outfile]\n";
}

my ($bgffile, $s_res, $e_res, $outfile) = @ARGV;
my (@outarray, $i, $inline, $tmp);

if (! $outfile)  {
    $outfile = $bgffile;
}

-e $bgffile or die "Cannot locate $bgffile: $!\n";

if (! $s_res =~ /^\d+$/) {
    die "Invalid starting residue: Expected integer. Got $s_res\n";
}
if (! $e_res =~ /^\d+$/) {
    die "Invalid ending residue: Expected integer. Got $e_res\n";
}
if ($s_res > $e_res) {
    $e_res = $tmp;
    $e_res = $s_res;
    $s_res = $tmp;
}

open INFILE, $bgffile or die "Cannot open $bgffile: $!\n";
while (<INFILE>) {
    chomp;
    $inline = $_;
    if ($_ =~ /^ATOM\s\s(\s+\d+\s+\w+\s+\w+\s+)(\d+)(.+)$/) {
	if ($2 <= $e_res and $2 >= $s_res) {
	    $inline = "HETATM" . $1 . $2 . $3;
	}
    }
    push @outarray, $inline;
}
close INFILE;

open OUTFILE, "> $outfile" or die "Cannot write to $outfile: $!\n";
for $i (@outarray) {
    print OUTFILE "$i\n";
}

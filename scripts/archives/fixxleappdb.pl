#!/usr/bin/perl -w
# fixes the names of the bases in the curve pdb file
# to make it compatible with xleap
# e.g changes G to DG, T to  DT etc.

# Variable Declaration Section
use strict;

my ($i, @outputarry, $flenm, $in_text, $base_nm);

if (!@ARGV) {
    die "usage: fixpdb.pl pdbfile.pdb\n";
}

$flenm = $ARGV[0];

-e $flenm or die "Cannot locate $flenm: $!\n";

print "Fixing Curve file $flenm to xleap format...";

open PDBFILE, $flenm or die "Cannot open $flenm, $!\n";

while (<PDBFILE>) {
    chomp;
    if ($_ =~ /^(ATOM\s+\d+\s+\w+\'?\s+)(\w+\s)(.+)$/) {
	$base_nm = $2;
	chop $base_nm;
	$in_text = $1 . "D" . $base_nm . $3;
    } else {
	$in_text = $_;
    }
    push @outputarry, $in_text;

}

close PDBFILE;

# Writes the changes
open OUTFILE, "> " . $flenm or die "Cannot write to $flenm. $!\n";

for $i (@outputarry) {
    print OUTFILE "$i\n";
}

print "Done\n";


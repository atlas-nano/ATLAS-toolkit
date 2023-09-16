#! /usr/bin/perl -w
# fixpdb.pl - fixes the names of the bases in the namot2 pdb file
# to make it compatible with biograf and MSCFF4.2
# e.g changes GUA to G, THY to  T etc.

# Variable Declaration Section

use strict;


if (!@ARGV) {
    die "usage: fixpdb.pl pdbfile [savename]\n";
}

my (@outputarry, $flenm, $in_text, $i, $saveName);

# Opens the file and gets the input

($flenm, $saveName) = @ARGV;

if (! -e $flenm || ! -r $flenm || ! -T $flenm) {
    die "ERROR: Cannot access $flenm: $!\n";
}

$saveName = $flenm if (! $saveName);

$|++;
print "Parsing $flenm...";
open PDBFILE, $flenm or die "Cannot open $flenm, $!\n";
while (<PDBFILE>) {
    $in_text = $_;
    $in_text =~ s/HETATM/ATOM  /g;
    if ($in_text =~ /O3T|O5T/i) {
    } elsif ($in_text =~ /^ATOM\s+\d+\s+[O|C|N].+[GUA|THY|ADE|CYT]/i) {
	push @outputarry, uc($in_text);
    } else {
	if ($in_text =~ /HE|WAT|Na\+|Cl\-/i) {
	    push @outputarry, "TER\n";
	}
    }
}
close PDBFILE;
print "Done\n";

print "Creating $saveName...";
open OUTFILE, "> $saveName" or die "ERROR: Cannot create file $saveName: $!\n";
for $i (@outputarry) {
    print OUTFILE $i;
}
close OUTFILE;

print "Done\n";

#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use FileFormats qw(GetPDBFileInfo);
use Getopt::Std qw(getopt);

sub init;
sub printSequence;
sub numerically { ($a<=>$b); }

my ($pdbFile, $saveName, $ATOMS);

$|++;
&init;
print "Parsing PDB file $pdbFile...";
($ATOMS, undef, undef) = GetPDBFileInfo($pdbFile);
print "Done\n";
&printSequence($ATOMS, $saveName);

sub init {
    my (%OPTS);
    getopt('ps',\%OPTS);
    
    die "usage: $0 -p pdb file -s [save name]\n" if (! defined($OPTS{p}));

    ($pdbFile, $saveName) = ($OPTS{p}, $OPTS{s});

    print "Initializing...";
    FileTester($pdbFile);
    print "Done\n";
}

sub printSequence {
    my ($atoms, $save) = @_;
    my ($OUTPUT, $i, $res, $resid);
    
    if (defined($save)) {
	print "Printing sequence to $save...";
	open $OUTPUT, "> $save" or die "ERROR: Cannot create $save: $!\n";
    } else {
	$OUTPUT = \*STDOUT;
    }

    $resid = 0;
    for $i (sort numerically keys %{ $atoms }) {
	next if ($atoms->{$i}{RESNUM} == $resid);
	$res = $atoms->{$i}{RESNAME};
	$resid = $atoms->{$i}{RESNUM};
	$res =~ s/[\d+|D]//g;
	printf $OUTPUT "%-4d $res\n", $resid;
    }
    
    if (defined($save)) {
	close $OUTPUT;
	print STDOUT "Done\n";
    }

}

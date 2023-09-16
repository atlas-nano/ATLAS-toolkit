#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use MolData;
use Getopt::Std qw(getopt);
use LAMMPS qw(ParseLAMMPSTrj CreateLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType);

sub init;

my ($molStructure, $fileName, $fileType, $readFunc, $saveName);
my ($start, $end);

$|++;
&init;

print "Gettting data from $fileType file $fileName...";
$molStructure->read($fileName, $fileType);
print "Done\nDeleting atoms...";
my $currMol = $molStructure->molecule(2);
$molStructure->deleteMol($currMol);
print "Done\nWriting to $fileType file $saveName...";
$molStructure->write($saveName, $fileType);
print "Done\n";



sub init {
    my (%OPTS);

    getopt('ftsl', \%OPTS);
    die "usage: $0 -f file name -t (file type)\n" if (! exists($OPTS{f}));

    print "Initialzing...";
    $molStructure =  MolData->new();
   ($fileName, $fileType, $saveName) = ($OPTS{f}, $OPTS{t}, $OPTS{s});
    $molStructure->testFile($fileName);
    if (! defined($fileType)) {
	$fileType = "bgf";
	if ($fileName =~ /\.(\w+)$/) {
	    $fileType = lc $1;
	}
    }

    $saveName = $molStructure->getFileName($fileName) if (! defined($saveName));
    print "Done\n";
}

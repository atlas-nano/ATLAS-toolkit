#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use MolData;
use Getopt::Std qw(getopt);

sub init;
sub stressCell;

my ($fileName, $fileType, $saveName);
my ($molStructure, $readFunc, $DIMS);

$|++;
&init;
print "Gettting data from $fileType file $fileName...";
$molStructure->read($fileName, $fileType);
$molStructure->moveMol("all", "com_center");
&stressCell($molStructure, $DIMS);
print "Done\nWriting $fileType file $saveName...";
$molStructure->write($saveName, $fileType);
print "Done\n";

sub stressCell {
    my ($mols, $dimOpts) = @_;
    my ($i);

    for $i (keys %{ $dimOpts }) {
	$mols->stressCell($i, $dimOpts->{$i});
    }
}

sub init {
    my (%OPTS, $dimStr);

    getopt('ftds', \%OPTS);
    die "usage: $0 -f file name -d \"cell length scaling x: xx.xx y: xx.xx ...\" -t (file type) -s (savename)\n" 
	if (! exists($OPTS{f}) or ! exists($OPTS{d}));

    print "Initialzing...";
    $molStructure =  MolData->new();
   ($fileName, $fileType, $dimStr, $saveName) = ($OPTS{f}, $OPTS{t}, $OPTS{d}, $OPTS{s});
    while ($dimStr =~ /(x|y|z):\s*(\d+\.?\d*)/gi) {
	$DIMS->{$1} = $2;
    }
    die "ERROR: Invalid cell length information found. Expected [x|y|z]: x.xx. Got \"$dimStr\"\n"
	if (! $DIMS);
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

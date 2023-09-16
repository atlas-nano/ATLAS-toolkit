#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use MolData;
use Getopt::Std qw(getopt);

sub init;
sub addBoxToStructure;

my ($fileName, $fileType, $readFunc, $saveName, $BOX);
my ($molStructure, $start, $end);

$|++;
$molStructure =  MolData->new();
&init;
print "Gettting data from $fileType file $fileName...";
$molStructure->read($fileName, $fileType);
print "Done\nAdding box...";
&addBoxToStructure($molStructure, $BOX);
print "Done\nWriting to $fileType file $saveName...";
$molStructure->write($saveName, $fileType);
print "Done\n";


sub addBoxToStructure {
    my ($struct, $box) = @_;
    my ($i, $cellLength);

    for $i ("a", "b", "c") {
	$cellLength = $struct->vdwbox->{$i}{max} - $struct->vdwbox->{$i}{min};
	if (! defined($box)) {
	    if (! exists($struct->cell->{$i})) {
		$struct->cell->{$i} = $cellLength;
		$struct->cell->{alpha} = 90;
		$struct->cell->{beta} = 90;
		$struct->cell->{gamma} = 90;
		$struct->cell->{image}{$i}{p} = 1;
		$struct->cell->{image}{$i}{n} = -1;
	    }
	} elsif ($box->{$i} > $cellLength) {
	    $struct->cell->{$i} = $box->{$i};
	    $struct->cell->{alpha} = $box->{alpha};
	    $struct->cell->{beta} = $box->{beta};
	    $struct->cell->{gamma} = $box->{gamma};
            $struct->cell->{image}{$i}{p} = 1;
            $struct->cell->{image}{$i}{n} = -1;
	} elsif (! exists($struct->cell->{$i})) {
	    $struct->cell->{$i} = $cellLength;
	    $struct->cell->{alpha} = 90;
	    $struct->cell->{beta} = 90;
	    $struct->cell->{gamma} = 90;
            $struct->cell->{image}{$i}{p} = 1;
            $struct->cell->{image}{$i}{n} = -1;

	}
    }
    $struct->cell->{"valid"} = 1;        
}

sub init {
    my (%OPTS, $boxOpt, @tmp);

    getopt('ftsb', \%OPTS);
    die "usage: $0 -f file name -b (\"box dims\") -t (file type) -s (save name)\n" if (! exists($OPTS{f}));

    print "Initialzing...";
   ($fileName, $fileType, $saveName, $boxOpt) = ($OPTS{f}, $OPTS{t}, $OPTS{s}, $OPTS{b});
    $molStructure->testFile($fileName);
    if (! defined($fileType)) {
	$fileType = "bgf";
	if ($fileName =~ /\.(\w+)$/) {
	    $fileType = lc $1;
	}
    }
    
    if (defined($boxOpt)) {
	push @tmp, $1 while ($boxOpt =~ /(\d+\.?\d*)/g);
	if ($#tmp > 4) {
	    $BOX = {
		"a"     => $tmp[0],
		"b"     => $tmp[1],
		"c"     => $tmp[2],
		"alpha" => $tmp[3],
		"beta"  => $tmp[4],
		"gamma" => $tmp[5]
		};
	}
    }
    $saveName = $molStructure->getFileName($fileName) if (! defined($saveName));
    print "Done\n";
}

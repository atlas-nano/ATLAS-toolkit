#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use FileFormats qw(GetBGFFileInfo addHeader createHeaders createBGF addBoxToHeader);
use BOX qw(GetBox);
use ManipAtoms qw(GetMols GetAtmData);
use General qw(FileTester CoM);

sub init;
sub getNewBox;
sub stressCell;

my ($bgfFile, $saveFile, $newCell);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $MOLS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, undef) = GetBGFFileInfo($bgfFile);
$BOX = GetBox($ATOMS, undef, undef);
$MOLS = GetMols($ATOMS, $BONDS);
&getNewBox($BOX, $newCell);
print "Done\nScaling cell by $newCell->{STRING}...";
&stressCell($ATOMS, $MOLS, $newCell);
print "Done\nCreating BGF file $saveFile...";
$HEADERS = createHeaders(undef, $saveFile);
&addBoxToHeader($HEADERS, $BOX);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub stressCell {
	my ($atoms, $mols, $cellScaling) = @_;
	my ($i, $j, $k, $cellExtrema, $molCoM, $currMol, $molScale, $offset);

	my @dims = ("XCOORD", "YCOORD", "ZCOORD");
	for $i (@dims) {
		$cellExtrema->{$i}{MIN} = 99999999999999999;
		$cellExtrema->{$i}{MAX} = -99999999999999999;
	}
	for $i (keys %{ $atoms }) {
		for $j (@dims) {
			$cellExtrema->{$j}{MIN} = $atoms->{$i}{$j} if ($cellExtrema->{$j}{MIN} > $atoms->{$i}{$j});
			$cellExtrema->{$j}{MAX} = $atoms->{$i}{$j} if ($cellExtrema->{$j}{MAX} < $atoms->{$i}{$j});
		}
	}

	for $i (keys %{ $mols }) {
		$currMol = GetAtmData($atoms, $mols->{$i}{MEMBERS});
		$molCoM = CoM($currMol);
		for $j (@dims) {
			$molScale = ($molCoM->{$j} - $cellExtrema->{$j}{MIN}) * $cellScaling->{$j};
			$offset = $molCoM->{$j} - $molScale - $cellExtrema->{$j}{MIN};
			for $k (keys %{ $mols->{$i}{MEMBERS} }) {
				$atoms->{$k}{$j} -= $offset;
			}
		}
	}
}

sub getNewBox {
	my ($oldBox, $newBox) = @_;
	my ($i);

	for $i ("X","Y","Z") {
		$oldBox->{$i}{hi} *= $newBox->{"${i}COORD"};
		$oldBox->{$i}{lo} *= $newBox->{"${i}COORD"};
		$oldBox->{$i}{len} *= $newBox->{"${i}COORD"};
	}
}

sub init {
	my (%OPTS, $cStr);
	getopt('bcs',\%OPTS);
	for ("b", "c") {
		die "usage: $0 -b bgf file -c \"x y z\" cell length scaling -s (save name)\n"
			if (! exists($OPTS{$_}));
	}
	print "Initializing...";
	($bgfFile, $cStr, $saveFile) = ($OPTS{b}, $OPTS{c}, $OPTS{s});
	FileTester($bgfFile);
	if ($cStr !~ /(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)/) {
		die "ERROR: Expected integers for x,y and z cell length. Got \"$cStr\"\n";
	} else {
		$newCell = (
					{
						"XCOORD" => $1,
						"YCOORD" => $2,
						"ZCOORD" => $3,
						"STRING" => "${1}x${2}x${3}",
					}
					);
	}

	if (! defined($saveFile)) {
		$saveFile = basename ($bgfFile);
		$saveFile =~ s/\.\w+$//;
		$saveFile .= "_compress_" . $newCell->{STRING} . ".bgf";
	}
	print "Done\n";
}

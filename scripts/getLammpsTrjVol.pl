#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester TrjSelections GetStats);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use FileFormats qw(GetBGFFileInfo);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub saveCellCoords;
sub writeData;

my ($lammpsTrj, $bgfFile, $SELECT, $saveFile); 
my ($BGF, $BONDS, $LAMMPSOPTS, $pStr, $outData, $cVol);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($BGF, $BONDS, undef) = GetBGFFileInfo($bgfFile, 1);
print "Done\n";
&GetLammpsByteOffset($SELECT, $lammpsTrj, scalar keys %{ $BGF });
&GetLammpsTrjType($SELECT, $lammpsTrj, "atom", \%{ $LAMMPSOPTS });
$pStr = "Parsing LAMMPS trajectory $lammpsTrj...";
ParseLAMMPSTrj($BGF, $lammpsTrj, $SELECT, "atom", \&saveCellCoords, $pStr, \*OUTDATA);
&writeData($outData, $saveFile);

sub writeData {
    my ($outText, $saveName) = @_;
    my ($avg, $stats, $list, $i, $j);

    @{ $list } = ("VOL", "X", "Y", "Z", "XY");
    open OUTDATA, "> $saveName" or die "ERROR: Cannot write to $saveName: $!\n";
    printf OUTDATA "%-8s%12s%12s%12s%12s%12s\n", "#Timestep", "Volume", "X", "Y", "Z", "XY";
    print OUTDATA $outText;
    for $i (@{ $list }) {
	$stats->{$i} = GetStats($cVol->{$i});
    }
    for $i ("AVG", "STDEV") {
        printf OUTDATA "%-8s", "#${i}";
	for $j (@{ $list }) {
	    printf OUTDATA "%12.3f", $stats->{$j}{$i};
	}
	print OUTDATA "\n";
    }
    close OUTDATA or die "ERROR: Cannot close $saveName: $!\n";
}

sub saveCellCoords {
    my ($DATA, $bgfInfo, $fileHandle) = @_;
    my ($BOX, $vol, $x, $y, $z);

    $BOX = ConvertLammpsBox($DATA->{"BOX BOUNDS"});
    for ("X", "Y", "Z") {
	push @{ $cVol->{$_} }, $BOX->{"${_}COORD"}{len};
    }
    ($x, $y, $z) = ($BOX->{XCOORD}{len}, $BOX->{YCOORD}{len}, $BOX->{ZCOORD}{len});
    $vol = $x * $y * $z;
    push @{ $cVol->{XY} }, ($x*$y);
    push @{ $cVol->{VOL} }, $vol;
    $outData .= sprintf("%-8d%12.3f%12.3f%12.3f%12.3f%12.3f\n", $DATA->{TIMESTEP}[0], $vol, $x, $y, $z, ($x * $y));
}


sub init {
    my (%OPTS, $trajSelect);

    getopt('blst',\%OPTS);
    for ("l", "b") {
	die "usage: $0 -b bgf file -l lammps trajectory file -t (trajectory selection) -s (save file)\n" 
	    if (! exists($OPTS{$_}));
    }

    print "Initializing...";
    ($bgfFile, $lammpsTrj, $trajSelect, $saveFile) = ($OPTS{b}, $OPTS{l}, $OPTS{t}, $OPTS{s});
    FileTester($lammpsTrj);
    FileTester($bgfFile);

    $trajSelect = "*" if (! defined($trajSelect));
    $SELECT = TrjSelections($trajSelect);    
    if (! defined($saveFile)) {
	$saveFile = basename($lammpsTrj);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_vol.dat";
    }
    print "Done\n";
}

#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester TrjSelections CoM GetSoluteAtoms);
use FileFormats qw(GetBGFFileInfo);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use ManipAtoms qw(UnwrapAtoms ImageAtoms SelectAtoms BuildAtomSelectionString GetAtmData GetMols);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub numerically { ($a<=$b); }
sub writeData;
sub calcAvg;
sub calcBins;
sub computeAvg;
sub adjustTrjVals;
sub sortAtomsByField;
sub init;
sub showUsage;

my ($lammpsTrj, $bgfFile, $selection, $reImageCoord, $savePrefix, $slabDim, $field, $binsize);
my ($BGF, $BONDS, $FIELDS, $atmCount, $MOLS, $SOLUTEATMS, $bounds);
my ($SELECT, $pStr, $count, $BINS, $STATS, $LAMMPSOPTS, $frameNum);
my ($imgPoint, $soluteSelect, $reScaleCoord, $com_center);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($BGF, $BONDS) = GetBGFFileInfo($bgfFile, 0);
die "ERROR: Field $field is not valid\n" 
    if (! exists($BGF->{1}{$field}));
$FIELDS = sortAtomsByField($BGF, $field);
$atmCount = scalar(keys %{ $BGF });
$MOLS = GetMols($BGF, $BONDS);
$SOLUTEATMS = SelectAtoms($soluteSelect, $BGF) if (defined($soluteSelect));
print "Done\n";
&GetLammpsByteOffset($SELECT, $lammpsTrj, scalar keys %{ $BGF });
&GetLammpsTrjType($SELECT, $lammpsTrj, "atom", \%{ $LAMMPSOPTS });
$pStr = "Parsing LAMMPS trajectory $lammpsTrj to get averages...";
$count = 0;
&ParseLAMMPSTrj($BGF, $lammpsTrj, $SELECT, "atom", \&calcAvg, $pStr, \*OUTDATA);
&computeAvg($STATS, $frameNum);
$pStr = "Parsing LAMMPS trajectory $lammpsTrj to get histograms...";
&ParseLAMMPSTrj($BGF, $lammpsTrj, $SELECT, "atom", \&calcBins, $pStr, \*OUTDATA);
print "Writing Data...";
&writeData($BINS, $bounds, $savePrefix);
print "Done\n";

sub writeData {
    my ($bins, $minMax, $prefix) = @_;
    my ($i, $j, $k, $l, $outname, $index, $val, $prec);

    if ($binsize !~ /\.\d+/) {
	$prec = 1;
    } else {
        $prec = $binsize;
	$prec =~ s/.*\.//;
	$prec = length($prec) + 1;
    }
    
    for $i (keys %{ $minMax }) {
	$minMax->{$i}{TOT} = ($minMax->{$i}{max} - $minMax->{$i}{min})/$binsize;
    }

    for $l (keys %{ $bins }) {
	$outname = "${prefix}.${l}.dat";
	open OUTFILE, "> $outname" or die "ERROR: Cannot write to $outname: $!\n";
	for $i (sort numerically keys %{ $bins->{$l} }) {
	    for $j (sort numerically keys %{ $bins->{$l}{$i} }) {
		for $k (sort numerically keys %{ $bins->{$l}{$i}{$j} }) {
		    $val = $bins->{$l}{$i}{$j}{$k};
		    printf OUTFILE "%.${prec}f %.${prec}f %.${prec}f %10d\n", ($i)*$binsize,($j)*$binsize,($k)*$binsize,$val;
		}
	    }
	    printf OUTFILE "\n";
	}
	close OUTFILE;
    }
}

sub calcAvg {
    my ($DATA, $bgfInfo, $fileHandle) = @_;
    my ($i, $j, $BOX);

    $frameNum++;
    $BOX = ConvertLammpsBox($DATA->{"BOX BOUNDS"});
    &adjustTrjVals($DATA->{ATOMS}, $BOX);

    for $i (keys %{ $DATA->{ATOMS} }) {
	for $j ("XCOORD", "YCOORD", "ZCOORD") {
	    $STATS->{$i}{$j}{TOT} += $DATA->{ATOMS}{$i}{$j};
	}
    }
}

sub calcBins {
    my ($DATA, $bgfInfo, $fileHandle) = @_;
    my ($i, $j, $BOX, $del, $index);

    $frameNum++;
    $BOX = ConvertLammpsBox($DATA->{"BOX BOUNDS"});
    &adjustTrjVals($DATA->{ATOMS}, $BOX);

    for $i (keys %{ $DATA->{ATOMS} }) {
	for $j ("XCOORD", "YCOORD", "ZCOORD") {
	    $del = $STATS->{$i}{$j}{AVG} - $DATA->{ATOMS}{$i}{$j};
	    $index->{$j} = sprintf("%.0f",$del/$binsize);
	    $bounds->{$j}{min} = $index->{$j} if ($index->{$j} < $bounds->{$j}{min});
	    $bounds->{$j}{max} = $index->{$j} if ($index->{$j} > $bounds->{$j}{max});
	}
	$BINS->{ $FIELDS->{$i} }{$index->{XCOORD}}{$index->{YCOORD}}{$index->{ZCOORD}}++;
	
    }
}

sub computeAvg {
    my ($data, $factor) = @_;
    my ($i);

    for $i (keys %{ $data }) {
	for (keys %{ $data->{$i} }) {
	    $data->{$i}{$_}{AVG} = $data->{$i}{$_}{TOT}/$factor;
	}
    }
}

sub adjustTrjVals {
    my ($atoms, $box) = @_;
    my (@tmp, $i, $j, $MOLECULE, $CENTER, $index);

    if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
	UnwrapAtoms($atoms, $box, $LAMMPSOPTS->{scaled});
    }
    if ($reImageCoord) {
        @tmp = ("XCOORD", "YCOORD", "ZCOORD");
        for $j (@tmp) {
            $CENTER->{$j} = $box->{$j}{CENTER} = $box->{$j}{len}/2;
	}
        if (defined($SOLUTEATMS)) {
	    $MOLECULE = GetAtmData($atoms, $SOLUTEATMS);
            $CENTER = CoM($MOLECULE);
	} elsif (defined($imgPoint)) {
	    $CENTER = $imgPoint;
	} elsif (defined($com_center)) {
	    $CENTER = CoM($atoms);
	}
	for $j (@tmp) {
	    $box->{$j}{hi} = $box->{$j}{len};
	    $box->{$j}{lo} = 0;
	    $CENTER->{$j} = $box->{$j}{CENTER} if(defined($slabDim) and $j eq $slabDim);
	    for $i (keys %{ $atoms }) {
		$atoms->{$i}{$j} += ($box->{$j}{CENTER} - $CENTER->{$j});
	    }
	}
        for $index (keys %{ $MOLS }) {
	    $MOLECULE = GetAtmData($atoms, $MOLS->{$index}{MEMBERS});
	    $CENTER = CoM($MOLECULE);
	    ImageAtoms($MOLECULE, $CENTER, $box);
        }
    }
}

sub sortAtomsByField {
    my ($atoms, $field) = @_;
    my ($i, $FIELDS);

    for $i (keys %{ $atoms }) {
	$FIELDS->{$i} = $atoms->{$i}{$field};
    }

    return $FIELDS;
}

sub init {
    my ($i, %OPTS, $usage, $atomSelection, @tmp, $coordOpt, $pointStr);

    getopt('lbtsmapdnfz',\%OPTS);
    $usage = &showUsage;

    ($lammpsTrj,$bgfFile,$selection,$savePrefix,$coordOpt,$atomSelection,$pointStr,$slabDim,$com_center,$field,$binsize) = 
       ($OPTS{l}, $OPTS{b}, $OPTS{t}, $OPTS{s}, $OPTS{m}, $OPTS{a}, $OPTS{p}, $OPTS{d}, $OPTS{n}, $OPTS{f}, $OPTS{z});

    for ($lammpsTrj, $bgfFile, $selection) {
	die "$usage\n" if (! defined($_));
    }
    print "Initializing...";

    FileTester($lammpsTrj);
    FileTester($bgfFile);
    
    $coordOpt = 0 if (! defined($coordOpt) or $coordOpt !~ /^(1|2|3)$/);
    if ($coordOpt == 0) {
	$reScaleCoord = 0;
	$reImageCoord = 0;
    } elsif($coordOpt==1) {
        $reScaleCoord = 0;
        $reImageCoord = 1;
    } elsif($coordOpt==2) {
        $reScaleCoord = 1;
        $reImageCoord = 1;
    } else {
        $reScaleCoord = 1;
        $reImageCoord = 0;
    }
    print "with coordinate imaging.." if ($reImageCoord);
    print "without coordinate imaging..." if (! $reImageCoord);

    if (! defined($savePrefix)) {
	$savePrefix = basename($lammpsTrj);
	$savePrefix =~ s/\.\w+$//;
    }

    $SELECT = TrjSelections($selection);
    
    if (defined($atomSelection)) {
	$soluteSelect = BuildAtomSelectionString($atomSelection);
    }
    if (defined($pointStr) and $pointStr =~ /(\-?\d+\.?\d*)\s+(\-?\d+\.?\d*)\s+(\-?\d+\.?\d*)/) {
	($imgPoint->{XCOORD}, $imgPoint->{YCOORD}, $imgPoint->{ZCOORD}) = ($1, $2, $3);
    }
    undef($soluteSelect) if defined($imgPoint);
    if(defined($slabDim) && $slabDim =~ /^(x|y|z)$/) {
	$slabDim = uc $1 . "COORD";
    }

    $com_center = 0 if (! defined($com_center) or $com_center !~ /1|yes/i);
    $com_center = 1 if ($com_center =~ /1|yes/i);

    $frameNum = 0;
    $field = "FFTYPE" if (! defined($field));
    $field = uc $field;
    $binsize = 0.01 if (! defined($binsize) or $binsize !~ /^\d+\.\d*/);
    $bounds->{XCOORD}{max} = $bounds->{YCOORD}{max} = $bounds->{ZCOORD}{max} = -99999999;
    $bounds->{XCOORD}{min} = $bounds->{YCOORD}{min} = $bounds->{ZCOORD}{min} = 999999999;

    print "Done\n";
}

sub showUsage {
    my ($usage) = <<DATA;
usage: $0 -b bgf file -l lammps trj -t trajectory selection -d [periodic dimension] -z [bin size] -f [atom group field] -n [com center] -m [reimage/rescale] -p [spacial point] -a [solute atoms] -s [save prefix]
options:
	-b bgf file (Required): the location of the bgf file
	-l lammps trj (Required): the location of the lammps trajectory file
	-t trajectory selection: The number of frames to use. Can be a single integer, or several integers in quotes
		To specify a range specify it as :Ita-b:c, which will take frames a - b, every c. Specify multiple 
		ranges or a combination of ranges and single frames by enclosing them in quotes. \"*\" for all frames
	-m [reimage/rescale coordinates]: (Optional) Whether to reimage and rescale the coordinates back into the 
		unit cell. 0:No reimage or rescale (Default), 1:reimage no rescale, 2:rescale and rescale 
		3:rescale no reimage
	-p [spacial point]: (Optional) Image the coordinates wrt to point in space - or -
	-n [com center]: (Optional) Image coordinated wrt center of mass - or -
	-a [solute atoms]: (Optional) Image the coordinates wrt to the center of mass of certain atoms
		Any valid field (fftype, resname, xcoord etc) expression (e.g. "resname eq 'WAT'")
	-d [periodic dimension]: (Optional) Dimensionality. Default xyz
	-s [save prefix]: (Optional). The name the ouput files
	-f [atom group field]: (Optional) Field to group atoms by. Any valid bgf field can be used. Default fftype
	-z [bin size]: Displacement bin size. Default 0.01 angstroms
DATA

    return $usage;
}

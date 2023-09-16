#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use warnings;
no warnings "recursion";
use Packages::General qw(FileTester TrjSelections Trim CenterOnMol CoP GetSoluteAtoms GetBondLength GetSelections);
use Packages::FileFormats qw(GetBGFFileInfo GetBondList createHeaders addHeader createBGF);
use Packages::LAMMPS qw(ParseLAMMPSTrj CreateLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType);
use Packages::AMBER qw(CreateAmberTrj);
use Packages::ManipAtoms qw(UnwrapAtoms ScaleAtoms ImageAtoms GetAtmList);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub saveCoords;
sub makeBox;
sub init;
sub createLAMMPSHeader;
sub updateBGF;
sub getAtoms;
sub doResFix;
sub showUsage;

my ($lammpsTrj, $bgfFile, $selection, $reImageCoord, $saveName, $outType);
my ($printAtoms, $SELECT, $LAMMPSOPTS, $pStr, $SOLUTEATMS, $count);
my ($BGF, $BONDS, $atmBondList, $isTrj);

$|++;
&init;

print "Parsing BGF file $bgfFile...";
($BGF, $BONDS) = GetBGFFileInfo($bgfFile, 0);
if (! keys %{ $atmBondList }) {
    $atmBondList = GetBondList($BGF, $BONDS);
    $SOLUTEATMS = GetSoluteAtoms($BGF, $atmBondList);
} else {
    $SOLUTEATMS = GetAtmList($atmBondList, $BGF);
    $atmBondList = GetBondList($BGF, $BONDS);
}
die "ERROR: No valid solute atoms found!\n" if (! keys %{ $SOLUTEATMS });
print "Done\n";
&GetLammpsByteOffset($SELECT, $lammpsTrj, scalar keys %{ $BGF });
&GetLammpsTrjType($SELECT, $lammpsTrj, "atom", \%{ $LAMMPSOPTS });
$pStr = "Parsing LAMMPS trajectory $lammpsTrj...";
$count = 0;
ParseLAMMPSTrj($atmBondList, $lammpsTrj, $SELECT, "atom", \&saveCoords, $pStr, \*OUTDATA);
close OUTDATA or die "ERROR: Cannot close $saveName: $!\n" if ($isTrj);

sub saveCoords {
    my ($DATA, $moleculeList, $fileHandle) = @_;
    my ($BOX, $mID, $HEADERS, @tmp, $CENTER, $MOLECULE, $i, $fName);

    $BOX = makeBox($DATA->{"BOX BOUNDS"});
    UnwrapAtoms($DATA->{ATOMS}, $BOX, $LAMMPSOPTS->{scaled});

    if ($reImageCoord) {
        for $i (1 .. $DATA->{"NUMBER OF ATOMS"}[0]) {
            next if (exists($DATA->{ATOMS}{$i}{VISITED}));
            doResFix($DATA->{ATOMS}, $BOX, $i);
        }

        $MOLECULE = getAtoms($DATA->{ATOMS}, $SOLUTEATMS);
        $CENTER = CoP($MOLECULE);
        CenterOnMol($DATA->{ATOMS}, $CENTER);

        for $mID (keys %{ $moleculeList }) {
	    $MOLECULE = getAtoms($DATA->{ATOMS}, $moleculeList->{$mID});
	    $CENTER = CoP($MOLECULE);
	    ImageAtoms($DATA->{ATOMS}, $moleculeList->{$mID}, $CENTER, $BOX);
        }
    }
    if ($outType eq "AMBER") {
	CreateAmberTrj($DATA->{ATOMS},$BOX, $fileHandle);
    } elsif ($outType eq "LAMMPS") {
	ScaleAtoms($DATA->{ATOMS}, $BOX);
	$HEADERS = createLAMMPSHeader($DATA, $BOX);
	CreateLAMMPSTrj($DATA->{ATOMS}, $HEADERS, $fileHandle);
    } else { #bgf
	$count++;
	$fName = $saveName;
	@tmp = keys %{ $SELECT };
	if (defined($SELECT) && $#tmp > 0) {
	    $fName =~ s/\.\w+$//;
	    $fName .= "_${count}.bgf";
	}
	updateBGF($BGF, $DATA->{ATOMS}, $BOX, $DATA->{TIMESTEP}[0]);
	createBGF($BGF, $BONDS, $fName);
    }
}

sub doResFix {
    my ($ATOMS, $BOX, $currAtm) = @_;
    my ($i, @dims, $j, $bondPos, $atomPos);

    $ATOMS->{$currAtm}{VISITED} = 1;
    @dims = ("X", "Y", "Z");
    for $j (@dims) {
	$atomPos->{$j} = $ATOMS->{$currAtm}{$j . "COORD"};
    }
    for $i ( @{ $BONDS->{$currAtm} }) {
	next if (exists($ATOMS->{$i}{VISITED}));
	for $j (@dims) {
	    $bondPos = $ATOMS->{$i}{$j . "COORD"};
	    if(abs($bondPos - $atomPos->{$j}) > 5) { #image problem
		if ($atomPos->{$j} > $bondPos) { # move in + direction
		    while (($atomPos->{$j} - $bondPos) > 5) {
			$bondPos += ($BOX->{$j}{hi} - $BOX->{$j}{lo});
		    }
		} else { # move in - direction
                    while (($bondPos - $atomPos->{$j}) > 5) {
                        $bondPos -= ($BOX->{$j}{hi} - $BOX->{$j}{lo});
                    }
		}
		$ATOMS->{$i}{$j . "COORD"} = $bondPos;
		$ATOMS->{$i}{$j . "INDEX"} = $ATOMS->{$currAtm}{$j . "INDEX"};
	    }
	}
	doResFix($ATOMS, $BOX, $i);
    }
}	     
	    
sub getAtoms {
    my ($allAtoms, $atomList) = @_;
    my (%ATOMS, $i);

    for $i (keys %{ $atomList }) {
	$ATOMS{$i} = $allAtoms->{$i};
    }

    return \%ATOMS;
}

sub updateBGF {
    my ($REF, $MOD, $BOX, $TSTEP) = @_;
    my (@tmp, $i, $HEADERS, $j);

    $BOX->{1}{DATA} = 90;
    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    for $i (2 .. 4) {
	$BOX->{$i}{DATA} = ($BOX->{$tmp[($i - 2)]}{hi} - $BOX->{$tmp[($i - 2)]}{lo});
    }

    $TSTEP = "ts_${TSTEP}";
    $HEADERS = createHeaders($BOX, $TSTEP);

    for $i (keys %{ $REF }) {
	for $j (keys %{ $MOD->{$i} }) {
	    $REF->{$i}{$j} = $MOD->{$i}{$j};
	}
    }
    addHeader($REF, $HEADERS);
}

	
sub makeBox {
    my ($box) = $_[0];
    my (%BOX, $i, @dim);

    @dim = ("X", "Y", "Z");
    for $i (0 .. $#{ $box }) {
	$BOX{$dim[$i]}{lo} = $box->[$i]{lo};
	$BOX{$dim[$i]}{hi} = $box->[$i]{hi};
	$BOX{$dim[$i] . "COORD"}{lo} = $box->[$i]{lo};
	$BOX{$dim[$i] . "COORD"}{hi} = $box->[$i]{hi};
	$BOX{$dim[$i] . "COORD"}{len} = ($box->[$i]{hi} - $box->[$i]{lo});
    }

    return \%BOX;
}
	    	
sub init {
    my ($i, %OPTS, $usage, $atomSelection, @tmp);

    getopt('lstboam',\%OPTS);
    $usage = &showUsage;

    ($lammpsTrj, $bgfFile, $selection, $outType, $saveName, $reImageCoord, $atomSelection) = 
    ($OPTS{l},$OPTS{b},$OPTS{t},$OPTS{o},$OPTS{s},$OPTS{m},$OPTS{a});

    for ($lammpsTrj, $bgfFile, $selection) {
	die "$usage\n" if (! defined($_));
    }
    print "Initializing...";

    FileTester($lammpsTrj);
    FileTester($bgfFile);
    
    $outType = "lammps" if (! defined($outType));
    if ($outType !~ /(lammps|amber|bgf)/i) {
	print "WARNING: Expected lammps, amber or bgf for out_type. Got $outType...\n" . 
	    "Setting to default value: lammps\n";
	$outType = "lammps";
    } else {
	$outType = $1;
    }
    print "$outType traj ";

    if (! defined($reImageCoord) or $reImageCoord !~ /^1$/) {
	$reImageCoord = 0;
    } else {
	$reImageCoord = 1;
    }
    print "with coordinate imaging.." if ($reImageCoord);
    print "without coordinate imaging..." if (! $reImageCoord);

    if (! defined($saveName)) {
	$saveName = basename($lammpsTrj);
	$saveName =~ s/\.\w+$/_mod\.lammpstrj/;
    }

    $SELECT = TrjSelections($selection);
    $outType = uc($outType);
    
    $isTrj = 0;
    $isTrj = 1 if ($outType =~ /LAMMPS|AMBER/i);
    if ($isTrj) {
	open OUTDATA, "> $saveName" or die "Cannot write to trajectory file $saveName: $!\n";
	$pStr = "Creating $outType trajectory $saveName...";
	print OUTDATA "TITLE: LAMMPS Trajectory file create by lammpsTrj2amber on " . scalar(localtime(time())) . "\n";
    } else {
	$pStr = "Creating $outType files for each frame...";
    }

    if (defined($atomSelection)) {
	@tmp = split /\s+/, $atomSelection;
	$atmBondList = GetSelections(\@tmp, 0);
    }    
    print "Done\n";
}

sub createLAMMPSHeader {
    my ($data, $box) = @_;
    my (%HEADERS, $i, $count);
    
    for $i ("TIMESTEP", "NUMBER OF ATOMS") {
	@{ $HEADERS{$i} } = @{ $data->{$i} };
    }
    
    $count = 0;
    for $i ("XCOORD", "YCOORD", "ZCOORD") {
	for ("hi", "lo") {
	    $HEADERS{"BOX BOUNDS"}[$count]{$_} = $box->{$i}{$_};
	}
	$count++;
    }
	
    return \%HEADERS;
}

sub showUsage {
    my ($usage) = "usage: $0 -b bgf file -l lammps trj -t trajectory selection ". 
	"-o [output type] -m [reimage coordinates] -a [solute atoms] -s [save name]\n" .
	"options:\n" . 
	"-b bgf file: the location of the bgf file. This is required\n" . 
	"-l lammps trj: the location of the lammps trajectory file. This is required.\n" .
	"-t trajectory selection: The number of frames to use. Can be a single integer, or several integers in quotes\n" .
	"\tTo specify a range specify it as :Ita-b:c, which will take frames a - b, every c. Specify multiple ranges\n" .
	"\tor a combination of ranges and single frames by enclosing them in quotes. Specify \"*\" for all frames\n" .
	"-o [output type]: (Optional). Either lammps (default) or amber or bgf. If bgf, will print file for every frame\n" .
	"-m [reimage coordinates]: (Optional) Whether to reimage the coordinates back into the unit cell. Default yes\n" .
	"-a [solute atoms]: (Optional) Specifies the range for the solute atoms. If not specified the code will attempt\n" .
	"\tto determine this automatically. Expected [:][I|N|T][a|r]x[-y:z] for a range(:) or single Index|Name|Type of\n" .
	"\tAtom|Residue x (to y every z)\n" .
	"-s [save name]: (Optional). The name the save the output trajectory\n";
    return $usage;
}

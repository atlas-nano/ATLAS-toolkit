#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use General qw(TrjSelections FileTester GetBondLength CoM);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use FileFormats qw(GetBGFFileInfo GetMOL2FileInfo GetPDBFileInfo);
use ManipAtoms qw(SelectAtoms BuildAtomSelectionString SplitAtomsByMol GetMols GetAtmData ImageAtoms);
use Getopt::Std qw(getopt);

sub init;
sub calcNumAtoms;

my ($readFunc, $dataFile, $atmSELECT, $fileType, $trjSELECT, $trjFile, $saveFile, $refATMS);
my ($ATOMS, $BONDS, $HEADERS, $calcMol, $MOLS, $cATOMS);
my ($field, $pStr, $LAMMPSOPTS, $getByteOffset, $getSnapshots, $trjType, $OUTFILE);

$|++;

&init;

print "Parsing $fileType file $dataFile...";
($ATOMS, $BONDS, $HEADERS) = $readFunc->($dataFile);
delete $ATOMS->{HEADER} if (exists($ATOMS->{HEADER}));
print "Done\n";
$pStr = "Selecting Relevant Atoms...";
$refATMS = SelectAtoms($refATMS, $ATOMS) if (defined($refATMS));
if (! defined($trjFile)) {
    print "$pStr";
    &GetMols($ATOMS, $BONDS);
    ($cATOMS, $BONDS) = SelectAtoms($atmSELECT, $ATOMS);
    $MOLS = SplitAtomsByMol($ATOMS, $cATOMS) if ($calcMol);
    &calcNumAtoms($cATOMS, $MOLS);
} else {
    $field = scalar keys %{ $ATOMS };
    $getByteOffset->($trjSELECT, $trjFile, $field);
    if ($trjType == 2) {
        &GetLammpsTrjType($trjSELECT, $trjFile, "coord", \%{ $LAMMPSOPTS });
        $field = "coord";
    }
    $getSnapshots->($ATOMS, $trjFile, $trjSELECT, $field, \&calcNumAtoms, $pStr, $OUTFILE);
    close $OUTFILE;
}

sub calcNumAtoms {
    my ($atoms, $box, $frameNum, $fileHandle) = @_;
	my (@tmp, $j, $CENTER, $MOLECULE, $i, $cMOLS);

    if (defined($trjFile)) {
		if ($trjType == 2) { #LAMMPS
			$fileHandle = $frameNum;
			$frameNum = $atoms->{TIMESTEP}[0];
			$box = ConvertLammpsBox($atoms->{"BOX BOUNDS"});
			%{ $box->{X} } = %{ $box->{XCOORD} };
			%{ $box->{Y} } = %{ $box->{YCOORD} };
			%{ $box->{Z} } = %{ $box->{ZCOORD} };
			$atoms = $atoms->{ATOMS};
			if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
				UnwrapAtoms($atoms,  $box, $LAMMPSOPTS->{scaled});
			}
		} elsif ($trjType == 1) { #AMBER
			$box = ConvertAmberBox(\%{ $box });
		}
		@tmp = ("XCOORD", "YCOORD", "ZCOORD");
		for $j (@tmp) {
			$CENTER->{$j} = $box->{$j}{CENTER} = $box->{$j}{len}/2;
		}
		if (defined($refATMS)) {
			$MOLECULE = GetAtmData($atoms, $refATMS);
			$CENTER = CoM($MOLECULE);
			for $j (@tmp) {
				$box->{$j}{hi} = $box->{$j}{len};
				$box->{$j}{lo} = 0;
				$CENTER->{$j} = $box->{$j}{CENTER};
				for $i (keys %{ $atoms }) {
					$atoms->{$i}{$j} += ($box->{$j}{CENTER} - $CENTER->{$j});
				}
			}
			for $i (keys %{ $MOLS }) {
				$MOLECULE = GetAtmData($atoms, $MOLS->{$i}{MEMBERS});
				$CENTER = CoM($MOLECULE);
				ImageAtoms($MOLECULE, $CENTER, $box);
			}
		}
		for $i (keys %{ $atoms }) {
			for $j ("XCOORD", "YCOORD", "ZCOORD") {
				$ATOMS->{$i}{$j} = $atoms->{$i}{$j};
			}
		}
	}

    if (!defined($trjFile)) {
		print "Done\nFound " . scalar(keys %{ $atoms }) . " atoms\n";
		print "Found " . scalar(keys %{ $MOLS }) . " molecules\n" if ($calcMol);
    } else {
		$cATOMS = SelectAtoms($atmSELECT, $ATOMS);
		printf $fileHandle "%-20d %10d", $frameNum, scalar(keys %{ $cATOMS });
		if ($calcMol) {
			$cMOLS = GetMols($ATOMS, $BONDS, $cATOMS);
			printf $fileHandle " %10d", scalar(keys %{ $cMOLS });
		}
		print $fileHandle "\n";
    }
}

sub init {
    my (%OPTS, $select, $trjSel, $list, $i, $refSelect);

    getopt('ftsmlwric',\%OPTS);
    die "usage: $0 -f data file -t (file type [bgf|pdb|mol2]) -i (reimage system atom selection) -l (traj file) -c (traj type=lammps(default) or amber) -r (traj range) -s [atom selection] -m (calc mols = no) -w (write traj data)\n" 
	if (! exists($OPTS{f}));

    print "Initializing...";
    ($dataFile, $select, $fileType, $calcMol, $trjFile, $trjType, $trjSel, $saveFile, $refSelect) = 
	($OPTS{f}, $OPTS{s}, $OPTS{t}, $OPTS{m}, $OPTS{l}, $OPTS{c}, $OPTS{r}, $OPTS{w}, $OPTS{i});
    FileTester($dataFile);
    $select = "index > 0" if (! defined($select));
    if (! defined($fileType)) {
		if ($dataFile =~ /.*\.(\w+)$/) {
			$fileType = $1;
		    if ($fileType !~ /^(bgf|mol2|pdb)$/i) {
			$fileType = "bgf";
			}
		} else {
		    $fileType = "bgf";
		}
    }
    $fileType = uc($fileType);
    
    if ($fileType =~ /bgf|mol2|pdb/i) {
		$readFunc = eval ('\&Get' . $fileType . 'FileInfo');
    } else {
		die "ERROR: Expected bgf|mol2|pdb while parsing filetype. Got $fileType\n";
    }
    
	$refATMS = BuildAtomSelectionString($refSelect) if (defined($refSelect));
    $atmSELECT = BuildAtomSelectionString($select);
    $calcMol = 0 if (! defined($calcMol) or $calcMol !~ /(1|yes)/i);
    $calcMol = 1 if ($calcMol =~ /(1|yes)/i);

    if (defined($trjFile) and -e $trjFile and -r $trjFile and -T $trjFile) {
        if (! defined($trjType)) {
            if ($trjFile =~ /\.lammps/) {
                $trjType = "lammps";
            } else {
                $trjType = "amber";
            }
        }

        if (lc($trjType) ne "lammps") {
            $trjType = 1;
            $getSnapshots = \&ParseAmberTrj;
            $getByteOffset = \&GetAmberByteOffset;
        } else {
            $trjType = 2;
            $getSnapshots = \&ParseLAMMPSTrj;
            $getByteOffset = \&GetLammpsByteOffset;
        }
        if (! defined($trjSel)) {
            $trjSel = "*";
        }
        $list = TrjSelections($trjSel);
        for $i (keys %{ $list }) {
            $trjSELECT->{$i} = $list->{$i};
        }
        die "ERROR: No valid frames selected with selection $trjSel!\n"
            if (! keys %{ $trjSELECT } and $trjSel ne "*");
        if (! defined($saveFile)) {
           $saveFile = basename($trjFile);
           $saveFile =~ s/\.\w+$//;
           $saveFile .= ".count.dat";
        }
	 	open $OUTFILE, "> $saveFile" or die "ERROR: Cannot write to $saveFile: $!\n";
		printf $OUTFILE "%-20s %10s","#tstep","numAtoms";
	        printf $OUTFILE " %10s", "numMols" if ($calcMol);
		print $OUTFILE "\n";
    }

    print "Done\n";
    
}

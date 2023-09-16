#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester LoadFFs TrjSelections CoM ReadFFs);
use FileFormats qw(GetBGFFileInfo)
use Getopt::Std qw(getopt);
use ManipAtoms qw(SelectAtoms BuildAtomSelectionString UnwrapAtoms ScaleAtoms GetMols GetAtmData ImageAtoms);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use File::Basename qw(basename);

sub init;
sub addVDWradii;
sub saveXYZR;
sub getSASAandVOL;
sub getSnaps;

my ($bgfFile, $ff, $atmSelection, $atmSELECT, $trjSelection, $trjSELECT);
my ($ATOMS, $BONDS, $FFILES, $FF, $ffList, $msms, $data, $saveFile, $MOLS, $sprefix);
my ($trjFile, $field, $pStr, $LAMMPSOPTS, $getByteOffset, $getSnapshots, $trjType, $OUTFILE, $reImageCoord);

$|++;
&init;
print "Parsing $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
$atmSELECT = SelectAtoms($atmSelection, $ATOMS);
die "ERROR: No valid atoms found!\n" if (! keys %{ $atmSELECT });
$MOLS = GetMols($ATOMS, $BONDS, $atmSELECT);
print "Done\n";
($FFILES, undef) = ReadFFs($ffList);
$FF = LoadFFs($FFILES);
print "Adding VDW radii...";
&addVDWradii($ATOMS, $FF);
if (! defined($trjFile)) {
	print "Done\nGetting surface area and enclosed volume...";
	&getSASAandVOL($ATOMS, $atmSELECT, \%{ $data });
	printf "Done\n%-20s %20s\n%20.5f %20.5f\n","SES-AREA","SES-VOLUME",$data->{sesa},$data->{sesv};
} else {
	$field = scalar keys %{ $ATOMS };
	$getByteOffset->($trjSELECT, $trjFile, $field);
	if ($trjType == 2) {
		&GetLammpsTrjType($trjSELECT, $trjFile, "coord", \%{ $LAMMPSOPTS });
		$field = "coord";
	}
	$pStr = "Getting Solvent Accessible Surface Area and Volume from $trjFile...";
	open $OUTFILE, "> $saveFile" or die "ERROR: Cannot write to $saveFile: $!\n";
	printf $OUTFILE "%-20s %20s %20s\n","#Tstep","SES-AREA","SES-VOLUME";
	$getSnapshots->($ATOMS, $trjFile, $trjSELECT, $field, \&getSnaps, $pStr, $OUTFILE);
	close $OUTFILE;
}

sub getSnaps {
    my ($atoms, $box, $frameNum, $fileHandle) = @_;
    my ($tot, $geomVal, $atmList, $data);
	my (@tmp, $MOLECULE, $CENTER, $j, $i);

    if ($trjType == 2) { #LAMMPS
		$fileHandle = $frameNum;
        $frameNum = $atoms->{TIMESTEP}[0];
        $box = ConvertLammpsBox($atoms->{"BOX BOUNDS"});
        $tot = $ATOMS->{"NUMBER OF ATOMS"}[0];
        $atoms = $atoms->{ATOMS};
        if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
            UnwrapAtoms($atoms,  $box, $LAMMPSOPTS->{scaled});
        }
    } else {
        $box = ConvertAmberBox(\%{ $box });
        $tot = scalar(keys %{ $atoms });
    }
	@tmp = ("XCOORD", "YCOORD", "ZCOORD");
	if ($reImageCoord) {
		$MOLECULE = GetAtmData($atoms, $atmSELECT);
		$CENTER = CoM($MOLECULE);
		for $j (@tmp) {
			$box->{$j}{hi} = $box->{$j}{len};
			$box->{$j}{lo} = 0;
			for $i (keys %{ $atmSELECT }) {
				$atoms->{$i}{$j} += ($box->{$j}{CENTER} - $CENTER->{$j});
				$atoms->{$i}{RADII} = $ATOMS->{$i}{RADII};
			}
		}
		for $i (keys %{ $MOLS }) {
			$MOLECULE = GetAtmData($atoms, $MOLS->{$i}{MEMBERS});
			$CENTER = CoM($MOLECULE);
			&ImageAtoms($MOLECULE, $CENTER, $box);
		}
	} else {
		for $j (@tmp) {
			for $i (keys %{ $atmSELECT }) {
				$atoms->{$i}{RADII} = $ATOMS->{$i}{RADII};
			}
		}
	}
	$data = ();
	&getSASAandVOL($atoms, $atmSELECT, \%{ $data });
	printf $fileHandle "%-20d %20.5f %20.5f\n",$frameNum, $data->{sesa}, $data->{sesv};	

}

sub getSASAandVOL {
	my ($atoms, $aSel, $data) = @_;
	my ($calc_Cmd, $valid, $sesa, $sesv);

	$valid = 0;

	&saveXYZR($atoms,$aSel,$sprefix);
	open MSMSCMD, "$msms -if $sprefix | " or die "ERROR: Cannot execute $msms -if __xyzr: $!\n";
	while(<MSMSCMD>) {
		chomp;
		if ($_ =~ /^\s*Comp. probe_radius SES_volume SES_area/) {
			$valid = 1;
		} elsif ($valid and $_ =~ /\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
			$data->{sesv} = $1;
			$data->{sesa} = $2;
			$valid = 2;
		}
	}
	close MSMSCMD;

	die "ERROR: Could not determine the surface area or volume!\n" if ($valid < 2);
	system("rm -fr $sprefix");
}

sub saveXYZR {
	my ($atoms, $select, $sname)  = @_;
	my ($i, $atom, $tot);

	$tot = scalar(keys %{ $atoms });
	open XYZRFILE, "> ${sname}" or die "ERROR: Cannot write to ${sname}: $!\n";
	for $i (1 .. $tot) {
		next if (!exists($select->{$i}));
		$atom = $atoms->{$i};
		printf XYZRFILE "%10.5f %10.5f %10.5f %10.5f\n",$atom->{XCOORD},$atom->{YCOORD},$atom->{ZCOORD},$atom->{RADII};
	}
	close XYZRFILE;
}

sub addVDWradii {
	my ($atoms, $ff) = @_;
	my ($i, $fftype, $radii);

	for $i (keys %{ $atoms }) {
		$fftype = $atoms->{$i}{"FFTYPE"};
		die "ERROR: Atom $atoms->{$i}{ATMNAME} does not have a force field type\n" if (! defined($fftype));
		die "ERROR: Atom $atoms->{$i}{ATMNAME} does not have a VDW entry\n" if (!exists($ff->{VDW}{$fftype}{$fftype}{1}{VALS}));
		$radii = $ff->{VDW}{$fftype}{$fftype}{1}{VALS}[1];
		$atoms->{$i}{RADII} = $radii/2;
	}
}

sub init {
    my (%OPTS, $aSel, $tSel, $i, $list, $probeR);
	

    getopt('bfamtswrp',\%OPTS);
    for ("b", "f") {
	die "usage: $0 -b bgf file -f force field -a (atom selection = all) -m (location of msms file) -t (traj file) -s (traj selection) -w (save file) -r (reimage coordinates = 0|no) -p (probe_radius)\n"
	    if (! exists($OPTS{$_}));
    }

    print "Initializing...";
    ($bgfFile, $ffList, $aSel, $trjFile, $tSel, $saveFile, $reImageCoord, $probeR) = 
		($OPTS{b}, $OPTS{f}, $OPTS{a}, $OPTS{t}, $OPTS{s}, $OPTS{w}, $OPTS{r}, $OPTS{v});
    FileTester($bgfFile);

	$reImageCoord = 0 if (! defined($reImageCoord) or $reImageCoord !~ /1|yes/i);
	$reImageCoord = 1 if ($reImageCoord =~ /1|yes/i);

    $aSel =  "index > 0" if (! defined($aSel));
    $atmSelection = BuildAtomSelectionString($aSel);

	$msms = "/home/tpascal/codes/msms/2.6.1/msms.x86_64Linux2.2.6.1.staticgcc";
	$msms = $OPTS{m} if (exists($OPTS{m}));
	if (defined($probeR) and $probeR =~ /(\d+\.?\d*)/) {
		$msms .= "-p $1"
	}
	die "ERROR: Cannot access $msms!\n" if (! -e $msms or ! -r $msms);
	$sprefix = basename($bgfFile);
	$sprefix =~ s/\.\w+$//;
	if (defined($trjFile) and -e $trjFile and -r $trjFile and -T $trjFile) {
		$trjType = "amber";
		$trjType = "lammps" if ($trjFile =~ /\.lammps/);
		$getSnapshots = \&ParseLAMMPSTrj;
		$getByteOffset = \&GetLammpsByteOffset;
		if (lc($trjType) ne "lammps") {
			$trjType = 1;
			$getSnapshots = \&ParseAmberTrj;
			$getByteOffset = \&GetAmberByteOffset;
		} else {
			$trjType = 2;
		}
		$tSel = "*" if (! defined($tSel));
		$list = TrjSelections($tSel);
		for $i (keys %{ $list }) {
			$trjSELECT->{$i} = $list->{$i};
		}
		die "ERROR: No valid frames selected with selection $tSel!\n"
			if (! keys %{ $trjSELECT } and $tSel ne "*");
		$sprefix = basename($trjFile);
		$sprefix =~ s/\.\w+$//;
		$saveFile = "${sprefix}.sesdata.dat" if (! defined($saveFile));
		$sprefix = "__${sprefix}";
	}

    print "Done\n";

}


#!/usr/bin/perl -w

use strict;
no warnings "recursion";
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use GRID;
use Getopt::Std qw(getopt);
use FileFormats qw(GetBGFFileInfo GetMOL2FileInfo GetPDBFileInfo GetXYZFileInfo);
use General qw(FileTester CoM);
use File::Basename qw(basename);
use BOX qw(GetBox);
use Storable qw(dclone);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use ManipAtoms qw(GetAtmData ImageAtoms GetMols);

sub usage;
sub calcVoids;
sub writeResults;

my ($inFile, $OUTFILE, $probe_r, $trjFile);
my ($readFunc, $getSnapshots, $getByteOffset); 
my ($trjType, $field, $pStr, $SELECT, $LAMMPSOPTS);
my ($ATOMS, $BONDS, $BOX, $MOLS, $HEADERS, $DATA);

$|++;
&init;
print "Getting atom information from $inFile...";
($ATOMS, $BONDS, $HEADERS) = $readFunc->($inFile,1);
$MOLS = GetMols($ATOMS, $BONDS);
$BOX = GetBox($ATOMS, undef, $HEADERS);
print "Done\n";
if ($trjType == 0) {
	&calcVoids($ATOMS, $BOX, 1, $OUTFILE);
} else {
	$field = scalar keys %{ $ATOMS };
	$getByteOffset->($SELECT, $trjFile, $field);
	if ($trjType == 2) {
		&GetLammpsTrjType($SELECT, $trjFile, "coord", \%{ $LAMMPSOPTS });
		$field = "coord";
	}
	$pStr = "Performing void analysis on $trjFile...";
	$getSnapshots->($ATOMS, $trjFile, $SELECT, $field, \&calcVoids, $pStr, $OUTFILE);
	close $OUTFILE;
}
print "Done\nSaving output...";
&writeResults($DATA, $OUTFILE);
print "Done\n";

sub calcVoids {
	my ($atoms, $box, $frame, $fileH) = @_;
	my ($i, $j, $grid, $tot, $gs, $voids, $CoM, $cMOL, $gm);

	if ($trjType == 2) { #LAMMPS
		$frame = $atoms->{TIMESTEP}[0];
		$box = ConvertLammpsBox($atoms->{"BOX BOUNDS"});
		$tot = $ATOMS->{"NUMBER OF ATOMS"}[0];
		$atoms = $atoms->{ATOMS};
		if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
			&UnwrapAtoms($atoms,  $box, $LAMMPSOPTS->{scaled});
		}
	} elsif ($trjType == 1) {
		$box = ConvertAmberBox(\%{ $box });
		$tot = scalar(keys %{ $atoms });
	}

	#initialize grid
	$gs = $probe_r*2; #cell length
	@{ $gm } = (int($box->{X}{len}/$gs)+1,int($box->{Y}{len}/$gs)+1,int($box->{Z}{len}/$gs)+1);
	$grid = new GRID (@{ $gm }, $gs); #initialize grid

	#reimage atoms into box based on molecule COM
	$CoM = CoM($atoms, $box);
	for $i (keys %{ $atoms }) {
		for $j ("X", "Y", "Z") {
			$atoms->{$i}{"${j}COORD"} += (($box->{$j}{hi}-$box->{$j}{lo})/2 - $CoM->{"${j}COORD"});
		}
	}
	for $j ("X", "Y", "Z") {
		$box->{"${j}COORD"} = $box->{$j};
	}
	for $i (keys %{ $MOLS }) {
		$cMOL = GetAtmData($atoms, $MOLS->{$i}{MEMBERS});
		$CoM = CoM($cMOL, $box);
		&ImageAtoms($cMOL, $CoM, $box);
	}

	
	#place atoms on grid
	for $i (keys %{ $atoms }) {
		$grid->store($atoms->{$i});
	}

	#determine void in structure
	$voids = $grid->find_voids();
}


sub init {
	my (%OPTS, $i, $tSel, $list, $saveFile, $itype);

	getopt('ijrltsw',\%OPTS);
	($inFile, $itype, $probe_r, $trjFile, $trjType, $tSel, $saveFile) = 
		($OPTS{i},$OPTS{j},$OPTS{r},$OPTS{l},$OPTS{t},$OPTS{s},$OPTS{w});

	for ($inFile, $probe_r) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	FileTester($inFile);
	die "ERROR: Expected number/decimal for probe radius. Got '$probe_r'\n"
		if($probe_r !~ /^\d+\.?\d*$/);

	die "ERROR: Cannot determine file type for $inFile!\n"
		if($inFile !~ /\.(\w+)$/ and ! defined($itype));
	$itype = $1 if (!defined($itype));
	die "ERROR: Invalid file type $itype. Expected bgf|mol2|pdb|xyz\n"
		if($itype !~ /(pdb|bgf|mol2|xyz)/i);
	$itype = uc $1;
	$readFunc = \&GetBGFFileInfo;
	$readFunc = \&GetMOL2FileInfo if ($itype =~ /MOL2/);
	$readFunc = \&GetPDBFileInfo  if ($itype =~ /PDB/);
	$readFunc = \&GetXYZFileInfo  if ($itype =~ /XYZ/);

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
		if (! defined($tSel)) {
			$tSel = "*";
		}
		$list = TrjSelections($tSel);
		for $i (keys %{ $list }) {
			$SELECT->{$i} = $list->{$i};
		}
		die "ERROR: No valid frames selected with selection $tSel!\n"
			if (! keys %{ $SELECT } and $tSel ne "*");
		if(!defined($saveFile)) {
			$saveFile = basename($trjFile);
			$saveFile =~ s/\.\w+$//;
			$saveFile = "${saveFile}.voids.dat";
		}
	} else {
		$trjType = 0;
	}
	if (defined($saveFile)) {
		open $OUTFILE, "> $saveFile" or die "ERROR: Cannot write to $saveFile: $!\n";
	} else {
		$OUTFILE = \*STDOUT;
	}
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -i structure_file -j (file_type) -r probe_radius -l (trajectory_file) -t (trajectory_type) -s (trajectory_selection) -w (savename)
Arguments:
  structure_file: structure file for analysis (required)
  file_type: type of structure file. Expected XYZ|MOL2|BGF|PDB (Optional)
  probe_radius: radius of test particle (Required)
  trajectory_file: location of trajectory file. (Optional)
  trajectory_type: type of trajectory file. Expected LAMMPS(.lammps)|AMBER (.crdbox|.mdcrd) (Optional)
  trajectory_selection: frame/frame range to use. Can be 'first|last|frame_number'. Specify range using ':Ita-b:c'
						to select frames a - b, every c. (Optional)
  savename: name of output file to save results. Default is none and write to screen. (optional) 						
DATA
die "\n";

}

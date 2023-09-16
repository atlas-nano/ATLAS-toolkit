#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(ParseStructFile);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use General qw(GetFileTypeStr FileTester GetBondLength GetAngle GetTorsion TrjSelections CoM);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use BOX qw(GetBox box2H);
use ManipAtoms qw(GetAtmData GetMols UnwrapAtoms ScaleAtoms SelectAtoms 
					BuildAtomSelectionString RemoveMolLongBonds ReimageMols);

sub numerically { $a<=>$b }

my ($structFile, $selection, $geomType, $trjFile, $SELECT, $saveFile, $binSize, $binFile, $comOpt, $sym_angle);
my ($field, $pStr, $LAMMPSOPTS, $getByteOffset, $getSnapshots, $trjType, $OUTFILE, $rec);
my ($atmList, $ATOMS, $BONDS, $MOLS, $HEADERS, $geomVal, $getGeomVal, $DATA, $BOX);

$|++;
&init;
print "Parsing structure file $structFile...";
($ATOMS, $BONDS, $HEADERS) = ParseStructFile($structFile,1);
&GetMols($ATOMS, $BONDS);
$BOX = GetBox($ATOMS, undef, $HEADERS);
($atmList,$geomType) = getAtoms($ATOMS, $BONDS, $selection, $BOX, $comOpt);
$MOLS = GetMols($ATOMS, $BONDS);
print "Done\n";
if (! defined($trjFile)) {
	@{ $rec } = map { \%{ $ATOMS->{$_} } } @{ $atmList->[0] };
	$geomVal = $getGeomVal->(@{ $rec },$BOX, 0);	
  	print "$geomType = $geomVal\n";
} else {
	$field = scalar keys %{ $ATOMS };
	$getByteOffset->($SELECT, $trjFile, $field);
	if ($trjType == 2) {
		&GetLammpsTrjType($SELECT, $trjFile, "coord", \%{ $LAMMPSOPTS });
		$field = "coord";
	}
	$pStr = "Calculating $geomType data from $trjFile...";
	open $OUTFILE, "> $saveFile" or die "ERROR: Cannot write to $saveFile: $!\n";
	printf $OUTFILE "%-10s%10s\n","Tstep",$geomType;
	$getSnapshots->($ATOMS, $trjFile, $SELECT, $field, \&getTrjGeom, $pStr, $OUTFILE);
	close $OUTFILE;
}
if (defined($binFile)) {
  print "Writing $geomType distribution data using $binSize sized bins to $binFile...";
  &writeDistribution($DATA, $binSize, $binFile);
  print "Done\n";
}

sub writeDistribution {
	my ($data, $bsize, $bfile) = @_;
	my ($i, $vals, $min, $max, $nbins, $tot, $counts, $index, $maxI, $minI);

	@{ $vals } = sort numerically keys %{ $data };
	$max = $vals->[$#{ $vals }];
	$min = $vals->[0];
	$minI = 999999999;
	$maxI = -1;
	$tot = 0;
	for $i (@{ $vals }) {
		$index = sprintf("%.0f",($i - $min)/$bsize);
		$counts->{$index} += $data->{$i};
		$maxI = $index if ($index > $maxI);
		$minI = $index if ($index < $minI);
		$tot += $data->{$i};
	}
	open BFILE, "> $bfile" or die "ERROR: Cannot write to $bfile: $!\n";
	printf BFILE "#%10s%10s\n",$geomType,"count";
	for $i ($minI .. $maxI) {
		printf BFILE "%-10.3f%10.5f\n",$i*$bsize+$min,$counts->{$i}/$tot;
	}
	close BFILE;
}

sub getTrjGeom {
	my ($atoms, $box, $frameNum, $fileHandle) = @_;
	my ($tot, $geomVal, $i, $CoM, $cATOMS, $rec, $j);

	if ($trjType == 2) { #LAMMPS
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

	&RemoveMolLongBonds($atoms, $BONDS, $MOLS, $box);
	for $i (keys %{ $atoms } ) {
		delete $atoms->{$i} if (! keys %{ $atoms->{$i} }); 
	}
	
	for $i (keys %{ $ATOMS }) {
		next if (ref($ATOMS->{$i}) ne "HASH" or ! exists($ATOMS->{$i}{aList}));
		$cATOMS = GetAtmData($atoms, $ATOMS->{$i}{aList});
		$CoM = CoM($cATOMS, $box);
		for ("XCOORD", "YCOORD", "ZCOORD") {
			$atoms->{$i}{$_} = $CoM->{$_};
		}
	}
	&box2H($box);
	#undef $box;
	for $i (@{ $atmList }) {
		@{ $rec } = map { \%{ $atoms->{$_} } } @{ $i };
		$geomVal = $getGeomVal->(@{ $rec },$box);
		if(defined($sym_angle)) {
			$geomVal -= 360 if ($geomVal>$sym_angle);
			$geomVal += 360 if ($geomVal<-${sym_angle});
		}
		$DATA->{$geomVal}++;
		printf $OUTFILE "%10d %10.3f\n",$frameNum,$geomVal;
	}
}

sub getAtoms {
	my ($atoms, $bonds, $atmSelect, $box, $comOpt)  = @_;
	my ($i, $j, $k, $l, @ATMLIST, $selectName, $geomStr, $selection, $aselect, $idx);
	my ($cATOMS, $CoM, $tmp, $natoms, $n, $tmp1);

	$idx = scalar(keys %{ $atoms })+1;
	$n = 0;
	for $i (@{ $atmSelect }) {
		$i =~ s/^['|"](.*)['|"]$/$1/;
		$selection = BuildAtomSelectionString($i);
		$aselect = SelectAtoms($selection, $ATOMS);
		@{ $tmp } = keys %{ $aselect };
		$natoms = $#{ $tmp } + 1;
		die "ERROR: Atom selection $i is not valid!\n" if ($natoms = 0);
		if ($natoms == 1) {
			$tmp1->[$n][0] = $tmp->[0];
			$selectName .= "$tmp->[0] ($atoms->{$tmp->[0]}{FFTYPE}) ";
		} else {
			if ($comOpt) {
				$atoms->{$idx}{aList} = $aselect;
				$cATOMS = GetAtmData($atoms, $atoms->{$idx}{aList});
				$CoM = CoM($cATOMS, $box);
				for ("XCOORD", "YCOORD", "ZCOORD") {
					$atoms->{$idx}{$_} = $CoM->{$_};
				}
				$ATMLIST[0][$n] = $idx++;
				$selectName .= "$idx ";
			} else {
				@{ $tmp1->[$n] } = keys %{ $aselect };
			}
			$n++;
		}
	}
	$geomStr = $selectName . $geomType;
	&getValidVal($bonds, $tmp1->[0], $tmp1, \@ATMLIST, 0, $n-1, undef) 
		if (! $comOpt);
	die "ERROR: No valid atoms found according to selection!\n" 
		if ($#ATMLIST == -1);

	return (\@ATMLIST, $geomStr);
}	

sub getValidVal {
	my ($bonds, $iList, $alist, $vList, $i, $n, $pAtoms) = @_;
	my ($j, $k, $rec);

	if ($i == $n) {
		@{ $rec } = @{ $pAtoms };
		push @{ $vList }, $rec;
		return;
	}
	for $j (@{ $iList }) {
		$pAtoms->[$i] = $j;
		for $k (@{ $alist->[$i+1]}) {
			if (validBond($bonds->{$j},$k, $pAtoms)) {
				$pAtoms->[$i+1] = $k;
				&getValidVal($bonds, [ $k ], $alist, $vList, $i+1, $n, $pAtoms);
			}
		}
	}
}

sub validBond {
	my ($bList, $bAtom, $pList) = @_;
	my ($i, $vBond);

	$vBond = 0;
	for $i (@{ $bList }) {
		if($i == $bAtom) {
			$vBond = 1;
			last;
		}
	}
	$vBond = checkAlreadyVisited($pList, $bAtom) if($vBond);
	return $vBond;
}

sub checkAlreadyVisited {
	my ($alist, $i) = @_;
	my ($j, $valid);

	$valid = 1;
	for $j (@{ $alist }) {
		if($j==$i) {
			$valid =0;
			last;
		}
	}
	return $valid;
}

sub init {
	my (%OPTS, $atm, $tSel, $list, $i);

	getopt('btawsfdco',\%OPTS);
	for ("b", "a") {
		die &usage if (! exists($OPTS{$_}));
	}
	print "Initializing...";
	($structFile, $atm, $geomType, $trjFile, $tSel, $saveFile, $binSize, $comOpt, $sym_angle) = 
		($OPTS{b}, $OPTS{a}, $OPTS{f}, $OPTS{t}, $OPTS{s}, $OPTS{w}, $OPTS{d}, $OPTS{c}, $OPTS{a});
	FileTester($structFile);
	while ($atm =~ /\(([^)]+)\)|(\S+)/g) {
		$i = $1;
		$i = $2 if (! defined($1));
		if ($i !~ /^\d+$/) {
			push @{ $selection }, $i;
		} else {
			push @{ $selection }, "index==$i";
		}
	}
	die "ERROR: Expected valid atom seection for atoms! Got \"$atmList\"\n" if (! $selection);
	die "ERROR: Expected more than 1 valid atom number!\n" if ($#{ $selection } == 0);
	if ($#{ $selection } > 3) {
		print "warning: max number of atoms should be 4..taking first 4...";
		while ($#{ $selection } > 3) {
			pop @{ $selection };
		}
	}
	
	undef $sym_angle if (defined($sym_angle) and $sym_angle !~ /^\-?\d+\.?\d*$/);
	if ($#{ $selection } == 1) { #bond
		print "will compute bond..." if (defined($geomType) and $geomType !~ /^bond$/i);
		$geomType = "Bond";
		$getGeomVal = \&GetBondLength;
		undef $sym_angle;
	} elsif ($#{ $selection } == 2) { 
		print "will compute angle..." if (defined($geomType) and $geomType !~ /^angle$/i);
		$geomType = "Angle";
		$getGeomVal = \&GetAngle;
	} else {
		print "will compute torsion..." if (defined($geomType) and $geomType !~ /^(torsion|dihedral)$/i);
		$geomType = "Torsion";
		$getGeomVal = \&GetTorsion;
	}

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
		if (! defined($saveFile)) {
			$structFile =~ /^(\S+)/;
		   $saveFile = basename($1);
		   $saveFile =~ s/\.\w+$//;
		   $saveFile .= ".${geomType}.dat";
		}
	} else {
		undef $trjFile;
	}

	if(defined($binSize) and $binSize =~ /(\d+\.?\d*)/) {
		$binSize = $1;
		$binFile = basename($saveFile);
		$binFile =~ s/\.\w+$//;
		$binFile .= ".hist.dat";
	}

	$comOpt = 0 if (! defined($comOpt) or $comOpt =~ /0|no/i);
	$comOpt = 1 if ($comOpt =~ /1|yes/i);
	print "Done\n";
}
sub usage {
	my ($fTypeStr) = GetFileTypeStr;
	return <<DATA;
This script will analyze a structure (trajectory) file for various valence values.	
usage: $0 -b struct_file -a \"atom(s) selection\" -c (com_opt) -t (traj_file) -s (traj_selection) -w (save_file) -f [measurement_type] -d (bin_size) -o (symmetric_angle_option)
Required Arguments:
	-b struct_file: Location of input file
$fTypeStr	
	-a atom_selection: any valid field expression. E.g. index==3 will select atom index 3
		while (resname eq 'WAT') will select all the "WAT" residues. 
		NOTE1: if only number specified, will assume index=={num} for each entry.
		NOTE2: You can perform multiple types of measurements concurrently, enclosing in {} brackets.
		E.g. -a "{1 2 3 4} {sqrt((xcoord-10)**2+(ycoord-20)**2)<5}" will perform a torsion calculation
		on atoms 1 2 3 4 as well as a calculation on the atoms within a circle of radius 5, centered at 10,20
		NOTE3: If more than 4 atoms specified in any selection, only the first 4 are retained
Optional Arguments:
	-c com_opt: Use the center of mass of the atoms in the group selection. E.g. -a "'resnum==3' 'resnum==4'" -c 1
		will calculate the distance between the center of mass of residue 3 and residue 4. Default 0
	-t traj_file: Location of LAMMPS of AMBER coordinate trajectory file. When activated, will save the results
		from each timestep to a filename specified below. Default none
	-s traj_selection: The number of frames to use. Can be a single integer, or several integers in quotes
		To specify a range specify it as :Ita-b:c, which will take frames a - b, every c. Specify multiple 
		ranges or a combination of ranges and single frames by enclosing them in quotes. \"*\" for all frames.
		Default "*"
	-w save_file: Name of file to save results when analyzing a trajectory. Default {trajectory_name}.{measure_type}.dat
	-f measurement_type: Either bond|angle|torsion. If not specfied, will determine from number of atoms in selection
	-d bin_size: Size of interval to create histogram. If activated, will generate a histogram {save_file}.hist.dat.
		Default is none
DATA

}

#!/usr/bin/perl

$DB::deep = 50000000; # or more if necessary
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Storable qw(dclone);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use FileFormats qw(ParseStructFile createBGF addHeader createPQR PrintCoord
				 DeleteAtoms InsertMol PrintCoord GetSystemCharge GetBGFAtoms);
use General qw(FileTester LoadFFs IsDecimal CoM GetBondLength ReadFFs GetFileTypeStr);
use ManipAtoms qw(SplitAtomsByMol GetMols GroupAtomsByField SelectAtoms AddMolsToSelection
				 CenterSystem MakeFieldSequential BuildAtomSelectionString GetAtmData);
use BOX qw(GetBox);
use warnings;
use POSIX;

use constant TOLERANCE => 0.0001;

sub numerically { ($a<=>$b); }


my ($structFile, $saveName, $ffFiles, $solvSelect, $ionType, $ionConc);
my ($ATOMS, $BONDS, $HEADERS, $PARMS, $SOLVENT, $IONS, $BOX);
my ($soluCharge, $solvEng, $sMOLS, $numSolv, $randOpt);

$|++;
my ($start) = time();
$IONS = &init;
$PARMS = LoadFFs($ffFiles);
print "Getting Ion parameters...";
&getIonParms($IONS, $PARMS);
print "Done\nParsing Structure file $structFile...";
($ATOMS, $BONDS, $HEADERS) = ParseStructFile($structFile, 1);
&GetMols($ATOMS,$BONDS);
print "Done\nGetting system dimensions...";
$BOX = GetBox($ATOMS, $PARMS, $HEADERS);
print "Done\nObtaining Solvent atoms...";
&checkAtomTypes($ATOMS, $PARMS);
$SOLVENT = SelectAtoms($solvSelect, $ATOMS);
&AddMolsToSelection($SOLVENT, $ATOMS);
$sMOLS = GetMols($ATOMS, $BONDS, $SOLVENT);
$soluCharge = &updateSolvent($ATOMS, $SOLVENT);
$numSolv = scalar(keys %{ $sMOLS });
print "Found " . scalar(keys %{ $SOLVENT }) . " atoms in $numSolv molecules..Done\n";
&getIonOpts($IONS, $BOX, $ionConc, $soluCharge, $numSolv);
if (!$randOpt) {
  print "Calculating Electrostatic Potential...";
  $solvEng = calcElecPot($ATOMS, $BONDS, $BOX, $HEADERS);
  print "Done\nDetermining Ion placement based on potential...";
  &determineIonPlacement($IONS, $solvEng);
  print "Done\n";
} else {
  &placeRand($IONS, $ATOMS, $sMOLS);
} 
&placeIons($ATOMS, $BONDS, \%{ $IONS->{1} }, $sMOLS, $BOX);
if (exists($IONS->{2})) {
	&placeIons($ATOMS, $BONDS, \%{ $IONS->{2} }, $sMOLS);
}
print "Creating BGF file $saveName...";
&GroupAtomsByField($ATOMS,$BONDS,"RESNUM");
&MakeFieldSequential($ATOMS,"RESNUM");
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

$soluCharge = GetSystemCharge($ATOMS);
printf "Total System Charge: %.3f\n", $soluCharge;

my ($end) = time();

print "Elapsed: " . ($end - $start) . " seconds\n";

sub placeRand {
	my ($ions, $atoms, $solv) = @_;
	my ($list, $i, $t, $c, $m);

	srand (time ^ ($$ + ($$ << 11)));
	@{ $list } = keys %{ $solv };
	for $t (keys %{ $ions }) {
		for $i (1 .. $ions->{$t}{total}*5) { #list five times as many spots, in case the selected mols are already removed due to overlap
			$c = int(rand() * $#{ $list }) + 1;
			$m = splice @{ $list }, $c, 1;
			push @{ $ions->{$t}{list} }, $m;
		}
	}
}

sub placeIons {
	my ($atoms, $bonds, $ions, $solvMols, $box) = @_;
	my ($currMol, $solvAtoms, $currAtom, $i, $index, $rindex, $count); 
	my ($resId, $resName, $currIon, $resname, $soluCoM, @tmp, $solvGrid);

	$index = 0;
	$rindex = 0;
	for $i (keys %{ $atoms }) {
		$index = $i if ($i > $index); #get the maximum number of atoms
		$rindex = $atoms->{$i}{RESNUM} if ($atoms->{$i}{IS_SOLVENT} and $atoms->{$i}{RESNUM} > $rindex); #max resid
		#$atoms->{$i}{RESNUM} += $count if ($atoms->{$i}{IS_SOLVENT}); #??
	}

	#first place the solvent atoms on a grid if using custom ion (needed to ensure correct overlap removal)
	$solvGrid = &placeSolvOnGrid(\%{ $atoms }, $solvMols, $box) 
		if(exists($ions->{custom}));
	
	$rindex++;
	$count = 0;
	for $currMol (@{ $ions->{list} }) {
		$solvAtoms = ();
		@tmp = keys %{ $solvMols->{$currMol}{MEMBERS} };
		$currAtom = shift @tmp;
		next if (! exists($atoms->{$currAtom}) || ! exists($atoms->{$currAtom}{MOLECULE}{MEMBERS}));
		for $i (keys %{ $atoms->{$currAtom}{MOLECULE}{MEMBERS} }) {
			$solvAtoms->{$i} = 1;
		}
		$currAtom = $atoms->{$currAtom};
		$resId = $currAtom->{"RESNUM"};
		$resname = $currAtom->{"RESNAME"};
		$rindex++;
		$soluCoM = CoM(GetAtmData($atoms,$solvAtoms));
		if(! exists($ions->{custom})) {
			$currIon = updateIonFields($ions, $soluCoM, $index, $rindex);
		} else {
			$currIon = updateIonCoords($ions, $soluCoM, $index, $rindex);
			&updateOverlap($atoms, $currIon->{atoms}, \%{ $solvAtoms }, $solvMols, $box, $solvGrid);
		}
		&DeleteAtoms($solvAtoms, $atoms, $bonds);
		&InsertMol($currIon->{atoms}, $currIon->{bonds}, $atoms, $bonds, $index);
		$index += scalar(keys%{ $currIon->{atoms} });
		print "Replacing $resname (molecule # $currMol) centered at " . PrintCoord($soluCoM) . " with " . $currIon->{NAME} . "\n";
		$count++;
		last if ($count== $ions->{total});
	}

	die "ERROR: Could not place $ions->{total} ions... rerun scripts and try again!\n"
		if ($count < $ions->{total});

}

sub placeSolvOnGrid {
	my ($atoms, $molSolv, $box) = @_;
	my ($i, $j, $gspacing, $pos, $mgrid);

	for $i ("X","Y","Z") {
		$gspacing->{$i} = 2.5/$box->{$i}{len};
	}
	#first place solvent molecules on grid
	for $j (keys %{ $molSolv }) {
		for $i (keys %{ $molSolv->{$j}{MEMBERS}}) {
			$pos = ();
			%{ $pos } = (XCOORD=>$atoms->{$i}{XCOORD},YCOORD=>$atoms->{$i}{YCOORD},ZCOORD=>$atoms->{$i}{ZCOORD});
			$mgrid->{int($atoms->{$i}{FA}/$gspacing->{X})}{int($atoms->{$i}{FB}/$gspacing->{Y})}{int($atoms->{$i}{FC}/$gspacing->{Z})}{$i} = \%{ $pos };
			$atoms->{$i}{GRID}{X} = int($atoms->{$i}{FA}/$gspacing->{X});
			$atoms->{$i}{GRID}{Y} = int($atoms->{$i}{FB}/$gspacing->{Y});
			$atoms->{$i}{GRID}{Z} = int($atoms->{$i}{FC}/$gspacing->{Z});
		}
	}

	return $mgrid;
}

sub updateOverlap {
	my ($atoms, $currIon, $solvAtoms, $molSolv, $box, $mgrid) = @_;
	my ($i, $j, $k, $l, $gspacing, $pos, $molData, $soluGIndex, $solvGNeigh);

	for $i ("X","Y","Z") {
		$gspacing->{$i} = 2.5/$box->{$i}{len};
	}

	#now get the overlapping solvent atoms for each solute atom and update hash
	for $i (keys %{ $currIon }) {
		@{ $soluGIndex } = (int($atoms->{$i}{FA}/$gspacing->{X}),int($atoms->{$i}{FB}/$gspacing->{Y}),int($atoms->{$i}{FC}/$gspacing->{Z}));
		$solvGNeigh = getNeighbors($mgrid,@{ $soluGIndex });
		for $j (@{ $solvGNeigh }) {
			foreach $k (keys %{ $j }) {
				next if (exists($solvAtoms->{$k}) or ! exists($atoms->{$k}) or ! scalar(keys %{ $atoms->{$k} }));
				if(GetBondLength($atoms->{$i},$atoms->{$k},$box) < 2.0) {
					for $l (keys %{ $atoms->{$k}{MOLECULE}{MEMBERS} }) {
						$solvAtoms->{$l} = 1;
					}
				}
			}
		}
	}
}

sub getNeighbors {
	my ($mgrid, $ix, $iy, $iz) = @_;
	my ($cgrid, $i, $j, $k);

	for $i (-1 .. 1) {
		next if !(exists($mgrid->{$ix+$i}));
		for $j (-1 .. 1) {
			next if !(exists($mgrid->{$ix+$i}{$iy+$j}));
			for $k (-1 .. 1) {
				next if !(exists($mgrid->{$ix+$i}{$iy+$j}{$iz+$k}));
				push @{ $cgrid }, $mgrid->{$ix+$i}{$iy+$j}{$iz+$k};
			}
		}
	}
	return $cgrid;
}

sub updateIonCoords {
	my ($ion, $soluCoM, $startIndex, $resId) = @_;
	my ($ionCoM, $i, $offset, $newIon, $j);

	$ionCoM = CoM($ion->{atoms});
	for $i ("XCOORD", "YCOORD", "ZCOORD") {
		$offset->{$i} = $soluCoM->{$i} - $ionCoM->{$i};
	}

	$newIon->{atoms} = dclone($ion->{atoms});
	$newIon->{bonds} = dclone($ion->{bonds});
	for $i (values %{ $newIon->{atoms} }) {
		$i->{RESNUM} = $resId;
		$i->{INDEX} += $startIndex;
		for $j ("XCOORD", "YCOORD", "ZCOORD") {
			$i->{$j} += $offset->{$j};
		}
	}
	$newIon->{NAME} = "custom ion";
	return $newIon;
}

sub updateIonFields {
	my ($ionParms, $ionPos, $ionIndex, $ionResID) = @_;
	my (%ION, $newIon);
	
	%ION = %{ $ionPos };
	$ION{INDEX} = $ionIndex;
	$ION{ATMNAME} = $ionParms->{atoms}{1}{ATMNAME};
	$ION{RESNAME} = $ionParms->{atoms}{1}{ATMNAME};
	$ION{RESNUM} = $ionResID;
	$ION{FFTYPE} = $ionParms->{atoms}{1}{FFTYPE};
	$ION{NUMBONDS} = 0;
	$ION{LONEPAIRS} = 0;
	$ION{CHARGE} = $ionParms->{charge};
	$ION{RADII} = $ionParms->{atoms}{1}{RADII};
	$ION{LABEL} = "HETATM";
	$newIon->{atoms}{1} = \%ION;
	$newIon->{bonds}{1} = ();
	$newIon->{NAME} = $ionParms->{atoms}{1}{ATMNAME};
	return $newIon;
}

sub determineIonPlacement {
	my ($ions, $solvEnergies) = @_;
	my (@tmp, %CHARGE, $i);

	$ions->{1}{placed} = 0;
	$ions->{2}{placed} = 0;

	for $i (keys %{ $solvEnergies }) {
		$CHARGE{$solvEnergies->{$i}} = $i;
	}

	if ($ions->{1}{charge} > 0) {
		@tmp = sort numerically keys %CHARGE;
	} else {
		@tmp = reverse sort numerically keys %CHARGE;
	}

	while(@tmp) {
		last if ($ions->{1}{placed} >= $ions->{1}{total});
		$i = shift @tmp;
		push @{ $ions->{1}{list} }, $CHARGE{$i};
		$ions->{1}{placed}++;
	}
	if (keys %{ $ions->{2}{atoms}{1} }) {
		while (@tmp) {
			last if ($ions->{2}{placed} >= $ions->{2}{total});
			$i = pop @tmp;
			push @{ $ions->{2}{list} }, $CHARGE{$i};
			$ions->{2}{placed}++;
		}
	} else {
		delete $ions->{2};
	}
}

sub computeElecEnergyPerAtom {
	my ($atoms, $pqrFile) = @_;
	my ($apbsCmd, $isValid, %ENERGIES, $molID);
	
	$apbsCmd = "/home/tpascal/codes/bin/coulomb -e ./_tmp.pqr";
	open COULCMD, "$apbsCmd |" or die "ERROR: Cannot run cmd $apbsCmd: $!\n";
	while (<COULCMD>) {
	chomp;
	if ($_ =~ /Atom\s+(\d+):\s+Energy\s+\=\s+(\-?\d+\.\d+E?.?\d*)/) {
		if (exists($atoms->{$1}{IS_SOLVENT})) {
			$molID = ${ $atoms->{$1}{MOLECULEID} };
			$ENERGIES{$molID} += $2;
			$atoms->{$1}{ENERGY} = $2;
		}
	}
	}
	die "ERROR: No valid energies obtained from $apbsCmd!\n" if (! %ENERGIES);

	system("rm -fr _tmp.pqr io.mc");
	return \%ENERGIES;
}

sub calcElecPot {
	my ($atoms, $bonds, $box, $headers) = @_;
	my ($pqrAtoms, $center, $sEng, $i);

	$pqrAtoms = $atoms;
	$center = CoM($atoms);
	for $i ("X", "Y", "Z") {
		$center->{"${i}COORD"} = -($box->{$i}{len}/2 - $box->{$i}{lo} - $center->{"${i}COORD"});
	}
	&CenterSystem($pqrAtoms, $center, ());
	&createPQR($pqrAtoms, $bonds, "_tmp.pqr",$headers);
	$sEng = computeElecEnergyPerAtom($atoms, "_tmp.pqr");
	return $sEng;
}

sub getIonParms {
	my ($ionInfo, $parms) = @_;
	my ($counter, $element, $radii, $i, $count, $valid);

	if(exists($ionInfo->{1}{custom})) {
		delete $ionInfo->{2};
		return 0;
	}	

	for $counter (keys %{ $parms->{"ATOMTYPES"} }) {
		$element = $parms->{"ATOMTYPES"}{$counter}{"ATOM"};
		if ($element =~ /^$ionInfo->{1}{atoms}{1}{"ATMNAME"}/) {
			$ionInfo->{1}{atoms}{1}{"FFTYPE"} = $counter;
			&addRadii(\%{ $ionInfo->{1}{atoms}{1} }, $parms);
		}
		if (exists($ionInfo->{2}) and $element =~ /^$ionInfo->{2}{atoms}{1}{"ATMNAME"}/) {
			$ionInfo->{2}{atoms}{1}{"FFTYPE"} = $counter;
			&addRadii(\%{ $ionInfo->{2}{atoms}{1} }, $parms);
		}
	}

	if (! exists($ionInfo->{1}{atoms}{1}{"FFTYPE"})) {
		die "ERROR: The " . $ionInfo->{1}{atoms}{1}{"ATMNAME"} .
		" ion does not have any parameters in the provided forcefield\n";
	} elsif (exists($ionInfo->{2}) and ! exists($ionInfo->{2}{atoms}{1}{"FFTYPE"})) {
		die "ERROR: The " . $ionInfo->{2}{atoms}{1}{"ATMNAME"} .
		" ion does not have any parameters in the provided forcefield\n";
	}
	delete $ionInfo->{2} if (! keys %{ $ionInfo->{2} });
}

sub getIonOpts {
	my ($ionInfo, $BOX, $num_ions, $charge, $watCount) = @_;
	my ($mol_vol, $dim, $molar);

	$mol_vol = 1;
	for $dim ("X", "Y", "Z") {
		$mol_vol = $mol_vol * ($BOX->{$dim}{"hi"} - $BOX->{$dim}{"lo"});
	}
	printf "CELL VOLUME %.3G A^3", $mol_vol;
	if ($watCount) {
		$mol_vol = 18 * $watCount/.6023; # vol in cm^3
		printf "\nADJUSTED CELL MOLAR VOLUME based on 1g/cm3 WATER density %.3G cm^3", $mol_vol;
	}
	$mol_vol *= 6.0221415E-7;  #Angstroms^3 to meter^3
	$mol_vol *= 1000; #meter^3 to liter
	printf "(%.3G liter)\n", $mol_vol;
	if(! $ionInfo->{1}{custom}) {
		$ionInfo->{1}{charge} = $ionInfo->{1}{"total"} = 0;
		$ionInfo->{2}{charge} = $ionInfo->{2}{"total"} = 0;
	}

	if ((abs($charge - sprintf("%.0f", $charge))) > TOLERANCE) {
		print "WARNING: System has non-integer net charge of $charge!\n";
	}
	$charge = sprintf("%.0f", $charge);

	if (! $ionInfo->{1}{custom} and $ionInfo->{1}{atoms}{1}{ATMNAME} =~ /Na|K|Li|Rb|Cs/i) { # +1 charge
		$ionInfo->{1}{charge} = "+1";
	} elsif (! $ionInfo->{1}{custom} and $ionInfo->{1}{atoms}{1}{ATMNAME} =~ /cl|f|br|i/i) {
		$ionInfo->{1}{charge} = -1;
	} elsif(! $ionInfo->{1}{custom}) {
		$ionInfo->{1}{charge} = "+2";
	}

	if (IsDecimal($num_ions) and $num_ions > 0) { #placing a concentration
		$ionInfo->{1}{"total"} = sprintf("%.0f",($mol_vol * $num_ions));
		$molar = $num_ions;
	} elsif ($num_ions > 0 and ! IsDecimal($num_ions)) {
		$ionInfo->{1}{"total"} = $num_ions;
		$molar = $ionInfo->{1}{"total"} / ($mol_vol);
	} elsif ($num_ions == 0 and ! IsDecimal($num_ions)) { # neutralize
		if ($ionInfo->{1}{charge} != 1) {
			if (($charge % 2) == 1) {
				print "WARNING: Charge is not even but divalent cation chosen! System will have net charge\n";
			}
	}
	die "ERROR: Ion and system have same charge but attempting to neutralize!\n" if(($charge/$ionInfo->{1}{charge})>0);
		$ionInfo->{1}{total} = sprintf("%.0f",-$charge/$ionInfo->{1}{charge});
		$molar = $ionInfo->{1}{total} / ($mol_vol);
	} else {
		die "ERROR: Invalid option for ion/salt type and concentration\n";
	}

	if (exists($ionInfo->{2}{atoms}{1}{ATMNAME})) {
		$ionInfo->{2}{charge} = -1;
		$ionInfo->{2}{total} = $ionInfo->{1}{total} * $ionInfo->{1}{charge};
		delete $ionInfo->{2} if ($ionInfo->{2}{total} < 0);
	}

	if (! exists($ionInfo->{2}{atoms}{1}{ATMNAME}) and ! exists($ionInfo->{1}{custom})) {
		printf "Placing " . $ionInfo->{1}{total} . " molecule(s) of " . $ionInfo->{1}{atoms}{1}{ATMNAME} .
			$ionInfo->{1}{charge} . " ion (%.3G molar)\n", $molar;
	} elsif (exists($ionInfo->{1}{custom})) {
	printf "Placing " . $ionInfo->{1}{total} . " molecule(s) of custom ion (" . $ionInfo->{1}{charge} . "e) " .
		"(%.3G molar)\n", $molar;
	} else {
		printf "Placing " . $ionInfo->{1}{total} . " molecule(s) of " . $ionInfo->{1}{atoms}{1}{ATMNAME} .
			$ionInfo->{2}{atoms}{1}{ATMNAME} . " salt (%.3G molar)\n", $molar;
	}
	delete $ionInfo->{2} if (! exists($ionInfo->{2}{atoms}{1}{ATMNAME}));
}

sub updateSolvent {
	my ($atoms, $solvList) = @_;
	my ($i, $charge);

	$charge = 0;
	for $i (keys %{ $atoms }) {
		if (exists($solvList->{$i})) {
			$atoms->{$i}{IS_SOLVENT} = 1;
		} else {
			$charge += $atoms->{$i}{CHARGE};
		}
	}

	return $charge;
}

sub checkAtomTypes {
	my ($atoms, $parms) = @_;
	my ($i, $ffType);

	for $i (keys %{ $atoms }) {
	$ffType = $atoms->{$i}{FFTYPE};
	die "ERROR: Force field type $ffType not found in forcefield(s). Aborting\n"
		if (! exists($parms->{ATOMTYPES}{$ffType}));
		&addRadii($atoms->{$i}, $parms);
	}
}

sub addRadii {
	my ($atom, $parms) = @_;
	my ($i, $ffType);

	$ffType = $atom->{FFTYPE};
	$atom->{RADII} = $parms->{VDW}{$ffType}{$ffType}{1}{VALS}[1];
}

sub init {
	my (%OPTS, $forceFields, $ffType, $select, @tmp, $ION);
	getopt('bfwsinr', \%OPTS);
	
	for ("i", "b", "f") {
		if (! exists($OPTS{$_})) {
			&usage;
			die "\n";
		}
	}
	print "Initializing...";
	($structFile, $forceFields, $saveName, $select, $ionType, $ionConc, $randOpt) = 
	($OPTS{b}, $OPTS{f}, $OPTS{w}, $OPTS{s}, $OPTS{i}, $OPTS{n}, $OPTS{r});
	
	if (! defined($saveName)) {
		$structFile =~ /^\s+(\S+)/;
		$saveName = basename($1);
		$saveName =~ s/\.\w+$//;
		$saveName .= "_ion.bgf";
	}

	($ffFiles, $ffType) = ReadFFs($forceFields);

	$select = "resname eq 'WAT' and molsize == 3" if (! defined($select));
	$solvSelect = BuildAtomSelectionString($select);

	die "ERROR: Available ion types are - Na, F, Cl, I, Br, Mg, K, Ca, Cu " .
		"and the corresponding salts (e.g. KCl)\n"
		if ($ionType !~ /F|Na|I|Br|Cl|Mg|K|Ca|Cu|Rb|Cs/i and ! -e $ionType);

	$ionConc = 0 if (! defined($ionConc));
	if ($ionType !~ /^(Na|Li|Mg|K|Ca|Cl|F|Cl|Br|I|Cu|Rb|Cs)(\w*)$/i and ! -e $ionType and $ionType !~ /\s+/) {
		print "ERROR: Invalid ion/salt type $ionType\n";
		&usage();
		die "\n";
	} elsif( ! -e $ionType and $ionType !~ /\s+/) {
		$ION->{1}{atoms}{1}{"ATMNAME"} = $1;
		if ($2 and $1 ne $2) {
			$ION->{2}{atoms}{1}{"ATMNAME"} = $2;
		}
	} else {
		$ION->{1} = getIonType($ionType);
	}
	if ($ionConc eq "0" and length($ionType) > 2 and ! -e $ionType) {
		die "ERROR: Cannot use a salt to neutralize the system\n";
	}
	$randOpt = 0 if (! defined($randOpt) or $randOpt !~ /1|yes/i);
	$randOpt = 1 if ($randOpt =~ /1|yes/i);
	print "Done\n";

	return $ION;
}

sub getIonType {
	my ($ionFile) = $_[0];
	my ($ionAtoms, $ionBonds, $charge, $i, $data);

	($ionAtoms, $ionBonds) = ParseStructFile($ionFile,0,1);
	$charge = 0.0;
	for $i (keys %{ $ionAtoms }) {
		$charge += $ionAtoms->{$i}{CHARGE};
	}
 
	$charge += 0;
	die "ERROR: Charge read in from $ionFile ($charge) is not interger!\n" if(!isint($charge));

	$data->{atoms} = $ionAtoms;
	$data->{bonds} = $ionBonds;
	$data->{charge} = $charge;
	$data->{total} = 0;
	$data->{custom} = 1;
	return $data;   
}

sub usage {
	my ($fTypeStr) = GetFileTypeStr;
	print <<USAGE
usage: $0 -b struct_file -f force field(s) -i ion type -s [solvent selection] -n [\# ions] -w [save_name] -r [rand_opt=0]
options:
	struct_file: the name of the Structure file.
$fTypeStr	
	force field: 1 (or more if in quotes) cerius2/mpsim/cmdf formatted forcefile for the bgf file
	ion_type: Na, F, Cl, Br, I, Mg, K, Ca and the combination salts (e.g. NaCl) or add custom by giving bgf file
	solvent selection: will select atoms based on selection string. Fields can be used and combined
			For e.g. -s "resname eq 'WAT'" will select all residues named 'WAT' and
			-s "resname eq 'WAT' and zcoord < 50 and zcoord > 25" will only select waters 
			between 25 and 50 angstroms in the z direction
	[num_ions] : number of ions
				 0 - neutralize the system (default) - ions only
				 X - X number of ions/salt
				 X.x - X.x molar ion/salt
	[rand opt] : random options. Will select solvent atoms at random instead of based on best electrostic
		 potential
	[save_name]: name of output BGF file (optional)
USAGE
}

sub isint{
	my $val = shift;
	$val =~ s/^\-?\d+/0/;
	$val  = sprintf("%.5g",$val);
	return 1 if($val<0.0001);
	return 0;
}

#!/usr/bin/perl -w

$DB::deep = 500000; # or more if necessary
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
no warnings 'recursion';

use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use Storable qw(dclone);

use FileFormats qw(ParseStructFile insertHeaderRemark addHeader createBGF);
use General qw(PrintProgress GetTime GetFileTypeStr GetTime CrossProduct GetBondLength 
				Normalize ReadFFs CoM ShuffleArray CombineMols Normalize);
use Superimpose; 
use Math::MatrixReal;
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);
use BOX qw(GetBox CreateGrid GetGridDims GetNeighbours PlaceAtomsOnGrid);

sub numerically { ($a<=>$b) }

my ($FILES, $selStr, $num, $randomize, $functAttachID, $qeqOpt, $exclude, $minOpt, $uFieldsOpt); 
my ($structure, $funct, $selection, $aList, $atomQ, $GRID, $bbox);

$|++;
&init;
print "Gettting structure data from $FILES->{STRUCT}{LOC}...";
($structure->{ATOMS}, $structure->{BONDS}, $structure->{HEADERS}) = ParseStructFile($FILES->{STRUCT}{LOC}, 1);
&GetMols($structure->{ATOMS}, $structure->{BONDS});
print "Done\nSelecting structure atoms according to \"$selStr\"...";
$aList = SelectAtoms($selection, $structure->{ATOMS});
die "ERROR: No valid atoms found\n" if (! keys %{ $aList });
&removeExclusions($structure, \%{ $aList }, $exclude) if ($exclude);
print "found " . scalar(keys %{ $aList }) . " atoms...Done\n"; 
print "Getting functional group data from file $FILES->{FUNCT}{LOC}...";
($funct->{ATOMS}, $funct->{BONDS}, $funct->{HEADERS}) = ParseStructFile($FILES->{FUNCT}{LOC}, 0);
($atomQ, $num, $aList) = getFunctOpts($funct->{ATOMS}, $aList, $num, $functAttachID, $qeqOpt, $randomize, $structure->{ATOMS});
print "...will add $num function groups..Done\nSetting up Grids...";
($GRID, $bbox) = setupGRIDS($structure, $funct);
print "Done\n";
&insertFuncGrp($structure, $funct, $aList, $functAttachID, $atomQ, $uFieldsOpt, $GRID, $bbox);
print "Creating BGF file $FILES->{SAVE}{LOC}...";
&insertHeaderRemark($structure->{HEADERS}, "REMARK selected atoms according to $selStr");
&addHeader($structure->{ATOMS},$structure->{HEADERS});
&createBGF($structure->{ATOMS}, $structure->{BONDS}, $FILES->{SAVE}{LOC});
print "Done\n";
&minimizeStructure($FILES->{SAVE}{LOC}, $minOpt, $qeqOpt, $aList) 
	if ($minOpt and $FILES->{SAVE}{TYPE} eq "bgf");

sub insertFuncGrp {
	my ($struct, $fgrp, $selAtoms, $fgrpID, $qOpt, $uOpt, $grid, $box) = @_;
	my ($i, $j, $gAtoms, $gBonds, $tot, $field, $pStr, $start, $nSAtoms, $tmp); 
	my ($p, $strLen, $count, $sV, $gV, $sVcom, $gVcom, $nGatoms, $sAtoms);

	$pStr = "Inserting functional group...";
	print "${pStr}Calculating time remaining...\r";
	$start = time();
	$tot = scalar(keys %{ $selAtoms });
	$nSAtoms = scalar(keys %{ $struct->{ATOMS} });
	$nGatoms = scalar(keys %{ $fgrp->{ATOMS} });
	$count=0;

	#first get the center of mass of the structure
	$sVcom = CoM($struct->{ATOMS});
	#now get the dipole vector of the connecting group atom to its center of mass
	$gVcom = CoM($fgrp->{ATOMS});
	for $i ("XCOORD", "YCOORD", "ZCOORD") {
		$gV->{1}{$i} = $fgrp->{ATOMS}{$fgrpID}{$i} - $gVcom->{$i};
		$gV->{0}{$i} = 0;
	}

	for $i (keys %{ $selAtoms }) {
		$struct->{ATOMS}{$i}{CHARGE} += $qOpt; #charge update

		if($uOpt) { #fields update
			for $j (values %{ $gAtoms }) {
				for $field ("RESNUM", "RESNAME", "LABEL", "CHAIN") {
					$j->{$field} = $struct->{ATOMS}{$i}{$field};
				}
			}
		}

		for $p ("1", "-1") { #try both phases
			#make local copies of atom arrays
			$gAtoms = dclone($fgrp->{ATOMS});
			$gBonds = dclone($fgrp->{BONDS});
			$sAtoms = dclone($struct->{ATOMS});

			#function group atom positions update
			for $j ("XCOORD", "YCOORD", "ZCOORD") {
				#functionalized atom position vector points away from CoM 
				$sV->{1}{$j} = ($sVcom->{$j}-$struct->{ATOMS}{$i}{$j})*$p;
				$sV->{0}{$j} = 0;
			}
			&positionGrpAtoms($struct->{ATOMS}{$i}, $gAtoms, $fgrpID, $sV, $gV, $fgrp->{ATOMS}{$fgrpID});

			#merge structure and functional group
			($tmp->{$p}{ATOMS}, $tmp->{$p}{BONDS}) = CombineMols($sAtoms, $gAtoms, $struct->{BONDS}, $gBonds, 0);
			$tmp->{$p}{overlaps} = findOverlaps($gAtoms, $grid, $box);
		}

		#select the rotation phase with least overlaps
		if($tmp->{"-1"}{overlaps} < $tmp->{1}{overlaps}) {
			$struct->{ATOMS} = $tmp->{"-1"}{ATOMS};
			$struct->{BONDS} = $tmp->{"-1"}{BONDS};
		} else {
			$struct->{ATOMS} = $tmp->{1}{ATOMS};
			$struct->{BONDS} = $tmp->{1}{BONDS};
		}
		$tmp = (); #free up some memory

		#now create bonds
		push @{ $struct->{BONDS}{$i} }, $nSAtoms + $fgrpID;
		push @{ $struct->{BONDS}{$nSAtoms + $fgrpID} }, $i;
		$nSAtoms += $nGatoms; #total number of atoms in structure (dynamically updated)
		#show progress
		$count++;
		$strLen = PrintProgress($count, $tot, $start, $pStr);
	}
	$tot = GetTime(time() - $start);
	printf "${pStr}%-${strLen}s\n", "${tot}s elapsed...Done";
}

sub positionGrpAtoms {
	my ($sAtom, $fgAtoms, $fgAtomID, $sV, $gV, $gVref) = @_;
	my ($i, $j, $offset, $vec, $rotMat);

	#first we place the center of mass of the functional group at the origin
	for $i ("XCOORD", "YCOORD", "ZCOORD") {
		$offset = $gVref->{$i};
		for $j (values %{ $fgAtoms }) {
			$j->{$i} -= $offset;
		}   
	}

	#then use Kabash to position functional group by rotating supplied vectors
	&Kabash::SuperimposeAtoms($sV, $gV, $fgAtoms); 

	#normalize the supplied vector
	for $i ("XCOORD", "YCOORD", "ZCOORD") {
		$vec->{$i} = $sV->{1}{$i}-$sV->{0}{$i};
	}
	&Normalize($vec);

	#now position to functioal group near the functionalized atom, assume a 1.5A bond length
	for $i ("XCOORD", "YCOORD", "ZCOORD") {
		$offset = $sAtom->{$i} - $fgAtoms->{$fgAtomID}{$i} - $vec->{$i}*1.5;
		for $j (values %{ $fgAtoms }) {
			$j->{$i} += $offset;
		}   
	}
}

sub findOverlaps {
	my ($atoms, $grid, $box) = @_;
	my ($i, $c, $j, $tmpgrid, $nCells, $N);

	#place the functional group on the grid for faster searching
	$tmpgrid = dclone($grid);
	&PlaceAtomsOnGrid($atoms, $tmpgrid, $box, 2.5, 1);

	$N = 0; #number of overlaps
	for $i (values %{ $atoms }) {
		$nCells = GetNeighbours($tmpgrid, $i->{CELL}, 1);
		for $c (@{ $nCells }) {
			next if (!exists($c->{ATOMS}));
			for $j (@{ $c->{ATOMS}}) {
				$N++ if(GetBondLength($i, $j, $box) < 2);
			}
		}
	}
	return $N;
}

sub setupGRIDS {
	my ($s, $f) = @_;
	my ($i, $bl, $bbox, $grids);

	#first remove the BONDATOM array from the atoms (makes it difficut to dclone)
	for $i (keys %{ $s->{ATOMS} }) {
		delete $s->{ATOMS}{$i}{BONDATOM};
	}

	#get largest size of functional group
	$bbox = GetBox($f->{ATOMS}, undef, undef);
	$bl = sqrt(($bbox->{X}{hi}-$bbox->{X}{hi})**2+($bbox->{Y}{hi}-$bbox->{Y}{hi})**2+($bbox->{Z}{hi}-$bbox->{Z}{hi})**2);
	$bbox = GetBox($s->{ATOMS}, undef, undef, $bl);
	($grids, undef, undef) = CreateGrid($s->{ATOMS}, 3, $bbox, 2.5, 0);

	return ($grids, $bbox);
}

sub getFunctOpts {
	my ($atoms, $selAtoms, $n, $groupAttachID, $chargeOpt, $randOpt, $struct) = @_;
	my ($i, $tmp);

	die "ERROR: Cannot find  functional group atom id $groupAttachID\n"
		if(! exists($atoms->{$groupAttachID}));

	if ($chargeOpt) {
		$chargeOpt = 0;
		for $i (keys %{ $atoms }) {
			$chargeOpt += $atoms->{$i}{CHARGE};
		}
		$chargeOpt *= -1;
	}
	if($n eq "all") {
		for $i (keys %{ $struct }) {
			$selAtoms->{$i} = 1;
		}
	} elsif (scalar(keys %{ $selAtoms }) > $n) {
		@{ $tmp } = keys %{ $selAtoms };
		&ShuffleArray($tmp);
		$selAtoms = ();
		for $i (0 .. ($n-1)) {
			$selAtoms->{ $tmp->[$i] } = 1;
		}
	}
	$n = scalar(keys%{ $selAtoms });
	return ($chargeOpt, $n, $selAtoms);
}

sub removeExclusions {
	my ($struct, $atomList, $level) = @_;
	my ($i, $exclude, $pStr, $tot, $tmp, $j, $include);

	for $i (1 .. $level) {
		$pStr .= "1 - " . ($i+1) . ", ";
	}
	chop $pStr;
	chop $pStr;
	print "excluding $pStr...";

	@{ $tmp } = sort numerically keys %{ $atomList }; #index for all atoms currently selected
	$exclude = (); #excluded atoms
	$include = (); #included atoms
	&getAtomIDs($struct, $level, \%{ $exclude }, \%{ $include }, $tmp->[0], $level);

	for $i (keys %{ $atomList }) {
		delete $atomList->{$i} if (! exists($include->{$i})); #remove all the excluded atoms
	}
}

#recursively search through the bond list and exclude/include atoms
sub getAtomIDs {
	my ($struct, $max_depth, $excluded, $included, $atomID, $curr_depth) = @_;
	my ($i, $bList);

	#first, test to see if current depth is zero, in which case exclude atom and reset depth_counter
	if ($curr_depth == $max_depth) {
		$included->{ $atomID } = 1;
		$curr_depth--;
	} elsif ($curr_depth == 0) {
		$excluded->{ $atomID } = 1;
		$curr_depth = $max_depth;
	} else {
		$excluded->{ $atomID } = 1;
		$curr_depth--;
	}
	$bList = (); # empty the bList hash (stores the bonds of the current atom)
	&getBlist($struct->{ATOMS}, $atomID, $struct->{BONDS}, \%{ $bList }, 1, 1); #get the bonds of current atom
	for $i (keys %{ $bList }) {
		next if (exists($excluded->{$i}) or exists($included->{$i}));
		&getAtomIDs($struct, $max_depth, \%{ $excluded }, \%{ $included }, $i, $curr_depth);
	}
}

sub getBlist {
	my ($atoms, $atomID, $bonds, $blist, $currDepth, $maxDepth) = @_;
	my ($i);

	for $i (@{ $bonds->{$atomID} }) {
		$blist->{ $i } = 1;
		if($currDepth < $maxDepth) {
			$currDepth++;
			&getBlist($atoms, $i, $bonds, \%{ $blist }, $currDepth, $maxDepth);
			$currDepth--;
		}
	}
}

sub minimizeStructure {
	my ($file, $ffList, $qOpt, $newAtoms) = @_;
	my ($execCmd, $prefix, $tot, $outStr, $inStr, $freeze, $i, $tmp, $pqeqOpt, $dumpStr);

	$pqeqOpt = "";
	$pqeqOpt = "-o 'pqeq' -q $Bin/../ff/pqeq1.par" if ($qOpt == 2);
	print "Minimizing structure...Creating LAMMPS files\r";
	$prefix = basename($file);
	$prefix =~ s/\.\w+$//;
	$tmp = "";
	for $i (keys %{ $newAtoms }) {
		$tmp .= "index==" . $i . " or ";
	}
	$tmp =~ s/ or $//;
	$freeze = "group		freeze molecule 1\nfix				restraint freeze spring/self 10.0";
	$execCmd = "$Bin/syncBonds.pl -b $file -s $file";
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
	$execCmd = "$Bin/modifyAtomData.pl -s $file -a \"index>0\" -f \"RESNUM:1\" -w ${prefix}_tmp.bgf";
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
	$execCmd = "$Bin/modifyAtomData.pl -s ${prefix}_tmp.bgf -a \"$tmp\" -f \"RESNUM:2\" -w ${prefix}_tmp.bgf";
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
	$execCmd = "$Bin/createLammpsInput.pl -b ${prefix}_tmp.bgf -f \"$ffList\" -s $prefix $pqeqOpt";
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
	
	print "Minimizing structure...Updating LAMMPS files\r";
	# get the atom types of the inputed functional groups

	#update input script
	open LMPINFILE, "in.${prefix}_singlepoint" or die "ERROR: Cannot read from to in.${prefix}_singlepoint: $!\n";
	while (<LMPINFILE>) {
		chomp;
		$inStr = $_;
		if (! $qOpt) {
			$inStr =~ s/kspace_style.*/kspace_style		none/g;
			$inStr =~ s/lj\/charmm\/coul\/long\/opt/lj\/charmm\/coul\/charmm/g;
			$inStr =~ s/coul\/long/coul\/charmm/g;
		}
		$outStr .= "$inStr\n";
	}
	close LMPINFILE; 

	$dumpStr = "dump			1 all custom 25 \${sname}_min.lammpstrj id type xu yu zu";
	$dumpStr = "" if($qOpt == 2); #PQEq fix
	open LMPINFILE, "> in.${prefix}_singlepoint" or die "ERROR: Cannot write to in.${prefix}_singlepoint: $!\n";
	print LMPINFILE <<DATA;
$outStr
$freeze
thermo		  10 
$dumpStr
min_style	   sd
#fix				boxrel all box/relax x 1.0 y 1.0 couple none vmax 0.001
minimize		1.0e-4 1.0e-4 500 5000
#unfix				boxrel
min_style	   cg
minimize		1.0e-4 1.0e-4 500 5000
min_style		hftn
min_modify		line quadratic
minimize		1.0e-4 1.0e-4 500 5000
undump		  1
DATA

	close LMPINFILE;

	# run LAMMPS
	print "Minimizing structure...Running LAMMPS minimization\r";
	#$execCmd = "/home/caltech/openmpi-1.3.3/bin/mpirun -np 1 /home/caltech/lammps-stable/src/lmp_parallel -in in.${prefix}_singlepoint -screen none -log ${prefix}_log.lammps";
	$execCmd = "/home/tpascal/codes/bin/lmp_expanse -in in.${prefix}_singlepoint -screen none -log ${prefix}_log.lammps";
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));

	#update coordinates
	print "Minimizing structure...updating coordinates\r";
	if ($qOpt == 2) {
		$tot = `grep -i -c timestep ${prefix}.lammps`;
		chop $tot;
		$execCmd = "$Bin/convertLammpsTrj.pl -b $file -o bgf -s ${prefix}.min.bgf -t ${tot} -l ${prefix}.lammps";
	} else {
		$tot = `grep -i -c timestep ${prefix}_min.lammpstrj`;
		chop $tot;
		$execCmd = "$Bin/convertLammpsTrj.pl -b $file -o bgf -s ${prefix}.min.bgf -t ${tot} -l ${prefix}_min.lammpstrj";
	}
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));

	system("rm -fr in.${prefix} in.${prefix}_singlepoint data.${prefix} tmp"); 
	system("rm -fr ${prefix}.lammps.pbs ${prefix}_log.lammps ${prefix}_tmp.bgf ${prefix}_sync.bgf");
	system("rm -fr ${prefix}_min.lammpstrj _byte_offset_${prefix}* ${prefix}.lammps");
	print "Minimizing structure...created ${prefix}.min.bgf...Done\n";
}

sub init {
	my (%OPTS, $fileType, $FFS, $usageStr, $ext);
	
	getopt('fsatgcmqxnru', \%OPTS);
	$usageStr = &usage();
	for ("f", "a", "g", "c") {
		die "$usageStr" if (! exists($OPTS{$_}));
	}

	print "Initializing...";
	($FILES->{STRUCT}{LOC}, $FILES->{FUNCT}{LOC}, $FILES->{SAVE}{LOC}) = ($OPTS{f}, $OPTS{g}, $OPTS{s});
	($fileType,  $selStr, $functAttachID) = ($OPTS{t}, $OPTS{a}, $OPTS{c});
	($minOpt, $qeqOpt, $exclude, $uFieldsOpt) = ($OPTS{m}, $OPTS{q}, $OPTS{x}, $OPTS{u});
	($num, $randomize) = ($OPTS{n}, $OPTS{r});
	
	die "ERROR: Expected integer for funct group attach atom id . Got \"$functAttachID\"\n"
		if ($functAttachID !~ /^\d+$/);

	if (! defined($FILES->{SAVE}{LOC})) {
		$FILES->{SAVE}{LOC} = basename($FILES->{STRUCT}{LOC});
		$ext = "bgf";
		if ($FILES->{SAVE}{LOC} =~ /^(.+)\.(\w+)$/) {
			$FILES->{SAVE}{LOC} = $1;
			$ext = $2;
		}
		$FILES->{SAVE}{LOC} =~ s/\.\w+$//;
		$FILES->{SAVE}{LOC} .= basename($FILES->{FUNCT}{LOC});
		$FILES->{SAVE}{LOC} =~ s/\.\w+$//;
		$FILES->{SAVE}{LOC} .= ".${ext}";
	}

	for ("STRUCT", "FUNCT", "SAVE") {
		if (defined($fileType)) {
			$FILES->{$_}{TYPE} = $fileType;
			next
		}
		$fileType = "bgf";
		if ($FILES->{$_}{LOC} =~ /\.(\w+)$/) {
			$fileType = lc $1;
		}
		$FILES->{$_}{TYPE} = $fileType;
	}

	$num = "all" if (! defined($num) || $num !~ /^\d+/);
	$num = "all" if (! $num);

	$selStr = "index>0" if (! defined($selStr));
	$selection = BuildAtomSelectionString($selStr);

	$randomize = 0 if (! defined($randomize) or $randomize !~ /1|yes/i);
	$randomize = 1 if ($randomize =~ /1|yes/);

	$qeqOpt = 0 if (! defined($qeqOpt) or $qeqOpt !~ /^(1|2|yes)/i);
	$qeqOpt = 1 if ($qeqOpt =~ /yes/i);
	$qeqOpt = 2 if ($qeqOpt =~ /2/);

	$exclude = 0 if (! defined($exclude) or $exclude !~ /^\d+$/);

	$uFieldsOpt = 1 if (! defined($uFieldsOpt) or $uFieldsOpt !~ /0|no/i);
	$uFieldsOpt = 0 if ($uFieldsOpt =~ /0|no/i);

	print "Done\n";
}

sub usage {
	my ($fTypeStr) = GetFileTypeStr;
	return <<DATA;
usage: $0 -f struct file -g funct grp file -a struct atom sel -c funct grp attach id -n (number atoms) 
		-r (randomize) -q (equalize charge=no) -x (exclusion options) -u (update_fields)
		-t (file type = bgf/mol2/pdb) -m (forcefield for minimization) -s (save name)
options:
		-f struct file: Required. The structure file
$fTypeStr		
		-g funct grp file: Required. The functional group file.
		-a struct atom sel: Required. Atom(s) to functionalize. Any of the fields in the
				structure file can be use (resid atmname fftype etc). Check file type for details
		-c funct grp attach id: Required. The atom # in the funct grp file to attach 
				to the functionalized atoms
		-q equalize charge: Optional. If set to 1/yes will adjust the charge of the functionalized
				atom to be the opposite of the attached functional group. Set to 2 to use PQEq in LAMMPS minimization.
				Default is no
		-o update_fields: Optional. If set to 0 then will not update the resname|chain|resnum of the functional
				group to match the 	molecule. Default 1.	
		-x exculde options: Optional. Specify the bonded neighbors to exclude. 1 exclude neighbors, 
				2 excludes nearest neighbors, 3 next nearest neighbors... Default is none/0
		-t file type: Optional. The type of structure/functional group files being used. Either
				BGF|MOL2|PDB|MSI. If not specified, will attempt to determine from the extension.
		-m force field for minimization: Optional. Specifies the CERIUS2/MPSIM force field(s) used
				to minimize the functionalized structure with LAMMPS.
		-s save name: Optional. The name of the output file. Will be {structure file}_{functional group}.xxx
				if not specified. The minimized structure will have _min appended to the save name.
DATA
}

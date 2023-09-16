#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";

use Getopt::Std qw(getopt);
use File::Basename qw(basename);

use FileFormats qw(GetBGFAtoms ParseStructFile createBGF createMOL2
							addHeader createHeaders addBoxToHeader insertHeaderRemark);
use General qw(FileTester CombineMols GetBondLength CoM);
use BOX qw(GetBox);
use ManipAtoms qw(GetMols GetAtmData);
use REPLICATE qw(ReplicateCell);

sub numerically { ($a<=> $b); }

my ($structFile, $dim, $num, $saveFile, $head, $tail, $chargeOpt, $alignOpt, $periodicOpt);
my ($SYS, $pStr, $mol2File, $bondL, $count, $ffStr); 
my ($ATOMS, $BONDS, $HEADERS, $BOX);

$|++;
&init;
print "Parsing structure file $structFile...";
($SYS->{ATOMS}, $SYS->{BONDS}, $SYS->{HEADERS}) = ParseStructFile($structFile, 1);
$SYS->{MOLS} = GetMols($SYS->{ATOMS}, $SYS->{BONDS});
die "ERROR: BGF file contains more than 1 molecule. Cannot continue!\n"
   if(scalar(keys %{ $SYS->{MOLS} }) > 1);
print "Done\n";
$bondL = splitSystem($SYS,$head,$tail,$num,$dim);
($ATOMS, $BONDS, $HEADERS, $BOX) = &genMonomers($SYS,$dim,$num,$bondL,$chargeOpt,$periodicOpt);
print "Creating BGF file $saveFile...";
&insertHeaderRemark($HEADERS, "REMARK polymer unit replicated by $num in $dim direction");
&insertHeaderRemark($HEADERS, "REMARK polymer head definition '$head'");
&insertHeaderRemark($HEADERS, "REMARK polymer tail definition '$tail'");
&addBoxToHeader($HEADERS, $BOX);
&createMOL2($ATOMS, $BONDS, $mol2File,1);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";
&minimizeStructure($saveFile, $ffStr, $dim, $alignOpt, $periodicOpt)
	if (defined($ffStr) and $ffStr ne "");

sub minimizeStructure {
	my ($file, $ffList, $d, $aOpt, $pOpt) = @_;
	my ($execCmd, $prefix, $tot, $outStr, $inStr, $i, $tmp, $s1, $s2);

	print "Minimizing structure...Creating LAMMPS files\r";
	$prefix = basename($file);
	$prefix =~ s/\.\w+$//;
	$tmp = "";
	$execCmd = "$Bin/syncBonds.pl -b $file -s $file";
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
	if ($aOpt) {
		$s1 = $s2 = "";
		for $i ("X", "Y", "Z") {
			if($i ne $d) {
				$s1 .= " 0.0 ";
				$s2 .= " NULL ";
			} else {
				$s1 .= " NULL ";
				$s2 .= " 0.0 ";
			}
		}
		$execCmd = "$Bin/createLammpsInput.pl -b ${file} -f \"$ffList\" -s $prefix -t 'min' -o finite ";
		die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
		$execCmd = "sed -i '/^timestep/a fix sf all setforce $s1' in.${prefix}";
		die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
		$execCmd = "sed -i '/^undump/i fix sf all setforce $s2\\nminimize             1.0e-4 1.0e-4 500 5000' in.${prefix}";
		die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
	} else {
		$execCmd = "$Bin/createLammpsInput.pl -b ${file} -f \"$ffList\" -s $prefix -t 'min heat'";
		$execCmd .= " -o finite" if (! $pOpt);
		die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
	}
	
	print "Minimizing structure...Running LAMMPS minimization\r";
	$execCmd = "sed -i 's/charmmfsw/charmm/' in.${prefix}";
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
	$execCmd = "module load cpu/0.15.4  gcc/9.2.0  openmpi/3.1.6 amdfftw/2.2 &> /dev/null 2>&1; mpirun -n 1 -mca btl vader,self /home/tpascal/codes/bin/lmp_expanse -in in.${prefix} -screen none -log ${prefix}.lammps.log";
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));

	#update coordinates
	print "Minimizing structure...updating coordinates\r";
	if ($aOpt) {
		$execCmd = "$Bin/convertLammpsTrj.pl -b $file -o bgf -s ${file} -t last -l ${prefix}.min.lammps";
	} else {
		$execCmd = "$Bin/convertLammpsTrj.pl -b $file -o bgf -s ${file} -t last -l ${prefix}.heat.lammpstrj";
	}
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));

	system("rm -fr in.${prefix} in.${prefix}_singlepoint data.${prefix} tmp"); 
	system("rm -fr ${prefix}.lammps.slurm ${prefix}.lammps.log log.cite");
	system("rm -fr ${prefix}.min.lammps _byte_offset_${prefix}.heat.lammpstrj ${prefix}.heat.lammpstrj");
	print "Minimizing structure...Done\n";
}

sub genMonomers {
	my ($sys, $dim, $num, $bl, $qAdjust, $pBonds) = @_;
	my ($i, $atoms, $bonds, $tmp, $box, $rep, $N); 
	my ($totQ, $absQ, $qDiff);

	$pStr =  "Generating $num monomers along $dim direction...";
	for $i ("X", "Y", "Z") {
		$rep->{$i} = 1;
		$rep->{$i} = $num if($dim eq $i);
		#$sys->{OFFSET}{$i}{len} = 1;
	}
	if($qAdjust) {
		$qDiff->{X} = $qDiff->{Y} = $qDiff->{Z} = 0;
		($totQ, $absQ) = getCharge($sys->{MONM}{ATOMS});
		#$qDiff->{$dim} = -$totQ/2;
		if (! $pBonds) {
			$qDiff->{$dim} = -$totQ/($count->{t}+$count->{h});
		} else {
			$N = scalar(keys %{ $sys->{MONM}{ATOMS} });
			for $i (values %{ $sys->{MONM}{ATOMS}}) {
				$i->{CHARGE} -= $totQ/$N;
				#$i->{CHARGE} -= $i->{CHARGE}/$absQ;
			}
		}
	}

	#make polymer
	($atoms, $bonds, $box) = 
		ReplicateCell($sys->{MONM}{ATOMS}, $sys->{MONM}{BONDS}, $sys->{MONM}{BOX}, $rep, 0, $pStr, 0, 0, $qDiff,0);

	if(!$pBonds) {
		#reconnect head and tail
		($atoms, $bonds) = addHeadAndTail($atoms, $bonds, $sys);

		$box = GetBox($sys->{ATOMS}, undef, $sys->{HEADERS});
		$box->{$dim}{len} *= $num;
		$box->{$dim}{hi} = $box->{$dim}{len};
	}
	return ($atoms, $bonds, $sys->{HEADERS}, $box);
}

sub splitSystem {
	my ($sys, $headStr, $tailStr, $n, $d) = @_;
	my ($i, $j, $head, $tail, $monm, $rec, $a_box); 
	my ($bl, $c, $mols, $tmp, $natoms, $CoM);

	while ($headStr =~ /(\d+)\s+(\d+)/g) {
		#store pair
		($rec->{atm1}, $rec->{atm2}) = ($1, $2);

		#get bond index
		$i = findBondIndex($sys->{BONDS}{$rec->{atm1}},$rec->{atm2});
		$j = findBondIndex($sys->{BONDS}{$rec->{atm2}},$rec->{atm1});

		#remove bonds from head atm hash to main(monomer) atm hash
		splice @{ $sys->{BONDS}{ $rec->{atm1} } }, $i, 1;
		splice @{ $sys->{BONDS}{ $rec->{atm2} } }, $j, 1;
		
		#save atoms
		push @{ $monm }, $rec->{atm1};
		push @{ $head }, $rec->{atm2};
	}

	$c = 0;
	while ($tailStr =~ /(\d+)\s+(\d+)/g) {
		#store pair
		($rec->{atm1}, $rec->{atm2}) = ($1, $2);

		if($c == 0) { #store monomer offset
			for $i ("X", "Y", "Z") {
				$j = $sys->{ATOMS}{ $rec->{atm1} }{"${i}COORD"} - $sys->{ATOMS}{ $head->[0] }{"${i}COORD"};
				$sys->{MONM}{BOX}{$i}{len} = $sys->{MONM}{BOX}{$i}{len}{hi} = abs($j);
				$sys->{MONM}{BOX}{$i}{lo} = 0;
			}
		}
		#store tail offset
		for $i ("X", "Y", "Z") {
			$sys->{TAIL}{OFFSET}{$c}{$i} = $sys->{ATOMS}{ $rec->{atm2} }{"${i}COORD"} - $sys->{ATOMS}{ $rec->{atm1} }{"${i}COORD"}; 
		}

		#get bond index
		$i = findBondIndex($sys->{BONDS}{$rec->{atm1}},$rec->{atm2});
		$j = findBondIndex($sys->{BONDS}{$rec->{atm2}},$rec->{atm1});

		#remove bonds from head atm hash to main(monomer) atm hash
		splice @{ $sys->{BONDS}{ $rec->{atm1} } }, $i, 1;
		splice @{ $sys->{BONDS}{ $rec->{atm2} } }, $j, 1;
		
		#save atoms
		push @{ $monm }, $rec->{atm1};
		push @{ $tail }, $rec->{atm2};

		#now record monomer head -> tail atm bond
		$i = $monm->[$c];
		$j = $rec->{atm1};
		&saveMbond($sys,$i,$j,$d,1); #connect monomer head to tail
		&saveMbond($sys,$j,$i,$d,-1); #connect monomer tail to head
		$c++;
	}

	$mols = GetMols($sys->{ATOMS}, $sys->{BONDS});
	$i = scalar(keys %{ $mols });
	die "ERROR: Only $i fragment while seperating head and tail units! Should have at least 3!\n"
		if ($i <3);
	$tmp = getMolAtmList($sys->{ATOMS}, $mols, $monm, "monomer");
	($sys->{MONM}{ATOMS}, $sys->{MONM}{BONDS}) = 
		GetBGFAtoms($tmp,$sys->{ATOMS}, $sys->{BONDS});
	$tmp = getMolAtmList($sys->{ATOMS}, $mols, $head, "head");
	($sys->{HEAD}{ATOMS}, $sys->{HEAD}{BONDS}) = 
		GetBGFAtoms($tmp,$sys->{ATOMS}, $sys->{BONDS});
	$tmp = getMolAtmList($sys->{ATOMS}, $mols, $tail, "tail");
	($sys->{TAIL}{ATOMS}, $sys->{TAIL}{BONDS}) = 
		GetBGFAtoms($tmp,$sys->{ATOMS}, $sys->{BONDS});
	$tmp = GetMols($sys->{TAIL}{ATOMS}, $sys->{TAIL}{BONDS});

	#record bonds for tail and head atoms
	@{ $natoms } = (scalar(keys %{ $sys->{HEAD}{ATOMS} }), scalar(keys %{ $sys->{MONM}{ATOMS} }), scalar(keys %{ $sys->{TAIL}{ATOMS} }));
	$j = $natoms->[0] + $natoms->[1] * ($n - 1); #atom offset for tail atoms
	#update bond pointers
	for $i (0 .. ($c - 1)) {
		#record future head bond
		$tmp = getOIndex($sys->{HEAD}{ATOMS}, $head->[$i], "head");
		$sys->{HEAD}{ATOMS}{$tmp}{BNDATMINDX} = getOIndex($sys->{MONM}{ATOMS}, $monm->[$i], "monomer") + $natoms->[0];
		$sys->{HEAD}{ATOMS}{$tmp}{HEADATM} = 1;
		#record future tail bond
		$tmp = getOIndex($sys->{TAIL}{ATOMS}, $tail->[$i], "tail");
		$sys->{TAIL}{ATOMS}{$tmp}{BNDATMINDX} = getOIndex($sys->{MONM}{ATOMS}, $monm->[$i+$c], "monomer") + $j;
		$sys->{TAIL}{ATOMS}{$tmp}{TAILATM} = $i;
	}

	#update periodic bonds
	for $i (values %{ $sys->{MONM}{ATOMS} }) {
		next if (! exists($i->{PBOND})); 
		for $n (keys %{ $i->{PBOND} }) {
			if($i->{PSIGN}{$n}<0) {
				delete $i->{PSIGN};
				delete $i->{PBOND};
			} else {
				$tmp = getOIndex($sys->{MONM}{ATOMS}, $i->{PBOND}{$n}, "monomer") + $j; #get index for future long periodic bond
				$i->{PBOND}{$n} = $tmp; #update the periodic bond to have the correct index
			}
		}
	}

	#record monomer offset vector angle
	$tmp = $sys->{MONM}{BOX}{$d}{len};
	$sys->{MONM}{BOX} = GetBox($sys->{MONM}{ATOMS}, undef, $sys->{HEADERS});
	$sys->{MONM}{BOX}{$d}{hi} = $sys->{MONM}{BOX}{$d}{len} = $tmp;
	$sys->{MONM}{BOX}{$d}{lo} = 0;

	for $i ("X", "Y", "Z") {
		next if ($i eq $d);
		$sys->{MONM}{BOX}{$i}{hi} += 30;
		$sys->{MONM}{BOX}{$i}{len} += 30;
		$sys->{MONM}{BOX}{$i}{lo} = 0;
	}
}

#add the head and tail to the structure
sub addHeadAndTail {
	my ($atoms, $bonds, $sys) = @_;
	my ($i, $j, $k, $l, $tmp);
		
	($tmp->{atoms}, $tmp->{bonds}) = CombineMols($sys->{HEAD}{ATOMS}, $atoms, $sys->{HEAD}{BONDS}, $bonds);
	($atoms, $bonds) = CombineMols($tmp->{atoms}, $sys->{TAIL}{ATOMS}, $tmp->{bonds}, $sys->{TAIL}{BONDS});
	&GetMols($atoms, $bonds);

	$tmp = $sys->{MONM}{BOX}{$dim}{len}; #half the box length for determining long bonds
	#head atoms
	for $i (keys %{ $atoms }) {
		next if (! exists($atoms->{$i}{HEADATM}) and ! exists($atoms->{$i}{TAILATM}));
		$j = $atoms->{$i}{BNDATMINDX}; #bond atom on polymer

		#align tail atoms
		&alignTailAtms($atoms, $i, $j, $sys->{TAIL}{OFFSET}) if(exists($atoms->{$i}{TAILATM}));
		
		#bond head/tail to polymer
		push @{ $bonds->{$i} }, $j;
		push @{ $bonds->{$j} }, $i;

		next if (! exists($atoms->{$i}{HEADATM}));
		
		#now delete polymer - polymer bond long bond for atom j
		next if (! exists($atoms->{$j}{PBOND}));
		for $k (values %{ $atoms->{$j}{PBOND} }) {
			$l = findBondIndex($bonds->{$j}, $k);
			splice @{ $bonds->{$j} }, $l, 1;
			delete $atoms->{$j}{MONMBND};
			$l = findBondIndex($bonds->{$k}, $j);
			splice @{ $bonds->{$k} }, $l, 1;
			delete $atoms->{$j}{MONMBND};
			delete $atoms->{$i}{HEADATM};
		}
	}
	return ($atoms, $bonds);
}

#This was created to handle the case where the new monomer bond already exists
sub saveMbond {
	my ($sys, $i, $j, $d, $s) = @_;
	my ($n, $t, $bexists);

	$bexists = 0;
	for $n (@{ $sys->{BONDS}{$i} }) {
		if ($n == $j) {
			$bexists = 1;
			last;
		}
	}
	$t = $#{ $sys->{BONDS}{$i}}; #total number of bonds to atom i
	if ($bexists and ! exists($sys->{ATOMS}{$i}{"DISP${d}"})) { #if the DISP array is not set
		for $n (0 .. $t) {
			$sys->{ATOMS}{$i}{"DISP${d}"}[$n] = 0; #set all "normal" bond DISP to 0
		}
	}
	push @{ $sys->{BONDS}{$i} }, $j; #create i -> j bond
	$sys->{ATOMS}{$i}{MONMBND} = 1; #update flag for use when replicating monomers
	$sys->{ATOMS}{$i}{"DISP${d}"}[$t+1] = -$s * 1; #set the appropriate DISP flag
	#record periodic bond, will be updated later to have the correct id
	$sys->{ATOMS}{$i}{PBOND}{$j} = $j;
	$sys->{ATOMS}{$i}{PSIGN}{$j} = $s;
}

sub alignTailAtms {
	my ($atoms, $i, $j, $offset) = @_;
	my ($k, $l, $molAtms, $d);

	$l = $atoms->{$i}{TAILATM}; #index for offset hash
	for $k ("X", "Y", "Z") {
		$d->{$k} = $atoms->{$j}{"${k}COORD"} - $atoms->{$i}{"${k}COORD"} + $offset->{$l}{$k};
	}

	for $l (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
		for $k ("X", "Y", "Z") {
			$atoms->{$l}{"${k}COORD"} += $d->{$k};
		}
	}
	delete $atoms->{$i}{TAILATM};
}

sub getOIndex {
	my ($atoms, $id, $iType) = @_;
	my ($i, $oindex);

	for $i (keys %{ $atoms }) {
		if($atoms->{$i}{OINDEX} == $id) {
			$oindex = $i;
			last;
		}
	}

	die "ERROR: Cannot find $iType atom # $id!\n"
		if(! defined($oindex));
	return $oindex;
}

sub getMolAtmList {
	my ($atoms, $mols, $atmArray, $tStr) = @_;
	my ($i, $j, $m, $atmList, $uMols);

	for $i (@{ $atmArray }) {
		$m = ${ $atoms->{$i}{MOLECULEID} };
		die "ERROR: Fragment $m does not exists while looking for atom # $i ($tStr)! Check head and tail definiton\n"
			if (! exists($mols->{$m}));
		$uMols->{$m} = 1;
		for $j (keys %{ $mols->{$m}{MEMBERS} }) {
			$atmList->{$j} = 1;
			delete $atoms->{$j}{BONDATOM};
		}
	}
	for $m (keys %{ $uMols }) {
		delete $mols->{$m};
	}

	return $atmList;
}

sub findBondIndex {
	my ($bonds, $atom) = @_;
	my ($i, $index);
	
	for $i (0 .. $#{ $bonds }) {
		if($bonds->[$i] == $atom) {
			$index = $i;
			last;
		}
	}

	return $index;
}

sub getCharge {
	my ($atoms) = $_[0];
	my ($i, $totQ, $absQ);

	for $i (keys %{ $atoms }) {
		$totQ += $atoms->{$i}{CHARGE};
		$absQ += abs($atoms->{$i}{CHARGE});
	}

	return ($totQ, $absQ);
}

sub init {
	my (%OPTS, $rStr, $tStr, $tmp);

	getopt('bdnshtqfap',\%OPTS);
	for ("b", "n", "h", "t") {
		die &showUsage if (! exists($OPTS{$_}));
	}
	print "Initializing...";
	($structFile, $dim,     $num,     $saveFile, $head,    $tail,    $chargeOpt, $ffStr,   $alignOpt, $periodicOpt) = 
	($OPTS{b}, $OPTS{d}, $OPTS{n}, $OPTS{s},  $OPTS{h}, $OPTS{t}, $OPTS{q},   $OPTS{f}, $OPTS{a}, $OPTS{p});
	@{ $tmp } = split /\s+/,$structFile;
	for (@{ $tmp }) {
		FileTester($_); 
	}
	die "ERROR: Expected x|y|z for growth direction, got \"$dim\"!\n"
		if($dim !~ /(x|y|z)/i);
	$dim = uc $1;
	die "ERROR: Expected integer for number monomers, got \"$num\"!\n"
		if($num !~ /(\d+)/);
	$num = $1;
	$tStr = "";
	$count->{h} = 0;
	while($head =~ /(\d+\s+\d+)/g) {
		$tStr .= "$1 ";
		$count->{h}++;
	}
	die "ERROR: Atom specification for head is invalid"
		if (! $count->{h});
	chop $tStr;
	$head = $tStr;
	$tStr = "";
	$count->{t} = 0;
	while($tail =~ /(\d+\s+\d+)/g) {
		$tStr .= "$1 ";
		$count->{t}++;
	}
	die "ERROR: Atom specification for tail is invalid"
		if (! $count->{t});
	chop $tStr;
	$tail = $tStr;
	undef $head if (defined($head) and $head !~ /\d+\s+\d+/);
	undef $tail if (defined($tail) and $tail !~ /\d+\s+\d+/);
	die "ERROR: $count->{h} atom pair selected for head but $count->{t} atom pairs selected for tail!\n"
		if($count->{h} != $count->{t});
	if (! defined($saveFile)) {
		$structFile =~ /^\s*(\S+)/;
		$saveFile = basename ($1);
		$saveFile =~ s/\.\w+$//;
		$saveFile .= "x${num}.bgf";
	}
	$mol2File = $saveFile;
	$mol2File =~ s/\.\w+$//;
	$mol2File .= ".mol2";
	srand(time() ^($$ + ($$ <<15)));
	$chargeOpt = 0 if (! exists($OPTS{q}) or $OPTS{q} !~ /1|yes/i);
	$chargeOpt = 1 if ($chargeOpt =~ /1|yes/i);
	$alignOpt = 0 if (! exists($OPTS{a}) or $OPTS{a} !~ /1|yes/i);
	$alignOpt = 1 if ($alignOpt =~ /1|yes/i);
	$periodicOpt = 0 if (! defined($periodicOpt) or $periodicOpt !~ /1|yes/i);
	$periodicOpt = 1 if ($periodicOpt =~ /1|yes/i);
	print "Done\n";
}

sub showUsage {
	return <<DATA;
usage: $0 -b monomer file -d growth direction -n #monomer 
			-t [tail atoms] -h [head atoms] -f (forcefield)
			-a (align) -q (q equalize) -p (make periodic) -s (save name)
	Required:
	-b monomer file: Location of monomer structure file in BGF format
	-d growth direction: X|Y|Z axis along with to grow polymer
	-n #monomer: number of monomers units in final polymer
	-t tail atoms: atom indices specifying start of the polymer structure. (see below)
	    Each pair of numbers represents the bond that will initially be broken 
		to create the polymer. Can use multiple pairs of integers to represent
		multiple polymerization points, but need same number of pairs in head and tail
	-h head atoms: atom indices specifying end of the polymer structure. Same as tail

	Optional:
	-f forcefield: Monomer forcefield location. When specified will optimize the polymer in LAMMPS
	-a polymer align: Expected: 1|yes|0|no. When specified, will constrain the minimization 
		so that the polymer is generally aliged along the growth direction. Default: 0
	-q q equalize: Expected: 1|yes|0|no. When specified will adjust the charge at the attachment
		points so that the total charge of the polymer is the same as n x monomer_charge. Default: 0
	-p make_periodic: Expected 0|no(default)|1|yes. make the system periodic by making an additional 
		bond between the head and tail.
	-s save name: Name of final polymer structure. Default: monomerx{N}.bgf

	Atom specification schematic
                                               N -->
                                R    ...........................    R
                                |    .    |               |    .    | 
                            R-- z -- X -- y -- R -- R --  b -- X -- a -- R
                                |    .    |               |    .    |
                                R    ...........................    R
	Head spec: b a 
	Tail spec: y z 
	Dotted box is repeated unit (with b <--> y bond)
DATA

}

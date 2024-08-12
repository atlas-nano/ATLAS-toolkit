#!/usr/bin/perl

$DB::deep = 50000; # or more if necessary
use strict;
use warnings;
use POSIX;

use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";

use Storable qw(dclone);
use File::Basename qw(basename);


use General qw(FileTester LoadElements Permutate GetSoluteAtoms FindElement
			CRadDegrees LoadFFs IsDecimal AddRingField AddElementField ReadFFs
			GetFileTypeStr GetFFTypeStr CoM LoadConverter);
use FileFormats qw(GetBondList WriteSWfile addHeader createBGF insertHeaderRemark 
			ParseStructFile);
use CERIUS2 qw(GenerateUFFParms getLammpsOpts findDuplicate);
use BOX qw(GetBox);
use ManipAtoms qw(GetMols BuildAtomSelectionString SelectAtoms 
                  SelectAtomsByField AddMolsToSelection);
use LAMMPS qw(CreateInputFile);
use Getopt::Std qw(getopt);
use Math::Polynomial::Solve qw(cubic_roots GetLeastSquaresPlane);

my ($inputFile, $FF, $suffix, %OPTS, $inputType, $ScaleTorsionList, $reaxFF, $opts);
my ($ATOMS, $CONS, $HEADERS, $FILES, $PARMS, %ERRORS, $i, $BOX, $ffType, $nocross, $sysOpts);
my ($BONDS, $TORSIONS, $INVERSIONS, $ANGLES, $PI_TORSIONS, $BI_TORSIONS); 
my ($CMAP, @atmIndices, $totAtms, $MOLS, $QEqtmp, $qxFile, $readFile);

$|++;
my ($start) = time();
$FILES = &init;
$PARMS = LoadFFs($FILES, 2);
print "Step 3: Parsing structure file $inputFile...";
($ATOMS, $CONS, $HEADERS) = ParseStructFile($inputFile, 1);
$BOX = GetBox($ATOMS, $PARMS, $HEADERS);
$MOLS = GetMols($ATOMS, $CONS);
$PARMS->{PARMS}{MOLS} = $MOLS;
&addCoreShellAtoms($ATOMS,$CONS, $PARMS, $sysOpts->{coreShell}) if (defined($sysOpts->{coreShell}));
&checkAtomTypes($ATOMS, $PARMS, $ffType, $sysOpts);
&GenerateUFFParms($ATOMS, $CONS, $PARMS) if exists($PARMS->{PARMS}{UFF});
&fixMesoDNAAtoms($ATOMS, $CONS) if ($sysOpts->{mesodna}); #mesodna fix
&findDuplicateParms($PARMS, "VDW BONDS ANGLES TORSIONS PI_TORSIONS");
&updateTorsionList($PARMS->{TORSIONS});
&setOpts($ATOMS, $PARMS, $opts, $QEqtmp, $MOLS);
#&AddRingField($ATOMS, $CONS, $ATOMS);
&addQEqData($ATOMS,$PARMS, $QEqtmp) if (defined($QEqtmp));
# remove all valence interactions for reax/3 body force fields or as specified
&deleteBonds($ATOMS, \%{ $CONS }, $sysOpts->{nobonds}, $ffType) if ($ffType == 5 or exists($sysOpts->{nobonds})); 
&setElectrodeOpts($ATOMS, $CONS, $PARMS, $sysOpts->{electrode}, $QEqtmp) if (defined($sysOpts->{electrode}));
print "Done\n";
&printErrors(\%ERRORS, 1) if (keys %ERRORS);
print "Step 4: Determining valence list from connectivities...";
@atmIndices = sort numerically keys %{ $ATOMS };
&removeCrossParms($PARMS) if ($nocross);
&getValenceParms($ATOMS, $CONS, $PARMS, \@atmIndices);
&findDuplicateParms($PARMS, "INVERSIONS");
&updateParmIndex($PARMS);
&determineTinkerOpts($PARMS, $suffix, $sysOpts->{amoeba_monopole});
&addPQEqOpt($ATOMS, $CONS, $PARMS) if (defined($QEqtmp) and $QEqtmp->{Opt} == 2);
&getTorsionScalingFactor($TORSIONS);
print "will use single inversion..." if ($PARMS->{PARMS}{single_inversion});
&getFEPopt($PARMS,\%{ $ATOMS },$CONS,$sysOpts->{fep}) if(defined($sysOpts->{fep}));
&getRigidOpt($PARMS, $ATOMS, $CONS, $sysOpts) if (defined($sysOpts->{rigid}));
&addLammpsParms($PARMS, $ATOMS, $MOLS, $BONDS, $ANGLES, $CMAP, $suffix, $sysOpts, $sysOpts->{coreShell});
print "Done\nStep 5: Writing data file data.${suffix}...\r";
open DATFILE, "> data.${suffix}" or die "ERROR: Cannot create data.${suffix}: $!\n";
&createDatFileHeader($PARMS, $BOX, \*DATFILE, $sysOpts);
&printAtoms($ATOMS, $BOX, \*DATFILE, $PARMS->{QEq});
for $i ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS", "CMAP", "BI_TORSIONS", "PI_TORSIONS") {
	&printValence(eval('$' . $i), \*DATFILE, $i);
}
&writeTinkerTypes($ATOMS, $PARMS->{ATOMTYPES}, \*DATFILE) if (exists($PARMS->{PARMS}{ALL_TINKER}) and $PARMS->{PARMS}{ALL_TINKER} == 1);
&writeIsotopeMass($ATOMS, \*DATFILE) if ($PARMS->{PARMS}{OPTIONS}{ISOTOPES});
close DATFILE;
&writeTinkerBiTorsionFile($PARMS, $suffix) if(exists($PARMS->{BI_TORSIONS}) and $PARMS->{BI_TORSIONS}{counter}>0);
if(exists($PARMS->{PARMS}{POLARIZATION}) and keys %{ $PARMS->{PARMS}{POLARIZATION} }) {
	&insertHeaderRemark($HEADERS, "REMARK added shell atoms");
	&addHeader($ATOMS,$HEADERS);
	&createBGF($ATOMS, $CONS, "${suffix}.core-shell.bgf");
}
print "Step 6: Creating $PARMS->{PARMS}{INPUTNAME} input files in.${suffix}...";
&getTIP4Popts($PARMS, $ATOMS, $BONDS) if (exists($PARMS->{PARMS}{is_tip4p}) and $PARMS->{PARMS}{is_tip4p});
&writeTypeCharges($PARMS, $ATOMS)
	if($sysOpts->{TypeCharges} and $PARMS->{QEq}{Opt}==0);
&writeCMAPfile($PARMS, $suffix) if(defined($CMAP) and $CMAP->{counter}>0);	
&CreateInputFile($PARMS);
&writeQEqParms($PARMS, $suffix) if ($PARMS->{QEq}{Opt}==1);
print "Created in.${suffix}_singlepoint...Done\n";
print "Step 7: Creating LAMMPS cluster script file ${suffix}.lammps.slurm...";
&createLammpsClusterScript($FILES);
print "Done\n";
my ($end) = time();
printf "Elapsed time: %.3f secs\n", ($end - $start);
&printErrors(\%ERRORS, 0) if (keys %ERRORS);

sub genInvList {
	my ($blist) = @_;
	my (@tmp,$i,$j,$k,$inv);

	@tmp = keys %{ $blist };
	for $i (@tmp) {
	for $j (@tmp) {
		next if ($i eq $j);
		for $k (@tmp) {
			next if ($j eq $k);
			$inv->{$i}{$j}{$k} = 1;
		}
	}
	}
	return $inv;

}

sub setElectrodeOpts {
	my ($atoms, $bonds, $parms, $eltrde, $qeqParm) = @_;
	my ($tmp, $i, $sel, $max_res);

	$sel = SelectAtoms($eltrde->{top}{aStr}, $atoms);
	$parms->{PARMS}{ELECTRODE}{top}{MOLID} = SelectAtomsByField($atoms, $bonds, "RESNUM", $sel);
	$sel = SelectAtoms($eltrde->{bot}{aStr}, $atoms);
	$parms->{PARMS}{ELECTRODE}{bot}{MOLID} = SelectAtomsByField($atoms, $bonds, "RESNUM", $sel);
	for $i (keys %{ $atoms }) {
		$max_res = $atoms->{$i}{RESNUM} if ($atoms->{$i}{RESNUM} > $max_res);
	}
	$parms->{PARMS}{ELECTRODE}{max_res} = $max_res;
	$parms->{PARMS}{ELECTRODE}{type} = $eltrde->{type};

	if ($eltrde->{type} eq "qeq") {
		@{ $tmp } = grep {!/TYPE|counter/i} keys %{ $parms->{VDW} };
		$sel = SelectAtoms($eltrde->{top}{aStr}, $atoms);
		$eltrde->{top}{FFTYPE} = SelectAtomsByField($atoms, $bonds, "FFTYPE", $sel);
		$sel = SelectAtoms($eltrde->{bot}{aStr}, $atoms);
		$eltrde->{bot}{FFTYPE} = SelectAtomsByField($atoms, $bonds, "FFTYPE", $sel);
		$parms->{PARMS}{ELECTRODE}{type} = $eltrde->{type};
		for $i ($tmp) {
			$parms->{PARMS}{ELECTRODE}{CHI}{ $parms->{ATOMTYPES}{$i}{INDEX} } = $qeqParm->{CHI};
			$parms->{PARMS}{ELECTRODE}{top}{FFTYPE}{$i} = $qeqParm->{CHI}
				if(exists($eltrde->{top}{FFTYPE}{$i}));
			$parms->{PARMS}{ELECTRODE}{bot}{FFTYPE}{$i} = $qeqParm->{CHI}
				if(exists($eltrde->{bot}{FFTYPE}{$i}));
		}
	}
}

sub addPQEqOpt {
	my ($atoms, $bonds, $parms) = @_;
	my ($i, $j, $k, $indx, $tmp, $curr, $qeqParm, @types); 
	my ($ldata, $lindx, $special_coul, $nspecial);


	$parms->{PARMS}{HYBRID_VDW} = 2;
	@types = grep {!/TYPE|counter/i} keys %{ $parms->{VDW} };
	#iterate the vdw types array count the number of 1-2/1-3/1-4 coul types
	for $i (@types) {
		next if (!exists($parms->{QEq}{$i}));
		for $j (keys %{ $parms->{VDW}{$i}{$i} }) {
			next if (! defined($parms->{VDW}{$i}{$i}{$j}{COUL}));
			$special_coul->{line} = $parms->{VDW}{$i}{$i}{$j}{COUL};
			$special_coul->{indx}{$special_coul->{line}}++;
		}
	}
	#now if more that one special coul flags then select the one with the most atom types and set the lammps flags
	$nspecial = scalar(keys %{ $special_coul->{indx} });
	if($nspecial > 1) {
		$special_coul->{max}{val} = 0;
		for $i (keys %{ $special_coul->{indx} }) {
			if($special_coul->{indx}{$i} > $special_coul->{max}{val}) {
				$special_coul->{max}{val} = $special_coul->{indx}{$i};
				$special_coul->{max}{type} = $i;
			}
		}
		$curr = $special_coul->{max}{type};
		$parms->{PARMS}{scale_cou_12} = substr($curr,0,1);
		$parms->{PARMS}{scale_cou_13} = substr($curr,1,1);
		$parms->{PARMS}{scale_cou_14} = substr($curr,3);
	}

	$parms->{VDW}{TYPE}{"coul/pqeqgauss"} = (
												{
													"OPTS"     => "0.00 12.50",
													"NAME"     => "PQEq",
													"PAIRMIX"  => 0,
													"ORDER"    => 1,
													"LMPNAME"  => "coul/pqeqgauss",
													"LMPONAME" => "coul/pqeqgauss",
													"SCC"      => 1,
													"FEP" => {
														'pair'      => 'coul/pqeqgauss',
														'pair_opts' => '',
														'addsoft'   => 0,
														'parms'     => ['lambda'],
													},
												}
											);
	for $i (2 .. 4) {
		$parms->{PARMS}{"scale_cou_1${i}"} = 1e-12 if (! $parms->{PARMS}{"scale_cou_1${i}"} );
	}

	for $i (@types) {
		next if (!exists($parms->{QEq}{$i}));
		$curr = $parms->{VDW}{$i}{$i};
		$qeqParm = $parms->{QEq}{$i};
		@{ $tmp } = sort numerically keys %{ $curr };
		$ldata = ();
		$lindx = 0;
		$indx = 0;
		if(@{ $tmp }) {
			$ldata = \%{ $curr->{$tmp->[$#{ $tmp }] }{DATA} };
			$lindx = $curr->{ $tmp->[$#{ $tmp }] }{INDEX};
			$indx = $tmp->[$#{ $tmp }] + 1;
		}
		$tmp = ();
		$curr->{$indx} = (
							{
								"ATOM"     => "$i $i ",
								"IGNORE"   => 0,
								"IT"       => "vdw",
								"KEY"      => "$i",
								"DERIVED"  => 1,
								"Lammps"   => \%{ $parms->{VDW}{TYPE}{"coul/pqeqgauss"} },
								"TYPE"     => "PQEq",
								"USED"     => 0,
								"DATA"     => $ldata,
								"INDEX"    => $lindx,
								"VALS"    => [$qeqParm->{CHI},$qeqParm->{ETA},$qeqParm->{Rc},$qeqParm->{P},$qeqParm->{Zcore},
											  $qeqParm->{Rs},$qeqParm->{Ks1},$qeqParm->{Ks2}],
							  }
		);

		if (exists($parms->{PARMS}{POLARIZATION}) and 
		(exists($parms->{PARMS}{POLARIZATION}{SHELLS}{$i}) or exists($parms->{PARMS}{POLARIZATION}{CORES}{$i}))) {
			$curr->{$indx}{VALS}[3] = 0; #set polarization flag to 0
			#now check for drude bond and update force constant to match pair
			$j = "${i}_shell";
			if($i =~ /_shell/) {
				$j = $i;
				$j =~ s/_shell//;
			}
			if(exists($parms->{BONDS}{$i}) and exists($parms->{BONDS}{$i}{$j})) {
				$parms->{BONDS}{$i}{$j}{1}{VALS}[0] = $curr->{$indx}{VALS}[6]; 
			} elsif (exists($parms->{BONDS}{$j}) and exists($parms->{BONDS}{$j}{$i})) {
				$parms->{BONDS}{$j}{$i}{1}{VALS}[0] = $curr->{$indx}{VALS}[6];
			}
		}
		undef $k;
		for $j (sort numerically keys %{ $parms->{VDW}{$i}{$i} }) {
			if(!exists($parms->{VDW}{$i}{$i}{$j}{DERIVED}) and defined($parms->{VDW}{$i}{$i}{$j}{COUL})) {
				$k = $j;
				last;
			}
		}
		if($nspecial > 1 and defined($k) and $parms->{VDW}{$i}{$i}{$k}{COUL} ne $special_coul->{max}{type}) {
			push @{ $curr->{$indx}{VALS} }, substr($parms->{VDW}{$i}{$i}{1}{COUL},0,1);
			push @{ $curr->{$indx}{VALS} }, substr($parms->{VDW}{$i}{$i}{1}{COUL},1,1);
			push @{ $curr->{$indx}{VALS} }, substr($parms->{VDW}{$i}{$i}{1}{COUL},3);
		}
	}

}

sub getPQEQchi {
	my ($parms, $fftype) = @_;
	my ($vdw, $i, $rVal);

	return if (! exists($parms->{VDW}{$fftype}));
	$vdw = $parms->{VDW}{$fftype}{$fftype};
	for $i (values %{ $vdw }) {
		if ($i->{TYPE} eq "PQEq") {
			$rVal = $i->{VALS}[0];
			last;
		}
	}
	return $rVal;
}

sub addQEqData {
	my ($atoms, $parms, $tmp) = @_;
	my ($i, $ele, $ftype, $qeqParms);

	$qeqParms = \%{ $parms->{QEq} };
	$qeqParms->{Opt} = $tmp->{Opt};
	$qeqParms->{File} = $tmp->{File};
	$qeqParms->{sys_charge} = $tmp->{sys_charge};
	return if ($qeqParms->{Opt} == 0);

	if($qeqParms->{Opt} == 1) {
		$qeqParms->{DATA} = parseQEqFile($qeqParms->{File});
	} elsif ($qeqParms->{Opt} == 2) {
		$qeqParms->{DATA} = parsePQEqFile($qeqParms->{File});
	}

	for $i (keys %{ $atoms }) {
		$ele = $atoms->{$i}{ELEMENT}{SYMBOL};
		$ftype = $atoms->{$i}{FFTYPE};
		next if exists($qeqParms->{$ftype});
		if(exists($qeqParms->{DATA}{$ftype})) {
			$qeqParms->{$ftype} = $qeqParms->{DATA}{$ftype};
			next;
		}elsif(exists($atoms->{$i}{CORE_PARENT}) and exists($qeqParms->{DATA}{$atoms->{$i}{CORE_PARENT}{TYPE}})) {
			#handle the case for shell atoms
			$qeqParms->{$ftype} = $qeqParms->{DATA}{$atoms->{$i}{CORE_PARENT}{TYPE}};
			next;
		}
		die "ERROR: No QEq paramaters found for atom $i (type: $ftype element: $ele)!\n"
			if (!exists($qeqParms->{DATA}{$ele}));
		$ele .= "w" 
			if(isWater($atoms, $i) and $ele =~ /^H|O/ and exists($qeqParms->{DATA}{"${ele}w"}));
		$qeqParms->{$ftype} = $qeqParms->{DATA}{$ele};
	}

}

sub isWater {
	my ($atoms, $i) = @_;
	my ($j, $eleList);

	return 0 if (${$atoms->{$i}{MOLSIZE}} != 3);
	@{ $eleList } = map { $atoms->{$_}{ELEMENT}{SYMBOL} } keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} };
	return 1 if ("@{ $eleList }" =~ /(H H O|O H H|H O H)/);
	return 0;
}

sub writeQEqParms {
	my ($parms,$prefix) = @_;
	my ($DATA, $count, $atype, $curr, $index);

	print "writing QEq data to ${prefix}.param.qeq...";
	open QEqFile, "> ${prefix}.param.qeq" or die "ERROR: Cannot write to param.qeq: $!\n";
	($DATA, $count) = sortParmByIndex($parms->{VDW}, 1);
	for $index (1 .. $count) {
		$atype = $DATA->{$index}{DATA}{LABEL};
		die "ERROR: Atomtype $atype does not have any QEq paramaters!\n"
			if (! exists($parms->{QEq}{$atype}));
		$curr = $parms->{QEq}{$atype};
		printf QEqFile "%-5d %8.4f %8.4f %8.4f %8.4f %8.4f\n",$index,$curr->{CHI},$curr->{ETA},$curr->{GAMMA},$curr->{ZETA},$curr->{CHRG};
	}
	close QEqFile;
}

sub parseQEqFile {
	my ($pfile) = $_[0];
	my ($parm);

	print "Parsing QEq param file $pfile...";
	open QEqPARM, $pfile or die "ERROR: Cannot access $pfile: $!\n";
	while(<QEqPARM>) {
		chomp;
		if ($_ =~ /^([a-zA-Z]+)\s+(\d+)\s+(\-?\d+)\s+(\-?\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
			$parm->{$1} = (
							{
								"n"     => $2,
								"qmin"  => $3,
								"qmax"  => $4,
								"CHI"   => $5*23.509,
								"ETA"   => $6*2*23.509,
								"GAMMA" => 0.1,
								"R"     => $7,
								"ZETA"  => (2*$2+1)/(4*$7*1.88973),
								"CHRG"  => 0,
							}
						);
		}
	}
	close QEqPARM;
	die "ERROR: No valid data read\n" if (! defined($parm));
	return $parm;
}

sub parsePQEqFile {
	my ($pfile) = $_[0];
	my ($parm);

	print "Parsing PQEq param file $pfile...";
	open QEqPARM, $pfile or die "ERROR: Cannot access $pfile: $!\n";
	while(<QEqPARM>) {
		chomp;
		if ($_ =~ /^\s*([a-zA-Z_0-9\+?\-?]+)\s+(\-?\d+\.?\d*)\s+(\-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*(\d*\.?\d*)/) {
			$parm->{$1} = (
							{
								"P"     => int($2),
								"CHI"   => $3,
								"ETA"   => $4,
								"Zcore" => $5,
								"Rc"    => $6,
								"Rs"    => $7,
								"Ks1"   => $8,
								"Ks2"   => 0,
							}
						);
			$parm->{$1}{Ks2} = $9 if(defined($9) and $9);			
		}
	}
	close QEqPARM;
	die "ERROR: No valid data read\n" if (! defined($parm));
	return $parm;
}

sub printValence {
	my ($valence, $datFile, $header) = @_;
	my ($i, $j, $count);
	
	return if (!$valence->{counter});
	$header = "dihedrals" if ($header eq "TORSIONS");
	$header = "impropers" if ($header eq "INVERSIONS");
	
	if ($header =~ /^(.)(.*)_(.)(.*)$/) {
		$header = uc($1) . lc($2) . uc($3) . lc($4);
	} elsif($header ne "CMAP") {
		$header = uc(substr($header, 0, 1)) . lc (substr($header,1,length($header)));
	}
	print "Step 5: Writing data file...$header\r";
	print $datFile "\n$header\n";
	for $i (1 .. $valence->{counter}) {
		printf $datFile "\n%8d %8d ", $i, $valence->{LIST}{$i}{DATA}{INDEX};
		$count = 0;
		for $j (@{ $valence->{LIST}{$i}{ATOMS} }) {
			printf $datFile "%8d ", $j;
			$count++;
		}
	}
	print $datFile "\n";
}
	
sub createDatFileHeader {
	my ($parms, $box, $datFile, $datWrite) = @_;
	my ($index, $curr, $parm, $DATA, $count, $j, $valence, $a, $b, $c, $parmHybrid, $crossterms);
	my ($cos_alpha, $sin_alpha, $cos_beta, $sin_beta, $cos_gamma, $sin_gamma, $tmp, $shouldWrite);

	$parms->{PARMS}{isTriclinic} = 0;
	print "Step 5: Writing data file...header\r";
	print $datFile "Created by $0 on " . scalar(localtime) . "\n\n";
	printf $datFile "%12d atoms\n", $totAtms;
	for $i ("UREY_BRADLEY", "BOND_ANGLE") {
		$crossterms->{$i}{str} = "";
		$crossterms->{$i}{header} = $i;
		$crossterms->{$i}{header} =~ /^(.)(.*)_(.)(.*)$/;
		$crossterms->{$i}{header} = "\n\n" . uc($1) . lc($2) . uc($3) . lc($4) . " Coeffs\n";
	}

	for $index ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS", "PI_TORSIONS", "BI_TORSIONS" ) {
		$curr = eval('$' . $index);
		next if ($curr->{counter} == 0);
		$parm = lc($index);
		$parm =~ s/_//;
		$parm = "dihedrals" if ($parm eq "torsions");
		$parm = "impropers" if ($parm eq "inversions");

		printf $datFile "%12d $parm\n", $curr->{counter};		
	}
	printf $datFile "%12d crossterms\n", $CMAP->{counter} 
		if (defined($CMAP) and $CMAP->{counter} > 0);

	printf $datFile "\n%12d atom types\n", $parms->{ATOMTYPES}{counter};
	for $index ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS", "PI_TORSIONS") {
		$parm = lc(substr($index,0,-1));
		$parm =~ s/_//;
		$parm = "dihedral" if ($parm eq "torsion");
		$parm = "improper" if ($parm eq "inversion");
		printf $datFile "%12d $parm types\n", $parms->{$index}{counter} 
			if (exists($parms->{$index}) and $parms->{$index}{counter} > 0);	
	}

	print "Step 5: Writing data file...box\r";
	$a = $box->{X}{len};
	$b = $box->{Y}{len};
	$c = $box->{Z}{len};
	$cos_alpha = cos(CRadDegrees($box->{X}{angle}, 0));
	$sin_alpha = sqrt(1-$cos_alpha*$cos_alpha);
	$cos_beta = cos(CRadDegrees($box->{Y}{angle}, 0));
	$sin_beta = sqrt(1-$cos_beta*$cos_beta);
	$cos_gamma = cos(CRadDegrees($box->{Z}{angle}, 0));
	$sin_gamma = sqrt(1-$cos_gamma*$cos_gamma);
	for $index ("X","Y","Z") {
		$parms->{PARMS}{isTriclinic} = 1 if ($box->{$index}{angle} != 90);
	}
	if($parms->{PARMS}{OPTIONS}{PERIODICITY} == 0) {
		for $index ("X", "Y", "Z") {
			$box->{$index}{lo} -= 100; $box->{$index}{hi} += 100;
		}
	}
	if($box->{X}{angle} == 90 and $box->{Y}{angle} == 90 and $box->{Z}{angle} == 90) {
		$box->{LAMMPS}{X}{max} = $box->{X}{hi}; $box->{LAMMPS}{X}{min} = $box->{X}{lo};
		$box->{LAMMPS}{Y}{max} = $box->{Y}{hi}; $box->{LAMMPS}{Y}{min} = $box->{Y}{lo};
		$box->{LAMMPS}{Z}{max} = $box->{Z}{hi}; $box->{LAMMPS}{Z}{min} = $box->{Z}{lo};
	} else {
		$box->{LAMMPS}{X}{max} = $a; $box->{LAMMPS}{X}{min} = $box->{X}{lo};
		$box->{LAMMPS}{Y}{max} = $sin_gamma*$b; $box->{LAMMPS}{Y}{min} = $box->{Y}{lo};
		$box->{LAMMPS}{Z}{max} = sqrt($c*$c*$sin_beta*$sin_beta - $c*($cos_alpha-$cos_gamma*$cos_beta)/$sin_gamma); $box->{LAMMPS}{Z}{min} = $box->{Z}{lo};
	}
	printf $datFile "\n		 %10.6f %10.6f xlo xhi", $box->{LAMMPS}{X}{min},$box->{LAMMPS}{X}{max};
	printf $datFile "\n		 %10.6f %10.6f ylo yhi", $box->{LAMMPS}{Y}{min},$box->{LAMMPS}{Y}{max};
	printf $datFile "\n		 %10.6f %10.6f zlo zhi", $box->{LAMMPS}{Z}{min},$box->{LAMMPS}{Z}{max};

	if ($parms->{PARMS}{isTriclinic}) { #write xy xz yz
		printf $datFile "\n		 %10.6f %10.6f %10.6f xy xz yz",
			$b*$cos_gamma,
			$c*$cos_beta,
			$c*($cos_alpha-$cos_gamma*$cos_beta)/$sin_gamma;
	}
	print "Step 5: Writing data file...coeffs\r";
	%{ $tmp } = ();
	for $i (keys %{ $parms->{ATOMTYPES} }) {
		next if ($i eq "counter" or !$parms->{ATOMTYPES}{$i}{USED});
		$tmp->{$parms->{ATOMTYPES}{$i}{TYPEID}} = $i;
	}
	print $datFile "\n\nMasses\n\n";
	$i = 1;
	for $index (sort numerically keys %{ $tmp }) {
		printf $datFile "%5d %8.4E", $i, $parms->{ATOMTYPES}{$tmp->{$index}}{MASS};
		printf $datFile "%-10s",' # ' . $parms->{ATOMTYPES}{$tmp->{$index}}{NAME} if ($datWrite->{Labels});
		print $datFile "\n";
		$i++;
	}
	($DATA, $count) = sortParmByIndex($parms->{VDW}, 1);
	if (! $datWrite->{InputCoeffs}) {
		$shouldWrite = 1; #should we write pair coeffs in datafile?
		if ($ffType =~ /5|8/ or scalar(keys %{ $parms->{VDW}{TYPE} }) > 1 or $datWrite->{InputCoeffs} or ! keys %{ $DATA }) {
			$shouldWrite = 0;
		} else {
			for $index (1 .. $parms->{ATOMTYPES}{counter}) {
				$shouldWrite = 0 if ($DATA->{$index}{Lammps}{LMPNAME} =~ /(sw|eam)/i); 
			}
		}
	}
	$shouldWrite = 0 if ($parms->{PARMS}{HYBRID_VDW} > 0);
	$shouldWrite = 0 if (exists($parms->{PARMS}{FEP}) and exists($parms->{PARMS}{FEP}{valid}) and $parms->{PARMS}{FEP}{valid} == 1);
	$shouldWrite = 0 if (exists($parms->{PARMS}{ALL_TINKER}) and $parms->{PARMS}{ALL_TINKER} == 1);

	if($shouldWrite) {
		print $datFile "\nPair Coeffs\n";
		for $index (1 .. $parms->{ATOMTYPES}{counter}) {
			printf $datFile "\n%5d ", $index;
			for $j (@{ $DATA->{$index}{VALS} }) {
				last if (! defined($j) or $j eq "" or $j =~ /\#/);
				if (IsDecimal($j)) {
					printf $datFile "%10.5f ", $j;
				} else {
					printf $datFile "%10d ", $j;
				}
			}
			print $datFile '# ' . $DATA->{$index}{KEY} if ($datWrite->{Labels});
		}
	} elsif ($parms->{PARMS}{HYBRID_VDW}  == 0) {
		for $index (1 .. $parms->{ATOMTYPES}{counter}) {
			next if (! exists($DATA->{$index}{KEY}));
			next if ($DATA->{$index}{TYPE} =~ /(SW|REAX|REXPON|EAM)/);
			$parms->{PARMS}{PAIR_COEFFS} .= sprintf("\npair_coeff %5d %5d ", $index, $index);
			$parms->{PARMS}{PAIR_COEFFS} .= sprintf("%15s ", $DATA->{$index}{Lammps}{LMPNAME})
				if (scalar(keys %{ $parms->{VDW}{TYPE} }) > 1);
			for $j (@{ $DATA->{$index}{VALS} }) {
				last if (! defined($j) or $j eq "" or $j =~ /\#/);
				if (IsDecimal($j)) {
					$parms->{PARMS}{PAIR_COEFFS} .= sprintf("%15.8f ", $j);
				} else {
					$parms->{PARMS}{PAIR_COEFFS} .= sprintf("%15d ", $j);
				}
			}
			$parms->{PARMS}{PAIR_COEFFS} .= sprintf('# %s',$DATA->{$index}{KEY}) if ($datWrite->{Labels});
		}
	}

	$valence = "";   
	for $curr ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS", "PI_TORSIONS") {
		next if (exists($parms->{$curr}) and $parms->{$curr}{counter} == 0);
		if ($curr =~ /^(.)(.*)_(.)(.*)$/) {
			$valence = uc($1) . lc($2) . uc($3) . lc(substr($4,0,-1));
		} else {
			$valence = uc(substr($curr, 0, 1)) . lc (substr($curr,1,-1));
		}
		$valence = "Dihedral" if ($valence eq "Torsion");
		$valence = "Improper" if ($valence eq "Inversion");
		($DATA, $count) = sortParmByIndex($parms->{$curr});
		if (scalar(keys %{ $parms->{$curr}{TYPE} }) == 1 and ! $datWrite->{InputCoeffs}) {
			print $datFile "\n\n${valence} Coeffs\n";
			for $index (1 .. $parms->{$curr}{counter}) {
				printf $datFile "\n%5d ", $index;
				if (! $parms->{PARMS}{single_inversion} and $curr eq "INVERSIONS" and 
					$DATA->{$index}{Lammps}{name} =~ /umbrella/i) {
					$DATA->{$index}{VALS}[0] /= 3;
				}
				for $j (@{ $DATA->{$index}{VALS} }) {
					#last if (! defined($j) or $j eq "" or $j =~ /\#/);
					printf $datFile "%8g ", $j;
				}
				print $datFile '# ' . $DATA->{$index}{KEY} if ($datWrite->{Labels});
				for $i (keys %{ $crossterms}) {
					if(exists($DATA->{$index}{$i}) and exists($DATA->{$index}{$i}{VALS})) {
						$crossterms->{$i}{str} .= sprintf("\n%5d ", $index);
						for $j (@{ $DATA->{$index}{$i}{VALS} }) {
							$crossterms->{$i}{str} .= sprintf("%8g ", $j);
						}
						$crossterms->{$i}{str} .= " # " . $DATA->{$index}{KEY} if ($datWrite->{Labels});
					}
				}
			}
			
		} else {
			$valence = lc $curr;
			chop $valence;
			$valence = "dihedral" if ($valence eq "torsion");
			$valence = "improper" if ($valence eq "inversion");
			$parms->{HYBRID_VALENCE} .= "\n";
			for $index (1 .. $parms->{$curr}{counter}) {
				$parms->{HYBRID_VALENCE} .= sprintf("${valence}_coeff\t%5d ",$index);
				$parms->{HYBRID_VALENCE} .= sprintf("%15s ",$DATA->{$index}{Lammps}{name}) 
					if(scalar(keys %{ $parms->{$curr}{TYPE} }) > 1);
				if (! $parms->{PARMS}{single_inversion} and $curr eq "INVERSIONS" and
					$DATA->{$index}{Lammps}{name} =~ /umbrella/i) {
						$DATA->{$index}{VALS}[0] /= 3;
				}
				for $j (@{ $DATA->{$index}{VALS} }) {
					#last if (! defined($j) or $j eq "" or $j =~ /\#/);
					if (IsDecimal($j) and $j !~ /\.(0)+$/) {
						$parms->{HYBRID_VALENCE} .= sprintf("%8.5f ", $j);
					} elsif($j =~ /[a-zA-Z]/) {
						$parms->{HYBRID_VALENCE} .= sprintf("%8s ", $j);
					} else {
						$parms->{HYBRID_VALENCE} .= sprintf("%8d ", $j);
					}
				}
				$parms->{HYBRID_VALENCE} .= sprintf('# ' . $DATA->{$index}{KEY}) if ($datWrite->{Labels});
				$parms->{HYBRID_VALENCE} .= "\n";

				for $i (keys %{ $crossterms}) {
					if(exists($DATA->{$index}{$i}) and exists($DATA->{$index}{$i}{VALS})) {
						$crossterms->{$i}{str} .= sprintf("%5d ", $index);
						for $j (@{ $DATA->{$index}{$i}{VALS} }) {
							$crossterms->{$i}{str} .= sprintf("%8.3f ", $j);
						}
						$crossterms->{$i}{str} .= " # " . $DATA->{$index}{KEY} if ($datWrite->{Labels});
						$crossterms->{$i}{str} .= "\n";
					}
				}
			}
			$parms->{HYBRID_VALENCE} .= "\n";
		}
	}
	for $i (keys %{ $crossterms}) {
		if(length($crossterms->{$i}{str}) > 0) {
			print $datFile $crossterms->{$i}{header};
			print $datFile $crossterms->{$i}{str};
		}
	}
	printf $datFile "\n\n";
	print "Step 5: Writing data file...All Done                                 \n";
}

sub printAtoms {
	my ($atoms, $box, $datFile, $qeq) = @_;
	my ($counter, $type_id, $atm_name, $fmt, $out_string, $index, $dim, %IMAGE, $molid);

	$fmt = "%8d %8d %8d %11.8f %10.5f %10.5f %10.5f %4d %4d %4d\n"; #full: atom-ID molecule-ID atom-type q x y z
	#$fmt = "%8d %8d %11.8f %10.5f %10.5f %8d %11.5f %11.5f %5d %5d %5d %4d %4d %4d\n" 
	#	if($qeq->{Opt} == 2000); #hybrid full pqeq: atom-ID atom-type x y z molecule-ID q q Rsx Rsy Rsz-defunct
	for $dim ("X", "Y", "Z") {
		$box->{$dim}{"LEN"} = $box->{$dim}{"hi"} - $box->{$dim}{"lo"};
	}

	print $datFile "Atoms\n\n";
	$index = 1;
	for $counter (sort numerically keys %{ $atoms } ) {
		%IMAGE = ();
		for $dim ("X", "Y", "Z") {
			#if ($atoms->{$counter}{$dim . "COORD"} > $box->{LAMMPS}{$dim}{"hi"}) {
				#$IMAGE{$dim} = int($atoms->{$counter}{$dim . "COORD"}/$box->{$dim}{"LEN"}) + 1;
				#$IMAGE{$dim} = getImage($atoms->{$counter}{$dim . "COORD"}, $box->{LAMMPS}{$dim}{"hi"}, $box->{$dim}{"LEN"}, 1);
				#$atoms->{$counter}{$dim . "COORD"} -= (($IMAGE{$dim} - 1) * $box->{$dim}{"LEN"});
				#$atoms->{$counter}{$dim . "COORD"} -= ($IMAGE{$dim} * $box->{$dim}{"LEN"});
			#} elsif ($atoms->{$counter}{$dim . "COORD"} < $box->{LAMMPS}{$dim}{"lo"}) {
				#$IMAGE{$dim} = -1 * int(abs($atoms->{$counter}{$dim . "COORD"})/$box->{$dim}{"LEN"}) - 1;
				#$IMAGE{$dim} = getImage($atoms->{$counter}{$dim . "COORD"}, $box->{$dim}{"lo"}, $box->{$dim}{"LEN"}, 0);
				#$atoms->{$counter}{$dim . "COORD"} -= ($IMAGE{$dim} * $box->{$dim}{"LEN"});
			#} else {
				$IMAGE{$dim} = 0;
			#}
		}
		#$IMAGE{X} = $IMAGE{Y} = $IMAGE{Z} = 0 if ($PARMS->{PARMS}{isTriclinic});

		$molid = ${$atoms->{$counter}{MOLECULEID}};
		#$molid = $atoms->{$counter}{RESNUM};
		if ($qeq->{Opt} != 2000 or 1 < 2) {
			$out_string = sprintf($fmt, $index, $molid, $atoms->{$counter}{"PARMS"}{"INDEX"}, 
						$atoms->{$counter}{"CHARGE"}, $atoms->{$counter}{"XCOORD"}, $atoms->{$counter}{"YCOORD"}, 
						$atoms->{$counter}{"ZCOORD"}, $IMAGE{"X"}, $IMAGE{"Y"}, $IMAGE{"Z"});
#		} else { #defunct
#			$out_string = sprintf($fmt, $index, $atoms->{$counter}{"PARMS"}{"INDEX"}, 
#						$atoms->{$counter}{"XCOORD"}, $atoms->{$counter}{"YCOORD"}, 
#						$atoms->{$counter}{"ZCOORD"}, $atoms->{$counter}{"CHARGE"},
#						$atoms->{$counter}{"CHARGE"}, $atoms->{$counter}{RSXS}{X},
#						$atoms->{$counter}{RSXS}{Y},$atoms->{$counter}{RSXS}{Y},
#						$IMAGE{"X"}, $IMAGE{"Y"}, $IMAGE{"Z"});
		}

		print $datFile $out_string;
		$index++;
	}

	print $datFile "\n";

}

sub getValenceParms {
	my ($atoms, $cons, $parms, $atmList) = @_;
	my ($i, $j, $k, $l, $m, $n, @currIndices, $PLIST, $tmp);

	for $i (@{ $atmList }) {
		if ($#{ $cons->{$i} } == 2) { #Inversion
			@currIndices = ($i,@{ $cons->{$i} });
			&searchForValence($atoms, $parms->{INVERSIONS}, \@currIndices, "INVERSIONS", 0, 1);
		}
		for $j (@{ $cons->{$i} }) {
			next if ($j == $i);
			if ($j > $i) {
				@currIndices = ($i, $j);
				&searchForValence($atoms, $parms->{BONDS}, \@currIndices, "BONDS", 1, 1);
			}
			for $k (@{ $cons->{$j} }) {
				next if (($k == $i) or ($k == $j));
				if ("${i}${j}${k}" > "${k}${j}${i}") {
					@currIndices = ($i, $j, $k);
					&searchForValence($atoms, $parms->{ANGLES}, \@currIndices, "ANGLES", 1, 1);
				}
				for $l (@{ $cons->{$k} }) {
					next if (($l == $i) or ($l == $j) or ($l == $k));
					if ("${i}${j}${k}${l}" > "${l}${k}${j}${i}") {
						@currIndices = ($i, $j, $k, $l);
						&searchForValence($atoms, $parms->{TORSIONS}, \@currIndices, "TORSIONS", 1, 1);
					}
					for $m (@{ $cons->{$l} }) {
						next if (($m == $i) or ($m == $j) or ($m == $k) or ($m == $l));
						if("${i}${j}${k}${l}${m}" > "${m}${l}${k}${j}${i}") { 
							@currIndices = ($i, $j, $k, $l, $m);
							&searchForValence($atoms, $parms->{CMAP}, \@currIndices, "CMAP", 2, 0); #CMAP 5 body term
							&searchForValence($atoms, $parms->{BI_TORSIONS}, \@currIndices, "BI_TORSIONS", 1, 0); #Tinker 5 body term
						}
						for $n (@{ $cons->{$m} }) {
							next if (($n == $i) or ($n == $j) or ($n == $k) or ($n == $l) or ($n == $m));
							if("${i}${j}${k}${l}${m}${n}" > "${n}${m}${l}${k}${j}${i}") { 
								@currIndices = ($i, $j, $k, $l, $m, $n);
								&searchForValence($atoms, $parms->{PI_TORSIONS}, \@currIndices, "PI_TORSIONS", 1, 0); #Tinker 6 body term
							}
						}
					}
				}
			}
		}
	}
}

sub searchForValence {
	my ($atoms, $parmList, $indices, $TYPE, $parmTypeFlag, $writeError) = @_;
	my (@currTypes, $result, $error_code, $curr, $tmp); 
	my ($count, $IndexList, $valType, $bestParm, $i);

	$writeError = 1 if (! defined($writeError));
	@currTypes = ();
	$result = ();
	$error_code = "";
	for $count (@{ $indices }) {
		push @currTypes, $atoms->{$count}{FFTYPE};
		$error_code .= $atoms->{$count}{FFTYPE} . "-";
	}
	chop $error_code;
	if ($parmTypeFlag == 1) { #all except inversions and cmap
		@{ $IndexList->[0] } = @currTypes;
		@currTypes = reverse @currTypes;
		@{ $IndexList->[1] } = @currTypes;
	} elsif($parmTypeFlag == 0) { #inversion
		$IndexList = [getPermutations(\@currTypes)];
	} else { #cmap
		@{ $IndexList->[0] } = reverse @currTypes;
	}

	for $curr (@{ $IndexList }) {
		($tmp, $count) = findValenceType($parmList, $curr, 0);
		if (keys %{ $tmp }) {
			for $i (keys %{ $tmp }) {
				push @{ $result }, $tmp->{$i} if (keys %{ $tmp->{$i} });
			}
		}
	}

	if ($#{ $result } > -1) {
		$valType = '$' . $TYPE;
		$bestParm = getLeastX($result);
		for $i (values %{ $bestParm }) {
			next if(exists($i->{ISCHILD}) or $i->{IGNORE});
			if ($parmTypeFlag > 0) { 
				&saveValence(\%{ eval($valType) }, $i, $indices, \%{ $parmList });
			} else { #have to figure out atom sequence for inversion only
				$IndexList = getCorrectInversionIndices($i, $indices, $atoms);
				for (@{ $IndexList }) {
					&saveValence(\%{ eval($valType) }, $i, $_, \%{ $parmList });
				}
			}
		}
		if ($TYPE eq "TORSIONS") {
			($indices->[1], $indices->[2]) = ($indices->[2], $indices->[1]) 
				if ($indices->[1] > $indices->[2]);
			$ScaleTorsionList->{$indices->[1]}{$indices->[2]} = 0 
				if (! exists($ScaleTorsionList->{$indices->[1]}{$indices->[2]}));
			$ScaleTorsionList->{$indices->[1]}{$indices->[2]}++;
		}
	}
	$ERRORS{$TYPE}{$error_code}++ if (($#{ $result } == -1 or ! keys %{ $bestParm }) and $writeError);		
}

sub findValenceType {
	my ($MPARM, $type_keys, $count) = @_;
	my ($curr_key, $results, $curr_parm, %DAT, @tmp);
	my ($i, $xCount, $minXCount, @junk, $j, @new_keys);

	$curr_parm = $MPARM;
	$curr_key = $type_keys->[0];
	for $j ($curr_key,"X", keys %{ $PARMS->{EQUIVALENCE}{$curr_key} }) {
		$curr_parm = $MPARM;
		next if (! exists($curr_parm->{$j}));
		$curr_parm = \%{ $curr_parm->{$j} };
		if ($#{ $type_keys } == 0) {
			$count++;
			$DAT{$count} = $curr_parm;
			return (\%DAT, $count);
		}
		@new_keys = @{ $type_keys };
		shift @new_keys;
		($results, $count) = findValenceType($curr_parm, \@new_keys, $count);
		@tmp = keys %{ $results };
		if ($#tmp > -1) { # multiple so get one with least amount of Xs
			for $i (@tmp) {
				$DAT{$i} = $results->{$i};
			}
		} 
	}
	return (\%DAT, $count);
}

sub findValenceType_old {
	my ($MPARM, $type_keys, $count) = @_;
	my ($curr_key, $results, $curr_parm, %DAT, @tmp);
	my ($i, $xCount, $minXCount, @junk, $j, @new_keys);

	$curr_parm = $MPARM;
	$curr_key = $type_keys->[0];
	$curr_key = $PARMS->{EQUIVALENCE}{$curr_key}
		if (! exists($curr_parm->{$curr_key}) and exists($PARMS->{EQUIVALENCE}{$curr_key}));
	if (exists($curr_parm->{$curr_key}) or exists($curr_parm->{X})) {
		for $j ($curr_key,"X") {
			$curr_parm = $MPARM;
			next if (! exists($curr_parm->{$j}));
			$curr_parm = \%{ $curr_parm->{$j} };
			if ($#{ $type_keys } == 0) {
				$count++;
				$DAT{$count} = $curr_parm;
				return (\%DAT, $count);
			}
			@new_keys = @{ $type_keys };
			shift @new_keys;
			($results, $count) = findValenceType($curr_parm, \@new_keys, $count);
			@tmp = keys %{ $results };
			if ($#tmp > -1) { # multiple so get one with least amount of Xs
				for $i (@tmp) {
					$DAT{$i} = $results->{$i};
				}
			} 
		}
		return (\%DAT, $count);
	} else {
		return (undef, $count);
	}
}

sub getLeastX {
	my ($parmList) = $_[0];
	my ($i, @tmp, $xCount, $minXCount, $retParm, @pKeys, $j, $index);

	$minXCount = 99999;
	for $index (@{ $parmList }) {
		@tmp = keys %{ $index };
		next if (! @tmp);
		$i = $index->{$tmp[0]};
		@pKeys = split /\s+/, substr($i->{KEY},0,-1);
		$xCount = 0;
		for $j (@pKeys) {
			$xCount++ if ($j eq "X");
		}
		if ($xCount < $minXCount) {
			$retParm = $index;
			$minXCount = $xCount;
		}
	}

	return $retParm;
}

sub saveValence {
	my ($TYPE, $VAL, $atomList, $parmType) = @_;
	my ($rec, $i, $type, $counter, $j);

	$counter = $TYPE->{counter};
	if (! $VAL->{USED}) {
		$VAL->{USED} = 1;
		$VAL->{INDEX} = ++$parmType->{counter};
	}
	$rec->{DATA} = $VAL;
	for $i (@{ $atomList }) {
		push @{ $rec->{ATOMS} }, $i;
		#$type = $ATOMS->{$i}{FFTYPE};
	}
	$counter++;
	$TYPE->{LIST}{ $counter } = $rec;
	if (exists($VAL->{NEXT}) and defined($VAL->{NEXT})) {
		for $i (@{ $VAL->{NEXT} }) {
			$rec = ();
			$i->{USED} = 1;
			$rec->{DATA} = $i;
			for $j (@{ $atomList }) {
				push @{ $rec->{ATOMS} }, $j;
			}
			$counter++;
			$TYPE->{LIST}{ $counter } = $rec;
		}
	}
	$TYPE->{counter} = $counter;
}

sub getCorrectInversionIndices {
	my ($inversionData, $aList, $atomData) = @_;
	my ($centralIndex, $atomList, $centralAtom, @inversionList, @indexCombo, $i, $j, @inversionKey, $sortedList);

	@inversionKey = split /\s+/, substr($inversionData->{KEY},0,-1);
	# for inversion of IJKL
	if ($inversionData->{TYPE} eq "UMBRELLA") { # umbrella torsion
		$centralIndex = 0; #central atom is index[1] J
	} elsif ($inversionData->{TYPE} eq "IT_JIKL") { # amber torsion
		$centralIndex = 2; #central atom is index[1] J
	} else { # charmm torsion
		$centralIndex = 0; # central atom is index[0] I
	}
	$centralAtom = $aList->[0];
	@{ $atomList } = @{ $aList };
	shift @{ $atomList }; # remove the central atom index from the atom list
	shift @inversionKey; # remove the central atom type form the key list
	@{ $sortedList } = sort numerically @{ $atomList }; # sort the remainding atom indices in decreasing order
	$atomList = ();

	for $i (0 .. $#inversionKey) { # search of first atom of current atom type
		next if ($inversionKey[$i] =~ /^X*$/i);
		for $j (0 .. $#{ $sortedList }) { # atom with lower id checked first
			if ($atomData->{ $sortedList->[$j] }{FFTYPE} eq $inversionKey[$i]) { # found atom type for current type
				$atomList->[$i] = $sortedList->[$j]; # save atom index in list
				splice @{ $sortedList }, $j, 1; # remove atom index from sorted list
				last; # exit sorted list loop
			}
		}
	}

	for $i (0 .. $#inversionKey) {
		if (! $atomList->[$i]) { # atom key was not found, must have been an X
			$atomList->[$i] = shift @{ $sortedList }; # place the smallest atom index
		}
	}

	@indexCombo = ([$atomList->[0], $atomList->[1], $atomList->[2]],
			[$atomList->[2], $atomList->[0], $atomList->[1]],
			[$atomList->[1], $atomList->[2], $atomList->[0]]); # list for all 3 inversion about central atom
	for $j (0 .. $#indexCombo) {
		$atomList = $indexCombo[$j];
		for $i (0 .. ($centralIndex - 1)) { # add all atoms before central
			$inversionList[$j][$i] = shift(@{ $atomList });
		}
		$inversionList[$j][$centralIndex] = $centralAtom; # add central atom
		while (@{ $atomList }) { # add all atoms after central
			push @{ $inversionList[$j] }, shift @{ $atomList };
		}
		last if ($PARMS->{PARMS}{single_inversion}); # end if only using single inversion
	}

	return \@inversionList;
}

sub sortParmByIndex {
	my ($parm, $jVal) = @_;
	my (%INDICES, $i, $PLIST, $count, $j, $tmp);

	$PLIST = ();
	$count = 0;
	$tmp = 0;
	&getValParent($parm, \%{ $PLIST }, \$tmp);
	$count = 0;
	for $i (keys %{ $PLIST }) {
		for $j (keys %{ $PLIST->{$i} }) {
			next if (defined($jVal) and ($j != $jVal));
			if (exists($PLIST->{$i}{$j}{INDEX}) and defined($PLIST->{$i}{$j}{INDEX}) and 
			   (! exists($PLIST->{$i}{$j}{IGNORE}) or ! $PLIST->{$i}{$j}{IGNORE}) and 
				$#{$PLIST->{$i}{$j}{VALS}} > -1 and @{ $PLIST->{$i}{$j}{VALS} }) {
				$count++;
				#$INDICES{$count} = \%{ $PLIST->{$i}{$j} };
				$INDICES{ $PLIST->{$i}{$j}{INDEX} } = \%{ $PLIST->{$i}{$j} };
			}
		}
	}
	return (\%INDICES, $count);
}

sub saveParmIndex {
	my ($parms, $i, $index) = @_;

	$parms->{ATOMTYPES}{$i}{INDEX} = $$index;
	$parms->{ATOMTYPES}{$i}{NAME} = $i;
	$parms->{ATOMTYPES}{$i}{TYPEID} = $$index;
	$$index++;

}

sub deleteParmIndex {
	my ($parms, $i) = @_;
	my ($tmp, $j);

	@{ $tmp } = keys %{ $parms->{ATOMTYPES} };
	for $j (@{ $tmp }) {
		delete $parms->{VDW}{$j}{$i};
	}
	delete $parms->{ATOMTYPES}{$i};
	delete $parms->{VDW}{$i};
	delete $parms->{QEq}{$i} if ($parms->{QEq}{Opt} > 0);
}

sub updateParmIndex {
	my ($parms) = $_[0];
	my ($i, $count, $PLIST, $index, @tmp, $j, $l, $k, $lammpsType, $dflag, $d);
	my ($scale_14_coul, $scale_14_vdw, $parmIgnore, $makeZero, $nNoMix, $vdwTypes);

	$parms->{PARMS}{USE_HBOND} = $makeZero = 0;
	$index = 1;
	map { $parms->{ATOMTYPES}{$_}{USED} == 1 ? saveParmIndex($parms, $_, \$index):deleteParmIndex($parms,$_) } 
		sort { $parms->{ATOMTYPES}{$a}{TYPEID} <=> $parms->{ATOMTYPES}{$b}{TYPEID} } keys %{ $parms->{ATOMTYPES} };

	&setVDWtypes($parms);
	map { $vdwTypes->{$_} = 1 } grep {!/TYPE|counter/i} keys %{ $parms->{VDW} };
	for $j (keys %{ $vdwTypes }) {
		for $i (keys %{ $vdwTypes }) {
			$parmIgnore = 0;
			for $l (keys %{ $parms->{VDW}{$j}{$i} }) {
				$k = $parms->{VDW}{$j}{$i}{$l};
				next if (! keys %{ $k } or $k->{IGNORE} == 1);
				if ($k->{TYPE} eq "DREIDHB") { # dreiding hb fix
					unshift @{ $k->{VALS} }, $j;
				}
				if ($k->{IT} eq "hbond") {
					if (! $parms->{ATOMTYPES}{ $k->{VALS}[0] }{USED}) {
						delete $parms->{VDW}{$j}{$i}{$l};
						next;
					}
					if (! $parms->{BONDS}{$i}{$k->{VALS}[0]}{1}{USED} and
						! $parms->{BONDS}{$k->{VALS}[0]}{$i}{1}{USED} and
						! $parms->{BONDS}{$j}{$k->{VALS}[0]}{1}{USED} and
						! $parms->{BONDS}{$k->{VALS}[0]}{$j}{1}{USED}) {
						delete $parms->{VDW}{$j}{$i}{$l};
						next;
					}
					$k->{VALS}[0] = $parms->{ATOMTYPES}{ $k->{VALS}[0] }{INDEX}; # replace the hydrogen atom type with the type index
					$dflag = "j";
					$dflag = "i" if ($k->{DONOR} eq $i and $parms->{ATOMTYPES}{$i}{INDEX}<=$parms->{ATOMTYPES}{$j}{INDEX});
					$dflag = "i" if ($k->{DONOR} eq $j and $parms->{ATOMTYPES}{$i}{INDEX}>$parms->{ATOMTYPES}{$j}{INDEX});
					unshift @{ $k->{VALS} }, $dflag;
					# format is i_type j_type hydrogen_type is_itype_donor d0 (alpha for morse) r0 cos_power 
					($k->{VALS}[0], $k->{VALS}[1]) = ($k->{VALS}[1], $k->{VALS}[0]);
					$parms->{PARMS}{HAS_HBONDS} = 1;
					$parms->{PARMS}{USE_HBOND} = 1;
					$parms->{PARMS}{USE_HBOND} = 2 if ($k->{TYPE} eq "MORSE_COSP");
				}
				#delete $parms->{VDW}{$i}{$j}{$l} if ($k->{IGNORE} or ! exists($k->{VALS}) or scalar(@{ $k->{VALS} }) == 0);
			}
			next if ($i ne $j);
			for $k (values %{ $parms->{VDW}{$j}{$i} }) {
				$k->{INDEX} = $parms->{ATOMTYPES}{$j}{INDEX};
				$k->{DATA} = \%{ $parms->{ATOMTYPES}{$j} };
			}
			$makeZero = 1
				if($i eq $j and scalar(keys %{ $parms->{VDW}{$j}{$i} }) == 1 and $parmIgnore == 1);
			if(! $makeZero and $parmIgnore == 1) {
				$nNoMix = 1;
				for $l (values %{ $parms->{VDW}{$j}{$i} }) {
					$nNoMix++ if ($l->{Lammps}{PAIRMIX} == 0);
				}
				$makeZero = 1 if (scalar(keys %{ $parms->{VDW}{$j}{$i} }) == $nNoMix);
			}
		}
	}
	for $i (keys %{ $parms->{ATOMTYPES} }) {
		delete $parms->{ATOMTYPES}{$i} if (! keys %{ $parms->{ATOMTYPES}{$i} } or ! $parms->{ATOMTYPES}{$i}{USED});
	}

	$parms->{ATOMTYPES}{counter} = $parms->{VDW}{counter} = $index - 1;
	&WriteSWfile($parms,$suffix) if (exists($parms->{PARMS}{SW}) and $parms->{PARMS}{SW});
	$parms->{PARMS}{makeZero} = $makeZero;

	#buck/coul/long fix
	if(scalar(keys %{ $parms->{VDW}{TYPE} }) == 1 and exists($parms->{VDW}{TYPE}{"buck/coul/long"})) {
		$sysOpts->{WriteOffDiagPairs} = 1;
		#now we add a soft potential to prevent exponential 6 catastrope
		for $i (keys %{ $vdwTypes }) {
			for $j (keys %{ $vdwTypes }) {
				next if (! exists($parms->{VDW}{$i}) or ! exists($parms->{VDW}{$i}{$j}) or ! keys %{ $parms->{VDW}{$i}{$j}});
				&addSoftPair($parms, $parms->{VDW}{$i}{$j}, $i, $j,[10]);
			}
		}
	}

	for $i ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS", "PI_TORSIONS", "BI_TORSIONS") {
		$index = 0;
		$count = 0;
		$PLIST = ();
		&getValParent($parms->{$i}, \%{ $PLIST }, \$count);
		for $j (sort numerically keys %{ $PLIST }) {
			for $l (sort alphabetically keys %{ $PLIST->{$j} }) {
				if (! $PLIST->{$j}{$l}{USED} or $PLIST->{$j}{$l}{IGNORE}) {
					delete $PLIST->{$j}{$l};
					next;
				}
				$index++;
				$PLIST->{$j}{$l}{INDEX} = $index if($i =~ /TORSIONS|INVERSIONS/);
				if($i eq "TORSIONS") {
					$scale_14_coul = substr($PLIST->{$j}{$l}{COUL},3);
					$scale_14_vdw = substr($PLIST->{$j}{$l}{VDW},3);
					$parms->{"PARMS"}{"scale_cou_14"} = $scale_14_coul;
					$parms->{"PARMS"}{"scale_vdw_14"} = $scale_14_vdw;
					if(($PLIST->{$j}{$l}{TYPE} eq "SHFT_DIHDR" or $PLIST->{$j}{$l}{TYPE} eq "DIHEDRAL") and 
						($parms->{QEq}{Opt} == 2 or ! isCharmm($PLIST->{$j}{$l}{VDW}, 0) or ! isCharmm($PLIST->{$j}{$l}{COUL}, 1))) { 
						$PLIST->{$j}{$l}{Lammps}{name} = "harmonic";
						$d = 1;
						$d = -1 if ($PLIST->{$j}{$l}{VALS}[2]==180);
						@tmp = ($PLIST->{$j}{$l}{VALS}[0],$d,$PLIST->{$j}{$l}{VALS}[1]);
						@{ $PLIST->{$j}{$l}{VALS} } = @tmp;
					}
				}
				$lammpsType = $PLIST->{$j}{$l}{Lammps}{name} . " " . $PLIST->{$j}{$l}{Lammps}{opts} if (exists($PLIST->{$j}{$l}{Lammps}));
				$lammpsType =~ s/\s+$//;
				$parms->{$i}{TYPE}{$lammpsType} = 1 if (exists($PLIST->{$j}{$l}{Lammps}));
			}
		}
		&removeEmptyParms($parms->{$i});
		$parms->{$i}{counter} = $index;
	}
	&removeEmptyParms($parms->{CMAP}) if (defined($CMAP));
}

sub addAllPair {
	my ($parms, $pName, $optStr) = @_;
	my (@tmp, $k, $curr, $allCoul);

	$parms->{PARMS}{$pName} = 1;
	if(exists($parms->{VDW}{-1}{-1})) {
		@tmp = sort numerically keys %{ $parms->{VDW}{-1}{-1} };
		$k = $tmp[$#tmp] + 1;
	} else {
		$k = 0;
	}

	$optStr = 0 if (! defined($optStr));
	$optStr = $1 if ($optStr =~ /(\d+\.?\d*)\s*$/);
	$curr = \%{ $parms->{VDW}{-1}{-1} };
	$parms->{VDW}{TYPE}{lc($pName)} = (
									{
										"OPTS"    => "$optStr",
										"NAME"    => $pName,
										"PAIRMIX" => 0,
										"ORDER"   => 0,
										"LMPNAME" => lc $pName,
										"LMPONAME" => lc $pName,
										"USED"     => 1,
										"SCC"      => 1,
									}
								);
	$curr->{$k} = (
					{
						"ATOM"     => "* * ",
						"IGNORE"   => 0,
						"IT"       => "vdw",
						"KEY"      => "*",
						"Lammps"   => \%{ $parms->{VDW}{TYPE}{lc $pName} },
						"TYPE"     => $pName,
						"USED"     => 1,
						"DATA"     => undef,
						"INDEX"    => $k,
						"VALS"     => [0],
					  }
					);
	$parms->{PARMS}{HYBRID_VDW} = 2;
	#$curr->{$pName}{LMPNAME} = "coul/long";
	#$curr->{$pName}{OPTS} = "10";
}

sub determineAllCoul {
	my ($parms) = $_[0];
	my ($i, $j, $k, $hasCoulAll);

	$i = $j = $hasCoulAll = 0;
	for $k (keys %{ $parms->{VDW}{TYPE} }) {
		$hasCoulAll = 1 if ($k =~ /^coul/);
		next if ($parms->{VDW}{TYPE}{$k}{SCC} == 1); #skip coul check
		$i++;
		if($parms->{VDW}{TYPE}{$k}{LMPONAME} =~ /\/coul\/long|\/coul\/cut/i) {
			$j++;
		}
	}
	&fixCoulLong($parms) if($hasCoulAll or ($i > 1 and $i == $j)); # and ! defined($sysOpts->{fep}));
}

sub fixCoulLong {
	my ($parms) = $_[0];
	my ($k, $l, $cType, $curr, $isCoulLong);

	$isCoulLong = 0;
	$cType = "long";
	$curr = \%{ $parms->{VDW}{TYPE} };
	for $k (keys %{ $curr }) {
		$cType = "cut" if ($curr->{$k}{LMPNAME} =~ /\/coul\/charmm|\/coul\/cut/);
		$cType = "long/cs" if(exists($parms->{PARMS}{POLARIZATION}) and $parms->{PARMS}{POLARIZATION}{type}==0);
		$curr->{$k}{LMPNAME} =~ s/\/coul\/(long|cut|charmm)[\/opt]*//;
		$curr->{$k}{LMPNAME} = "lj/gromacs" if ($curr->{$k}{LMPNAME} eq "lj/charmm");
		$curr->{$k}{LMPNAME} = "lj/cut/soft" if ($curr->{$k}{LMPNAME} eq "lj/cutsoft"); #FEP fix
		$isCoulLong = 1 if ($curr->{$k}{LMPNAME} =~ /^coul/);
		$l = $k;
	}

	&addAllPair($parms,"coul/${cType}", $curr->{$l}{OPTS}) if (! $isCoulLong);
	$parms->{PARMS}{allCoul} = 1;
}

sub isCharmm {
	my ($inStr, $is_coul) = @_;
	my ($returnVal);

	$returnVal = 0;
	$returnVal = 1 if ($inStr =~ /0010.00000/);
	return $returnVal;
}

sub removeCrossParms {
	my ($parms) = $_[0];
	my ($i, $j, $l, $count, $PLIST);

	for $i ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS") {
		$count = 0;
		$PLIST = ();
		getValParent($parms->{$i}, \%{ $PLIST }, \$count);
		for $j (sort numerically keys %{ $PLIST }) {
			for $l (sort alphabetically keys %{ $PLIST->{$j} }) {
				if (exists($PLIST->{$j}{$l}{CTYPE}) and $PLIST->{$j}{$l}{CTYPE} ne "") {
					delete $PLIST->{$j}{$l};
				}
			}
		}
	}
}

sub removeEmptyParms {
	my ($parmData) = $_[0];
	my ($i);

	for $i (keys %{ $parmData }) {
		next if (($i eq "counter") or ($i eq "TYPE"));
		last if (exists($parmData->{$i}{INDEX}));
		if (! keys %{ $parmData->{$i} }) {
			delete $parmData->{$i};
		} else {
			&removeEmptyParms($parmData->{$i});
			delete $parmData->{$i} if (! keys %{ $parmData->{$i} });
		}
	}
}

sub getPermutations {
	my ($inArray) = @_;
	my (@PERMS, $firstAtm, $i);

	$firstAtm = shift @{ $inArray };
	@PERMS = Permutate($inArray, []);
	for $i (0 .. $#PERMS) {
		unshift @{ $PERMS[$i] }, $firstAtm;
	}
	return @PERMS;
}

sub updateTorsionList {
	my ($torsionList) = $_[0];
	my ($i, $j, $k, @tmp, $l);
	my ($TLIST, $count);


	$count = 0;
	$TLIST = ();
	&getValParent($torsionList, \%{ $TLIST }, \$count);

	for $i (keys %{ $TLIST }) {
		for $j (keys %{ $TLIST->{$i} }) {
			if ($j eq "") {
				delete $TLIST->{$i}{$j};
				$count--;
				next;
			}
			if ($TLIST->{$i}{$j}{NUM} > 1) { #multiple torsions
				@tmp = @{ $TLIST->{$i}{$j}{VALS} };
				$TLIST->{$i}{$j}{VALS} = ();
				for $l (0 .. ($TLIST->{$i}{$j}{PER} - 1)) {
					$TLIST->{$i}{$j}{VALS}[$l] = $tmp[$l];
				}
				if ($TLIST->{$i}{$j}{Lammps}{name} eq "charmm") {
					push @{ $TLIST->{$i}{$j}{VALS} }, getTorsionWeightFactor($TLIST->{$i}{$j});
				}
				$l = $TLIST->{$i}{$j}{PER};
				$k = 1;
				while ($k < $TLIST->{$i}{$j}{NUM}) {
					%{ $TLIST->{$i}{"${j}${k}"} } = %{ $TLIST->{$i}{$j} };
					$TLIST->{$i}{"${j}${k}"}{ISCHILD} = 1;
					$TLIST->{$i}{"${j}${k}"}{INDEX} += $k;
					$TLIST->{$i}{"${j}${k}"}{KEY} = "$k" . $TLIST->{$i}{"${j}${k}"}{KEY};
					$TLIST->{$i}{"${j}${k}"}{VALS} = ();
					for (0 .. ($TLIST->{$i}{$j}{PER} - 1)) {
						$TLIST->{$i}{"${j}${k}"}{VALS}[$_] = $tmp[$_ + $l];
					}
					if ($TLIST->{$i}{$j}{Lammps}{name} eq "charmm") {
						push @{ $TLIST->{$i}{"${j}${k}"}{VALS} }, 0;
					}
					$l += $TLIST->{$i}{$j}{PER};
					push @{ $TLIST->{$i}{$j}{NEXT} }, \%{ $TLIST->{$i}{"${j}${k}"} };
					$k++;
				}
			} else {
				if ($TLIST->{$i}{$j}{Lammps}{name} eq "charmm") {
					push @{ $TLIST->{$i}{$j}{VALS} }, getTorsionWeightFactor($TLIST->{$i}{$j});
				}
			}
		}
	}
}

sub getTorsionWeightFactor {
	my ($torsion) = $_[0];
	my ($typePath, $atomType, $atomTypeIndex, $j, $rs);

	if ($torsion->{FFTYPEID} =~ /^[2|4]$/) { #charmm/dreiding
		@{ $typePath } = split /\s+/, substr($torsion->{KEY},0,-1);
		for $i (0 .. $#{ $typePath }) {
			next if ($typePath->[$i] eq "X");
			$atomType = $typePath->[$i];
			$atomTypeIndex = $i;
			last;
		}
		if (! defined($atomType)) {
			return 1;
		}
		
		$rs = getRingSize($atomType, $atomTypeIndex, $typePath); #get the rings
		return 0 if (! defined($rs));
		if($torsion->{FFTYPEID} =~ /2/) {
			if(exists($rs->{4}) or exists($rs->{5})) {
				#in 4 or 5 membered rings the 1-4 dihedral also is counted as a 1-2 or 1-3 interaction 
				#when going around the ring in the opposite direction and thus the weighting factor is 0.0
				return 0;
			} elsif(exists($rs->{6})) {
				#In 6-membered rings, the same 1-4 interaction would be computed twice (once for the 
				#clockwise 1-4 pair in dihedral 1-2-3-4 and once in the counterclockwise dihedral 1-6-5-4) 
				#and thus the weighting factor has to be 0.5
				return 0.5;
			} else {
				return 1;
			}
		} else { #dreiding
			if($rs == 6) {
				return $PARMS->{PARMS}{"EXO_CYCLIC"}*$torsion->{"1_4scale"};
			} else {
				return $torsion->{"1_4scale"};
			}
		}
	} else {
		return 0;
	}
}

sub getTorsionScalingFactor {
	my ($torsions) = $_[0];
	my ($i, $j, $count, $atom1, $atom4, $index);

	$index = 0;
	for $i (1 .. $torsions->{counter}) {
		next if (! exists($torsions->{LIST}{$i}{DATA}{do_scale}) or ! $torsions->{LIST}{$i}{DATA}{do_scale});
		if (! $index) {
			print "scaling torsions...";
			$index = 1;
		}
		$atom1 = $torsions->{LIST}{$i}{ATOMS}[1];
		$atom4 = $torsions->{LIST}{$i}{ATOMS}[2];
		($atom1, $atom4) = ($atom4, $atom1) if ($atom1 > $atom4);
		$count = 1;
		$count = $ScaleTorsionList->{$atom1}{$atom4} if (exists($ScaleTorsionList->{$atom1}{$atom4}));
		#for $j (1 .. $torsions->{counter}) {
			#$count++ if (($torsions->{LIST}{$j}{ATOMS}[1] == $atom1 and
			#$torsions->{LIST}{$j}{ATOMS}[2] == $atom4) or
			#($torsions->{LIST}{$j}{ATOMS}[2] == $atom1 and
			#$torsions->{LIST}{$j}{ATOMS}[1] == $atom4));
		#}

		next if ($count == 1 or (exists($torsions->{LIST}{$i}{DATA}{CTYPE}) and $torsions->{LIST}{$i}{DATA}{CTYPE} ne ""));
		if (! exists($torsions->{LIST}{$i}{DATA}{scaled})) {
			$torsions->{LIST}{$i}{DATA}{VALS}[0] /= $count;
			$torsions->{LIST}{$i}{DATA}{scaled} = $count;
		} elsif ($torsions->{LIST}{$i}{DATA}{scaled} != $count) { #create new torsion
			$torsions->{LIST}{$i}{DATA} = getNewParm($torsions->{LIST}{$i}{DATA}, "TORSIONS", $count, $i);
		}
	}
}

sub getNewParm {
	my ($parmData, $parmType, $count, $newKey) = @_;
	my ($newParm, $i, $DATA, $j, @pkeys, $index, $oldParm, $retParm);
	
	($DATA, $j) = sortParmByIndex($PARMS->{$parmType});
	for $i (1 .. $PARMS->{$parmType}{counter}) {
		if ($DATA->{$i}{KEY} eq $parmData->{KEY} and
			exists($DATA->{$i}{scaled}) and
			$DATA->{$i}{scaled} == $count) {
			return $DATA->{$i}; # found similar parm so return it
		}
	}
	
	#else create a new parm
	@pkeys = split /\s+/, substr($parmData->{KEY},0,-1); 
	$index = $PARMS->{$parmType}{counter};
	$oldParm = $PARMS->{$parmType}{shift @pkeys};
	for $i (0 .. ($#pkeys -1)) {
		$oldParm = \%{ $oldParm->{$pkeys[$i]} };
	}
	$newKey = $pkeys[$#pkeys] . $newKey;
	$oldParm->{$newKey} = dclone($oldParm->{$pkeys[$#pkeys]});
	for $i (keys %{ $oldParm->{$newKey} }) {
		$newParm = $oldParm->{$newKey}{$i};
		$newParm->{INDEX} = ++$index;
		$newParm->{do_scale} = 1;
		$newParm->{VALS}[0] *= $parmData->{scaled}/$count;
		$newParm->{scaled} = $count;
		$retParm = $newParm if ($newParm->{KEY} eq $parmData->{KEY});
	}
	$PARMS->{$parmType}{counter} = $index;

	return $retParm; #return the new parm
}

sub getRingSize {
	my ($atomType, $atomPos, $typePath) = @_;
	my ($i, $blist, $rings, $j, $atomID, $tmp);

	for $i (0 .. $atomPos) {
		shift @{ $typePath };
	}

	$atomPos = 0;
	$rings = ();
	for $i (values %{ $ATOMS }) {
		if($i->{FFTYPE} eq $atomType) {
			$atomID = $i->{INDEX};
			$j = 6; #max size for ring
			&getAtmBondList($atomID, $atomID, 0, $ATOMS, $CONS, $typePath, \$atomPos, \$j, \%{ $rings });
		}
	}
	return $rings;
}

sub getAtmBondList {
	my ($sID, $aID, $parentID, $atoms, $bonds, $typePath, $pos, $depth, $rings) = @_;

	if ($$depth == 0) {
		$rings->{7} = 1; #maybe larger than 6
		return;
	}

	for $i (@{ $bonds->{$aID}}) {
		next if ($i == $parentID);
		next if(($$pos <= $#{ $typePath }) and ($typePath->[$$pos] ne "X") and ($typePath->[$$pos] ne $atoms->{$i}{FFTYPE}));
		if($i == $sID) {
			$rings->{6 - $$depth + 1} = 1;
			return;
		} else {
			$$depth--;
			$$pos++;
			&getAtmBondList($sID, $i, $aID, $atoms, $bonds, $typePath, $pos, $depth, \%{ $rings });
			$$pos--;
			$$depth++;
		}
	}
	return;
}

sub is5MemberRing_depr {
	#deprecated
	my ($atomType, $atomPos, $types, $count) = @_;
	my ($atom1List, $l, $atom4List, $i, $j, $is5Member, $k);

	if ($atomPos == 0) { #already the first atom
		$atom1List = getAtmList($atomType, 0);
		$atom4List = getAtmList($atomType, 4);
	} elsif ($atomPos == 3) { #already the last atom
		$atom4List = getAtmList($atomType, 0);
		$atom1List = getAtmList($atomType, 4);
	} else {
		$atom1List = getAtmList($atomType, $atomPos);
		$atom4List = getAtmList($atomType, (3 - $atomPos));
	}

	$is5Member = 1;
	MAIN: for $i (values %{ $atom1List }) {
		next if ($ATOMS->{$i}{FFTYPE} ne $types->[0] and $types->[0] ne "X");
		for $j (values %{ $atom4List }) {
			next if (($j == $i) or ($ATOMS->{$j}{FFTYPE} ne $types->[3] and $types->[3] ne "X"));
			for $k (@{ $BONDS->{$i} }) {
				for $l (@{ $BONDS->{$k} }) {
					if ($l == $i) {
						$is5Member = 0;
						last MAIN;
					}
				}
			}
			last MAIN;
		}
	}

	return $is5Member;
}

sub getAtmList {
	my ($atomType, $bondNum) = @_;
	my ($atom, $i, $BONDLIST, $count);

	for $i (keys %{ $ATOMS }) {
		next if ($ATOMS->{$i}{FFTYPE} ne $atomType);
		$atom = $i;
		last;
	}
	return undef if (! defined($atom));
	$BONDLIST = ();
	$count = 0;
	&getBondList(\%{ $BONDLIST }, $atom, \$bondNum, \$count);
	return $BONDLIST;
}

sub getBondList {
	my ($VALS, $atom, $bondC, $count) = @_;
	my ($i);

	if (${ $bondC } == 0) {
		${ $count }++;
		$VALS->{${ $count }} = $atom;
		} else {
		${ $bondC }--;
		for $i (@{ $CONS->{$atom} }) {
			&getBondList($VALS, $i, $bondC, $count);
		}
	}
}
sub getValParent {
	my ($valList, $VList, $counter) = @_;
	my ($rec, $i, $j, $valid);

	for $i (sort alphabetically keys %{ $valList }) {
		next if ($i =~ /counter|type/i or (exists($valList->{$i}{IGNORE}) and $valList->{$i}{IGNORE}));
		if (keys %{ $valList->{$i} }) {
			$valid = 0;
			for $j (sort alphabetically keys %{ $valList->{$i} }) {
				next if ($i =~ /counter|type/i);
				if (exists($valList->{$i}{$j}{VALS})) {
					$valid = 1;
					last;
				}
			}
			if (! $valid) {
				&getValParent($valList->{$i}, $VList, $counter);
			} else {
				${ $counter }++;
				$VList->{${ $counter }} = \%{ $valList->{$i} };
			}
		}
	}
}

sub printErrors {
	my ($errorlist, $fatal) = @_;
	my ($i, $j);

	for $i (keys %{ $errorlist }) {
		if ($fatal) {
			print "\n---===$i ERRORS===----\n";
		} else {
			print "\n---===MISSING $i TERMS===----\n";
		}
		for $j (keys %{ $errorlist->{$i} }) {
			print "$j: (occurred $errorlist->{$i}{$j} times)\n";
		}
	}

	die "The script cannot contine\n" if ($fatal);
}

sub addLammpsParms {
	my ($parms, $atoms, $mols, $bonds, $angles, $cmap, $sprefix, $datWrite, $csOpt) = @_;
	my ($i, $j, $k, $parmHybrid, $offDiag, $diag, $usedType, @tmp); 
	my ($mbff, $isZero, $lmpType, $hybridOpt, $nvdwtypes, $allIJpairs);

	@tmp = keys %{ $inputType->{names} };
	$parms->{PARMS}{SOLUTE} = GetSoluteAtoms($atoms, $mols);
	$parms->{PARMS}{SUFFIX} = $suffix;
	$parms->{PARMS}{FFTYPE} = $ffType;
	$parms->{PARMS}{NUM_FILES} = $#{ $FILES } + 1;
	$parms->{PARMS}{INPUTNAME} = "@tmp";
	$parms->{PARMS}{INPUTLOC} = $inputType->{loc};
	$parms->{PARMS}{NUM_ATOMS} = scalar(keys %{ $atoms });

	#shake bond and angle options
	&getShakeOpts($atoms, $bonds, $angles, $parms);

	#core-shell polarization options
	&getPolarizationInputCmds($parms, $csOpt) if (defined($csOpt));

	#morse with long range electrostatics fix
	&fixMorseCoul($parms);

	#mbff check
	($allIJpairs, $usedType) = checkIJPairs($parms, $diag);

	#coul options
	if(scalar(keys %{ $usedType }) == 1 and $parms->{PARMS}{HAS_CHARGE} == 1) {
		@tmp = keys %{ $usedType };
		$lmpType = $tmp[0];
		if(! $allIJpairs) {
			&fixCoulLong($parms);
		}
	}
	&determineAllCoul($parms, $usedType);

	#add zero
	&addAllPair($parms, "zero") if ($parms->{PARMS}{makeZero} == 1 and ! $parms->{PARMS}{allCoul}); 

	$parmHybrid = getHybridOpts($parms->{VDW});
	if(! exists($parms->{PARMS}{HYBRID_VDW})) {
		$parms->{PARMS}{HYBRID_VDW} = 0 if ($parmHybrid eq "");
		$parms->{PARMS}{HYBRID_VDW} = 1 if ($parmHybrid =~ /hybrid/);
		$parms->{PARMS}{HYBRID_VDW} = 2 if ($parmHybrid =~ /overlay/);
	}

	&doPairMix($parms) if ($parms->{PARMS}{HYBRID_VDW} > 0 or exists($sysOpts->{WriteOffDiagPairs}));
	
	$isZero = $nvdwtypes = 0;
	for $i (keys %{ $parms->{VDW}{TYPE} }) {
		if (! exists($parms->{VDW}{TYPE}{$i}{USED}) or ! $parms->{VDW}{TYPE}{$i}{USED}) {
			delete $parms->{VDW}{TYPE}{$i};
		} else { 
			$isZero = 1 if ($i =~ /zero/);
			$nvdwtypes++;
		}
	}
	#no need for zero if only 1 other atom type
	if($isZero and $nvdwtypes == 2) {
		delete $parms->{VDW}{TYPE}{zero};
		$parms->{PARMS}{HYBRID_VDW} = 0;
		delete $parms->{VDW}{-1};
	}

	#get diag and offdiag components
	($diag, $offDiag) = getVDWlist($parms, $parmHybrid, $datWrite);

	#coul/pqeqgauss fix
	if($parms->{QEq}{Opt} == 2) {
		$offDiag->{"coul/pqeqgauss"}{-1}{-1} = "pair_coeff * * coul/pqeqgauss 0.00000 0.00000 0.0000 0 0.00000 0.0000 0.00000 0 #dummy\n";
		$parms->{PARMS}{allPairs} = 1;
	}
	#mb fix
	for $i (qw/SW REAX REXPON EAM EAM\/FS/) {
		next if (! defined($parms->{VDW}{TYPE}{lc($i)}));
		$parms->{PARMS}{allPairs} = 1;
		$parms->{PARMS}{HYBRID_VDW} = 2;
		$parmHybrid = "hybrid/overlay ";
		if($nvdwtypes == 1) {
			$parmHybrid = "";
			$parms->{PARMS}{HYBRID_VDW} = 0;
		}
		$hybridOpt = "a";
		$hybridOpt = lc $i if ($parmHybrid ne "");
		$mbff = $reaxFF;
		$mbff = "${sprefix}." . lc($i) if ($i !~ /reax|rexpon/i);
		$mbff = $parms->{VDW}{TYPE}{lc($i)}{FILE} if (exists($parms->{VDW}{TYPE}{lc($i)}{FILE}));
		$offDiag->{$hybridOpt}{-1}{-1} = "pair_coeff       * * $hybridOpt $mbff @{ $diag->{lc($i)} }\n" if($parmHybrid);
		$offDiag->{$hybridOpt}{-1}{-1} = "pair_coeff       * * $mbff @{ $diag->{lc($i)} }\n" if(!$parmHybrid);
		$offDiag->{$hybridOpt}{-1}{-1} =~ s/NULL/X/g if ($i =~ /rexpon/i);
	}
	#no need for zero if allPairs used
	if($isZero and exists($parms->{PARMS}{allPairs})) {
		delete $parms->{VDW}{TYPE}{zero};
	}

	#save all pairs data
	if(exists($parms->{VDW}{-1}) and exists($parms->{VDW}{-1}{-1})) {
		for $i (values %{ $parms->{VDW}{-1}{-1} }) {
			next if (! exists($parms->{VDW}{TYPE}{ $i->{TYPE} }));
			$hybridOpt = "a";
			$hybridOpt = lc $i->{TYPE} if ($parmHybrid ne "");
			$offDiag->{$hybridOpt}{-1}{-1} .= "\npair_coeff       * * $i->{TYPE}\n";
		}
	}

	#save off diagnoal elements for later use
	for $k (sort { $parms->{VDW}{TYPE}{$a}{ORDER} <=> $parms->{VDW}{TYPE}{$b}{ORDER} } keys %{ $parms->{VDW}{TYPE} }) {
		$parms->{PARMS}{OFF_DIAG} .= "#$k\n" if ($parmHybrid =~ /overlay/);
		#$k = $hybridOpt if(scalar(keys %{ $parms->{VDW}{TYPE} }) > 1);
		next if (! defined($k) or !exists($offDiag->{$k}));
		next if (! exists($offDiag->{$k}));
		for $i (sort numerically keys %{ $offDiag->{$k} }) {
			for $j (sort numerically keys %{ $offDiag->{$k}{$i} }) {
				$parms->{PARMS}{OFF_DIAG} .= $offDiag->{$k}{$i}{$j};
			}
		}
	}
	$parms->{PARMS}{OFF_DIAG} .= "\n\n";

	#cmap fix
	delete $parms->{PARMS}{CMAP} if (! defined($cmap) or ! keys %{ $cmap } or $cmap->{counter} == 0);

	#include file
	$parms->{PARMS}{include_file} = $sysOpts->{include_file} if (exists($sysOpts->{include_file}));
}

sub fixMissingEquivVDW {
	my ($parms) = $_[0];
	my ($vdwTypes, $i, $j, $k, $l, $e, $eType, $cType, $tmp);

	#get all the types
	map { $vdwTypes->{$_} = 1 } grep {!/TYPE|counter/i} keys %{ $parms->{VDW} };
	for $i (keys %{ $vdwTypes }) {
		undef $e;
		for $j (keys %{ $parms->{EQUIVALENCE}}) {
			$e = $j if(exists($parms->{EQUIVALENCE}{$j}{$i}));
		}
		next if (! defined($e));
		#now go through the equivalence types and add any that are missing
		for $k (values %{ $parms->{VDW}{$e}{$e} }) {
			$eType = $k->{TYPE}; #equivalence vdwType
			#search in current atomtype for eType
			undef $cType;
			for $l (keys %{ $parms->{VDW}{$i}{$i}}) {
				if (! $l) {
					delete $parms->{VDW}{$i}{$i}{$l};
					next;
				}
				$cType = $parms->{VDW}{$i}{$i}{$l} if($parms->{VDW}{$i}{$i}{$l}{TYPE} eq $eType);
			}
			next if (defined($cType) and exists($cType->{VALS}) and scalar(@{ $cType->{VALS} }) > 0);
			if(defined($cType)) {
				#copy the equivalence vals
				@{ $cType->{VALS} } = @{ $k->{VALS} };
			} else {
				#no entry in atomtype for eType so copy entire hash
				@{ $tmp } = sort numerically keys %{ $parms->{VDW}{$i}{$i} };
				$l = 0;
				$l = $tmp->[$#{ $tmp }] + 1 if ($#{ $tmp } > -1);
				$l++;
				%{ $parms->{VDW}{$i}{$i}{$l} } = %{ $k };
				$cType = \%{ $parms->{VDW}{$i}{$i}{$l} };
				$cType->{ATOM} = $cType->{KEY} = "$i $i ";
			}
		}
	}
}

sub getShakeOpts {
	my ($atoms, $bonds, $angles, $parms) = @_;
	my ($i, $shakeOpts, @tmp);

	for $i (1 .. $bonds->{counter}) {
		if (! exists($parms->{PARMS}{POLARIZATION})) {
			if ($atoms->{ $bonds->{LIST}{$i}{ATOMS}[0] }{ELEMENT}{SYMBOL} eq "H") {
				$shakeOpts->{MASS} = $atoms->{ $bonds->{LIST}{$i}{ATOMS}[0] }{ELEMENT}{MASS};
				last;
			} elsif ($atoms->{ $bonds->{LIST}{$i}{ATOMS}[1] }{ELEMENT}{SYMBOL} eq "H") {
				$shakeOpts->{MASS} = $atoms->{ $bonds->{LIST}{$i}{ATOMS}[1] }{ELEMENT}{MASS};
				last;
			}
		} else {
			if (($atoms->{ $bonds->{LIST}{$i}{ATOMS}[0] }{ELEMENT}{SYMBOL} eq "H" or 
			    $atoms->{ $bonds->{LIST}{$i}{ATOMS}[1] }{ELEMENT}{SYMBOL} eq "H") and
				$atoms->{ $bonds->{LIST}{$i}{ATOMS}[0] }{FFTYPE} !~ /_shell/i and 
				$atoms->{ $bonds->{LIST}{$i}{ATOMS}[1] }{FFTYPE} !~ /_shell/i) {
				$shakeOpts->{BOND}{ $bonds->{LIST}{$i}{DATA}{INDEX} } = 1;
			}
		}
	}
	if(exists($shakeOpts->{MASS})) {
		$parms->{PARMS}{SHAKE_MASS} = $shakeOpts->{MASS};
	}elsif(exists($shakeOpts->{BOND})) {
		@tmp = sort numerically keys %{ $shakeOpts->{BOND} };
		$parms->{PARMS}{SHAKE_BOND} = "@tmp";
	}
	for $i (1 ... $angles->{counter}) {
		if ($#{ $ANGLES->{LIST}{$i}{ATOMS} } == 2) {
			if ($atoms->{ $angles->{LIST}{$i}{ATOMS}[0] }{ELEMENT}{SYMBOL} eq "H" and
				$atoms->{ $angles->{LIST}{$i}{ATOMS}[1] }{ELEMENT}{SYMBOL} eq "O" and
				$atoms->{ $angles->{LIST}{$i}{ATOMS}[2] }{ELEMENT}{SYMBOL} eq "H") {
				$parms->{PARMS}{SHAKE_ANGLE} = " a $angles->{LIST}{$i}{DATA}{INDEX}";
				last;
			}
		}
	}	

}

sub fixMorseCoul {
	my ($parms) = $_[0];

	if(exists($parms->{VDW}{TYPE}{"morse/opt"}) and $parms->{QEq}{Opt} == 0 and $parms->{PARMS}{HAS_CHARGE} == 1) {
		if ($parms->{PARMS}{OPTIONS}{PERIODICITY} > 0) {
			&addAllPair($parms,"coul/long", $parms->{VDW}{TYPE}{"morse/opt"}{OPTS});
			$parms->{PARMS}{allCoul} = 1;
		} else {
			&addAllPair($parms,"coul/cut", $parms->{VDW}{TYPE}{"morse/opt"}{OPTS});
			$parms->{PARMS}{allCoul} = 1;
		}
	}
}

sub getVDWlist {
	my ($parms, $parmHybrid, $datWrite) = @_;
	my ($i, $j, $k, $l, $m, @tmp, $hybridOpt, $vdwTypes);
	my ($index1, $index2, $atmN, $offDiag, $diag);

	for $i (keys %{ $parms->{VDW}{TYPE} }) {
		for $j (1 .. $parms->{ATOMTYPES}{counter}) {
			$diag->{lc $i}[$j-1] = "NULL";
		}
	}

	#now get all the types
	map { $vdwTypes->{$_} = 1 } grep {!/TYPE|counter/i} keys %{ $parms->{ATOMTYPES} };

	#loop over all types and get the offdiag components
	for $i (keys %{ $vdwTypes }) {
		next if (! defined($parms->{ATOMTYPES}{$i}) or ! exists($parms->{ATOMTYPES}{$i}{INDEX}));
		for $k (keys %{ $vdwTypes }) { 
			next if (! defined($parms->{ATOMTYPES}{$k}) or ! exists($parms->{ATOMTYPES}{$k}{INDEX}));
			@tmp = sort numerically keys  %{ $parms->{VDW}{$i}{$k} };
			if (!exists($parms->{ATOMTYPES}{$i}) and ! exists($parms->{ATOMTYPES}{$k})) {
				$index1 = $index2 = "*";
			} elsif (! exists($parms->{ATOMTYPES}{$i})) {
				$index1 = "*";
				$index2 = $parms->{ATOMTYPES}{$k}{INDEX};
			} elsif (! exists($parms->{ATOMTYPES}{$k})) {
				$index1 = "*";
				$index2 = $parms->{ATOMTYPES}{$i}{INDEX};
			} elsif ($parms->{ATOMTYPES}{$i}{INDEX} < $parms->{ATOMTYPES}{$k}{INDEX}) {
				$index1 = $parms->{ATOMTYPES}{$i}{INDEX};
				$index2 = $parms->{ATOMTYPES}{$k}{INDEX};
			} else {
				$index2 = $parms->{ATOMTYPES}{$i}{INDEX};
				$index1 = $parms->{ATOMTYPES}{$k}{INDEX};
			}
			for $l (@tmp) {
				$m = $parms->{VDW}{$i}{$k}{$l};
				if(($parms->{PARMS}{HYBRID_VDW} > 0 and exists($m->{IGNORE}) and $m->{IGNORE}==2) or 
				   (exists($m->{IGNORE}) and $m->{IGNORE}==3 and $#tmp > 0) or 
					(! defined($m->{VALS}) or $#{ $m->{VALS} } == -1)) {
					$m->{Lammps}{USED}--;
					delete $parms->{VDW}{$i}{$k}{$l};
					next;
				} 
				next if (exists($m->{IGNORE}) and $m->{IGNORE}==1);
				$atmN = $parms->{VDW}{$i}{$k}{$l}{ATOM};
				$atmN =~ s/\s.*//;
				#eam and eam/fs fix so that the listed fftype is the element
				$atmN = $parms->{ATOMTYPES}{$i}{ATOM} if ($m->{Lammps}{LMPNAME} =~ /eam/); 
				$diag->{lc $m->{Lammps}{LMPNAME} }[$index1-1] = $atmN
					if($index1 ne '*' and $index1==$index2);
				next if ($m->{TYPE} =~ /SW|REXPON|REAX|EAM/i and $m->{TYPE} !~ /REXPON_/i);
				if (defined($parmHybrid) and $parmHybrid ne "") {
					$hybridOpt = $m->{Lammps}{LMPONAME};
					$offDiag->{$hybridOpt}{$index1}{$index2} .= 
						sprintf("%-15s %-4s %-4s $m->{Lammps}{LMPNAME} ","pair_coeff",$index1,$index2);
				} elsif ($i ne $k) {
					$hybridOpt = "a";
					$hybridOpt = $m->{Lammps}{LMPONAME};
					$offDiag->{$hybridOpt}{$index1}{$index2} .= sprintf("%-15s %-4s %-4s ","pair_coeff",$index1,$index2);
				} else {
					next;
				}
				for $j (@{ $m->{VALS} }) {
					last if ($m->{Lammps}{LMPNAME} =~ /^coul\/(long|cut)/);
					if (IsDecimal($j)) {
						$offDiag->{$hybridOpt}{$index1}{$index2} .= sprintf("%25.5f ", $j);
					} elsif($j =~ /[a-zA-Z]/) {
						$offDiag->{$hybridOpt}{$index1}{$index2} .= sprintf("%25s ", $j);
					} else {
						$offDiag->{$hybridOpt}{$index1}{$index2} .= sprintf("%25d ", $j);
					}
				}
				$offDiag->{$hybridOpt}{$index1}{$index2} .= "# $m->{KEY}" if ($datWrite->{Labels});
				$offDiag->{$hybridOpt}{$index1}{$index2} .= "\n";
				delete $parms->{VDW}{$i}{$k}{$l} if ($i eq $k and $parmHybrid); #delete the diagonal element since it will be save in the off diagonal
			}
		}
	}
	return ($diag, $offDiag);
}

sub checkIJPairs {
	my ($parms, $diag) = @_;
	my ($allIJpairs, $i, $j, $k, $atmN, $usedType);

	$allIJpairs = 1;
	for $i (grep {!/TYPE|counter|-1/i} keys %{ $parms->{ATOMTYPES} }) {
		for $j (keys %{ $parms->{VDW}{$i} }) {
			for $k (values %{ $parms->{VDW}{$i}{$j} }) {
				$allIJpairs = 0 if (($i eq $j) and 
									scalar(keys %{ $parms->{VDW}{$i}{$j} }) == 1 and
									(exists($k->{IGNORE}) and $k->{IGNORE} == 1) or
									(!exists($k->{VALS}) or ! @{ $k->{VALS}}));
				if (! exists($k->{IGNORE}) or ! $k->{IGNORE}) {
					$k->{Lammps} = \%{ $parms->{VDW}{TYPE}{ $k->{Lammps}{LMPONAME} }};					
					$k->{Lammps}{USED}++;
					$usedType->{$k->{Lammps}{LMPONAME}} = $k->{Lammps};
					if($i eq $j and exists($parms->{ATOMTYPES}{$i}{INDEX})) {
						$atmN = $k->{ATOM};
						$atmN =~ s/\s.*//;
						#eam and eam/fs fix so that the listed fftype is the element
						$atmN = $parms->{ATOMTYPES}{$i}{ATOM} if ($k->{Lammps}{LMPNAME} =~ /eam/); 
						$diag->{lc $k->{Lammps}{LMPNAME} }[ $parms->{ATOMTYPES}{$i}{INDEX}-1 ] = $atmN;
					}
				}
			}
		}
	}
	#remove unused types
	for $i (keys %{ $parms->{VDW}{TYPE}}) {
		delete $parms->{VDW}{TYPE}{$i} if (! exists($usedType->{$i}) and ($i ne "coul/long" and $i ne "coul/cut"));
	}
	return ($allIJpairs, $usedType);
}

sub getReaxRexPoNpairLine {
	my ($parms, $parmHybrid) = @_;
	my (@tmp, $i, $j);

	@tmp = ();
	if (1>2) {
	for $i (keys %{ $parms->{ATOMTYPES} }) {
		next if ($i eq "counter");
		$j = $parms->{ATOMTYPES}{$i}{INDEX};
			$tmp[$j - 1] = $parms->{ATOMTYPES}{$i}{REAX};
		}
		$j = "";
		for $i (@tmp) {
			$j .= "$i ";
		}
		$parms->{PARMS}{OFF_DIAG} .= sprintf("%-15s %-4s %-4s%s $j\n\n", "pair_coeff", "*", "*", $reaxFF);
	} elsif ($ffType =~ /8/) {
		@tmp = ();
		for $i (keys %{ $parms->{ATOMTYPES} }) {
			next if ($i eq "counter");
			$j = $parms->{ATOMTYPES}{$i}{INDEX};
			$tmp[$j - 1] = $parms->{ATOMTYPES}{$i}{REAX};
		}
		$j = "";
		for $i (@tmp) {
			$j .= "$i ";
		}
		$parms->{PARMS}{OFF_DIAG} .= sprintf("%-15s %-4s %-4s%s $j\n\n", "pair_coeff", "*", "*", $reaxFF);
	}

}

sub checkAtomTypes {
	my ($atoms, $parms, $ffTypeID, $wOpts) = @_;
	my ($i, $atomTypes, $cffType, $sameParms, $atmName, $PLIST);
	my ($eleList, $cTmp, $typeIDlist, $typeIDcounter, $count);
	
	$atomTypes = $parms->{ATOMTYPES};
	$eleList = &LoadElements;
	&AddElementField($atoms, $eleList, $parms);

	$typeIDcounter = 0;

	#first we collapse all the parameters into single key->value hash lists
	for $i (qw/VDW BONDS ANGLES TORSIONS INVERSIONS/) {
		$PLIST->{$i} = ();
		$count =  0;
		&getValParent($parms->{$i}, \%{ $PLIST->{$i}{LIST} }, \$count);
		$PLIST->{$i}{count} = $count;
	}

	for $i (sort numerically keys %{ $atoms }) {
		$totAtms++;
		$cffType = $atoms->{$i}{FFTYPE};
		if(exists($parms->{OVERWRITE_FFTYPES}) and exists($parms->{OVERWRITE_FFTYPES}{$cffType})) {
			$cffType = $parms->{OVERWRITE_FFTYPES}{$cffType};
			$atoms->{$i}{FFTYPE} = $cffType;
		}
		$atmName = $atoms->{$i}{ATMNAME};
		if($wOpts->{compress_parms} == 1 and $atomTypes->{$cffType}{USED} == 0) {
			if(! exists($sameParms->{$cffType})) {
				$cTmp = checkSameParm($parms, $cffType, $sameParms, $PLIST);
				$cTmp = $sameParms->{$cTmp} if (exists($sameParms->{$cTmp}));
				if ($cTmp ne $cffType) {
					$sameParms->{$cffType} = $cTmp;
					$parms->{EQUIVALENCE}{$cTmp}{$cffType} = 1;
					$cffType = $atoms->{$i}{FFTYPE} = $cTmp;
				}
			} else {
				$cffType = $atoms->{$i}{FFTYPE} = $sameParms->{$cffType};
			}
		} 
		if (exists($parms->{EQUIVALENCE}) and ! exists($atomTypes->{$cffType})) {
			&checkAtomTypeEquiv($parms, $cffType) 
		}
		if (exists($atomTypes->{$cffType})) {
			$atomTypes->{$cffType}{USED} = 1;
			$atomTypes->{$cffType}{ELENAME} = $atoms->{$i}{ELEMENT}{SYMBOL};
			$atomTypes->{$cffType}{ELENUM} = $atoms->{$i}{ELEMENT}{NUM};
			$atoms->{$i}{PARMS} = \%{ $atomTypes->{$cffType} };
			$atoms->{$i}{CHARGE} = $atomTypes->{$cffType}{CHARGE} if ($atomTypes->{$cffType}{USE_CHARGE});
			if(!exists($typeIDlist->{$cffType})) {
				$typeIDlist->{$cffType} = ++$typeIDcounter;
				$atomTypes->{$cffType}{TYPEID} = $typeIDcounter;
			}
		}elsif ($ffTypeID =~ /5|8/) {
			($cTmp, $atoms->{$i}{ELENUM}) = getElement($eleList, $cffType, $atmName);
			if (! defined($cTmp) or ! exists($atomTypes->{$cTmp})) {
				$cTmp = uc $cffType;
				$cTmp = substr($cTmp,0,1);
				if (! defined($cTmp) or ! exists($atomTypes->{$cTmp})) {
					$ERRORS{ELEMENT}{$atoms->{$i}{FFTYPE}}++;
					next;
				}
			}
			%{ $atomTypes->{$cffType} } = %{ $atomTypes->{$cTmp} };
			$atomTypes->{$cffType}{USED} = 1;
			$atomTypes->{$cffType}{REAX} = $cffType;
			$atoms->{$i}{ELENAME} = $cTmp;
			$atoms->{$i}{FFTYPE} = $cffType;
			$atoms->{$i}{PARMS} = \%{ $atomTypes->{$cffType} };
		} else {
			$ERRORS{ATOMTYPES}{$cffType}++;
		}
	}
	if($wOpts->{sort_fftypes} == 1) {
		$typeIDcounter = 0;
		for $i (sort alphabetically keys %{ $typeIDlist }) {
			$typeIDlist->{$i} = ++$typeIDcounter;
			$atomTypes->{$i}{TYPEID} = $typeIDcounter;
		}
	}

	#fix missing VDW that have entry in equivalence array
	&fixMissingEquivVDW($parms);
}

#checks for the same parameters as the current one to reduce redudancies
sub checkSameParm {
	my ($parms, $cffType, $sameParms, $PLIST) = @_;
	my ($i, $j, $k, $l1, $tlist, $ct, $pt, $sP, $foundSP);

	$sP = $cffType;

	#now iterate through all other fftypes and determine if same as current
	@{ $tlist }= grep {!/$cffType/i} keys %{ $parms->{ATOMTYPES} }; #get list of all other fftypes
	$ct = $parms->{ATOMTYPES}{$cffType};
	for $i (@{ $tlist }) { #loop through them...
		$pt = $parms->{ATOMTYPES}{$i};
		next if($pt->{CHARGE} != $ct->{CHARGE} or $pt->{MASS} != $ct->{MASS});
		$foundSP = 1; #set the found same parm flag
		TLOOP: for $j (qw/VDW BONDS ANGLES TORSIONS INVERSIONS/) { #loop over all the types of parms
			$l1 = findTypeParms($PLIST->{$j}{LIST}, $cffType); #find all the parameters with keys corresponding to the current type
			for $k (@{ $l1 }) {
				if(! matchTypeParm($PLIST->{$j}{LIST}, $k, $cffType, $i, $sameParms)) {
					$foundSP = 1;
					last TLOOP;
				}
			}
		}
		if ($foundSP) {
			$sP = $i;
			last;
		}
	}

	return $sP;
}

sub matchTypeParm {
	my ($plist, $cParm, $cffType, $sffType, $sameParms) = @_;
	my ($i, $sKey, $found);

	$found = 1;
	$sKey = $cParm->{KEY};
	$sKey =~ s/$cffType/$sffType/g;
	for $i (keys %{ $plist }) {
		next if ($plist->{$i}{1}{KEY} !~ /^$sKey/);
		if(@{ $cParm->{VALS} } != @{ $plist->{$i}{1}{VALS} }) {
			$found = 0;
			last;
		}
	}
	return $found;
}

sub findTypeParms {
	my ($plist, $cffType) = @_;
	my ($i, $pL, $tmp);

	for $i (keys %{ $plist }) {
		%{ $tmp } = map { $_ => undef } split /\s+/, $plist->{$i}{1}{KEY};
		push @{ $pL }, $plist->{$i}{1} if (exists($tmp->{$cffType}));
	}
	return $pL;
}
sub checkAtomTypeEquiv {
	my ($parms, $cffType) = @_;
	my ($atomTypes, $vdw, $equiv, $j, $cTmp, $k, $tmp);

	$atomTypes = $parms->{ATOMTYPES};
	$vdw = $parms->{VDW};
	for $equiv (keys %{ $parms->{EQUIVALENCE}{$cffType} }) {
		next if (! exists($atomTypes->{$equiv}));
		#check for and save equivalence atomtype data if not present
		if (! exists($atomTypes->{$cffType}{MASS})) {
			$atomTypes->{$cffType} = dclone($atomTypes->{$equiv});
			$atomTypes->{$cffType}{LABEL} = $cffType;
		}
		#now we check for the vdw interactions. Here we simply add the equivalence data 
		#to the end of the hash and allow findDuplicateParms to deal with any duplicates
		for $j (keys %{ $vdw->{$equiv} }) {
			$cTmp = 0;
			if(exists($vdw->{$cffType}) and exists($vdw->{$cffType}{$j})) { 
				@{ $tmp } = sort numerically keys %{ $vdw->{$cffType}{$j} };
				$cTmp = pop @{ $tmp };
			}
			$cTmp++;
			for $k (keys %{ $vdw->{$equiv}{$j} }) {
				$vdw->{$cffType}{$j}{$k+$cTmp} = dclone($vdw->{$equiv}{$j}{$k});
				$vdw->{$cffType}{$j}{$k+$cTmp}{KEY} =~ s/$equiv /$cffType /;
			}		
		}
		for $j (keys %{ $vdw } ) {
			next if($j eq $equiv or !exists($vdw->{$j}{$equiv}));
			$cTmp = 0;
			if(exists($vdw->{$j}{$cffType})) {
				@{ $tmp } = sort numerically keys %{ $vdw->{$j}{$cffType} };
				$cTmp = pop @{ $tmp };
			}
			$cTmp++;
			for $k (keys %{ $vdw->{$j}{$equiv} }) {
				$vdw->{$j}{$cffType}{$k+$cTmp} = dclone($vdw->{$j}{$equiv}{$k});
				$vdw->{$j}{$cffType}{$k+$cTmp}{KEY} =~ s/$equiv /$cffType /;
			}	
		}
	}	
	#undef $equiv;
	#$equiv = $parms->{EQUIVALENCE}{$cffType}
	#	if(exists($parms->{EQUIVALENCE}) and exists($parms->{EQUIVALENCE}{$cffType}) and exists($atomTypes->{ $parms->{EQUIVALENCE}{$cffType} }));
	#if(!exists($atomTypes->{$cffType}) and defined($equiv)) {
	#	$atomTypes->{$cffType} = dclone($atomTypes->{$equiv});
	#	$atomTypes->{$cffType}{LABEL} = $cffType;
	#	$parms->{VDW}{$cffType}{$cffType} = dclone($parms->{VDW}{$equiv}{$equiv}) 
	#		if (! exists($parms->{VDW}{$cffType}{$cffType}) and exists($parms->{VDW}{$equiv}) and exists($parms->{VDW}{$equiv}{$equiv}));
	#	for $i (values %{ $parms->{VDW}{$cffType}{$cffType} }) {
	#		$i->{KEY} = "$cffType $cffType ";
	#	}
	#	for $i (keys %{ $parms->{VDW}{$equiv} }) {
	#		next if(exists($parms->{VDW}{$cffType}{$i}));
	#		$parms->{VDW}{$cffType}{$i} = dclone($parms->{VDW}{$equiv}{$i});
	#		for $cTmp (values %{ $parms->{VDW}{$cffType}{$i} }) {
	#			$cTmp->{KEY} =~ s/$equiv /$cffType /;
	#		}
	#	}
	#	for $i (keys %{ $parms->{ATOMTYPES} }) {
	#		next if (! exists($parms->{VDW}{$i}{$equiv}));
	#		next if ($i eq $equiv);
	#		$parms->{VDW}{$i}{$cffType} = dclone($parms->{VDW}{$i}{$equiv});
	#		for $cTmp (values %{ $parms->{VDW}{$i}{$cffType} }) {
	#			$cTmp->{KEY} =~ s/$equiv /$cffType /;
	#		}
	#	}
	#}
}

sub getElement {
	my ($ELEMENTS, $fftype, $atmName) = @_;
	my ($eleNum);
	
	$fftype =~ s/_.*//;
	$fftype =~ s/\d+.*//;
	$eleNum = FindElement($ELEMENTS, $fftype);
	$eleNum = FindElement($ELEMENTS, $atmName) if (! defined($eleNum));
	return ($ELEMENTS->{$eleNum}{SYMBOL},$ELEMENTS->{$eleNum}) if (defined($eleNum));
}

sub getElement_old {
	my ($eleList, $ffType, $atmName) = @_;
	my ($i, $element);

	for $i ($ffType, $atmName) {
		$element = $i;
		next if ($element !~ /^([A-Za-z])([a-z]?)/);
		$element = $1;
		$element .= $2 if ($2);
		return ($element, $eleList->{uc($element)}) if (exists($eleList->{uc($element)}));
		$element = $1;
		return ($element, $eleList->{uc($element)}) if (exists($eleList->{uc($element)}));
	}

}

sub getHybridOpts {
	my ($vdws) = $_[0];
	my ($parmHybrid, @tmp, $i, $k, $l, $pairs);

	return "" if (scalar(keys %{ $vdws->{TYPE} } ) == 1);
	#hybrid so determine if overlay
	$parmHybrid = "hybrid ";
	@tmp = grep {!/TYPE|counter/i} keys %{ $vdws };
	iLoop: for $i (@tmp) {
		for $k (keys %{ $vdws->{$i} }) {
			next if ($i gt $k);
			for $l (keys %{ $vdws->{$i}{$k} }) {
				$pairs->{$i}{$k} = 0 if (! exists($pairs->{$i}) or ! exists($pairs->{$i}{$k}));
				$pairs->{$i}{$k}++ if ((!exists($vdws->{$i}{$k}{$l}{IGNORE}) or $vdws->{$i}{$k}{$l}{IGNORE} != 1) or
				                       (!exists($vdws->{$i}{$l}{$k}{IGNORE}) or $vdws->{$i}{$l}{$k}{IGNORE} != 1));
				if ($pairs->{$i}{$k} > 1) {
					$parmHybrid = "hybrid/overlay ";
					last iLoop;
				}
			}
		}
	}
	return $parmHybrid;
}

sub doPairMix {
	my ($parms) = $_[0];
	my ($i, $j, $k, @tmp, $iName, $jType, $kType, $count, $ignore);
	
	@tmp = grep {!/^(TYPE|counter|\*)/i} keys %{ $parms->{ATOMTYPES} };
	for $i (keys %{ $parms->{VDW}{TYPE} }) {
		next if (! $parms->{VDW}{TYPE}{$i}{PAIRMIX});
		$iName = $parms->{VDW}{TYPE}{$i}{NAME};
		for $j (@tmp) {
			$jType = findVDWEntry($parms->{VDW}{$j}{$j}, $iName);
			next if (! $jType or !exists($parms->{ATOMTYPES}{$j}) or 
				!exists($parms->{VDW}{$j}{$j}{$jType}{VALS}) or scalar(@{ $parms->{VDW}{$j}{$j}{$jType}{VALS} }) == 0);
			for $k (@tmp) {
				next if ($k eq $j or !exists($parms->{ATOMTYPES}{$k}) or !exists($parms->{ATOMTYPES}{$k}{INDEX}) 
						or !exists($parms->{ATOMTYPES}{$j}) or !exists($parms->{ATOMTYPES}{$j}{INDEX}) 
						or $parms->{ATOMTYPES}{$k}{INDEX} > $parms->{ATOMTYPES}{$j}{INDEX});
				$kType = findVDWEntry($parms->{VDW}{$k}{$k}, $iName);
				next if (! $kType or !exists($parms->{VDW}{$k}{$k}{$kType}{VALS}) or scalar(@{ $parms->{VDW}{$k}{$k}{$kType}{VALS} }) == 0);
				next if ((findVDWEntry($parms->{VDW}{$j}{$k}, $iName) and exists($parms->{VDW}{$j}{$k}{$jType}{VALS}) and scalar(@{ $parms->{VDW}{$j}{$k}{$jType}{VALS} }) > 0) or 
						 (findVDWEntry($parms->{VDW}{$k}{$j}, $iName) and exists($parms->{VDW}{$k}{$j}{$jType}{VALS}) and scalar(@{ $parms->{VDW}{$k}{$j}{$jType}{VALS} }) > 0));
				$count = scalar(keys %{ $parms->{VDW}{$j}{$k} }) + 1;
				%{ $parms->{VDW}{$j}{$k}{$count} } = %{ $parms->{VDW}{$j}{$j}{$jType} };
				$parms->{VDW}{$j}{$k}{$count}{ATOM} = $parms->{VDW}{$j}{$k}{$count}{KEY} = "${j} ${k} ";
				$ignore = 0;
				$ignore = 1 if ($parms->{VDW}{$j}{$j}{$jType} == 1 and $parms->{VDW}{$k}{$k}{$kType}{IGNORE} == 1);
				$ignore = 3 if ($parms->{VDW}{$j}{$j}{$jType} > 1 and $parms->{VDW}{$k}{$k}{$kType}{IGNORE} > 1);
				$parms->{VDW}{$j}{$k}{$count}{IGNORE} = $ignore;
				$parms->{VDW}{$j}{$k}{$count}{VALS} = getPairMix($parms->{VDW}{$j}{$j}{$jType}{VALS},$parms->{VDW}{$k}{$k}{$kType}{VALS}, $iName);
			}
		}
	}

}

sub getPairMix {
	my ($iVals, $jVals, $pairType) = @_;
	my (@VALS, $mixType);

	$mixType = $PARMS->{PARMS}{mix_rule};

	$VALS[0] = sqrt($iVals->[0]*$jVals->[0]);
	if ($mixType eq "geometric") {
		$VALS[1] = sqrt($iVals->[1]*$jVals->[1]);
	} else {
		$VALS[1] = 0.5*($iVals->[1]+$jVals->[1]);

	}
  
	if ($pairType ne "LJ_12_10" and $#{ $iVals } == 2) {
		$VALS[2] = 0.5*($iVals->[2]+$jVals->[2]);
	}
	
	if ($pairType eq "LJ_M_N") {
		$VALS[2] = $iVals->[2];
		$VALS[3] = $iVals->[3];
	}
	return \@VALS;
}

sub findVDWEntry {
	my ($vdw, $lammpsName) = @_;
	my ($i, $returnVal);

	$returnVal = 0;
	for $i (keys %{ $vdw }) {
		if ($vdw->{$i}{TYPE} eq $lammpsName and (scalar(@{ $vdw->{$i}{VALS} }) > 0 or $vdw->{$i}{IGNORE})) {
			$returnVal = $i;
		} elsif ($vdw->{$i}{TYPE} eq $lammpsName) {
			delete $vdw->{$i};
		}
	}

	return $returnVal;
}

sub getImage {
	my ($atomPos, $boxDim, $boxLen, $isHi) = @_;
	my ($imageVal) = 0;

	if ($isHi) {
		while ($boxDim < $atomPos) {
			$imageVal++;
			$atomPos -= $boxLen;
		}
	} else {
		while ($atomPos < $boxDim) {
			$imageVal--;
			$atomPos += $boxLen;
		}
	}
	return $imageVal;
}

sub findDuplicateParms {
	my ($parms, $pL) = @_;
	my ($PLIST, $parmList, $i, $searchStr, $count, $j, $k, $SLIST, @tmp);

	$pL = "VDW BONDS ANGLES TORSIONS INVERSIONS" if (! defined($pL));
	@{ $parmList } = split /\s+/,$pL;
	for $i (@{ $parmList }) {
		$SLIST = ();
		$count = 0;
		$PLIST = ();
		getValParent($parms->{$i}, \%{ $PLIST }, \$count);
		for $j (keys %{ $PLIST }) {
			for $k (keys %{ $PLIST->{$j} }) {
				if (! $k) {
					delete $PLIST->{$j}{$k};
					next;
				}
				next if ($PLIST->{$j}{$k}{KEY} eq "*" or $PLIST->{$j}{$k}{TYPE} eq "SW" or (exists($PLIST->{$j}{$k}{IGNORE}) and $PLIST->{$j}{$k}{IGNORE} ));
				@tmp = split /\s+/,$PLIST->{$j}{$k}{KEY};
				$searchStr = "@tmp $PLIST->{$j}{$k}{TYPE} @{$PLIST->{$j}{$k}{VALS}}";
				if (exists($SLIST->{$searchStr})) {
					delete $PLIST->{$j}{$k};
				} else {
					@tmp = reverse @tmp;
					$searchStr = "@tmp $PLIST->{$j}{$k}{TYPE} @{$PLIST->{$j}{$k}{VALS}}";
					if (exists($SLIST->{$searchStr})) {
						delete $PLIST->{$j}{$k};
					} else {
						$SLIST->{$searchStr} = 1;
					}
				}
			}
		}
	}
}

sub createLammpsClusterScript {
	my ($ffs) = $_[0];
	my ($ffList, $scriptFile, $scriptCmd, $reaxff);
   
	if ($#{ $ffs } == 0) {
		$ffList = $ffs->[0]{FF};
		$reaxff = $ffs->[0]{FF} if ($ffs->[0]{FFTYPE} =~ /REAX|REXPON/i);
	} else {
		$ffList = '"';
		for (@{ $ffs }) {
			$ffList .= "$_->{FF} ";
		}
		chomp $ffList;
		$ffList .= '"';
	}

	$scriptFile = "$Bin/createClusterScript.pl";
	if (! -e $scriptFile) {
		print "Failed!\n";
		return;
	}
	$scriptFile .= " -r $reaxFF" if (defined($reaxff));
	$scriptCmd = "$scriptFile -p $suffix -b \"$inputFile\" -i in.${suffix} -d data.${suffix} -s ${suffix}.lammps.slurm -f \"$ffList\"> _createscript";
	system($scriptCmd);
	system("rm -fr _createscript");
}

sub getTIP4Popts {
	my ($parms, $atoms, $bonds) = @_;
	my ($tip4Str, $oType, $hType, $bType, $aType);

	# first get the type of the ow/hw atoms
	for $i (keys %{ $parms->{ATOMTYPES} }) {
		if ($i =~ /OW/) {
			$oType = $parms->{ATOMTYPES}{$i}{INDEX};
		} elsif ($i =~ /HW/i) {
			$hType = $parms->{ATOMTYPES}{$i}{INDEX};
		}
	}

	#get ow - hw bond type
	for $i (1 ... $bonds->{counter}) {
		if ($bonds->{LIST}{$i}{DATA}{KEY} =~ /(H|D|T)W OW/ or $bonds->{LIST}{$i}{DATA}{KEY} =~ /OW (H|T|D)W/) {
			$bType = $bonds->{LIST}{$i}{DATA}{INDEX};
			last;
		}
	}
	die "ERROR: OW - HW bond type not found for TIP4P water!\n" if (! defined($bType));
	#angle type is already recorded for shake
	$aType = $parms->{PARMS}{SHAKE_ANGLE};
	$aType =~ /(\d+)/;
	$aType = $1;

	#$tip4Str = "lj/cut/tip4p/long/opt $oType $hType $bType $aType $parms->{PARMS}{tip4_om_dist} 10";
	$tip4Str = "$oType $hType $bType $aType $parms->{PARMS}{tip4_om_dist} 10";
	$parms->{VDW}{TYPE}{"lj/cut/tip4p/long/opt"}{OPTS} = $tip4Str;
	if(exists($parms->{VDW}{TYPE}{"lj/cut/tip4p/long/soft"})) {
		$parms->{VDW}{TYPE}{"lj/cut/tip4p/long/soft"}{OPTS} = "$oType $hType $bType $aType $parms->{PARMS}{tip4_om_dist} 1 0.5 10 10";
	}
	#delete $parms->{VDW}{TYPE}{"lj/cut/tip4p/long/opt"};
	#$parms->{VDW}{TYPE}{$tip4Str} = "";
}

sub getInputTypes {
	my ($isStr) = $_[0];

	my ($searchDir) = "$Bin/dat/LAMMPS/";
	my (@list) = `ls ${searchDir}/in.lammps.*`;
	my ($i, $typeStr);

	for $i (@list) {
		chop $i;
		if ($i =~ /in.lammps.(\S+)/) {
			if ($isStr) {
				$typeStr .= "$1 ";
			} else {
				$typeStr->{$1} = "$i";
			}
		}
	}
	return $typeStr;
}

sub setOpts {
	my ($bgf, $parms, $opts, $qeqOpts, $mols) = @_;
	my ($i, $j, $hasCharge, $tmp, $count, $add);

	$parms->{PARMS}{OPTIONS}{PERIODICITY} = 3;
	$parms->{PARMS}{OPTIONS}{PERIODICITY} = 0 if ($opts =~ /finite/i);
	$parms->{PARMS}{OPTIONS}{PERIODICITY} = 2 if ($opts =~ /2d|slab|dimension 2/i);
	$parms->{PARMS}{OPTIONS}{SHAKE} = 1;
	$parms->{PARMS}{OPTIONS}{SHAKE} = 0 if ($opts =~ /no shake/i);
	$parms->{PARMS}{OPTIONS}{SHAKE_WHAT} = $1 if ($opts =~ /shake (solute|solvent)/i);
	$parms->{PARMS}{OPTIONS}{VDW_EWALD} = 0;
	$parms->{PARMS}{OPTIONS}{HYBRID_OVERLAY} = 1;
	$parms->{PARMS}{OPTIONS}{HYBRID_OVERLAY} = 0 if ($opts =~ /no hybrid overla(p|y)/i);
	$parms->{PARMS}{OPTIONS}{ISOTOPES} = 1 if ($opts =~ /isotope/);
	$parms->{PARMS}{OPTIONS}{MD_TIMESTEP} = 1 if (! exists($parms->{PARMS}{OPTIONS}{MD_TIMESTEP}));
	$parms->{PARMS}{TYPE_CHARGES} = "";
	$nocross = 0;
	$nocross = 1 if ($opts =~ /no cross/i);
	$count = 0;

	delete $parms->{PARMS}{is_tip4p};
	if(exists($parms->{PARMS}{tip4_om_dist})) {
		for $i (keys %{ $mols }) {
			$tmp = ();
			@{ $tmp } = keys %{ $mols->{$i}{MEMBERS} };
			if(isWater($bgf, $tmp->[0])) {
				$parms->{PARMS}{is_tip4p} = 1;
				last;
			}
		}
	}

	$tmp = ();
	&getValParent($parms->{VDW}, \%{ $tmp }, \$count);
	$hasCharge = 0;
	for $i (values %{ $bgf }) {
		if ($i->{CHARGE} > 0 or $i->{CHARGE} < 0) {
			$hasCharge = 1;
			last;
		}
	}
	$parms->{PARMS}{HAS_CHARGE} = $hasCharge;
	if ($opts =~ /vdw ewald|ewald vdw/i) {
		$add = "long long";
		$add = "long off" if (! $hasCharge);
		$parms->{PARMS}{OPTIONS}{VDW_EWALD} = 1;
	}
	for $i (1 .. $count) {
		for $j (values %{ $tmp->{$i} }) {
			next if (!exists($j->{Lammps}));
			if($parms->{PARMS}{OPTIONS}{VDW_EWALD} == 1) {
				$j->{Lammps}{name} =~ s/lj\/charmm\/coul\/long\S*/lj\/coul $add/g;
				$j->{Lammps}{name} =~ s/lj\/cut\/coul\/long\S*/lj\/coul $add/g;
				$j->{Lammps}{name} =~ s/buck\/coul\S*/buck\/coul $add/g;
				$j->{Lammps}{opts} = " $parms->{PARMS}{cut_vdw} $parms->{PARMS}{cut_coul}";
			} elsif ($qeqOpts->{Opt} == 2) {
				$j->{Lammps}{name} = "lj/gromacs" if($j->{Lammps}{name} =~ /charmm/);
				$j->{Lammps}{name} =~ s/\/coul.*//;
			} elsif (! $hasCharge) {
				$j->{Lammps}{name} =~ s/\/coul\S+//;
				$j->{Lammps}{name} = "lj/charmm/coul/charmm" if ($j->{Lammps}{name} eq "lj/charmm");
			} elsif ($parms->{PARMS}{OPTIONS}{PERIODICITY} == 0) {
				$j->{Lammps}{name} =~ s/lj\/charmm\/coul\/long\/opt/lj\/charmm\/coul\/charmm/g;
				$j->{Lammps}{name} =~ s/\/coul\/long/\/coul\/cut/g;
			} elsif ($j->{TYPE} eq "LJ_6_12" and exists($parms->{PARMS}{is_tip4p}) and $parms->{PARMS}{is_tip4p}) { # tip4p fix
				$j->{Lammps}{name} = "lj/cut/tip4p/long/opt";
				$j->{Lammps}{fep}{pair} = "lj/cut/tip4p/long/soft";
				$j->{Lammps}{opts} = "";
			}
		}
	}
	$parms->{PARMS}{FEP}{valid} = 0;
	$parms->{PARMS}{FEP}{var_str} = "variable             lambda equal ramp(1.0,0.0)\n";
	$parms->{PARMS}{FEP}{var_str} .= "variable             dlambda equal -0.05\n";

}

sub fixMesoDNAAtoms {
	my ($atoms, $bonds) = @_;
	my ($i, $j, $bases, $molidmax, $molid);

	$molidmax = 0;
	for $i (keys %{ $atoms }) {
		$molidmax = ${ $atoms->{$i}{MOLECULEID} } 
			if (${ $atoms->{$i}{MOLECULEID} } > $molidmax);
		next if ($atoms->{$i}{NUMBONDS} != 1);
		for $j (@{ $bonds->{$i} }) {
			$bases->{$atoms->{$j}{INDEX}} = 1;
		}
	}

	$PARMS->{MESODNA} = "group	nuc molecule > $molidmax\ngroup	backbone subtract solute nuc\ngroup	movable subtract all nuc\n";
	$molidmax++;
	for $i (sort numerically keys %{ $bases }) {
		$molid = ();
		$molid->{index} = $molidmax;
		$atoms->{$i}{MOLECULEID} = \$molid->{index};
		for $j (@{ $bonds->{$i} }) {
			next if ($atoms->{$j}{NUMBONDS} > 1);
			$atoms->{$j}{MOLECULEID} = \$molid->{index};
		}
		$molidmax++;
	}
}

sub determineTinkerOpts {
	my ($parms, $sfix, $monopole_flag) = @_;
	my ($i, $j, $k, $l, $amoeba_flag, $amoeba_ffs, $ff, $tinker_flag, $a_types, $CNV, $rec, $curr);
	my ($key_file);

	$CNV = LoadConverter();
	@{ $a_types } = grep {!/^counter|TYPE/i} keys %{ $parms->{ATOMTYPES} };
	#determine if amoeba
	$amoeba_flag = 0;
	for $i (@{ $a_types }) {
		if(exists($parms->{ATOMTYPES}{$i}{AMOEBA_FLAG}) and $parms->{ATOMTYPES}{$i}{AMOEBA_FLAG} == 1) {
			$amoeba_flag = 1;
			$parms->{PARMS}{AMOEBA_FLAG} = 1;
			last;
		}
	}
	return if (! $amoeba_flag);
	#now let's ensure that all the atomtypes are amoeba...
	if($amoeba_flag) {
		for $i (@{ $a_types }) {
			if(! exists($parms->{ATOMTYPES}{$i}{AMOEBA_FLAG}) or $parms->{ATOMTYPES}{$i}{AMOEBA_FLAG} != 1) {
				die "ERROR: Some atoms are AMOEBA polarizable atom types but other aren't.\n" . 
				"Cannot (yet) mix polarizable/non-polarizable types forcefields!\n";
			}
			$amoeba_ffs->{ $parms->{ATOMTYPES}{$i}{TINKER_PRM} } = 1;
			$ff = $parms->{ATOMTYPES}{$i}{TINKER_PRM};
		}
	} else {
		#non-polarizable, so determine if we need extra tinker options (for bi- and pi-torsions)
		$tinker_flag = 0;
		for $i (@{ $a_types }) {
			if (!exists($parms->{ATOMTYPES}{$i}{TINKER_ID}) or $parms->{ATOMTYPES}{$i}{TINKER_ID} == 0) {
				$tinker_flag = 1;
			}
			next if (! exists($parms->{VDW}{$i}));
			for $j (values %{ $parms->{VDW}{$i} }) {
				for $k (values %{ $j }) {
					next if ($k->{TYPE} ne "AMOEBA");
					$k->{TYPE} = "LJ_6_12";
					$k->{Lammps} = getLammpsOpts("LJ_6_12","vdw", $CNV, 0, $parms->{PARMS});
				}
			}
		}
		delete $parms->{VDW}{TYPE}{amoeba};
		&setVDWtypes($parms);
		$parms->{PARMS}{ALL_TINKER} = 1 if($tinker_flag and (
			( exists($parms->{BI_TORSIONS}) and keys %{ $parms->{BI_TORSIONS} })
		 or ( exists($parms->{PI_TORSIONS}) and keys %{ $parms->{PI_TORSIONS} })));
	}
	#now make sure that every Tinker angle has associated urey-bradley and bond-angle entries
	for $i (@{ $a_types }) {
		next if (! exists($parms->{ANGLES}{$i}));
		for $j (values %{ $parms->{ANGLES}{$i} }) {
			for $k (values %{ $j }) {
				for $l (values %{ $k }) {
					next if (! exists($l->{TINKER}) or $l->{TINKER} == 0);
					@{ $l->{UREY_BRADLEY}{VALS} } = (0,0) if (! exists($l->{UREY_BRADLEY}));
					@{ $l->{BOND_ANGLE}{VALS} } = (0,0,0,0) if (! exists($l->{BOND_ANGLE}));
				}
			} 
		}
	}

	#now determine if multiple prm TINKER forcefield (prm) files were specified, and if so, consolidate and write new ff
	if($amoeba_flag and ! $monopole_flag) { 
		if(scalar(keys %{ $amoeba_ffs }) > 1) {
			$ff = &writeAmoebaFF($parms, $amoeba_ffs, $sfix);
		}
		$parms->{VDW} = ();
		&addAllPair($parms, "amoeba");
		$parms->{VDW}{"-1"}{"-1"}{0}{VALS} = ($ff);
		$parms->{PARMS}{OFF_DIAG} = "pair_coeff * * $ff";
		$key_file = $ff;
		$key_file =~ s/\.\w+$/\.key/;
		$parms->{PARMS}{OFF_DIAG} .= " $key_file" if (-e $key_file); 
		$parms->{PARMS}{ALL_TINKER} = 1;
	}

	&convertAmoeba($parms) if ($monopole_flag);
}

sub convertAmoeba {
	my ($parms) = $_[0];
	my ($i, $j, $k, $l, $m, $a_types, $CNV);

	$CNV = LoadConverter();
	@{ $a_types } = grep {!/^counter|TYPE/i} keys %{ $parms->{ATOMTYPES} };

	#here we switch from AMOEBA multipole to a monopole description.
	#first remove bi and pi torsions
	delete $parms->{BI_TORSIONS};
	delete $parms->{PI_TORSIONS};
	#now switch from AMOEBA vdw to LJ_6_12
	for $i (@{ $a_types }) {
		next if (! exists($parms->{VDW}{$i}));
		for $j (values %{ $parms->{VDW}{$i} }) {
			for $k (values %{ $j }) {
				next if ($k->{TYPE} ne "AMOEBA");
				$k->{TYPE} = "LJ_M_N";
				$k->{Lammps} = getLammpsOpts("LJ_M_N","vdw", $CNV, 0, $parms->{PARMS});
				#remove extra
				for $l (2 .. $#{ $k->{VALS} }) {
					pop @{ $k->{VALS} };
				}
				push @{ $k->{VALS} }, (14,7);
			}
		}
	}
	delete $parms->{VDW}{TYPE}{amoeba};
	&setVDWtypes($parms);

	#finally switch from improper AMOEBA to improper harmonic
	for $i (@{ $a_types }) {
		next if (! exists($parms->{INVERSIONS}{$i}));
		for $j (values %{ $parms->{INVERSIONS}{$i} }) {
			for $k (values %{ $j }) {
				for $l (values %{ $k }) {
					for $m (values %{ $l }) {
						next if ($m->{TYPE} ne "AMOEBA");
						$m->{TYPE} = "HARMONIC";
						$m->{Lammps} = getLammpsOpts("HARMONIC","inversion", $CNV, 0, $parms->{PARMS});
						push @{ $m->{VALS} }, 0;
					}
				}
			}
		}
	}
	delete $parms->{INVERSIONS}{TYPE}{amoeba};
	$parms->{INVERSIONS}{TYPE}{harmonic} = 1;
}

sub writeTinkerTypes {
	my ($atoms, $parms, $datFile) = @_;
	my ($i);

	print $datFile "\nTinker Types\n\n";
	for $i (sort numerically keys %{ $atoms } ) {
		printf $datFile "%5d %8d\n", $i, $atoms->{$i}{PARMS}{TINKER_ID};
	}
}	

sub writeTinkerBiTorsionFile {
	my ($parms, $suffix, $header_flag, $BITORFILE) = @_;
	my ($i, $j, $k, $l, $m, $n, $c, $sfile);

	if(! defined($header_flag) or $header_flag > 0) {
		$sfile = "${suffix}.bitorsion.dat";
		open $BITORFILE, "> $sfile" or die "ERROR: Cannot create $sfile: $!\n";
		print "writing Tinker bitorsion data to ${sfile}...";
		printf $BITORFILE "Tinker BiTorsion parameter file for fix bitorsion\n\n";
		printf $BITORFILE "%d bitorsion types\n", $parms->{BI_TORSIONS}{counter};
	}

	for $i (grep {!/^counter|TYPE/i} keys %{ $parms->{BI_TORSIONS} }) {
		for $j (keys %{ $parms->{BI_TORSIONS}{$i} }) {
			for $k (keys %{ $parms->{BI_TORSIONS}{$i}{$j} }) {
				for $l (keys %{ $parms->{BI_TORSIONS}{$i}{$j}{$k} }) {
					for $m (keys %{ $parms->{BI_TORSIONS}{$i}{$j}{$k}{$l} }) {
						$c = $parms->{BI_TORSIONS}{$i}{$j}{$k}{$l}{$m}{1};
						next if (!exists($c->{VALS}));
						printf $BITORFILE "%-7s %7d %7d %7d %7d %7d %7d %7d\n","tortors",getAtmTypeID($i),getAtmTypeID($j),
						getAtmTypeID($k),getAtmTypeID($l),getAtmTypeID($m),$c->{NX},$c->{NY};
						for $n (@{ $c->{VALS }}) {
							printf $BITORFILE "%9.3f %9.3f %9.3f\n", $n->[0],$n->[1],$n->[2];
						}
						printf $BITORFILE "\n";
					}
				}
			}
		}
	}
	close $BITORFILE;
}

sub writeIsotopeMass {
	my ($atoms, $datFile) = @_;
	my ($i);

	print $datFile "\nIsotopes\n\n";
	for $i (sort numerically keys %{ $atoms } ) {
		printf $datFile "%5d %8.4f\n", $i, $atoms->{$i}{PARMS}{MASS};
	}
}

sub writeTypeCharges {
	my ($parms, $atoms) = @_;
	my ($i, @types);

	$parms->{PARMS}{TYPE_CHARGES} = "variable scaleQ index 1\n";
	map { 
		$parms->{ATOMTYPES}{$_}{USE_CHARGE} == 1 ? 
			$parms->{PARMS}{TYPE_CHARGES} .= getTypeIDStr($parms->{ATOMTYPES}{$_}) : 
			$parms->{PARMS}{TYPE_CHARGES} .= setTypeQ($parms->{ATOMTYPES}{$_}, $atoms) 
		} 
		sort { 
			$parms->{ATOMTYPES}{$a}{TYPEID} <=> $parms->{ATOMTYPES}{$b}{TYPEID} 
		} grep {
			!/TYPE|counter/i
		} keys %{ $parms->{ATOMTYPES} };
}

sub getTypeIDStr {
	my ($parm) = $_[0];
	my ($str, $q, $tID); 
	
	$q = $parm->{CHARGE};
	$tID = $parm->{TYPEID};
	$str = "variable q${tID} equal ${q}*\${scaleQ}\n" .
			"set type ${tID} charge \${q${tID}}\n";

	return $str;
}

sub setTypeQ {
	my ($parm, $atoms) = @_;
	my ($ffType, $i, $qList, $typeQ, $typeID, $tmp, $str);

	$ffType = $parm->{NAME};
	for $i (values %{ $atoms }) {
		if ($i->{FFTYPE} eq $ffType) {
			$qList->{$i->{CHARGE}}++;
		}
	}

	@{ $tmp } = keys %{ $qList };
	$typeQ = $tmp->[0];
	$typeID = $parm->{TYPEID};
	$str = "variable q${typeID} equal ${typeQ}*\${scaleQ}\n" .
			"set type ${typeID} charge \${q${typeID}}\n";

	$ERRORS{ATOM_CHARGE}{"MULTIPLE CHARGES for $ffType DETECTED: @{$tmp}"}++ if($#{$tmp} > 0);
	
	return $str;
}

sub getRigidOpt {
	my ($parms, $atoms, $bonds, $sOpt) = @_;
	my ($rigidSelect, $mols, $i, $rigidStr, $nMatched);

	$rigidStr = $sOpt->{rigid};
	$rigidSelect = SelectAtoms($rigidStr,$atoms);
	&AddMolsToSelection($rigidSelect, $atoms);
	if (scalar(keys %{ $rigidSelect } == scalar(keys %{ $atoms }))) {
		$parms->{PARMS}{RIGID}{MOLS}{all} = 1;
		return;
	}
	for $i (keys %{ $rigidSelect }) {
		$parms->{PARMS}{RIGID}{MOLS}{${$atoms->{$i}{MOLECULEID}}}{$i} = 1;
	}
	#if polarizable drude or thole
	$nMatched = 0;
	if(exists($sOpt->{coreShell}) and $sOpt->{coreShell}{type}>0) {
		#figure out of the coreshell atoms are the same as the rigid atoms
		for $i (keys %{ $sOpt->{coreShell}{atoms} }) {
			$nMatched++ if (exists($rigidSelect->{$i}));
		}
		$parms->{PARMS}{POLARIZATION}{rigidOpt} = "FLEX";
		if($nMatched == scalar(keys %{ $sOpt->{coreShell}{atoms} })) {
			$parms->{PARMS}{POLARIZATION}{rigidOpt} = "RIGID";
			$parms->{PARMS}{POLARIZATION}{rigidOpt} = "ATOMS" if($nMatched == scalar(keys %{ $atoms }));
		} elsif ($nMatched > 0) {
			$parms->{PARMS}{POLARIZATION}{rigidOpt} = "RIGID|FLEX";
		}
	}
}

sub getFEPopt {
	my ($parms, $atoms, $bonds, $fepStr) = @_;
	my ($fepSelect, $fepAtomTypes, $fepPairType, $nonfepAtomTypes, $fepParentAtomTypes, $valid); 
	my ($newPair, $newType, $i, $j, $k, $pType, $cParm, $first, $idx, $q, $qIdx, $a1, $a2, $a3);

	print "fep...";
	$fepSelect = SelectAtoms($fepStr,$atoms);
	
	if (! keys %{ $fepSelect } ) {
		print "error occurred...skipping...";
		return;
	}
	print "found " . scalar(keys %{ $fepSelect }) . " atoms...writing data...";
	#now we have to get the atom types for the FEP
	for $i (sort numerically keys %{ $fepSelect }) {
		$a1 = $atoms->{$i}{PARMS}{LABEL};
		$a2 = "${a1}fep";
		$a3 = "${a2}" . $atoms->{$i}{CHARGE};
		if (! defined($fepAtomTypes->{$a1})) {
			$fepAtomTypes->{$a1} = $atoms->{$i}{CHARGE}; #fftype for FEP atom
			$fepParentAtomTypes->{$a1}{$a1} = 1;
		} elsif(exists($fepAtomTypes->{$a2}) and $fepAtomTypes->{$a2} == $atoms->{$i}{CHARGE}) {
			$atoms->{$i}{FFTYPE} = $a2;
		} elsif(exists($fepAtomTypes->{$a3}) and $fepAtomTypes->{$a3} == $atoms->{$i}{CHARGE}) {
			$atoms->{$i}{FFTYPE} = $a3;
		} else { 
			($newPair, $newType) = duplicatePair($parms, $a1, $atoms->{$i}{CHARGE}, "fep");
			$fepParentAtomTypes->{ $atoms->{$i}{PARMS}{LABEL} }{ $newType } = 1; #parent of new derived FEP fftype
			$atoms->{$i}{PARMS} = \%{ $newPair };
			$atoms->{$i}{FFTYPE} = $newType;
			$fepAtomTypes->{ $newType } = $atoms->{$i}{CHARGE}; 
		}
	}
	#now we need to handle the case where there are atoms (molecules) in the structure with the same fftype as the FEP atoms
	for $i (sort numerically keys %{ $atoms }) {
		next if (exists($fepSelect->{$i}));
		if (defined($fepAtomTypes->{ $atoms->{$i}{PARMS}{LABEL} })) { 
			$first = 1;
			for $j (keys %{ $fepSelect }) {
				next if($atoms->{$j}{PARMS}{LABEL} ne $atoms->{$i}{PARMS}{LABEL});
				#create a duplicate fftype for the fep atom 
				($newPair, $newType) = duplicatePair($parms, $atoms->{$j}{PARMS}{LABEL}, $atoms->{$j}{CHARGE}, "fep") 
					if ($first==1);
				$first++;
				#and update	
				$fepParentAtomTypes->{ $atoms->{$i}{PARMS}{LABEL} }{ $newType } = 1;
				$atoms->{$j}{PARMS} = \%{ $newPair };
				$atoms->{$j}{FFTYPE} = $newType;
				$fepAtomTypes->{ $newType } = $atoms->{$j}{CHARGE};
			}
			#remove the old atomlabel entry from fepAtomTypes and fepParentAtomTypes
			delete $fepAtomTypes->{ $atoms->{$i}{PARMS}{LABEL} };
			delete $fepParentAtomTypes->{ $atoms->{$i}{PARMS}{LABEL} }{ $atoms->{$i}{PARMS}{LABEL} };
		}
		$nonfepAtomTypes->{ $atoms->{$i}{PARMS}{LABEL} } = $atoms->{$i}{CHARGE};
	}

	&doPairMix($parms); #generate off diagonal entry for fep - nonfep fftypes

	#now create fep charge strings
	$parms->{PARMS}{FEP}{valid}=1;
	$parms->{QEq}{Opt} = 0 if (! exists($parms->{QEq})); #probably already defined, but just in case...
	for $i (sort {$parms->{ATOMTYPES}{$a}{INDEX} <=> $parms->{ATOMTYPES}{$b}{INDEX}} keys %{ $fepAtomTypes }) {
		$idx = ${parms}->{ATOMTYPES}{$i}{INDEX}; #type index
		$q =  $fepAtomTypes->{$i}; #atom charge
		$qIdx = getAtomIndexFromTypeID($atoms, $fepSelect, $idx); #atom index
		if ($parms->{QEq}{Opt}==2 or $q != 0) {
			$parms->{PARMS}{FEP}{var_str} .= "variable             v${idx} equal q[$qIdx] #q = ${q}\n";
		}
		if($q != 0 and $parms->{QEq}{Opt}==0) { #nonQEq case
			$parms->{PARMS}{FEP}{var_str} .= "variable             q${idx} equal \${v${idx}}*v_lambda\n";
			$parms->{PARMS}{FEP}{var_str} .= "variable             dq${idx} equal \${v${idx}}*v_dlambda\n";
			$parms->{PARMS}{FEP}{atom_fix_str} .= "                      atom charge $idx v_q${idx} \&\n";
			$parms->{PARMS}{FEP}{atom_compute_str} .= "                      atom charge $idx v_dq${idx} \&\n";
		}elsif($parms->{QEq}{Opt}==2) { #PQEq
			$parms->{PARMS}{FEP}{var_str} .= "variable             q${idx} equal f_v${idx}_avg\n";
			$parms->{PARMS}{FEP}{var_str} .= "variable             dq${idx} equal v_q${idx}*v_dlambda\n";
		}
	}
	#finally update the VDW hash for the FEP pairs
	&determineAllCoul($parms);
	$valid = 0;
	$valid = 1 if (keys %{ $fepAtomTypes } and keys %{ $nonfepAtomTypes } and keys %{ $fepParentAtomTypes} );
	&createFEPVDWPairStr($parms, $fepAtomTypes, $nonfepAtomTypes, $fepParentAtomTypes) if ($valid);
}

sub getAtomIndexFromTypeID {
	my ($atoms, $atomSelect, $typeIndx) = @_;
	my ($atomID, $i);

	for $i (keys %{ $atomSelect }) {
		if($atoms->{$i}{PARMS}{TYPEID} == $typeIndx) {
			$atomID = $i;
			last;
		}
	}

	return $atomID;
}

sub duplicatePair {
	my ($parms, $atmLabel, $atmCharge, $pairLabel) = @_;
	my ($newPair, $newPairLabel);

	$newPairLabel = "${atmLabel}${pairLabel}";
	$newPairLabel .= "${atmCharge}" if (defined($parms->{ATOMTYPES}{$newPairLabel}));

	$newPair = dclone($parms->{ATOMTYPES}{$atmLabel});
	$newPair->{LABEL} = $newPairLabel;
	$newPair->{NAME} = $newPairLabel;
	$parms->{ATOMTYPES}{counter}++;
	$newPair->{TYPEID} = $newPair->{INDEX} = $parms->{ATOMTYPES}{counter};
	$parms->{ATOMTYPES}{$newPairLabel} = \%{ $newPair };
	$parms->{VDW}{counter}++;

	return ($newPair, $newPairLabel);
}

sub createFEPVDWPairStr {
	my ($parms, $fepAtomTypes, $nonfepAtomTypes, $fepParentAtomTypes) = @_;
	my ($i, $j, $k, $l, $parent, $pOrder, $curr, $nKey);
	my ($nIdx, $addSoftFlag, $tmp, $childParentMap);

	$childParentMap = getChildFromParent($fepParentAtomTypes);
	for $i (keys %{ $fepParentAtomTypes }) { #enumerate over all the parent fep atomtypes
		for $j (keys %{ $parms->{ATOMTYPES} }) { #enumerate over all the atom types
			next if ($j =~ /counter/);
			$l = $j;
			$l = $childParentMap->{$j} if (exists($childParentMap->{$j}));
			if(exists($parms->{VDW}{$i}) and exists($parms->{VDW}{$i}{$l}) and keys %{ $parms->{VDW}{$i}{$l}}) {
				$parent = \%{ $parms->{VDW}{$i}{$l} };
				$pOrder = 1;
				delete $parms->{VDW}{$j}{$i} if ($j ne $i);
			} elsif (exists($parms->{VDW}{$l}) and exists($parms->{VDW}{$l}{$i}) and keys %{ $parms->{VDW}{$l}{$i}}) {
				$parent = \%{ $parms->{VDW}{$l}{$i} };
				$pOrder = 2;
				delete $parms->{VDW}{$i}{$j} if ($j ne $i);
			} else {
				#should never get here...
				undef $parent;
			}
			next if (! defined($parent));
			for $k (keys %{ $fepParentAtomTypes->{$i} }) { #enumerate over all the children of parent atomtype i
				$parms->{VDW}{counter}++;
				if(exists($parms->{VDW}{$j}) and exists($parms->{VDW}{$j}{$k})) { #check for reverse vdw entry for j k
					$parms->{VDW}{$j}{$k} = dclone($parent);
					$curr = \%{ $parms->{VDW}{$j}{$k} };
				} else {
					$parms->{VDW}{$k}{$j} = dclone($parent);
					$curr = \%{ $parms->{VDW}{$k}{$j}}
				}
				$nKey = "$k $j";
				$nKey = "$j $k" if ($pOrder == 2); # swap the keys order based on the original parent - j
				$nIdx = "$parms->{ATOMTYPES}{$k}{INDEX} $parms->{ATOMTYPES}{$j}{INDEX}";
				$nIdx = "$parms->{ATOMTYPES}{$j}{INDEX} $parms->{ATOMTYPES}{$k}{INDEX}" 
					if($parms->{ATOMTYPES}{$k}{INDEX}>$parms->{ATOMTYPES}{$j}{INDEX});

				#now record the nofep - fep pair
				next if (!exists($nonfepAtomTypes->{$j}));
				for $l (values %{ $curr }) {
					next if (! exists($l->{Lammps}{FEP})); #skip all no fep entries
					$addSoftFlag = $l->{Lammps}{FEP}{addsoft};
					#now replace the LAMMPS vdw options for the FEP ones
					$l->{KEY} = $l->{ATOM} = $nKey; #update key and atom strings
					$l->{INDEX} = $parms->{VDW}{counter}; #update vdw index counter
					%{ $tmp }  = %{ $l->{Lammps}{FEP} };
					if (! exists($parms->{VDW}{TYPE}{ $tmp->{pair} })) {
						#add fep pair style
						$parms->{VDW}{TYPE}{ $tmp->{pair} } = dclone($parms->{VDW}{TYPE}{ $l->{Lammps}{LMPNAME} });
						delete $parms->{VDW}{TYPE}{ $tmp->{pair} }{FEP};
						$parms->{VDW}{TYPE}{ $tmp->{pair} }{LMPNAME} = $parms->{VDW}{TYPE}{ $tmp->{pair} }{LMPONAME} = $tmp->{pair};
						$parms->{VDW}{TYPE}{ $tmp->{pair} }{OPTS} = $tmp->{pair_opts};
						$parms->{VDW}{TYPE}{ $tmp->{pair} }{SCC} = 1;
						$parms->{VDW}{TYPE}{ $tmp->{pair} }{USED} = 1;
						%{ $parms->{VDW}{TYPE}{ $tmp->{pair} }{FEP} } = %{ $tmp };
					}
					$l->{Lammps} = \%{ $parms->{VDW}{TYPE}{ $tmp->{pair} } }; 
					push @{ $l->{VALS} }, 1 if (! exists($l->{FEP_MOD}));
					$l->{FEP_MOD} = 1;
					$l->{Lammps}{USED} = 1;
				}
				&addSoftPair($parms, $curr, $j, $k) if($addSoftFlag);
				&recordFEPpair($parms, $curr, $nIdx);
			}
		}
	}
	#cleanup
	for $i (keys %{ $parms->{VDW}{TYPE} }) {
		next if ($i =~ /^(zero|coul)/i);
		delete $parms->{VDW}{TYPE}{$i} if(exists($parms->{VDW}{TYPE}{$i}{USED}) and ! $parms->{VDW}{TYPE}{$i}{USED});
	}
	if(defined($parms->{PARMS}{FEP}{var_list})) {
		for $i (1 .. 5) {
			chop $parms->{PARMS}{FEP}{var_list};
		}
	}
	#&findDuplicateParms($parms, "VDW"); # remove any duplicates...
}

sub createFEPVDWPairStr_old {
	my ($parms, $fepAtomTypes, $nofepAtomTypes, $fepParentAtomTypes) = @_;
	my ($parent, $pOrder, $curr, $nParm, $i, $j, $k, $l, $m, $nKey, $aS, $nIdx, $pTypes, $fpTypes, $tmp);

	for $i (keys %{ $parms->{VDW}{TYPE} }) {
		next if (!exists($parms->{VDW}{TYPE}{$i}{FEP}));
		$nParm = $parms->{VDW}{TYPE}{$i}{FEP}{pair};
		next if (exists($parms->{VDW}{TYPE}{$nParm}));
		$parms->{VDW}{TYPE}{$nParm} = dclone($parms->{VDW}{TYPE}{$i});
		delete $parms->{VDW}{TYPE}{$nParm}{FEP};
		$parms->{VDW}{TYPE}{$nParm}{OPTS} = $parms->{VDW}{TYPE}{$i}{FEP}{pair_opts};
		$parms->{VDW}{TYPE}{$nParm}{LMPONAME} = $parms->{VDW}{TYPE}{$nParm}{LMPNAME} = $nParm;
		$parms->{VDW}{TYPE}{$nParm}{USED} = 0;
	}
	#first fill pTypes array: all fftypes
	@{ $pTypes } = sort { $parms->{ATOMTYPES}{$a}{INDEX} <=> $parms->{ATOMTYPES}{$b}{INDEX} } grep {!/TYPE|counter|-1/i} keys %{ $parms->{VDW} };
	#next fill fpTypes array: all fep fftypes
	@{ $fpTypes } = sort {$parms->{ATOMTYPES}{$a}{INDEX} <=> $parms->{ATOMTYPES}{$b}{INDEX}} keys %{ $fepParentAtomTypes };
	for $i (@{ $fpTypes }) { #we search the vdw array for all the parents and write paramters for the FEP atoms
		for $j (@{ $pTypes }) { 
			undef $parent;
			# for all i(parent),j configurations in the vdw array
			if(exists($parms->{VDW}{$i}) and exists($parms->{VDW}{$i}{$j})) { #set the parent pointer
				$parent = $parms->{VDW}{$i}{$j};
				$pOrder = 1; #parent order, used for setting nKey variable
			} elsif(exists($parms->{VDW}) and exists($parms->{VDW}{$j}{$i})) {
				$parent = $parms->{VDW}{$j}{$i};
				$pOrder = 2;
			} 
			next if	(! defined($parent) or ! keys %{ $parent });
			#using the parent, record data for all FEP children
			for $k (keys %{ $fepParentAtomTypes->{$i} }) { 
				#create entry for the FEP child using the parent's i,j data
				if (! exists($parms->{VDW}{$k}{$j}) and (exists($parms->{VDW}{$j}) and ! exists($parms->{VDW}{$j}{$k}))) {
					#new fep entry from parent
					$parms->{VDW}{$k}{$j} = dclone($parent); 
					$curr = $parms->{VDW}{$k}{$j};
				} elsif (exists($parms->{VDW}{$j}) and exists($parms->{VDW}{$j}) and exists($parms->{VDW}{$j}{$k}) and keys %{ $parms->{VDW}{$j}{$k} }) {
					#point to vdw entry j k
					$curr = $parms->{VDW}{$j}{$k};
				} else {
					#point to vdw entry k j (create if necessary)
					$curr = $parms->{VDW}{$k}{$j};
				}
				$nKey = "$k $j";
				$nKey = "$j $k" if ($pOrder == 2); # swap the keys order based on the original parent - j
				$nIdx = "$parms->{ATOMTYPES}{$k}{INDEX} $parms->{ATOMTYPES}{$j}{INDEX}";
				$nIdx = "$parms->{ATOMTYPES}{$j}{INDEX} $parms->{ATOMTYPES}{$k}{INDEX}" 
					if($parms->{ATOMTYPES}{$k}{INDEX}>$parms->{ATOMTYPES}{$j}{INDEX});

				#i don't know you need the next few lines, since the parent - parent ids should already be present
				#if(exists($fepParentAtomTypes->{ $j })) { 
					#here we make an additional entry for fepParentAtom fepParentAtom, which is not taken care off otherwise
				#	for $m (keys %{ $fepParentAtomTypes->{ $j } }) {
				#		$parms->{VDW}{$k}{$m} = dclone($parent);
				#		for $l (values %{ $parms->{VDW}{$k}{$m} } ) {
				#			$parms->{VDW}{counter}++;
				#			$l->{KEY} = $l->{ATOM} = "$k $m";
				#			$l->{INDEX} = $parms->{VDW}{counter};
				#			delete $l->{Lammps}{FEP};
				#		}
				#	}
				#}
				#next if (exists($fepAtomTypes->{$j})); #this is the case where the jtype is also an FEP type, so don't record the fepPair data

				#now we record the fepPair data
				$aS = 0;
				$parms->{VDW}{counter}++;
				for $l (values %{ $curr }) {
					$l->{KEY} = $l->{ATOM} = $nKey;
					$l->{INDEX} = $parms->{VDW}{counter};
					if (exists($l->{Lammps}{FEP}{addsoft}) and $l->{Lammps}{FEP}{addsoft} == 1 and ! exists($fepParentAtomTypes->{$j})) {
						#add a soft potential to the list (for any style without an explicit soft potential, such as buck or yukawa etc)
						$aS=1;
					} else {
						#update the pair style if both atoms are not fepParentAtomsTypes
						%{ $tmp }  = %{ $l->{Lammps}{FEP} };
						$l->{Lammps} = \%{ $parms->{VDW}{TYPE}{ $l->{Lammps}{FEP}{pair} } }; 
						%{ $l->{Lammps}{FEP} } = %{ $tmp };
						push @{ $l->{VALS} }, 1 if (! exists($l->{FEP_MOD}));
						$l->{FEP_MOD} = 1;
						$l->{Lammps}{USED} = 1;
					}
				}
				&addSoftPair($parms, $curr, $j, $k) if($aS);
				&recordFEPpair($parms, $curr, $nIdx) if (!exists($nofepAtomTypes->{$j}));
			}
		}
	}
	#cleanup
	for $i (keys %{ $parms->{VDW}{TYPE} }) {
		next if ($i =~ /^(zero|coul)/i);
		delete $parms->{VDW}{TYPE}{$i} if(exists($parms->{VDW}{TYPE}{$i}{USED}) and ! $parms->{VDW}{TYPE}{$i}{USED});
	}

	if(defined($parms->{PARMS}{FEP}{var_list})) {
		for $i (1 .. 5) {
			chop $parms->{PARMS}{FEP}{var_list};
		}
	}
	&findDuplicateParms($parms, "VDW"); # remove any duplicates...
}

sub addSoftPair {
	my ($parms, $curr, $itype, $jtype, $vals) = @_;
	my (@tmp, $k);

	@{ $vals } = () if (! defined($vals));
	$curr = \%{ $parms->{VDW}{$itype}{$jtype} };
	$k = 0;
	if(keys %{ $curr }) {
		@tmp = (sort numerically keys %{ $curr } );
		$k = pop @tmp;
		$k++;
	}

	if (! exists($parms->{VDW}{TYPE}{soft})) {
		$parms->{VDW}{TYPE}{soft} = (
									{
										"OPTS"     => "1.0",
										"NAME"     => "soft",
										"PAIRMIX"  => 0,
										"ORDER"    => 1,
										"LMPNAME"  => "soft",
										"LMPONAME" => "soft",
										"USED"     => 1,
										"SCC"      => 1,
									}
								);
	}
	$parms->{VDW}{counter}++;
	$curr->{$k} = (
					{
						"ATOM"     => "$itype $jtype",
						"IGNORE"   => 0,
						"IT"       => "vdw",
						"KEY"      => "$itype $jtype",
						"Lammps"   => \%{ $parms->{VDW}{TYPE}{soft} },
						"TYPE"     => "soft",
						"USED"     => 1,
						"DATA"     => undef,
						"INDEX"    => $k,
						"VALS"     => [@{$vals}],
						"INDEX"    => $parms->{VDW}{counter},
					  }
					);
}

sub recordFEPpair {
	my ($parms, $curr, $nIdx) = @_;
	my ($i, $j, $k, $jVal, $idxStr);

	$idxStr = $nIdx;
	$idxStr =~ s/ /_/g;
	for $i (keys %{ $curr }) {
		next if (! exists($curr->{$i}{Lammps}{FEP}));
		for $k (0 .. $#{ $curr->{$i}{Lammps}{FEP}{parms} }) {
			$j = $curr->{$i}{Lammps}{FEP}{parms}[$k];
			if($j !~ /lambda/i) {
				$jVal = $curr->{$i}{VALS}[ $curr->{$i}{Lammps}{FEP}{idx}[$k] ];
				next if ($jVal == 0);
				$parms->{PARMS}{FEP}{var_str} .= "variable             ${j}var${idxStr} equal ramp(${jVal},0.0)\n";
				$parms->{PARMS}{FEP}{var_str} .= "variable             d${j}var${idxStr} equal ${jVal}*v_dlambda\n";
				$parms->{PARMS}{FEP}{fix_adapt_str} .= "                      pair $curr->{$i}{Lammps}{LMPNAME} ${j} $nIdx v_${j}var${idxStr} \&\n";
				$parms->{PARMS}{FEP}{compute_fep_str} .= "                      pair $curr->{$i}{Lammps}{LMPNAME} ${j} $nIdx  v_d${j}var${idxStr} \&\n";
				$parms->{PARMS}{FEP}{var_list} .= "                                            $curr->{$i}{Lammps}{LMPNAME} $nIdx: ${j} = \${${j}var${idxStr}}\\n &\n";
			} else {
				$parms->{PARMS}{FEP}{fix_adapt_str} .= "                      pair $curr->{$i}{Lammps}{LMPNAME} ${j} $nIdx v_${j} \&\n";
				$parms->{PARMS}{FEP}{compute_fep_str} .= "                      pair $curr->{$i}{Lammps}{LMPNAME} ${j} $nIdx v_dlambda \&\n";
			}
		}
	}
}

sub getChildFromParent {
	my ($parentList) = $_[0];
	my ($i, $j, $children);

	for $i (keys %{ $parentList }) {
		for $j (keys %{ $parentList->{$i} }) {
			$children->{$j}=$i;
		}
	}

	return $children;
}
sub addCoreShellAtoms {
	my ($atoms, $bonds, $parms, $csOpt) = @_;
	my ($atomSel, $ffTyAtoms, $allTypes, $i); 

	$atomSel = SelectAtoms($csOpt->{aStr}, $atoms);
	$ffTyAtoms = SelectAtomsByField($atoms, $bonds, "FFTYPE", $atomSel);
	for $i (keys %{ $atomSel }) {
		$csOpt->{atoms}{$i} = 1;
	}
	$allTypes = SelectAtomsByField($atoms, $bonds, "FFTYPE", $atoms);
	if (! defined($qxFile) and ($csOpt->{type} == 1 or $csOpt->{type} == 2)) {
		$qxFile = getQXfileFromFF($FILES);
		$qxFile = "$Bin/../ff/schwerdtfeger.dff" if (! defined($qxFile));
	}
	&parsePolarizationFile($qxFile, $parms->{ATOMTYPES}) if (defined($qxFile));
	&saveCoreShellData($parms, $ffTyAtoms, $allTypes, $atoms, $bonds, $csOpt);
}

sub getQXfileFromFF {
#TODO: need to enumerate FF files and search for xxx.drude.dat file. If found add to list	
	my ($ff_list) = @_;
	my ($i, $pfix, $d_files);

	$d_files = "";
	for $i (@{ $ff_list }) {
		$pfix = $i->{FF};
		$pfix =~ s/\.\w+$/\.drude\.dat/;
		if (-e $pfix and -r $pfix and -T $pfix) {
			$d_files .= "$pfix ";
		}
	}
	return $d_files;
}

sub parsePolarizationFile {
	my ($qFile, $aTypes) = @_;
	my ($polarInfo, $pFile);

	while ($qxFile =~ /(\S+)/g) {
		$pFile = $1;
		open POLARFILE, $pFile or die "ERROR: Cannot access $pFile: $!\n";
		while (<POLARFILE>) {
			chomp;
			# type m_D/u q_D/e    k_D   alpha/A3 thole
			if ($_ =~ /^\s*(\S+)\s+(\d+\.?\d*)\s+(\-?\d+\.?\d*)\s+(\d+\.?\d*)\s*(\d*\.?\d*)\s*(\d*\.?\d*)/) {
				next if ($1 eq "#");
				$polarInfo->{$1} = (
									{
										m     => $2,
										q     => $3,
										k     => $4/4.184,
									}
				);
				if($5) {
					$polarInfo->{$1}{alpha} = $5
				} else {
					$polarInfo->{$1}{alpha} = $polarInfo->{$1}{q}**2/3.0114;
				}
				$polarInfo->{$1}{thole} = $6 if ($6);
			}
		}
		close POLARFILE;
	}
	if(defined($polarInfo)) {
		for $i (keys %{ $aTypes }) {
			if(exists($polarInfo->{$i})) {
				$aTypes->{$i}{polarInfo} = \%{ $polarInfo->{$i} };
			} elsif (exists($polarInfo->{ $aTypes->{$i}{ATOM}})) {
				$aTypes->{$i}{polarInfo} = \%{ $polarInfo->{ $aTypes->{$i}{ATOM}} };
			}
		}
	}
}

sub getPolarizationInputCmds {
	my ($parms, $csOpt) = @_;
	my ($i, $j, $cType, $tmp);

	$parms->{PARMS}{POLARIZATION}{type} = $csOpt->{type};
	if ($csOpt->{type} == 0) {
		###adiabatic core-shell###
		#first change the pair/style to a core-shell compatible one
		for $i (keys %{ $parms->{PARMS}{POLARIZATION}{SHELLS}}) {
			$cType = $parms->{PARMS}{POLARIZATION}{SHELLS}{$i}{core};
			for $j (values %{ $parms->{VDW}{$i}{'*'} }, values %{ $parms->{VDW}{$cType}{$i}}) {
				$j->{Lammps}{LMPNAME} =~ s/coul.long\S*/coul\/long\/cs/i;
				if ($j->{Lammps}{LMPNAME} =~ /charmm/i) {
					$j->{Lammps}{LMPNAME} =~ s/charmm/cut/;
					$j->{Lammps}{opts} = 10;
				}
			}
		}
		#add the core and shell types
		$parms->{PARMS}{POLARIZATION}{input} = "group SHELLS type";
		for $i (keys %{ $parms->{PARMS}{POLARIZATION}{SHELLS}}) {
			$parms->{PARMS}{POLARIZATION}{input} .= " $parms->{ATOMTYPES}{$i}{INDEX}";
		}
		$parms->{PARMS}{POLARIZATION}{input} .= "\ngroup CORES type";
		for $i (keys %{ $parms->{PARMS}{POLARIZATION}{CORES}}) {
			$parms->{PARMS}{POLARIZATION}{input} .= " $parms->{ATOMTYPES}{$i}{INDEX}";
		}
		#now add some thermo and pressure modifys
		$parms->{PARMS}{POLARIZATION}{input} .= <<DATA;

comm_modify vel yes
compute TATOMS all temp/cs CORES SHELLS
compute thermo_press_lmp all pressure thermo_temp 
thermo_modify temp TATOMS press thermo_press_lmp
DATA
	} else {
		###Drude induced dipoles w(type==1)/wo(type==2) thole core-shell screening###
		#first make fix drude command
		$parms->{PARMS}{POLARIZATION}{input} = "fix DRUDE all drude ";
		for $i (sort { $parms->{ATOMTYPES}{$a}{TYPEID}<=>$parms->{ATOMTYPES}{$b}{TYPEID} } grep {!/TYPE|counter|-1/i} keys %{ $parms->{ATOMTYPES} }) {
			if(exists($parms->{PARMS}{POLARIZATION}{CORES}{$i})) {
				$cType = "C";
				$tmp->{coreTypes} .= " $parms->{ATOMTYPES}{$i}{TYPEID}";
			} elsif (exists($parms->{PARMS}{POLARIZATION}{SHELLS}{$i})) {
				$cType = "D";
				$tmp->{shellTypes} .= " $parms->{ATOMTYPES}{$i}{TYPEID}";
			} else {
				$cType = "N";
				$tmp->{nonpTypes} .= " $parms->{ATOMTYPES}{$i}{TYPEID}";
			}
			$parms->{PARMS}{POLARIZATION}{input} .= " $cType";
		}
		$parms->{PARMS}{POLARIZATION}{input} .= "\ngroup ATOMS type ";
		if(exists($tmp->{nonpTypes})) {
			$parms->{PARMS}{POLARIZATION}{input} .= "$tmp->{nonpTypes}";
		}
		$parms->{PARMS}{POLARIZATION}{input} .= <<DATA;
$tmp->{coreTypes}
group CORES type $tmp->{coreTypes}
group SHELLS type $tmp->{shellTypes}

comm_modify vel yes
compute TATOMS ATOMS temp
compute TDRUDES all temp/drude
variable Tshell equal c_TDRUDES[2]
variable Tcore equal c_TDRUDES[1]
variable drude_T equal 1.0 #tempeature of Drude particle thermostat

fix MOMENTUM all momentum 100 linear 1 1 1 

compute thermo_press_lmp all pressure thermo_temp 
thermo_modify temp TATOMS press thermo_press_lmp
thermo_style custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press vol v_Tshell v_Tcore
thermo_modify line multi format float %14.6f flush yes


DATA
	}
}

sub saveCoreShellData {
	my ($parms, $ffTyAtoms, $allTypes, $atoms, $bonds, $csOpt) = @_;
	my ($i, $cType, $sType, $sAtm, $j, $curr, $aType);

	&setVDWtypes($parms);
	for $i (keys %{ $ffTyAtoms }) {
		$cType = $i;
		if(scalar(keys %{ $ffTyAtoms->{$cType} }) != scalar(keys %{ $allTypes->{$cType} })) {
			$cType = updateCoreAtomType($atoms, $parms, $ffTyAtoms->{$cType}, $cType);
		}
		$sType = "${cType}_shell";
		$aType = $parms->{ATOMTYPES}{$cType}{polarInfo};
		if (! exists($parms->{ATOMTYPES}{$sType})) {
			#first create a new shell atom fftype
			$sAtm = (
						{
							m     => $aType->{m},
							q     => $aType->{q},
							k1    => $aType->{k},
							k2    => 0,
							alpha => $aType->{alpha}, #alpha = q^2/3.0114
						}
			);
			$sAtm->{thole} = $aType->{thole} if (exists($aType->{thole})); 
			$sAtm->{m} = 0 if ($csOpt->{type} == 0); #massless core-shell particle
			if(exists($parms->{ATOMTYPES}{$sType}{polarInfo})) {
				for $j (keys %{ $parms->{ATOMTYPES}{$sType}{polarInfo} }) {
					$sAtm->{$j} = $parms->{ATOMTYPES}{$sType}{polarInfo}{$j};
				}
			}			
			&createShellAtomType($parms->{ATOMTYPES}, $cType, $sType, $sAtm);
		}
		#create shell VDW parms
		&createShellVDW($parms, $cType, $sType, $csOpt, $sAtm);
		#add bond type and zero other valence
		&createShellValence($parms, $cType, $sType, $sAtm->{k1}, $sAtm->{k2});
		for $j (sort numerically keys %{ $ffTyAtoms->{$i}}) {
			#create the shell atom and core - shell bond
			&createShellAtom($atoms, $bonds, $parms, $sType, $j, $csOpt->{type});
		}
		#record core-shell data for future use
		$parms->{PARMS}{POLARIZATION}{CORES}{$cType} = $parms->{ATOMTYPES}{$cType};
		$parms->{PARMS}{POLARIZATION}{SHELLS}{$sType} = $parms->{ATOMTYPES}{$sType};
		$parms->{PARMS}{POLARIZATION}{SHELLS}{$sType}{core} = $cType;
	}

}
sub createShellAtom {
	my ($atoms, $bonds, $parms, $sType, $cID, $pType) = @_;
	my ($sID, $i, $j, $CoM, $molAtoms);

	$sID =  scalar(keys %{ $atoms }) + 1;
	%{ $atoms->{$sID} } = %{ $atoms->{$cID} };
	$atoms->{$sID}{MOLECULEID} = \${ $atoms->{$cID}{MOLECULEID} };
	$atoms->{$sID}{FFTYPE} = $sType;
	$atoms->{$sID}{INDEX} = $sID;
	$atoms->{$sID}{CHARGE} = $parms->{ATOMTYPES}{$sType}{polarInfo}{q};
	$atoms->{$sID}{MOLECULE}{MOLSIZE}++;
	$atoms->{$sID}{MOLECULE}{MEMBERS}{$atoms->{$sID}{MOLECULE}{MOLSIZE}} = $sID;
	$atoms->{$sID}{CORE_PARENT}{TYPE} = $atoms->{$cID}{FFTYPE};
	$atoms->{$sID}{CORE_PARENT}{ID} = $cID;
	push @{ $bonds->{$cID} }, $sID;
	push @{ $bonds->{$sID} }, $cID;
	$atoms->{$cID}{CHARGE} -= $parms->{ATOMTYPES}{$sType}{polarInfo}{q};
	if($pType > 0) {
		#get the normal to the least squares plane of the molecule
		$molAtoms = ();
		for $j (keys %{ $MOLS->{ ${ $atoms->{$cID}{MOLECULEID} } }{MEMBERS} } ) {
			next if ($j == $atoms->{$sID}{MOLECULE}{MOLSIZE});
			$molAtoms->{ATOMS}{$j} = $atoms->{$j};
		}
		$molAtoms->{NORMAL}{XCOORD} = $molAtoms->{NORMAL}{YCOORD} = $molAtoms->{NORMAL}{ZCOORD} = 1**(1/3);
#			if(scalar(keys %{ $molAtoms->{ATOMS} })>3) {
		#			$CoM = CoM($molAtoms->{ATOMS});
		#			&GetLeastSquaresPlane($molAtoms, $CoM);
			#place the shell particle 0.1A atoms away from the atom, along the plane normal
			#	for $i ("XCOORD", "YCOORD", "ZCOORD") {
			#	$atoms->{$sID}{$i} += $molAtoms->{NORMAL}{$i}*.1;
			#}
			#		} else {
			for $i ("XCOORD", "YCOORD", "ZCOORD") {
				$atoms->{$sID}{$i} += 0.1 * (rand() - 0.5);
			}
			#		}
	}
}

sub createShellVDW {
	my ($parms, $cType, $sType, $csOpt, $sAtm) = @_;
	my ($rec, $coul, $i, $lammpsType, $iType, $aType, $thole);

	if (! exists($parms->{ATOMTYPES}{$cType}{polarInfo})) {
		$ERRORS{"POLARIZATION TYPE"}{$cType}++;
		return;
	}
	$aType = $parms->{ATOMTYPES}{$cType}{polarInfo};

	for $i (keys %{ $parms->{VDW}{TYPE}}) {
		if(!$parms->{VDW}{TYPE}{$i}{SCC}) {
			$lammpsType = $i;
			last;
		}
	}
	$coul = "0111.00000";
	$rec = (
			{
				CORE    => $cType,
				ATOM    => "$sType * ",
				IGNORE  => 0,
				IT      => "vdw",
				KEY     => "$sType * ",
				TYPE    => "LJ_6_12",
				USED    => 0,
				VALS    => [0,0.1],
				COUL    => $coul,
			}
	);
	if (defined($lammpsType)) {
		$rec->{Lammps} = \%{ $parms->{VDW}{TYPE}{$lammpsType} };
		$rec->{TYPE} = $parms->{VDW}{TYPE}{$lammpsType}{NAME};
	}
	%{ $parms->{VDW}{$sType}{"*"}{1} } = %{ $rec }; #all types
	if (!exists($parms->{VDW}{$sType}{$sType})) {	
		$iType = \%{ $parms->{VDW}{$sType}{$sType} };
		%{ $iType->{1} } = %{ $rec };
		$iType->{1}{ATOM} = $iType->{1}{KEY} = "$sType $sType ";	
	}
	if ($csOpt->{type} == 2) { 
		#first create Thole type if not present
		if (! exists($parms->{VDW}{TYPE}{thole})) {
			$parms->{VDW}{TYPE}{thole} = (
											{
												"OPTS"     => "2.6 10",
												"NAME"     => "THOLE",
												"PAIRMIX"  => 1,
												"ORDER"    => 1,
												"LMPNAME"  => "thole",
												"LMPONAME" => "thole",
												"SCC"      => 1,
											}
										);
		}
		#now create new entry for thole type
		$thole = 1.22;
		$thole = $aType->{thole} if(exists($aType->{thole}));
		#first for core
		$rec = (
				{
					ATOM    => "$cType $cType ",
					IGNORE  => 0,
					IT      => "vdw",
					KEY     => "$cType $cType ",
					TYPE    => "THOLE",
					USED    => 0,
					VALS    => [$aType->{alpha},$thole],
					COUL    => $coul,
					Lammps  => $parms->{VDW}{TYPE}{thole},
				}
		);
		$rec->{IGNORE} = 1 if (! exists($aType->{thole}));
		$iType = \%{ $parms->{VDW}{$cType}{$cType} };
		$i = scalar(keys %{ $iType }) + 1;
		%{ $iType->{$i} } = %{ $rec };
		#and for shell
		$rec->{ATOM} = $rec->{KEY} = "$sType $sType ";
		$rec->{CORE} = $cType;
		$iType = \%{ $parms->{VDW}{$sType}{$sType} };
		$i = scalar(keys %{ $iType }) + 1;
		%{ $iType->{$i} } = %{ $rec };
	} elsif($csOpt->{type} == 3) {
		$parms->{VDW}{$sType}{$sType}{1}{IGNORE} = 1;
		$parms->{VDW}{$sType}{$sType}{1}{Lammps}{pmix} = 0;
	}
}

sub createShellValence {
	my ($parms, $t1, $t2, $fc1, $fc2) = @_;
	my ($curr, $i, $aSet, $tmp);

	#bonds
	@{ $tmp } = ($t2, $t1);
	$curr = addDummyValenceType($parms, "BONDS", $tmp, "HARMONIC");
	$curr->{1}{Lammps} = (
						{
							name    => "harmonic",
							opts    => "",
							missing => 0,
						}
	);
	$fc1 = 500 if ($fc1<=0);
	$fc2 = 500 if ($fc2<=0);
	$curr->{1}{VALS} = [$fc1/2,0];
	$curr->{1}{IGNORE} = 0;
	#angles
	push @{ $tmp },  'X';
	&addDummyValenceType($parms, "ANGLES", $tmp, "THETA_HARM");
	#torsions
	push @{ $tmp },  "X";
	&addDummyValenceType($parms, "TORSIONS", $tmp, "HARMONIC");
	#invesions
	@{ $tmp } = ($t1, $t2, "X", "X");
	&addDummyValenceType($parms, "INVERSIONS", $tmp, "IT_IJKL");
}

sub addDummyValenceType {
	my ($parms, $vType, $aTypes, $name) = @_;
	my ($i, $curr, $alreadyPresent);

	$alreadyPresent = 1;
	$curr = \%{ $parms->{$vType} };
	for $i (@{ $aTypes }) {
		if (exists($curr->{$i})) {
			$curr = \%{ $curr->{$i} };
		} else {
			$alreadyPresent = 0;
			last;
		}
	}
	return if ($alreadyPresent);
	$curr = \%{ $parms->{$vType} };
	for $i (reverse @{ $aTypes }) {
		if (exists($curr->{$i})) {
			$curr = \%{ $curr->{$i} };
		} else {
			$alreadyPresent = 0;
			last;
		}
	}
	return if ($alreadyPresent);

	$curr = \%{ $parms->{$vType} };
	for $i (@{ $aTypes }) {
		$curr = \%{ $curr->{$i} };
	}
	$curr->{1} = (
					{
						INDEX   => 1,
						TYPE    => "$name",
						VALS    => [],
						USED    => 0,
						KEY     => "@{ $aTypes } ",
						IGNORE  => 1,
						CTYPE   => "",
						Lammps  => {},	
					}
	);
	
	return $curr;
}

sub createShellAtomType {
	my ($aTypes, $cType, $nType, $polarInfo) = @_;
	
	$aTypes->{$nType} = dclone($aTypes->{$cType});
	$aTypes->{$nType}{MASS} = $polarInfo->{m};
	$aTypes->{$nType}{LABEL} = $nType;
	$aTypes->{$nType}{CHARGE} = $polarInfo->{q};
	$aTypes->{$nType}{polarInfo} = \%{ $polarInfo };
	$aTypes->{$nType}{USED} = 1;
	$aTypes->{$nType}{FFTYPE} = $nType;
	#adjust charge and mass on core to compensate for shell
	$aTypes->{$cType}{MASS} -= $polarInfo->{m};
	$aTypes->{$cType}{CHARGE} -= $polarInfo->{q};
}


sub updateCoreAtomType {
	my ($atoms, $parms, $aList, $ffType) = @_;
	my ($cType, $i);

	$cType = $ffType;
	while(exists($parms->{ATOMTYPES}{$cType})) {
		$cType .= "core";
	}
	$parms->{ATOMTYPES}{$cType} = dclone($parms->{ATOMTYPES}{$ffType});
	$parms->{ATOMTYPES}{$cType}{LABEL} = $cType;
	$parms->{EQUIVALENCE}{$cType}{$ffType} = 1;
	for $i (keys %{ $aList }) {
		$atoms->{$i}{FFTYPE} = $cType;
	}
	return $cType;
}

sub getOptSelStr {
	my ($optStr, $sstr, $retStr) = @_;
	
	$optStr =~ s/'/"/g;
	while ($optStr =~ /$sstr\S*:\s*(\"[^"]+\")/gi) {
		&getStrOpt($1,\@{ $retStr });
	}
	return \@{ $retStr };
}

sub getStrOpt {
	my ($i, $retStr) = @_;
	my ($j);

	$i =~ s/'/"/g;
	while ($i =~ /"([^"]+)"/g) {
		$j = $1;
		$j = "index>0" if($j =~ /all/i);
		push @{ $retStr }, $j;
	}
}

sub setVDWtypes {
	my ($parms) = $_[0];
	my ($j, $i, $l, $k, $lammpsType);

	for $j (grep {!/TYPE|Counter/i} keys %{ $parms->{VDW} }) {
		for $i (grep {!/TYPE|Counter/i} keys %{ $parms->{VDW}{$j} }) {
			delete $parms->{VDW}{$j}{$i} if (! keys %{ $parms->{VDW}{$j}{$i} });
			for $l (keys %{ $parms->{VDW}{$j}{$i} }) {
				$k = $parms->{VDW}{$j}{$i}{$l};
				$lammpsType = $k->{Lammps}{name};
				next if (! defined($lammpsType) or ! $lammpsType);
				$k->{Lammps}{scc} = 0 if (! defined($k->{Lammps}{scc}));
				if (! exists($parms->{VDW}{TYPE}{$lammpsType})) {
					$parms->{VDW}{TYPE}{$lammpsType} = (
													{
														"OPTS"     => $k->{Lammps}{opts},
														"NAME"     => $k->{TYPE},
														"PAIRMIX"  => $k->{Lammps}{pmix},
														"ORDER"    => 2,
														"LMPNAME"  => $lammpsType,
														"LMPONAME" => $lammpsType,
														"SCC"      => $k->{Lammps}{scc},
													}
												   );
   					$parms->{VDW}{TYPE}{$lammpsType}{FILE} = $k->{VALS}[0] if ($k->{TYPE} =~ /eam/i);
					$parms->{VDW}{TYPE}{$lammpsType}{FEP} = $k->{Lammps}{fep} if (exists($k->{Lammps}{fep}));
				}
				delete $k->{Lammps};
				$k->{Lammps} = \%{ $parms->{VDW}{TYPE}{$lammpsType} };	
			}
		}
	}
}

sub writeCMAPfile {
	my ($parms, $suffix) = @_;

	$parms->{PARMS}{CMAP}{FILE} = "${suffix}.cmap";
	open CMAPFILE, "> $parms->{PARMS}{CMAP}{FILE}" or die "ERROR: Cannot create CHARMM CMAP file $parms->{PARMS}{CMAP}{FILE}: $!\n";
	print CMAPFILE $parms->{PARMS}{CMAP}{FILE_DATA};
	close CMAPFILE;
}

sub writeAmoebaFF {
	my ($parms, $ffs, $sfx) = @_;
	my ($amoeba_file, $indx, $i, $j, $k, $l, $m, $n, $curr, @aTypes);
	my ($p);

	$amoeba_file = basename($sfx) . ".amoeba.prm";
	@aTypes = sort {$parms->{ATOMTYPES}{$a}{INDEX} <=> $parms->{ATOMTYPES}{$b}{INDEX}} grep {!/^counter/} keys %{ $parms->{ATOMTYPES} };
	push @aTypes, "X";
	for $i (@aTypes) { #1 body
		if($i ne "X") {
			#atoms
			$curr = $parms->{ATOMTYPES}{$i};
			$p->{atoms} .= sprintf("%4s %10d %5d    %-5s \"%-25s\" %10d%10.3f%5d\n","atom",$curr->{INDEX},$curr->{INDEX},
			$curr->{LABEL},$curr->{TINKER_DESCRP},$curr->{ELENUM},$curr->{MASS},$curr->{TINKER_NBOND});
			$curr->{TINKER_ID} = $curr->{INDEX};
			#multipoles
			$p->{multipole} .= sprintf("%-10s %5d","multipole",$curr->{INDEX});
			$k = 0;
			for $j (0 .. $#{ $curr->{MULTIPOLE}{ATMTYPES} }) {
				$p->{multipole} .= sprintf(" %5d",$curr->{MULTIPOLE}{AXIS}[$j] * getAtmTypeID( $curr->{MULTIPOLE}{ATMTYPES}[$j] ));
				$k++;
			}
			for $l (($k+1) .. 4) {
				$p->{multipole} .= sprintf("%5s","");
			}
			$p->{multipole} .= sprintf(" %10.5f\n", $curr->{MULTIPOLE}{MONOPOLE});
			$k = 0;
			for $j (@{ $curr->{MULTIPOLE}{DIPOLE} }, @{ $curr->{MULTIPOLE}{QUADRUPOLE} }) {
				$p->{multipole} .= sprintf("%39s"," ") if($k == 0 or $k == 3 or $k == 4 or $k == 6);
				$p->{multipole} .= sprintf("%10.5f",$j);
				$p->{multipole} .= sprintf("\n") if($k == 2 or $k == 3 or $k == 5);
				$k++;
			}
			$p->{multipole} .= sprintf("\n");
			#atomic polarizabilities
			$p->{polarize} .= sprintf("%-10s %5d %20.4f %10.4f","polarize",$curr->{INDEX},$curr->{POLARIZE}{p1},$curr->{POLARIZE}{p2});
			for $j (@{ $curr->{POLARIZE}{ATOMS} }) {
				$p->{polarize} .= sprintf(" %5d",getAtmTypeID( $j ));
			}
			$p->{polarize} .= "\n";
		}
		for $j (@aTypes) { #2 body
			if(exists($parms->{VDW}) and exists($parms->{VDW}{$i}) and 
			exists($parms->{VDW}{$i}{$j}) and keys %{ $parms->{VDW}{$i}{$j} }) { #vdw
				$indx = findDuplicate($parms->{VDW}{$i}{$j}, "AMOEBA ");
				$curr = $parms->{VDW}{$i}{$j}{$indx};
				if($i eq $j) {
					$p->{vdws} .= sprintf("%-10s %5d","vdw",getAtmTypeID($i));
					$p->{vdws} .= sprintf(" %20.4f %10.4f",$curr->{VALS}[1],$curr->{VALS}[0]);
					$p->{vdws} .= sprintf " %10.4f",$curr->{VALS}[2] if ($#{ $curr->{VALS} } == 2);
					$p->{vdws} .= "\n";
				} else {
					$p->{vdwspr} .= sprintf("%-10s %5d %5d","vdwpr",getAtmTypeID($i),getAtmTypeID($j));
					$p->{vdwspr} .= sprintf(" %20.4f %10.4f",$curr->{VALS}[1],$curr->{VALS}[0]);
					$p->{vdwspr} .= sprintf " %10.4f",$curr->{VALS}[2] if ($#{ $curr->{VALS} } == 2);
					$p->{vdwspr} .= "\n";
				}
			}
			if(exists($parms->{BONDS}) and exists($parms->{BONDS}{$i}) and exists($parms->{BONDS}{$i}{$j})) { #bonds
				$indx = findDuplicate($parms->{BONDS}{$i}{$j}, "CLASS2");
				$curr = $parms->{BONDS}{$i}{$j}{$indx};
				$p->{bonds} .= sprintf("%-10s %5d %5d %10.2f %10.4f\n","bond",getAtmTypeID($i),
				getAtmTypeID($j),$curr->{VALS}[1],$curr->{VALS}[0]);
			}
			if(exists($parms->{PI_TORSIONS}) and exists($parms->{PI_TORSIONS}{$i}) and 
			exists($parms->{PI_TORSIONS}{$i}{$j})) { #pi-torsions
				$indx = findDuplicate($parms->{PI_TORSIONS}{$i}{$j}, "PITORSION");
				$curr = $parms->{PI_TORSIONS}{$i}{$j}{$indx};
				$p->{pitors} .= sprintf("%-10s %5d %5d %10.2f\n","pitors",getAtmTypeID($i),
				getAtmTypeID($j),$curr->{VALS}[0]);
			}			
			for $k (@aTypes) { #3 body
				if(exists($parms->{ANGLES}) and exists($parms->{ANGLES}{$i}) and 
				exists($parms->{ANGLES}{$i}{$j}) and exists($parms->{ANGLES}{$i}{$j}{$k})) { #angles
					$indx = findDuplicate($parms->{ANGLES}{$i}{$j}{$k},"AMOEBA");
					$curr = $parms->{ANGLES}{$i}{$j}{$k}{$indx};
					if($curr->{VALS}[0] == 0) {
						$p->{angles} .= sprintf("%-10s ", "angle");
					} else {
						$p->{angles} .= sprintf("%-10s ", "anglep");
					}
					$p->{angles} .= sprintf("%5d %5d %5d %10.2f %10.4f",getAtmTypeID($i),getAtmTypeID($j),
					getAtmTypeID($k), $curr->{VALS}[3],$curr->{VALS}[2]);
					$p->{angles} .= sprintf(" %10.4f", $curr->{VALS}[10]) if ($#{ $curr->{VALS} } > 10);
					$p->{angles} .= sprintf(" %10.4f", $curr->{VALS}[18]) if ($#{ $curr->{VALS} } > 18);
					$p->{angles} .= "\n";
					$p->{bndang} .= sprintf("%-10s %5d %5d %5d %10.2f %10.4f\n","strbnd",getAtmTypeID($i),
					getAtmTypeID($j),getAtmTypeID($k), $curr->{BOND_ANGLE}{VALS}[1],$curr->{BOND_ANGLE}{VALS}[0])
						if($curr->{BOND_ANGLE}{VALS}[0] > 0 and $curr->{BOND_ANGLE}{VALS}[1] > 0);
					$p->{ub} .= sprintf("%-10s %5d %5d %5d %10.2f %10.4f\n","ureybrad",getAtmTypeID($i),
					getAtmTypeID($j),getAtmTypeID($k), $curr->{UREY_BRADLEY}{VALS}[0],$curr->{UREY_BRADLEY}{VALS}[1])
						if($curr->{UREY_BRADLEY}{VALS}[0] > 0 and $curr->{UREY_BRADLEY}{VALS}[1] > 0);
				}
				for $l (@aTypes ) { #4 body
					if(exists($parms->{INVERSIONS}) and exists($parms->{INVERSIONS}{$j}) and exists($parms->{INVERSIONS}{$j}{$i}) and 
					exists($parms->{INVERSIONS}{$j}{$i}{$k}) and exists($parms->{INVERSIONS}{$j}{$i}{$k}{$l})) { #impropers
						$indx = findDuplicate($parms->{INVERSIONS}{$j}{$i}{$k}{$l}, "AMOEBA");
						$curr = $parms->{INVERSIONS}{$j}{$i}{$k}{$l}{$indx};
						$p->{inversions} .= sprintf("%-10s %5d %5d %5d %5d %10.2f\n","opbend",getAtmTypeID($i),
						getAtmTypeID($j),getAtmTypeID($k),getAtmTypeID($l),$curr->{VALS}[0]);
					}
					if(exists($parms->{TORSIONS}) and exists($parms->{TORSIONS}{$i}) and exists($parms->{TORSIONS}{$i}{$j}) and 
					exists($parms->{TORSIONS}{$i}{$j}{$k}) and exists($parms->{TORSIONS}{$i}{$j}{$k}{$l})) { #dihedrals
						$indx = findDuplicate($parms->{TORSIONS}{$i}{$j}{$k}{$l}, "FOURIER");
						$curr = $parms->{TORSIONS}{$i}{$j}{$k}{$l}{$indx};
						$p->{torsions} .= sprintf("%-10s %d %d %d %d","torsion",getAtmTypeID($i),getAtmTypeID($j),
						getAtmTypeID($k),getAtmTypeID($l));
						$m = 0;
						while($m < $#{ $curr->{VALS}} ) {
							$p->{torsions} .= sprintf(" %5.3f %5.1f %2d",$curr->{VALS}[$m+1],$curr->{VALS}[$m+3],$curr->{VALS}[$m+2]);
							$m+=3;
						}
						$p->{torsions} .= "\n";
					}
					#5 body torsion torsion is dealt with below
				}
			}
		}
	}

	open AMOEBA, "> $amoeba_file" or die "ERROR: Cannot create $amoeba_file: $!\n";
	print AMOEBA writeTinkerFFheader("Force Field Definition"); 
	print AMOEBA "forcefield           AMOEBA-BIO-2018\n";
	for $i (sort alphabetically keys %{ $parms->{TINKER}{HEADERS}}) {
		printf AMOEBA "%-20s %s\n",$i,$parms->{TINKER}{HEADERS}{$i}
	}
	print AMOEBA writeTinkerFFheader("Literature References"); 
	print AMOEBA "Walker, B., Liu, C., Wait, E., Ren, P., J. Comput. Chem. 2022, 1. https://doi.org/10.1002/jcc.26954\n";
	print AMOEBA "Wu, J.C.; Chattree, G.; Ren, P.Y.; Automation of AMOEBA polarizable force field\n";
	print AMOEBA "parameterization for small molecules. Theor Chem Acc.\n";
	print AMOEBA writeTinkerFFheader("Atom Type Definitions");
	print AMOEBA $p->{atoms};
	if(exists($p->{vdws})) {
		print AMOEBA writeTinkerFFheader("Van der Waals Parameters");
		print AMOEBA $p->{vdws};
	}
	if(exists($p->{vdwspr})) {
		print AMOEBA writeTinkerFFheader("Van der Waals Pair Parameters");
		print AMOEBA $p->{vdwspr};
	}
	if(exists($p->{bonds})) {
		print AMOEBA writeTinkerFFheader("Bond Stretching Parameters");
		print AMOEBA $p->{bonds};
	}
	if(exists($p->{angles})) {
		print AMOEBA writeTinkerFFheader("Angle Bending Parameters");
		print AMOEBA $p->{angles};
	}
	if(exists($p->{bndang})) {
		print AMOEBA writeTinkerFFheader("Stretch-Bend Parameters");
		print AMOEBA $p->{bndang};
	}
	if(exists($p->{ub})) {	
		print AMOEBA writeTinkerFFheader("Urey-Bradley Parameters");
		print AMOEBA $p->{ub};
	}
	if(exists($p->{inversions})) {	
		print AMOEBA writeTinkerFFheader("Out-of-Plane Bend Parameters");
		print AMOEBA $p->{inversions};
	}
	if(exists($p->{torsions})) {	
		print AMOEBA writeTinkerFFheader("Torsional Parameters");
		print AMOEBA $p->{torsions};
	}
	if(exists($p->{pitors})) {	
		print AMOEBA writeTinkerFFheader("Pi-Torsion Parameters");
		print AMOEBA $p->{pitors};
	}
	if(exists($p->{strtors})) {	
		print AMOEBA writeTinkerFFheader("Stretch-Torsion Parameters");
		print AMOEBA $p->{strtors};
	}
	if(exists($p->{angtors})) {	
		print AMOEBA writeTinkerFFheader("Angle-Torsion Parameters");
		print AMOEBA $p->{angtors};
	}
	if(exists($parms->{BI_TORSIONS}) and $parms->{BI_TORSIONS}{counter} > 0) {
		print AMOEBA writeTinkerFFheader("Torsion-Torsion Parameters");
		&writeTinkerBiTorsionFile($parms,"",0,\*AMOEBA);
	}
	print AMOEBA writeTinkerFFheader("Atomic Multipole Parameters");
	print AMOEBA $p->{multipole};
	print AMOEBA writeTinkerFFheader("Dipole Polarizability Parameters");
	print AMOEBA $p->{polarize};

	close AMOEBA;
	return $amoeba_file;
}

sub getAtmTypeID {
	my ($itype) = $_[0];

	return 0 if ($itype eq "X");
	return $PARMS->{ATOMTYPES}{$itype}{INDEX} if (exists($PARMS->{ATOMTYPES}{$itype}));
	die "ERROR: Atomtype $itype not in forcefield!\n" if(!exists($PARMS->{ATOMTYPES}{$itype}));
}

sub writeTinkerFFheader {
	my ($header) = $_[0];
	my ($retStr, $n, $i);

	$n = length($header) + 8;
	$retStr = "\n\n    ";
	for $i (1 .. $n) {
		$retStr .= "#"
	}
	$retStr .= "\n";
	$retStr .= "    ##";
	for $i (1 .. ($n-4)) {
		$retStr .= " ";
	}
	$retStr .= "##\n    ##  $header  ##\n    ";
	$retStr .= "##";
	for $i (1 .. ($n-4)) {
		$retStr .= " ";
	}
	$retStr .= "##\n    ";	
	for $i (1 .. $n) {
		$retStr .= "#"
	}
	$retStr .= "\n\n\n";
	return $retStr;
}
sub deleteBonds {
	my ($atoms, $bonds, $alist, $ffType) = @_;
	my ($atomSel, $i, $j, $k);

	$alist->{aStr} = BuildAtomSelectionString("index>0") if ($ffType == 5);
	$atomSel = SelectAtoms($alist->{aStr}, $atoms);

	for $i (keys %{ $atomSel }) {
		#first delete all entries for atom $i in other atom bondlist
		for $j (@{ $bonds->{$i} }) {
			for $k (0 .. $#{ $bonds->{ $j }}) {
				if($bonds->{$j}[$k] == $i) {
					splice @{ $bonds->{$j} },$k, 1;
				}
			}
		}
		$bonds->{$i} = ();
	}
}

sub init {
	my ($i, $j, $inputFFType, $inputList, $inputStr, $count); 
	my ($FFILES, $rxFile, $inputsysOptss, $tmp);
	
	getopt('bfstyorqi',\%OPTS);
	srand (time ^ ($$ + ($$ << 11)));

	die &usage 
		if ((! exists($OPTS{b}) and ! exists($OPTS{m})) or ! exists($OPTS{f}));
	#($bgfFile, $msiFile, $FF, $suffix, $inputStr, $reaxFF, $inputFFType, $sysOptsLabels, $opts, $pqeqFile, $inputsysOptss) = 
	#($OPTS{b},$OPTS{m},$OPTS{f},$OPTS{s},$OPTS{t},$OPTS{r},$OPTS{y},$OPTS{l},$OPTS{o}, $pqeqFile);
	($inputFile, $FF,     $suffix, $inputStr, $inputFFType, $opts,    $qxFile,  $rxFile,   $inputsysOptss) = 
	($OPTS{b},   $OPTS{f},$OPTS{s},$OPTS{t},  $OPTS{y},     $OPTS{o}, $OPTS{q},  $OPTS{r}, $OPTS{i});
	
	print "Step 1: Initializing...";
	@{ $tmp } = split /\s+/,$inputFile;
	for (@{ $tmp }) {
		FileTester($_); 
	}

	$suffix = "lammps" if (! defined($suffix));
	#reax/rexpon option
	if ($FF =~ /reax/i and !defined($rxFile)) {
		print "using default Reax force field...";
		$reaxFF = "$Bin/../ff/Reax.ff";
		#$inputStr = "reax";
	} elsif ($FF =~ /rexpon/i and ($FF =~ /rexpon_wat/i or $FF !~ /rexpon_/i) and !defined($rxFile)) {
		print "using default RexPoN force field...";
		$reaxFF = "${Bin}/../ff/ffield.RexPoN";
		$FF .= " rexpon_wat";
		$opts = "" if (! defined($opts));
		$opts = "pqeq ${opts}" 
			if($opts !~ /pqeq/);
	} elsif ($FF =~ /reax|rexpon/i and $FF !~ /rexpon_/i) {
		$reaxFF = $rxFile;
		FileTester($reaxFF);
	} else {
		$reaxFF = $rxFile;
	}
	$qxFile = "${Bin}/../ff/pqeqRexPoN.par" 
		if ($FF =~ /(rexpon_wat|rexpon |rexpon\$)/i) and ! defined($qxFile);	
	$FF =~ s/REXPON_WAT/rexpon_wat/g;	
	($FFILES, undef) = ReadFFs($FF, $reaxFF);
	for (@{ $FFILES }) {
		$ffType .= $_->{FFTYPEID};
	}
	$opts .= ' polarizable:"fftype=~/OW/" thole rigid: "resname =~ /WAT/"' if ($FF =~ /swm4-ndp/);
	$inputStr = "mesodna" if ($ffType =~ /3/);
	#qeq option
	$QEqtmp->{Opt} = 0;
	$opts = "" if (! defined($opts));
	$sysOpts = (
				{
					Labels          => 1,
					InputCoeffs     => 0,
					TypeCharges     => 0,
					sort_fftypes    => 0,
					compress_parms  => 0,
					amoeba_monopole => 0,
				}
	);
	$opts = "write_inputfile_coeff write_inputfile_type_charge ${opts}" 
		if(defined($inputsysOptss) and $inputsysOptss =~ /1|yes/i);
	if ($opts) {
		if($opts =~ /qeq/i) {
			$QEqtmp->{File} = "$Bin/../ff/qeq.par";
			$QEqtmp->{Opt} = 1;
			if ($opts =~ /pqeq/i) {
				$QEqtmp->{File} = "$Bin/../ff/pqeq2.par";
				$QEqtmp->{Opt} = 2;
			}
			if ($opts =~ /charge (\-?\d+\.*\d*)/) {
				$QEqtmp->{sys_charge} = $1;
			} else {
				$QEqtmp->{sys_charge} = 0;
			}
		}
		$sysOpts->{amoeba_monopole} = $sysOpts->{WriteOffDiagPairs} = 1
			if ($opts =~ /amoeba_monopole/i);
		$sysOpts->{InputCoeffs} = 1
			if ($opts =~ /write_inputfile_coeff/i);
		$sysOpts->{Labels} = 0 
			if ($opts =~ /no_labels/i);
		$sysOpts->{TypeCharges} = 1
			if ($opts =~ /write_inputfile_type_charge/i);
		$sysOpts->{sort_fftypes} = 1
			if ($opts =~ /sort_fftypes/i);
		$sysOpts->{WriteOffDiagPairs} = 1
			if ($opts =~ /write_inputfile_pair_off_diag/i);
		if ($opts =~ /nobonds/i) {
			$tmp = ();
			&getOptSelStr($opts,"nobonds",\@{ $tmp });
			if ($#{ $tmp } > -1) {
				$sysOpts->{nobonds}{aStr} = BuildAtomSelectionString($tmp->[0]);
			} else {
				$sysOpts->{nobonds}{aStr} = BuildAtomSelectionString("index>0");
			}
		}	
		$QEqtmp->{File} = $1 if(defined($qxFile) and $qxFile =~ /(\S+)/);
		#Core-shell
		if ($opts =~ /polar/i) {
			$tmp = ();
			&getOptSelStr($opts,"polar",\@{ $tmp });
			if ($#{ $tmp } > -1) {
				for $i (@{ $tmp }) {
					$sysOpts->{coreShell}->{aStr} .= BuildAtomSelectionString($i) . " or "; 
				}
				$sysOpts->{coreShell}->{aStr} = substr($sysOpts->{coreShell}->{aStr},0,-4);
				$sysOpts->{coreShell}->{type} = 0; #adiabatic core-shell
				$sysOpts->{coreShell}->{type} = 1 if ($opts =~ /drude/i); #drude induced dipole
				$sysOpts->{coreShell}->{type} = 2 if ($opts =~ /thole/i); #drude induced dipole with thole screening
				$sysOpts->{coreShell}->{type} = 3 if ($opts =~ /pqeq/i);  #drude induced dipole with pqeq screeneing
			}
		}
		#use original PQEq method for adiabatic core-shell with QEq
		undef $sysOpts->{coreShell} 
			if (defined($sysOpts->{coreShell}) and $sysOpts->{coreShell}->{type}==0 && defined($QEqtmp) and $QEqtmp->{Opt} == 2);

		#PQEq and FIX_CONP electrode opts
		if ($opts =~ /electrode/i and $opts =~ /(fixed|qeq|conp)/i) {
			$tmp = ();
			&getOptSelStr($opts,"electrode",\@{ $tmp });
			if ($#{ $tmp } > -1 ) {
				$sysOpts->{electrode}->{type} = lc $1;
				$sysOpts->{electrode}->{top}{aStr} = $sysOpts->{electrode}->{bot}{aStr} = BuildAtomSelectionString($tmp->[0]);
				$sysOpts->{electrode}->{bot}{aStr} = BuildAtomSelectionString($tmp->[1]) if ($#{$tmp} > 0 );
				$opts .= " 2d"
			}
		}

		#FEP
		if ($opts =~ /fep/ig) {
			$tmp = ();
			&getOptSelStr($opts,"fep",\@{ $tmp });
			for $i (@{ $tmp }) {
				$sysOpts->{fep} .= BuildAtomSelectionString($i) . " or "; 
			}
			$sysOpts->{fep} = substr($sysOpts->{fep},0,-4) if ($sysOpts->{fep} ne "");
		}

		#Rigid bodies
		if ($opts =~ /rigid/ig) {
			$tmp = ();
			&getOptSelStr($opts,"rigid",\@{ $tmp });
			for $i (@{ $tmp }) {
				$sysOpts->{rigid} .= BuildAtomSelectionString($i) . " or "; 
			}
			$sysOpts->{rigid} = substr($sysOpts->{rigid},0,-4);			
		}

		#include file
		if ($opts =~ /include_file/ig) {
			$tmp = ();
			&getOptSelStr($opts, "include_file",\@{ $tmp });
			$sysOpts->{include_file} = "@{ $tmp }" if ($#{ $tmp } > -1);
		}

		#compress parms
		$sysOpts->{compress_parms} = 1 if ($opts =~ /compress_parms/i);
	}

	$inputStr = "full" if (! defined($inputStr)); 
	$inputStr .= " fep " if (defined($sysOpts->{fep}) and $inputStr !~ / fep/i);
	$inputStr =~ s/ fep//g if ($inputStr =~ / fep/ and ! defined($sysOpts->{fep}));
	$inputList = getInputTypes(0);
	while($inputStr =~ /(\S+)/ig) {
		if(-e $1) {
			$inputType->{names}{custom} = 1;
			push @{ $inputType->{loc} }, $1
		} elsif (exists($inputList->{$1})) {
			$inputType->{names}{$1} = 1;
			push @{ $inputType->{loc} },$inputList->{$1};
		}
	}
	$BONDS->{counter} = $ANGLES->{counter} = $TORSIONS->{counter} = $INVERSIONS->{counter} = $CMAP->{counter} = 0;
	$BI_TORSIONS->{counter} = $PI_TORSIONS->{counter} = 0;
	$totAtms = 0;
	$sysOpts->{mesodna} = 1 if ($ffType =~ /3/);
	$ffType = $inputFFType if (defined($inputFFType));
	print "Done\n";
	return $FFILES;
}

sub usage {
	my ($list) = getInputTypes(1);
	my ($fTypeStr) = &GetFileTypeStr;
	my ($ffTypeStr) = &GetFFTypeStr;
return <<DATA;
This script will generate LAMMPS data and input files from a bgf structure
usage: $0 -b structureFile -f \"ff1 ff2...\" -s [suffix] -t [sim template] -q [qeq/pqeq file] -r [reax/rexpon file] -i [inputFile_coeffs] -o [other_options] 
Arguments:
	-b structureFile: 
$fTypeStr
	-f \"forcefield1 forcefield2...\": 
$ffTypeStr
	-r [reax/rexpon file]: (optional) Specifies either the ReaxFF  or RexPoN forcefield file
	-q [qeq/pqeq/polarization file]: (optional) Specifies either the QEq/PQEq/Drude parameters set.
	-s [suffix]: (optional) When specified, the program will generate in.[suffix]
		and data.[suffix] as the files. If not specified, the output will be
		in.lammps, data.lammps
	-i [write inputfile options]: (optional). If set to 1, then assumes -o "write_inputfile_coeffs write_inputfile_type_charges". Default 0		

	-o [options]: (optional). Controls various options in input file. Valid entries include:
		"amoeba_monopole" - switch from full multipole to monopole only amoeba
		"2D or dimension 2" - for 2D simulations
		"compress_parms" - compress parameters to eliminate multiples with the same values.
		"finite" - for 0D (isolated) simulations
		"ewald vdw" - calculate long range vdw using ewald summations. Only valid for lj and exp6 potentials.
		"no shake|shake solute/solvent" - shake constraints on the hydrogens are turned on by default. 
			This turns it off or applies it to the solute or solvent only
		"no_labels" - Specfies whether to supress writing atom type based label for the coefficients. Default yes
		"nobonds:" 'atom selection' - Delete all bonds involving selected atoms. If no atom selection give, then will delete ALL bonds!!
		"sort_fftypes" - Sort the fftypes alphabetically and renumber accordingly. Default 0: fftypes sorted 
			based on encountered order in structure file
		"write_inputfile_coeffs" - Specifies whether to write the coefficients in the input file, instead of the data file. Default 0
		"write_inputfile_type_charges" - Specifies whether to write the charge of the atom types from the forcefield (if provided)
			in the forcefield. Default 0
		"write_inputfile_pair_off_diag" - Specifies whether to write the off diagonal components to the input file. By default, this is done
		"fep:" 'atom selection'. Write simulation parameters to perform a FEP simulation on atom selection(s) 
		"qeq" - dynamically determine the charge during dynamics using the QEQ charge equilibration scheme.
			See the -x option to specify the parameter set, else the default set will be used
		"pqeq" - dynamically determine the partial atomic charge during dynamics using the PQEq scheme.		
			See the -x option to specify the parameter set, else the default set will be used
			NOTE: The charges are repesented as gaussian distributions (not point charges) and associated with
				the coul/pqeq pair style.
		"charge x" - overall charge on system with pqeq. Default 0		
		"electrode:" 'atom selection_1' ('atom_selection_2')" - for QEq simulations, this invokes the ECHEMDID method, 
			based on the Chi parameter of the specified forcefield. The format is top electrode, bottom electrode. 
			If only is specified then the fftype is for both. Same for conp simulations (see above)
			"conp" - Activate the CONP options for electrochemical cell simulations. Must be use in conjunction with the 'electrode' flag
			"fixedQ" - Activate the fixed charge electrode option. Must be used in conjunction with the 'electrode' flag
		"polarizable:" 'atom selection(s)' [adiabatic/drude/thole/pqeq]" - Turns on shell polarization options. 
			The atom selection(s) should be enclosed in quotes and is based on the usual selection criteria. 
			The shell polarization options are: 
				adiabatic: adiabatic core/shell model, 
				drude: drude induced dipole, 
				thole: drude induced dipole with thole short range screening, 
				pqeq: represent the electrostatics between the atoms using overlap of gaussian types orbitals, 
					as opposed to point dipoles in the other options.
		"rigid: 'atom selection(s)'" - Treat the specified atoms (and their associated molecules) as rigid bodies during dynamics
		"include_file:" - Include a file with LAMMPS code after the data_read line. Use for further customization

	-t [control file template]: (optional). Specifies the type of input file to create. See $Bin/dat/LAMMPS for a list
		Current options include "$list"
		or you can specify your own input file

Report any bugs to tpascal\@ucsd.edu
DATA

}

sub numerically {
	($a<=>$b);
}

sub alphabetically {
	($a cmp $b);
}

#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Math::Amoeba qw(MinimiseND);
use File::Basename qw(basename dirname);
use General qw(ReadFFs FileTester LoadFFs LoadElements FindElement CenterText);
use FileFormats qw(GetBGFFileInfo);
use Getopt::Std qw(getopt);
use ManipAtoms qw(BuildAtomSelectionString SelectAtoms);
use LAMMPS qw(ReadDataFile);
use lmp_lib;
use Cwd qw(cwd abs_path);

sub numerically { ($a<=>$b) };

my ($paramFile, $saveName, $isHybrid, $isPQEQ, $outFFstr, $ELEMENTS);
my ($PARMS, $DATA, $ffData,$count,$DATFILE, $valTypes);
my ($lammpsBinary, $oldParms, $hybridTypes, $typeMap);

$|++;
&init;
print "Parsing paramater file $paramFile...";
&parseParamFile(\%{ $PARMS },$paramFile);
print "Done\n";
&getFileData(\%{ $DATA }, $PARMS);
print "Creating LAMMPS data files...";
&createLammpsFiles($DATA, $PARMS, \%{ $ffData });
&offDiag($ffData, $hybridTypes) if (exists($PARMS->{VDW_OPT}));
print "Done\nCreating objective function...";
&getBounds($ffData,$PARMS);
&updateAtomTypeData(\%{ $ffData });
die "ERROR: Nothing to optimize!\n"
	if(!%{ $ffData } and (! exists($PARMS->{CHECK_ONLY}) or $PARMS->{CHECK_ONLY} != 1));
print "Done\nStarting SIMPLEX optimizer...\n";
&getOptimalParams($PARMS, $saveName);
print "Creating optimized CERIUS2 FF $saveName...";
&createFF($DATA, $PARMS, $ffData, "$saveName");
&cleanup($DATA);
print "Done\n";

sub parseParamFile {
	my ($parms, $infile) = @_;
	my ($key, $val, $rec, $i, $j, $nopt, $valid);

	$nopt = 0;
	open PARMFILE, $infile or die "ERROR: Cannot open $infile: $!\n";
	while(<PARMFILE>) {
		chomp;
		next if ($_ =~ /^#/);
		if ($_ =~ /^(\S+):\s*(.*)$/) {
			($key,$val) = (uc $1,$2);
			if($key !~ /_OPT|VARY|PQEQ_FRAG_CHARGE|COUPLE/ or $key =~ /MOVABLE|LMP_ARGS/) {
				$parms->{$key} = $val;
			} elsif ($key =~ /PQEQ_FRAG_CHARGE/) {
				while($val =~ /(\S+)/g) {
					push @{ $parms->{$key} }, $1;
				}
			} else {
				$rec = \%{ $parms->{$key} };
				while($val =~ /(\S+)/g) {
					$rec = \%{ $rec->{$1} };
				}
				$rec->{VAL} = 1;
				$nopt++;
			}
		}
	}
	close PARMFILE;
	foreach $i (qw/BGF_FILES FF_FILES/) {
		die "ERROR: Expected keyword $i not found in $infile!\n"
			if (!exists($parms->{$i}));
	}

	$valid = 0;
	foreach $i (qw/ENG PRESSURE QMOMENT IP-EA/) {
		$valid =1 if(exists($parms->{"${i}_FILES"}));
	}
	$valid = 1 if ($valid == 0 and exists($parms->{ATOM_CHARGES}));
	die "ERROR: Expected either [ENG/QMOMENT/IP-EA/PRESSURE]_FILES in $infile!\n"
		if (!$valid);
		
	$parms->{MODE} = "ENG";

	$parms->{GROUP1} = BuildAtomSelectionString($parms->{MOVABLE}) if (exists($parms->{MOVABLE}));

	while($parms->{FF_FILES} =~ /(\S+)/g) {
		$parms->{FF_FILES} =~ s/$1// if !(-e $1 and -r $1 and -T $1);
	}
	($parms->{FFs},undef) = ReadFFs($parms->{FF_FILES},undef);

	foreach $j (qw/ENG BGF QMOMENT PRESSURE IP-EA/) {
		next if (! exists($parms->{"${j}_FILES"}));
		@{ $parms->{"${j}s"} } = `ls $parms->{"${j}_FILES"}`;
		if ($#{ $parms->{"${j}s"} } == -1) {
			die "ERROR: No valid $j files found\n";
		}
		for $i (0 .. $#{ $parms->{"${j}s"} }) {
			$parms->{"${j}s"}[$i] =~ s/\n//;
		}
	}
	die "ERROR: Nothing to optimize!\n"
		if(!$nopt);
	die "ERROR: Expected IP/EA FIT since IP-EA_FILES flag is activated!\n"
		if(exists($parms->{"IP-EA_FILES"}) and (! exists($parms->{IP_FIT}) and ! exists($parms->{EA_FIT})));
	$parms->{MODE}= "IP-EA" if (exists($parms->{"IP-EAs"}));		
	if(exists($parms->{QMOMENTs})) {
		$parms->{MODE} = "QMOMENT";
		die "ERROR: Expected dipole/quadrupole for QMOMENT_FIT keyword. Got '$parms->{QMOMENT_FIT}'!\n"
			if(!exists($parms->{QMOMENT_FIT}) or $parms->{QMOMENT_FIT} !~ /(dipo|quad)/i);
		if($parms->{QMOMENT_FIT} =~ /dip/i) {
			$parms->{QMOMENT_FIT} = 1;
		} else {
			$parms->{QMOMENT_FIT} = 2;
		}
	}
	if(exists($parms->{PRESSUREs})) {
		$parms->{MODE} = "PRESSURE";
		die "ERROR: Expected tot/diagonal/full for PRESSURE_FIT keyword. Got '$parms->{PRESSURE_FIT}'!\n"
			if (! exists($parms->{PRESSURE_FIT}) or $parms->{PRESSURE_FIT} !~ /(tot|diag|full)/i);
		if($parms->{PRESSURE_FIT} =~ /tot/i) {
			$parms->{PRESSURE_FIT} = 1;
		} elsif($parms->{PRESSURE_FIT} =~ /diag/i) {
			$parms->{PRESSURE_FIT} = 2;
		} else{
			$parms->{PRESSURE_FIT} = 3;
		}
	}
	if(exists($parms->{ATOM_CHARGES})) {
		$parms->{MODE} = "ATOMQ";
	}
	&determineIfSW($PARMS, $infile);
}

sub determineIfSW {
	my ($parms, $inFile) = @_;
	my ($i, $j, $k, $l, $prefix);

	return if (! exists($parms->{"3BODY_OPT"}));
	$prefix = basename($inFile);
	$prefix =~ s/\.\w+$//;

	ALOOP: for $i (keys %{ $parms->{"3BODY_OPT"} }) {
		for $j (keys %{ $parms->{"3BODY_OPT"}{$i} }) {
			for $k (keys %{ $parms->{"3BODY_OPT"}{$i}{$j} }) {
				for $l (keys %{ $parms->{"3BODY_OPT"}{$i}{$j}{$k} }) {
					if($l =~ /^sw$/i) {
						$parms->{SW} = "pair_coeff * * sw ${prefix}.sw";
						last ALOOP;
					}
				}
			}
		}
	}
}

sub getFileData {
	my ($data, $parms) = @_;

	if($parms->{MODE} eq "ENG") {
		&getEngFileData(\%{ $data }, $parms);
	} elsif ($parms->{MODE} eq "QMOMENT") {
		&getChargeFileData(\%{ $data }, $parms);
	} elsif ($parms->{MODE} eq "PRESSURE") {
		&getPressureFileData(\%{ $data }, $parms);
	} elsif ($parms->{MODE} eq "IP-EA") {
		&getIPEAFileData(\%{ $data }, $parms);
	} elsif ($parms->{MODE} eq "ATOMQ") {
		&setAtomList($parms);
		&setFileData(\%{ $data }, $parms);
	}
}

sub setAtomList {
	my ($parms) = $_[0];
	my ($i, $j, $tmp);

	if($parms->{ATOM_CHARGES} =~ /all/i) {
		$parms->{ATOMQ}{sstr} = "index>0";
		return;
	}

	@{ $tmp } = split /\s+/, $parms->{ATOM_CHARGES};
	for $i (@{ $tmp }) {
		if($i =~ /(\d+)\-(\d+)/) {
			for $j ($1 .. $2) {
				$parms->{ATOMQ}{sstr} .= "index==$j or ";
			}
		} elsif ($i =~ /(\d+)/) {
			$parms->{ATOMQ}{sstr} .= "index==$1 or ";
		}
	}
	$parms->{ATOMQ}{sstr} = BuildAtomSelectionString(substr($parms->{ATOMQ}{sstr},0,-3));
}

sub setFileData {
	my ($data, $parms) = @_;
	my ($i, $j, $index, $tmp);

	$index = 0;
	for $i (@{ $parms->{BGFs}}) {
		$data->{$index}{DISP} = $index;
		($data->{$index}{ATOMS}, undef) = GetBGFFileInfo($i,0,0);
		$tmp = SelectAtoms($parms->{ATOMQ}{sstr}, $data->{$index}{ATOMS});
		@{ $data->{$index}{aList} } = sort numerically keys %{ $tmp }; 
		for $j (0 .. $#{ $data->{$index}{aList}}) {
			$data->{$index}{aList}[$j]--;
		}
		$data->{$index++}{FILE} = $i;
	}
}

sub getEngFileData {
	my ($data, $parms) = @_;
	my ($i, $bgfLoc, $findCmd, $engFile, $dataTmp, $index); 
	my ($rec, $j, $factor, $emin, $offset, $isNB, $tmp, $nprefix, $nsuffix);
	my ($bgfs, $engs) = ($parms->{BGFs}, $parms->{ENGs});

	$isNB = 0;
	$isNB = 1 if (exists($parms->{VDW_OPT}) or exists($parms->{QEQ_OPT}) or exists($parms->{"3BODY_OPT"}));

	$nprefix = 2;
	$nsuffix = 2;
	$nprefix = $parms->{ENG_BGF_PREFIX_NMATCH} if (exists($parms->{ENG_BGF_PREFIX_NMATCH}));
	$nsuffix = $parms->{ENG_BGF_SUFFIX_NMATCH} if (exists($parms->{ENG_BGF_SUFFIX_NMATCH}));
	print "will match the bgf files using $nprefix prefix fields and $nsuffix suffix fields...";
	$index = $offset = 0;
	$factor = 1;
	$factor = $parms->{ENG_FACTOR} if (exists($parms->{ENG_FACTOR}));
	$parms->{ENG_OFFSET} = 0 if (!exists($parms->{ENG_OFFSET}));
	$parms->{DIST_OFFSET} = 0 if (!exists($parms->{DIST_OFFSET}));
	for $engFile (@{ $engs }) {
		print "Parsing Energy function file $engFile...";
		$dataTmp = getObjectiveFunctionData($parms, $engFile, $factor, \${emin});
		&assignBGFS($dataTmp, $bgfs, $engFile, $nprefix, $nsuffix);
		for $j (keys %{ $dataTmp }) {
			$rec = $dataTmp->{$j};
			#$rec->{DISP}=$index;
			$data->{$index} = $rec;
			$data->{$index}{DISP} += $offset;
			$index++;
		}
		@{ $tmp } = sort numerically keys %{ $dataTmp };
		$offset += pop(@{ $tmp }) + 1;
		print "Done\n";
	}
	if(defined($parms->{ENG_OFFSET}) and $isNB) {
		#$parms->{ENG_OFFSET} = 0;
	}
	&printEngVals($data);
}

sub getChargeFileData {
	my ($data, $parms) = @_;
	my ($i, $index, $factor, $chargeFile, $dataTmp, $rec);

	$index = 0;

	for $chargeFile (@{ $parms->{QMOMENTs} }) {
		print "Parsing charge moments data file $chargeFile...";
		$dataTmp = parseChargeFile($chargeFile, $parms->{QMOMENT_FACTOR});
		&assignBGFS($dataTmp, $parms->{BGFs}, $chargeFile);
		for $i (keys %{ $dataTmp }) {
			$rec = $dataTmp->{$i};
			$data->{$index} = $rec;
			$index++;
		}
		print "Done\n";
	}
	&printChargeVals($data, $parms);
}

sub parseChargeFile {
	my ($chargeFile, $parmFactors) = @_;
	my ($DATA, $i, $isValid, $factors, $index, $vals);

	$isValid = 0;
	@{ $factors } = (1, 1, 1); #dipole quadrupole octapole
	$i = 0;
	if (defined($parmFactors)) {
		while($parmFactors =~ /(\d+\.?\d*)/g) {
			$factors->[$i++]=$1;
		}
	}
	open CHRGFILE, $chargeFile or die "ERROR: Cannot open $chargeFile: $!\n";
	while (<CHRGFILE>) {
		chomp;
		if ($_ =~ /^(\d+\.?\d*)\s+(\-?\d+\.?\d*\s*.*)$/) {
			$index = $1;
			$vals = $2;
			$i = 0;
			while ($vals =~ /([X|Y|Z]+)\s*:\s*(\-?\d+\.?\d*)/g) {
				$DATA->{$index}{lc $1} = $2;
				$isValid = 1;
			}
			$DATA->{$index}{WEIGHT}=1;
			$DATA->{$index}{DISP} = $index;
		}
	}
	close CHRGFILE;

	die "ERROR: No valid charge data found while searching $chargeFile!\n"
		if (! $isValid);
	
	return $DATA;	
}

sub getPressureFileData {
	my ($data, $parms) = @_;
	my ($i, $index, $factor, $pressureFile, $dataTmp, $rec);

	$index = 0;

	for $pressureFile (@{ $parms->{PRESSUREs} }) {
		print "Parsing Pressure data file $pressureFile...";
		$dataTmp = parsePressureFile($pressureFile);
		&assignBGFS($dataTmp, $parms->{BGFs}, $pressureFile);
		for $i (keys %{ $dataTmp }) {
			$rec = $dataTmp->{$i};
			$data->{$index} = $rec;
			$index++;
		}
		print "Done\n";
	}
	&printPressureVals($data, $parms);
}

sub parsePressureFile {
	my ($pressureFile) = $_[0];
	my ($DATA, $tmp, $i, $isValid, $factors, $index, $vals);

	$isValid = 0;
	@{ $tmp } = ("tot", "Pxx", "Pyy", "Pzz", "Pxy", "Pyz");

	open PRESSFILE, $pressureFile or die "ERROR: Cannot open $pressureFile: $!\n";
	while (<PRESSFILE>) {
		chomp;
		if ($_ =~ /^(\d+\.?\d*)\s+(\-?\d+\.?\d*\s*.*)$/) {
			$index = $1;
			$vals = $2;
			$i = 0;
			while($vals =~ /(\-?\d+\.?\d*)/g && $i <= $#{ $tmp }) {
				$DATA->{$index}{$tmp->[$i]} = $1;
				$i++;
			}
			$isValid = 1;
			$DATA->{$index}{WEIGHT}=1;
			$DATA->{$index}{DISP} = $index;
		}
	}
	close PRESSFILE;

	die "ERROR: No valid pressure data found while searching $pressureFile!\n"
		if (! $isValid);
	
	return $DATA;	
}

sub getIPEAFileData {
	my ($data, $parms) = @_;
	my ($i, $index, $factor, $ipeaFile, $dataTmp, $rec);

	$index = 0;

	for $ipeaFile (@{ $parms->{"IP-EAs"} }) {
		print "Parsing IP/EA data file $ipeaFile...";
		$dataTmp = parseIPEAFile($ipeaFile, $parms->{"IP-EA_FACTOR"});
		&assignBGFS($dataTmp, $parms->{BGFs}, $ipeaFile);
		for $i (keys %{ $dataTmp }) {
			$rec = $dataTmp->{$i};
			$data->{$index} = $rec;
			$index++;
		}
		print "Done\n";
	}
	&printIPEAVals($data, $parms);
}

sub parseIPEAFile {
	my ($ipeaFile, $parmFactors) = @_;
	my ($DATA, $tmp, $i, $isValid, $factors, $index, $vals);

	$isValid = 0;
	@{ $tmp } = ("IP", "EA");
	@{ $factors } = (1, 1);
	$i = 0;
	if (defined($parmFactors)) {
		while($parmFactors =~ /(\d+\.?\d*)/g) {
			$factors->[$i++]=$1;
		}
	}
	open IPEAFILE, $ipeaFile or die "ERROR: Cannot open $ipeaFile: $!\n";
	while (<IPEAFILE>) {
		chomp;
		if ($_ =~ /^(\d+\.?\d*)\s+(\-?\d+\.?\d*\s*.*)$/) {
			$index = $1;
			$vals = $2;
			$i = 0;
			while($vals =~ /(\-?\d+\.?\d*)/g && $i <= $#{ $tmp }) {
				$DATA->{$index}{$tmp->[$i]} = $1*$factors->[$i];
				$i++;
			}
			$isValid = 1;
			$DATA->{$index}{WEIGHT}=1;
			$DATA->{$index}{DISP} = $index;
			$DATA->{$index}{IP} *= -1;
		}
	}
	close IPEAFILE;

	die "ERROR: No valid IP/EA data found while searching $ipeaFile!\n"
		if (! $isValid);
	
	return $DATA;	
}

sub printEngVals {
	my ($data) = $_[0];
	my ($i, $tmp);

	print "===========================================================\n";
	printf "%-10s %20s %6s %s\n", "INDEX", "ENERGY(kcal/mol)", "WEIGHT", "BGFFILE";
	print "===========================================================\n";
	@{ $tmp } = sort { $data->{$a}{DISP} <=> $data->{$b}{DISP} }keys %{ $data };
	for $i (@{ $tmp }) {
		printf "%-10s %20.5f %6.3f %s\n",$data->{$i}{DISP},$data->{$i}{ENG},$data->{$i}{WEIGHT},$data->{$i}{FILE};
	}
}

sub printChargeVals {
	my ($data, $parms) = @_;
	my ($i, $tmp);

	print "===========================================================\n";
	if($parms->{QMOMENT_FIT} == 1) {
		printf "%-10s %20s %6s %s\n", "INDEX", "DIPOLE MOMENT", "WEIGHT", "BGFFILE";
	} elsif($parms->{QMOMENT_FIT} == 2) {
		printf "%-10s %21s %6s %s\n", "INDEX", "Traceless", "WEIGHT", "BGFFILE";
		printf "%-10s %21s %6s %s\n", "", "QUADRUOPLE MOMENT", "", "";
		printf "%-10s %5s/%5s/%5s/%5s/%5s/%5s %6s %s\n", "", "Qxx", "Qxy", "Qxz", "Qyy", "Qyz", "Qzz", "", "";
	}
	print "===========================================================\n";
	@{ $tmp } = sort { $data->{$a}{DISP} <=> $data->{$b}{DISP} }keys %{ $data };
	for $i (@{ $tmp }) {
		if($parms->{QMOMENT_FIT} == 1) {
			printf "%-10s %20.5f %6.3f %s\n",$data->{$i}{DISP},$data->{$i}{'DIPOLE'},$data->{$i}{WEIGHT},$data->{$i}{FILE};
		} elsif($parms->{QMOMENT_FIT} == 2) {
			printf "%-10s %5.3f/%5.3f/%5.3f/%5.3f/%5.3f/%5.3f %6.3f %s\n",
				$data->{$i}{DISP},$data->{$i}{'Qxx'},$data->{$i}{'Qxy'},$data->{$i}{'Qxz'},
				$data->{$i}{'Qyy'},$data->{$i}{'Qyz'},$data->{$i}{'Qzz'},$data->{$i}{WEIGHT},$data->{$i}{FILE};
		}
	}

}

sub printPressureVals {
	my ($data, $parms) = @_;
	my ($i, $tmp);

	print "=======================================================================================\n";
	if($parms->{PRESSURE_FIT} == 1) {
		printf "%-10s %20s %6s %s\n", "INDEX", "Ptot", "WEIGHT", "BGFFILE";
	} elsif($parms->{PRESSURE_FIT} == 2) {
		printf "%-10s %20s %20s %20s %6s %s\n", "INDEX", "Pxx", "Pyy", "Pzz", "WEIGHT", "BGFFILE";
	} else {
		printf "%-10s %20s %20s %20s %20s %20s %20s %6s %s\n", "INDEX", "Pxx", "Pyy", "Pzz", "Pxy", "Pxz", "Pyz", "WEIGHT", "BGFFILE";
	}
	print "=======================================================================================\n";
	@{ $tmp } = sort { $data->{$a}{DISP} <=> $data->{$b}{DISP} }keys %{ $data };
	for $i (@{ $tmp }) {
		if($parms->{PRESSURE_FIT} == 1) {
			printf "%-10s %20.5f %6.3f %s\n",$data->{$i}{DISP},$data->{$i}{Ptot},$data->{$i}{WEIGHT},$data->{$i}{FILE};
		} elsif($parms->{PRESSURE_FIT} == 2) {
			printf "%-10s %20.5f %20.5f %20.5f %6.3f %s\n",$data->{$i}{DISP},$data->{$i}{Pxx},$data->{$i}{Pyy},$data->{$i}{Pzz},$data->{$i}{WEIGHT},$data->{$i}{FILE};
		} else {
			printf "%-10s %20.5f %20.5f %20.5f %20.5f %20.5f %20.5f %6.3f %s\n",$data->{$i}{DISP},$data->{$i}{Pxx},$data->{$i}{Pyy},$data->{$i}{Pzz},$data->{$i}{Pxy},$data->{$i}{Pxz},$data->{$i}{Pyz},$data->{$i}{WEIGHT},$data->{$i}{FILE};
		}
	}

}
sub printIPEAVals {
	my ($data, $parms) = @_;
	my ($i, $tmp);

	print "===========================================================\n";
	printf "%-10s ", "INDEX";
	if($parms->{IP_FIT} == 1) {
		printf "%20s ", "IP(kcal/mol)";
	} 
	if($parms->{EA_FIT} == 1) {
		printf "%20s ", "EA(kcal/mol)";
	}
	printf "%6s %s\n","WEIGHT", "BGFFILE";
	print "===========================================================\n";
	@{ $tmp } = sort { $data->{$a}{DISP} <=> $data->{$b}{DISP} }keys %{ $data };
	for $i (@{ $tmp }) {
		printf "%-10s ", $data->{$i}{DISP};
		if($parms->{IP_FIT} == 1) {
			printf "%20.5f ",$data->{$i}{'IP'};
		}
		if($parms->{EA_FIT} == 1) {
			printf "%20.5f ",$data->{$i}{'EA'};
		}
		printf "%6.3f %s\n",$data->{$i}{WEIGHT},$data->{$i}{FILE};
	}

}

sub offDiag {
	my ($parms) = $_[0];
	my ($i, $j, @vals, $aa, $bb, $t, $curr);

	for $i (keys %{ $parms->{VDWS} }) {
		next if (! exists($parms->{VDWS}{$i}{$i}));
		for $j (keys %{ $parms->{VDWS} }) {
			next if (! exists($parms->{VDWS}{$j}{$j}));
			next if ($i eq $j || $i lt $j);
			for $t (keys %{ $parms->{VDWS}{$i}{$i} }) {
				next if (!exists($parms->{VDWS}{$j}{$j}{$t}));
				next if ((exists($parms->{VDWS}{$i}{$j}) and exists($parms->{VDWS}{$i}{$j}{$t})) or 
						 (exists($parms->{VDWS}{$j}{$i}) and exists($parms->{VDWS}{$j}{$i}{$t})));
				$aa = $parms->{VDWS}{$i}{$i}{$t};
				$bb = $parms->{VDWS}{$j}{$j}{$t};
				@vals = ();
				if ($t =~ /morse/i) {
					$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
					$vals[1] = sqrt($aa->{VALS}[1]*$bb->{VALS}[1]);
					$vals[2] = 0.5*($aa->{VALS}[2]+$bb->{VALS}[2]);
				} elsif ($t =~ /charmm|gromacs|lj\/cut/) {
					$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
					$vals[1] = sqrt($aa->{VALS}[1]*$bb->{VALS}[1]);
					if ($#{ $aa->{VALS} } > 1) {
						$vals[2] = sqrt($aa->{VALS}[2]*$bb->{VALS}[2]);
						$vals[3] = sqrt($aa->{VALS}[3]*$bb->{VALS}[3]);
					}
				} elsif ($t =~ /buck/) {
					$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
					$vals[1] = 0.5*($aa->{VALS}[1]+$bb->{VALS}[1]);
					$vals[2] = sqrt($aa->{VALS}[2]*$bb->{VALS}[2]);
				}
				next if (! @vals);
				@{ $parms->{VDWS}{$i}{$j}{$t}{VALS} } = @vals;
			}
		}
	}
}

sub getNewPairCoeffs {
	my ($data, $ff_type) = @_;
	my ($i, $j, $diag, $offDiag, $curr);

	for $i (keys %{ $data->{VDW} }) {
		for $j (keys %{ $data->{VDW}{$i} }) {
			$curr = $data->{VDW}{$i}{$j}{1};
			if ($ff_type eq "MPSIM") {
				$diag .= sprintf("%-4s   %3d", $i, $curr->{NUM}) if ($i eq $j);
				$offDiag .= sprintf("%-4s -%-4s   %3d", $i, $j, $curr->{NUM}) if ($i ne $j);
				for (@{ $curr->{VALS} }) {
					$offDiag .= sprintf("%10.5f",$_) if ($i ne $j);
					$diag .= sprintf("%10.5f",$_)  if ($i eq $j);
				}
				$diag .= "\n" if ($i eq $j);
				$offDiag .= "\n" if ($i ne $j);
			} elsif ($ff_type eq "CERIUS2") {
				$diag .= sprintf(" %-12s%-13s",$i,$curr->{TYPE}) if ($i eq $j);
				$offDiag .= sprintf(" %-9s%-9s%-13s",$i,$j,$curr->{TYPE}) if ($i ne $j);
				for (@{ $curr->{VALS} }) {
					$diag .= sprintf("%8.4g",$_) if ($i eq $j);
					$offDiag .= sprintf("%8.4g",$_)  if ($i ne $j);
				}
				$diag .= "\n" if ($i eq $j);
				$offDiag .= "\n" if ($i ne $j);
			}
		}
	}

	return ($diag, $offDiag);
}

sub createFF {
	my ($data, $parms, $ff, $file) = @_;
	my ($i, $j, $k, $l, $m, $n, $offset, $tlist, $ffMod, $c2Name, $c2Header); 
	my ($fftype, $curr, $headerStr, $header, $lmpName, $tmp);
	my ($iType, $jType, $ffFile, $idx, $val, $scale, $tstr);

	$file = basename($file);
	$fftype = "CERIUS2";
	if(exists($parms->{VDW_OPT}) and $parms->{VDW_OPT}) {
		if($fftype eq "CERIUS2") {
			for $i (keys %{ $parms->{VDW_OPT} }) {
				$iType = getTypeIndex($i);
				for $j (keys %{ $parms->{VDW_OPT}{$i} }) {
					$jType = getTypeIndex($j);
					($iType,$jType) = ($jType,$iType) if ($iType > $jType);
					if($i eq $j) {
						$ffMod->{DIAGONAL_VDW}{N} = 1;
						$curr = \%{ $ffMod->{DIAGONAL_VDW}{DATA}{$i} };
					} else {
						$ffMod->{OFF_DIAGONAL_VDW}{N} = 2;
						$curr = \%{ $ffMod->{OFF_DIAGONAL_VDW}{DATA}{$i}{$j} };
					}
					$curr->{VALS} = "";
					for $k (keys %{ $parms->{VDW_OPT}{$i}{$j} }) {
						next if ($k eq "coul/pqeqgauss");
						($lmpName, $c2Name) = getTypeName("VDW",$k);
						die "ERROR: Cannot find C2 Name for $k!\n" if (! defined($c2Name));
						$curr->{VALS} = sprintf(" %-8s%-8s %10s    ",$i,$j,$c2Name) if ($i ne $j);
						$curr->{VALS} = sprintf(" %-8s %10s    ",$i,$c2Name) if ($i eq $j);
						for $l (0 .. $#{ $ff->{VDWS}{$iType}{$jType}{$k}{VALS} }) {
							$m = $typeMap->{VDW}{$lmpName}{c2INDX}[$l];
							next if ($m == -1); #dependence check
							$idx = $ff->{VDWS}{$iType}{$jType}{$k}{VALSIDXPTR}[$m];
							if(defined($idx)) {
								$scale = eval($typeMap->{VDW}{$lmpName}{c2FACT}[$l]);
								$val = $parms->{OPT}{final}[$idx]*$scale;
							} else {
								$scale = eval($typeMap->{VDW}{$lmpName}{c2FACU}[$l]);
								$val = $ff->{VDWS}{$iType}{$jType}{$k}{VALS}[$l]*$scale;
							}
							$curr->{VALS} .= sprintf(" %23.15E",$val);
						}
						$curr->{VALS} .= "\n";
						$curr->{VALS} =~ s/E\+00/   /g;
						$curr->{USED} = 0;
					}
				}
			}
		}
	}
	for $i (qw/3BODY BOND ANGLE TORSION INVERSION/) {
		next if (! exists($ff->{"${i}S"}));
		($header, $n) = getC2Header("${i}S");
		if($fftype eq "CERIUS2") {
			for $j (keys %{ $ff->{"${i}S"} }) {
				for $k (keys %{ $ff->{"${i}S"}{$j} }) {
					($lmpName, $c2Name, $c2Header) = getTypeName($i,$k);
					$c2Header = $header 
						if(!defined($c2Header));
					$ffMod->{$header}{N} = $n;
					$curr = \%{ $ffMod->{$c2Header}{DATA} };
					$tstr = " ";
					$tmp = ();
					&flattenHash($ff->{"${i}S"}{$j}{$k}{KEY}, \@{ $tmp });
					for $l (@{ $tmp }) {
						$curr = \%{ $curr->{$l} };
						$tstr .= sprintf("%-7s ",$l);
					}
					$tstr .= sprintf("%14s",$c2Name);
					$curr->{VALS} = $tstr;
					for $l (0 .. $#{ $ff->{"${i}S"}{$j}{$k}{VALS} }) {
						$m = $typeMap->{$i}{$k}{c2INDX}[$l];
						next if($m == -1); #dependence check
						$idx = $ff->{"${i}S"}{$j}{$k}{VALSIDXPTR}[$m];
						if(defined($idx)) {
							$scale = eval($typeMap->{$i}{$k}{c2FACT}[$l]);
							$val = $parms->{OPT}{final}[$idx]*$scale;
						} else {
							$scale = eval($typeMap->{$i}{$k}{c2FACU}[$l]);
							$val = $ff->{"${i}S"}{$j}{$k}{VALS}[$m]*$scale;
						}
						$curr->{VALS} .= sprintf(" %23.15E",$val);
					}
					$curr->{VALS} .= "\n";
					$curr->{VALS} =~ s/E\+00/   /g;
					$curr->{USED} = 0;
				}
			}
		}
	}

	&writePQEQFile($ff, $parms->{OPT}, $file) if($isPQEQ);

	undef $header;
	$curr = ();
	@{ $curr } = keys %{ $ffMod };
	$headerStr = "@{ $curr }";
	$headerStr =~ s/\s+/|/g;
	@{ $ffFile } = split /\s+/,$parms->{FF_FILES};
	open FF, "$ffFile->[0]" or die "ERROR: Cannot access $ffFile->[0]: $!\n";
	while(<FF>) {
		chomp;
		if($_ =~ /^($headerStr)/) {
			if(defined($header)) {
				&getUnUsedParms($ffMod->{$header});
			}
			$header = $1;
			$i = 0;
			$curr = \%{ $ffMod->{$header}{DATA} };
			$outFFstr .= "$_\n";
		} elsif($_ =~ /^\s*END/) {
			if(defined($header)) {
				&getUnUsedParms($ffMod->{$header});
			}
			undef $header;
			$outFFstr .= "$_\n";
		} elsif(defined($header)) {
			$i = 0;
			$curr = \%{ $ffMod->{$header}{DATA} };
			while($_ =~ / (\S+) /g) {
				last if ($i == $ffMod->{$header}{N});
				if (!exists($curr->{$1})) {
					undef $curr;
					last;
				}
				$i++;
				$curr = \%{ $curr->{$1} };
			}
			if(defined($curr) and exists($curr->{VALS})) {
				$outFFstr .= $curr->{VALS};
				delete $curr->{VALS};
			} else {
				$outFFstr .= "$_\n";
			}
		} else {
			$outFFstr .= "$_\n";
		}
	}
	close FF;

	open OUTFF, ">$file" or die "ERROR: Cannot create $file: $!\n";
	print OUTFF $outFFstr;
	close OUTFF;
	system("$Bin/pruneCerius2FF.pl -b $data->{0}{FILE} -f $file -s $file > /dev/null");

	chdir "_lmptmp";
	for $i (keys %{ $data }) {
		&saveCoords($data->{$i}{FULL}{LMP},"opt") if ($isPQEQ);
	}
	chdir "../";
}

sub writePQEQFile {
	my ($ff, $opt, $fName) = @_;
	my ($i, $k, $l, $m, $ele, $shouldWrite, $iType, $val);

	$shouldWrite = 0;
	for $i (keys %{ $ff->{VDWS} }) {
		for $k (keys %{ $ff->{VDWS}{$i}{$i} }) {
			next if ($k ne "coul/pqeqgauss");
			if(exists($ff->{VDWS}{$i}{$i}{$k}{USED})) {
				$shouldWrite = 1;
				last;
			}
		}
	}

	return if (!$shouldWrite);

	$fName =~ s/\.\w+$//;
	$fName .= ".pqeq.dat";

	open PQEQFILE, "> $fName" or die "ERROR: Cannot write too $fName: $!\n";
	print PQEQFILE <<DATA; 
##########################################################################################
#E    P          Xo        Jo          Z          Rc         Rs       Ks           Ks1
##########################################################################################
DATA

	for $i (keys %{ $ff->{VDWS} }) {
		@{ $iType } = keys %{ $ff->{ATOMTYPES}{$i}{KEY} };
		$ele = getElement($iType->[0]);
		for $k (keys %{ $ff->{VDWS}{$i}{$i} }) {
			next if ($k ne "coul/pqeqgauss");
			printf PQEQFILE "%2s ",$ele;
			for $l (0 .. $#{ $ff->{VDWS}{$i}{$i}{$k}{VALS} }) {
				$m = $typeMap->{VDW}{"coul/pqeqgauss"}{c2INDX}[$l];
				next if ($m == -1); #dependence check
				$val = $ff->{VDWS}{$i}{$i}{$k}{VALS}[$m];
				if(exists($ff->{VDWS}{$i}{$i}{$k}{VALSIDXPTR}) and defined($ff->{VDWS}{$i}{$i}{$k}{VALSIDXPTR}[$m])) {
					$val = $ff->{VDWS}{$i}{$i}{$k}{VALSIDXPTR}[$m];
					$val = $opt->{final}[$val];
				}
				printf PQEQFILE "%10.5f ", $val;
			}
			printf PQEQFILE "\n";
		}
	}
	close PQEQFILE;
}

sub getTypeName {
	my ($type, $sstr) = @_;
	my ($i, $j, $c2Name, $lmpName, $c2Header);

	for $i (keys %{ $typeMap->{$type} }) {
		if($sstr =~ /$i/) {
			$c2Name = $typeMap->{$type}{$i}{c2Name};
			$c2Header = $typeMap->{$type}{$i}{c2Head}
				if(exists($typeMap->{$type}{$i}{c2Head}));
			$lmpName = $i;
			last;
		}
	}

	return ($lmpName, $c2Name, $c2Header);
}

sub getUnUsedParms {
	my ($ff) = $_[0];
	my ($i);

	return if (ref($ff) ne 'HASH');
	for $i (keys %{ $ff }) {
		if(ref($ff->{$i}) eq 'HASH' and exists($ff->{$i}{VALS}) and $ff->{$i}{VALS}) {
			$outFFstr .= $ff->{$i}{VALS};
			delete $ff->{$i}{VALS};
		} else {
			&getUnUsedParms($ff->{$i});
		}
	}
}

sub writeQEqParms {
	my ($parms) = $_[0];
	my ($outStr, $i, $val);

	$outStr .= "QEq\n";
	for $i (keys %{ $parms } ) {
		$val = $parms->{$i}{1}{VALS};
		$outStr .= sprintf(" %-5s %8.5f %8.5f %8.5f %8.5f %8.5f\n", $i, @{ $val });
	}
	return $outStr;
}

sub updateAtomTypeData {
	my ($ff) = $_[0];
	my ($i, $j, $k);

	for $i (keys %{ $ff->{VDWS} }) {
		for $j (keys %{ $ff->{VDWS}{$i} }) {
			for $k (keys %{ $ff->{VDWS}{$i}{$j} }) {
				delete $ff->{VDWS}{$i}{$j}{$k} if (!exists($ff->{VDWS}{$i}{$j}{$k}{USED}) and $k ne "coul/pqeqgauss");
			}
			delete $ff->{VDWS}{$i}{$j} if (!%{ $ff->{VDWS}{$i}{$j} });
		}
		delete $ff->{VDWS}{$i} if (!%{ $ff->{VDWS}{$i} });
	}
	delete $ff->{VDWS} if (!%{ $ff->{VDWS} });

	for $i ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS") {
		next if (!exists($ff->{$i}));
		for $j (keys %{ $ff->{$i} }) {
			for $k (keys %{ $ff->{$i}{$j} }) {
				delete $ff->{$i}{$j}{$k} if (!exists($ff->{$i}{$j}{$k}{USED}));
			}
			delete $ff->{$i}{$j} if (! %{ $ff->{$i}{$j} });
		}
		delete $ff->{$i} if (!%{ $ff->{$i} });
	}
}

sub updateAtomTypeData_old {
	my ($data, $ffData) = @_;
	my ($i, @tmp, $datfile, $start, $j, $ffTypeList);
	
	$ffTypeList = ();

	@tmp = keys %{ $data };
	$datfile = $data->{ $tmp[0] }{FULL}{LMPDAT};
	$start = 0;
	open DAT, $datfile or die "ERROR: Cannot open $datfile: $!\n";
	while (<DAT>) {
		chomp;
		$start = 1 if ($_ =~ /^Masses/);
		last if ($_ =~ /Coeff/);
		if ($start and $_ =~ /(\d+)\s+\d+\.\d+\s+\#\s+(\S+)/) {
			$ffTypeList->{$2} = $1;
		}
	}
	close DAT;
	die "ERROR: Cannot parse LAMMPS data file $datfile!\n" if (! defined($ffTypeList));

	for $i (keys %{ $ffData->{ATOMTYPES} }) {
		if (! exists($ffTypeList->{$i})) {
			delete $ffData->{ATOMTYPES}{$i};
			delete $ffData->{VDW}{$i};
			for $j (keys %{ $ffData->{VDW} }) {
				delete $ffData->{VDW}{$j}{$i} if exists($ffData->{VDW}{$j}{$i});
			}
		} else {
			$ffData->{ATOMTYPES}{$i}{INDEX} = $ffTypeList->{$i};
		}
	}
	for $i (keys %{ $ffData }) {
		delete $ffData->{$i} if ($i !~ /ATOMTYPES|VDW/);
	}
		
}	

sub isPairHybrid {
	my ($parms, $inFile) = @_;

	open INFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
	while(<INFILE>) {
		chomp;
		if ($_ =~ /^pair_coeff.* sw /) {
			$parms->{SW} = $_;
			$parms->{SW} =~ s/ sw (\S+) / sw $parms->{SW_file} /;
			last;
		}
	}
	close INFILE;

	return $isHybrid;
}

sub parseSWfile {
	my ($swParms, $swFile) = @_;
	my ($tmp, $i);

	$i = 1;
	open SWFILE, $swFile or die "ERROR: Cannot open $swFile: $!\n";
	while (<SWFILE>) {
		chomp;
		@{ $tmp } = split /\s+/, $_;
		if (scalar(@{ $tmp }) < 14) {
			die "ERROR: Invalid SW file $swFile!\n";
		}
		${swParms}->{$i}{sw}{KEY}{shift @{ $tmp } }{shift @{ $tmp } }{shift @{ $tmp } } = ();
		@{ ${swParms}->{$i}{sw}{VALS} } = @{ $tmp };
		$i++;
	}
	close SWFILE;
}

sub createLammpsFiles {
	my ($data, $parms, $ff) = @_;
	my ($curr, $fflist, $ffs, $bgfFile, $saveName, $currCmd, $currDir); 
	my ($splitBGF, $tmp, $c, $swfile, $isNB, $pcl, $chgLine);

	$currDir = cwd();
	$isNB = 0;
	$isNB = 1 if (exists($parms->{VDW_OPT}) or exists($parms->{QEQ_OPT}) or exists($parms->{"3BODY_OPT"}));

	$splitBGF = "$Bin/getBGFAtoms.pl";
	$fflist = $parms->{FFs};
	system("mkdir -p _lmptmp");
	for $curr (@{ $fflist }) {
		$saveName = abs_path($curr->{FF});
		$ffs .= "$saveName ";
	}

	$c = 0;
	for $curr (keys %{ $data }) {
		$c++;
		$bgfFile = abs_path($data->{$curr}{FILE});
		$bgfFile = "../$data->{$curr}{FILE}" if (! defined($bgfFile));
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$//;
		chdir "_lmptmp";
		#first create lammps input and data file for full system
		system("rm -fr retry.out");
		&createSys($bgfFile, "$ffs", $saveName, $parms);
		if(exists($parms->{SW}) and $c == 1) {
			$swfile = $parms->{SW};
			$swfile =~ s/^.* //;
			$parms->{SW_file} = $swfile;
			&isPairHybrid($parms, "in.${saveName}_singlepoint");
			&parseSWfile(\%{ $ff->{"3BODYS"} }, "${saveName}.sw");
		}
		$data->{$curr}{FULL}{LMPDAT} = "_lmptmp/data.${saveName}";
		$data->{$curr}{FULL}{LMPTRJ} = abs_path("${saveName}.init.lammps");
#		$data->{$curr}{FULL}{LMP} = init lmp_lib "-log none -nc";
		$data->{$curr}{FULL}{LMP} = init lmp_lib "-screen none -log none -nc";
		$data->{$curr}{FULL}{LMP}->file("./in.${saveName}_singlepoint");
		if($parms->{MODE} eq "QMOMENT") {
			$tmp = <<'DATA';
compute cc1 all chunk/atom molecule
compute mydChunk all dipole/chunk cc1
compute myqChunk all quadrupole/chunk cc1
print "compute initiated"
#dipole
variable qx equal sum(c_mydChunk[1])
variable qy equal sum(c_mydChunk[2])
variable qz equal sum(c_mydChunk[3])
variable qt equal sum(c_mydChunk[4])
#quadrupole
variable Qxx equal sum(c_myqChunk[1])
variable Qxy equal sum(c_myqChunk[2])
variable Qxz equal sum(c_myqChunk[3])
variable Qyx equal sum(c_myqChunk[4])
variable Qyy equal sum(c_myqChunk[5])
variable Qyz equal sum(c_myqChunk[6])
variable Qzx equal sum(c_myqChunk[7])
variable Qzy equal sum(c_myqChunk[8])
variable Qzz equal sum(c_myqChunk[9])
variable Qxx_diag equal sum(c_myqChunk[10]) #diagonalized Qxx
variable Qyy_diag equal sum(c_myqChunk[11]) #diagonalized Qyy
variable Qzz_diag equal sum(c_myqChunk[12]) #diagonalized Qzz
#variable Qxx_yy equal v_qxx_diag-v_qyy_diag
#variable Q2zz_xx_yy equal 2*v_qzz_diag-v_qxx_diag-v_qyy_diag
variable Qxx_yy equal v_Qxx-v_qyy
variable Q2zz_xx_yy equal 2*v_qzz-v_Qxx-v_qyy
fix 1 all ave/time 1 1 1 c_mydChunk[*] c_myqChunk[*] file tmp.out mode vector
run 10
DATA

			$data->{$curr}{FULL}{LMP}->commands_string($tmp);
			if($parms->{QMOMENT_FIT} == 1) {
				$data->{$curr}{INIT} = sprintf("%20.5f ",$data->{$curr}{FULL}{LMP}->extract_variable("qt",0));
			} else {
				$data->{$curr}{INIT} = sprintf("%9.5f/%-9.5f/%9.5f/%-9.5f/%9.5f/%-9.5f ",
										$data->{$curr}{FULL}{LMP}->extract_variable("Qxx",0),
										$data->{$curr}{FULL}{LMP}->extract_variable("Qxy",0),
										$data->{$curr}{FULL}{LMP}->extract_variable("Qxz",0),
										$data->{$curr}{FULL}{LMP}->extract_variable("Qyy",0),
										$data->{$curr}{FULL}{LMP}->extract_variable("Qyz",0),
										$data->{$curr}{FULL}{LMP}->extract_variable("Qzz",0)); 
				}
		} elsif($parms->{MODE} eq "IP-EA") {
			$data->{$curr}{FULL}{LMP}->run(10);
			$data->{$curr}{FULL}{initENG} = $data->{$curr}{INIT} = getLMPeng($data->{$curr}{FULL},());
			#$data->{$curr}{INIT} = $data->{$curr}{FULL}{initENG};
			#IP
			if(exists($parms->{"IP_FIT"})) {
				$data->{$curr}{LMPIP}{LMP} = init lmp_lib "-screen none -log none -nc";
				$data->{$curr}{LMPIP}{LMP}->file("./in.${saveName}_singlepoint");
				$data->{$curr}{LMPIP}{LMP}->command("fix     pqeq all pqeq method 2 nevery 1 charge 1.0 tolerance 1.0e-6 damp 1.0");
				$data->{$curr}{LMPIP}{LMP}->run(10);
				$data->{$curr}{LMPIP}{initENG} = $data->{$curr}{INIT} - getLMPeng($data->{$curr}{LMPIP},());
			}	
			#EA
			if(exists($parms->{"EA_FIT"})) {
				$data->{$curr}{LMPEA}{LMP} = init lmp_lib "-screen none -log none -nc";
				$data->{$curr}{LMPEA}{LMP}->file("./in.${saveName}_singlepoint");
				$data->{$curr}{LMPEA}{LMP}->command("fix     pqeq all pqeq method 2 nevery 1 charge -1.0 tolerance 1.0e-6 damp 1.0");
				$data->{$curr}{LMPEA}{LMP}->run(10);
				$data->{$curr}{LMPEA}{initENG} = $data->{$curr}{INIT} - getLMPeng($data->{$curr}{LMPEA},());
			}
		} elsif($parms->{MODE} eq "ATOMQ") {
			@{ $data->{$curr}{REFQ} } = map { $data->{$curr}{ATOMS}{$_}{CHARGE} } sort numerically keys %{ $data->{$curr}{ATOMS} }; 
			@{ $data->{$curr}{FFTYPE} } = map { $data->{$curr}{ATOMS}{$_}{FFTYPE} } sort numerically keys %{ $data->{$curr}{ATOMS} }; 
			$data->{$curr}{INITQ} = getLMPAtomVar($data->{$curr}{FULL},undef,"q");
			#pop @{ $data->{$curr}{INITQ} };
		} elsif($parms->{MODE} eq "PRESSURE") {
			$tmp = <<'DATA';
thermo_style custom etotal pxx pyy pzz pxy pxz pyz press
run 0
DATA
 			$data->{$curr}{FULL}{LMP}->commands_string($tmp);
		} else {
			$data->{$curr}{FULL}{LMP}->run(10);
			$data->{$curr}{FULL}{initENG} = $data->{$curr}{INIT} = getLMPeng($data->{$curr}{FULL},());
		}
		$isHybrid = `grep -c 'pair_style\\s*hybrid' in.${saveName}_singlepoint`;
		$isPQEQ = (system("grep '^fix.*pqeq' in.${saveName}_singlepoint > /dev/null")) ? 0: 1;
		&saveCoords($data->{$curr}{FULL}{LMP},"init") if ($isPQEQ);

		&getFFParms($ff, $saveName);

		if (!$isNB) {
			chdir "../";
			next;
		}
		next if ($parms->{MODE} eq "QMOMENT" or !exists($parms->{MOVABLE}));

		#now isolated fragment 1
		$data->{$curr}{FRAG1}{LMP} = init lmp_lib "-screen none -log none -nc";
		$data->{$curr}{FRAG1}{LMP}->file("./in.${saveName}_singlepoint");
		$data->{$curr}{FRAG1}{LMP}->command("group frag1 $parms->{MOVABLE}");
		$data->{$curr}{FRAG1}{LMP}->command("delete_atoms group frag1");
		if(defined($parms->{PQEQ_FRAG_CHARGE}) and $#{$parms->{PQEQ_FRAG_CHARGE}} > 0) {
			if ($c == 1) {
				$pcl=getPQEQchargeLine("./in.${saveName}_singlepoint",$parms->{PQEQ_FRAG_CHARGE}[0]);
				push @{ $chgLine }, $pcl;
			}
			$data->{$curr}{FRAG1}{LMP}->command("$chgLine->[0]");
		}
		$data->{$curr}{FRAG1}{LMP}->run(10);
		$data->{$curr}{FRAG1}{initENG} = getLMPeng($data->{$curr}{FRAG1},());

		#now isolated fragment 2
		$data->{$curr}{FRAG2}{LMP} = init lmp_lib "-screen none -log none -nc";
		$data->{$curr}{FRAG2}{LMP}->file("./in.${saveName}_singlepoint");
		$data->{$curr}{FRAG2}{LMP}->command("group frag1 $parms->{MOVABLE}");
		$data->{$curr}{FRAG2}{LMP}->command("group frag2 subtract all frag1");
		$data->{$curr}{FRAG2}{LMP}->command("delete_atoms group frag2");
		if(defined($parms->{PQEQ_FRAG_CHARGE}) and $#{$parms->{PQEQ_FRAG_CHARGE}} > 0) {
			if ($c == 1) {
				$pcl=getPQEQchargeLine("./in.${saveName}_singlepoint",$parms->{PQEQ_FRAG_CHARGE}[1]);
				push @{ $chgLine }, $pcl;
			}
			$data->{$curr}{FRAG2}{LMP}->command("$chgLine->[1]");
		}
		$data->{$curr}{FRAG2}{LMP}->run(10);
		$data->{$curr}{FRAG2}{initENG} = getLMPeng($data->{$curr}{FRAG2},());

		$data->{$curr}{INIT} = $data->{$curr}{FULL}{initENG} - $data->{$curr}{FRAG1}{initENG} - $data->{$curr}{FRAG2}{initENG};

		chdir "../";
	}
	chdir $currDir;
}

sub getPQEQchargeLine {
	my ($inFile, $fcharge) = @_;
	my ($chgLine);

	open INFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
	while(<INFILE>) {
		chomp;
		if($_ =~ /^fix.*pqeq\s*method/) {
			$chgLine = $_;
			$chgLine =~ s/charge\s+(\S+)/charge $fcharge/;
			last;
		}
	}
	close INFILE;
	die "ERROR: Cannot find fix pqeq line in $inFile!"
		if (! defined($chgLine));
	return $chgLine;
}

sub saveCoords {
	my ($struct, $suf) = @_;

	$struct->command("undump 1");
	$struct->command("dump 1 all custom 1 \${sname}.${suf}.lammps id type xu yu zu q");
	$struct->command("dump_modify 1 first yes");
	$struct->run(0);
	#$struct->command("undump 1");
}

sub getFFParms {
	my ($ff, $fle) = @_;
	my ($datFile, $inpFile, $i, $headerStr, $header, $tmp, $key, $alist, $curr, $vals);

	$datFile = "data.${fle}";
	$inpFile = "in.${fle}_singlepoint";

	for $i ($datFile, $inpFile) {
		die "ERROR: Cannot access $i: $!!\n"
			if(! -e $i or ! -r $i or ! -T $i);
	}

	#parse input file
	$headerStr = "pair|bond|angle|torsion|inversion";
	open INPFILE, $inpFile or die "ERROR: Cannot access $inpFile: $!\n";
	while (<INPFILE>) {
		chomp;
		if ($_ =~ /^\s*($headerStr)_style\s+(.*)$/) {
			$header = uc $1;
			$tmp = $2;
			$header = "ATOMTYPE" if($header =~ /masses/i);
			$header = "VDW" if($header =~ /pair/i);
			$hybridTypes->{$header} = 1 if ($_ =~ /hybrid/);
			if (!exists($hybridTypes->{$header})) {
				@{ $vals } = split /\s+/, $tmp;
				$valTypes->{$header} = shift @{ $vals };
			}
		} 

		if ($_ =~ /^\s*($headerStr)_coeff\s+(\d+)\s+(.*)/i) {
			$header = uc $1;
			($key, $tmp) = ($2, $3);
			$header = "ATOMTYPE" if($header =~ /masses/i);
			$header = "VDW" if($header =~ /pair/i);
			($tmp, $alist) = split /\#/,$tmp;
			$curr = \%{ $ff->{"${header}S"}{$key} };
			@{ $vals } = split /\s+/,$tmp;
			$curr = \%{ $curr->{shift @{ $vals }} } if ($header eq "VDW");
			if ($hybridTypes->{$header}) {
				$curr = \%{ $curr->{shift @{ $vals }} };
			} else {
				$curr = \%{ $curr->{$valTypes->{$header}} };
			}
			@{ $curr->{VALS} } = @{ $vals };
			if(defined($alist)) {
				$alist =~ s/^\s*//;
				$alist =~ s/\s*$//;
				$tmp = \%{ $curr->{KEY} };
				while($alist =~ /(\S+)/g) {
					$tmp = \%{ $tmp->{$1} };
				}
			}
		}
	}
	close INPFILE;

	#parse datafile and get forcefield parameters
	undef $header;
	$headerStr = "masses|pair coeffs|bond coeffs|angle coeffs|torsion coeffs|inversion coeffs";
	open DATFILE, $datFile or die "ERROR: Cannot acccess $datFile: $!\n";
	while(<DATFILE>) {
		chomp;
		last if ($_ =~ /^Atoms/);
		if ($_ =~ /^($headerStr)/i) {
			$header = uc $1;
			$header =~ s/ coeffs//i;
			$header = "ATOMTYPE" if($header =~ /masses/i);
			$header = "VDW" if($header =~ /pair/i);
		} elsif(defined($header) and $_ =~ /^\s*(\d+)\s+(.*)$/) {
			($key, $tmp) = ($1,$2);
			($tmp, $alist) = split /\#/,$tmp;
			if (!defined($alist) and $header eq "VDW") {
				@{ $vals } = keys %{ $ff->{ATOMTYPES}{$key}{KEY} };
				$alist = "$vals->[0] $vals->[0]";
			}
			if($header eq "ATOMTYPE") {
				$curr = \%{ $ff->{"${header}S"}{$key} };
			} elsif ($header eq "VDW") {
				$curr = \%{ $ff->{"${header}S"}{$key}{$key}{ $valTypes->{$header} } };
			} else {
				$curr = \%{ $ff->{"${header}S"}{$key}{ $valTypes->{$header} } };
			}
			@{ $curr->{VALS} } = split /\s+/,$tmp;
			$alist =~ s/^\s*//;
			$alist =~ s/\s*$//;
			$tmp = \%{ $curr->{KEY} };
			while($alist =~ /(\S+)/g) {
				$tmp = \%{ $tmp->{$1} };
			}
		}
	}
	close DATFILE;

}

sub createSys {
	my ($bgfFile, $ffs, $saveName, $parms) = @_;
	my ($currCmd, $execStr);

	$execStr = "$Bin/createLammpsInput.pl ";
	$currCmd = "${execStr} -f \"$ffs\" -b $bgfFile -s $saveName $parms->{LMP_ARGS} > /dev/null";
	#print "$currCmd\n";
	if (system($currCmd)) {
		print "\nWARNING: Cannot execute '$currCmd'\nWill retry...\n";
		sleep 1;
		$currCmd = "${execStr} -f \"$ffs\" -b $bgfFile -s $saveName $parms->{LMP_ARGS}";
		system($currCmd);
		#die "\n" if(system($currCmd));
	}
}

sub createLammpsFiles_old {
	my ($data, $parms, $lmpOpts) = @_;
	my ($curr, $fflist, $ffs, $execStr, $bgfFile, $saveName, $tmp, $lmp_data, $new_data);

	$fflist = $parms->{FFs};
	system("mkdir -p _lmptmp");
	$execStr = "$Bin/createLammpsInput.pl ";
	for $curr (@{ $fflist }) {
		$saveName = abs_path($curr->{FF});
		$ffs .= "$saveName ";
	}

	@{ $tmp } = keys %{ $data };
	$curr = $tmp->[0];
	$bgfFile = abs_path($data->{$curr}{FILE});
	$saveName = "_tmp";
	chdir "_lmptmp";
	if (system("${execStr} -f \"$ffs\" -b $bgfFile -s $saveName $parms->{LMP_ARGS}> /dev/null")) {
		die "ERROR: Cannnot execute \"${execStr} -f \"$ffs\" -b $bgfFile\"\n";
	}
	chdir "../";
	$lmp_data = parseDataFile("_lmptmp/data._tmp", $parms);
	for $curr (keys %{ $data }) {
		$bgfFile = abs_path($data->{$curr}{FILE});
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$//;
		chdir "_lmptmp";
		$new_data = updateAtomCoord($lmp_data,$bgfFile);
		open LMPDATA, "> data.${saveName}" or die "ERROR: Cannot create data.${saveName}: $!\n";
		print LMPDATA $new_data;
		close LMPDATA;
		$data->{$curr}{LMPDAT} = "_lmptmp/data.${saveName}";
		chdir "../";
	}
}

sub updateAtomCoord {
	my ($lmp_data, $bgfFile) = @_;
	my ($i, $atoms, $outStr, $atomCoord, $tmp, $j);

	($atoms, undef, undef) = GetBGFFileInfo($bgfFile, 0, 0);
	$j = 0;
	foreach $i (@{ $lmp_data->{atomdata} }) {
		$i =~ s/^\s+//;
		@{ $tmp } = split /\s+/, $i;
		if ($#{ $tmp } < 7) {
			$atomCoord .= "$i\n";
			next;
		}
		$j++;
		die "ERROR: Inconsistent #atoms in BGF and data file!\n" if (!exists($atoms->{$j}));
		$atomCoord .= sprintf("%10d%10d%10d%10.5f%10.5f%10.5f%10.5f%10d%10d%10d\n",
			$tmp->[0],$tmp->[1],$tmp->[2],
			$atoms->{$j}{CHARGE},$atoms->{$j}{XCOORD},$atoms->{$j}{YCOORD},$atoms->{$j}{ZCOORD},
			$tmp->[7],$tmp->[8],$tmp->[9]);
	}
	$outStr  = "";
	for $i (@{ $lmp_data->{header} }) {
		$outStr .= "$i\n";
	}
	$outStr .= $atomCoord;
	for $i (@{ $lmp_data->{footer} }) {
		$outStr .= "$i\n";
	}
	return $outStr;
}

sub parseDataFile {
	my ($data_file, $parms) = @_;
	my ($lmp_data, $header, $masses);

	$header = "header";
	$masses = 0;
	open LMPDATA, $data_file or die "ERROR: Cannot open $data_file: $!\n";
	while(<LMPDATA>) {
		chomp;
		$header = "atomdata" if ($header eq "header" and $_ =~ /^Atoms/);
		$header = "footer" if ($header eq "atomdata" and $_ =~ /^Bonds/);
		push @{ $lmp_data->{$header} }, $_;
		if ($_ =~ /^(\S+)/) {
			$masses = 0;
			$masses = 1 if ($1 eq "Masses");
		}
		if ($masses and $_ =~ /^\s+(\d+)\s+\d+\.\d+\s+\#\s+(\S+)/) {
			$parms->{TYPE_MAP}{$2}=$1;
		}
	}
	close LMPDATA;

	return $lmp_data;
}


sub getOptimalParams {
	my ($parms, $outfile) = @_;
	my ($outprefix, $tol, $iter, $opt, $error, $i, $tmp, $fftype, $q, $cmds, $fname);

	$tol = 5e-07;
	$iter = 250;
	$iter = $parms->{ITER} if (defined($parms->{ITER}));

	$outprefix = $outfile;
	$outprefix =~ s/\.\w+$//;
	printf "%-10s %12s Operation\n","Iteration","Residual";
	open $DATFILE, "> ${outprefix}.residuals.dat" or die "ERROR: Cannot write to ${outprefix}.residuals.dat: $!\n";
	@{ $oldParms } = @{ $parms->{OPT}{pGuess} };
	$count = 0;
	if(exists($parms->{CHECK_ONLY}) and $parms->{CHECK_ONLY} == 1) {
		@{ $opt } = @{ $parms->{OPT}{pGuess} };
		$error = 0;
		print "Skipping...";
	} else {
		#### SIMPLEX DOWNHILL ALGORITHM #####
		($opt, $error) = MinimiseND($parms->{OPT}{pGuess}, $parms->{OPT}{pScale}, \&calcResiduals, $tol, $iter);
	}

	@{ $parms->{OPT}{final} } = @{ $parms->{OPT}{pCurr} };
	$cmds = updateCoeffs($ffData, \@{ $parms->{OPT}{final} });
	&updateQEqParms($parms, \@_, \@{ $cmds }) if(exists($PARMS->{QEQ_OPT}));
	&updateSWparms($parms, $ffData, \@{ $opt}) if (exists($PARMS->{SW}));
	push @{ $cmds }, "run 10";
	&calcError($parms, $DATA, $cmds);	

	print $DATFILE "#ERROR: $error\n";
	close $DATFILE;
	open OUTFILE, "> ${outprefix}.summary.dat" or die "ERROR: Cannot write to ${outprefix}.summary.dat: $!\n";
	print OUTFILE "\n=========================\n";
	print OUTFILE "SIMPLEX Optimizer Summary\n";
	print OUTFILE "=========================\n";
	printf OUTFILE "%-60s %10s %10s\n","Parm", "Initial","Final";
	for $i (1 .. 82) { print OUTFILE "-"; }
	printf OUTFILE "\n";
	for $i (0 .. $#{ $opt }) {
		printf OUTFILE "%-60s %10.5f %10.5f\n", $parms->{OPT}{key}[$i], $parms->{OPT}{init}[$i], $parms->{OPT}{final}[$i]; 
	}
	close OUTFILE;
	&printFile2Screen("${outprefix}.summary.dat");

	print "\n#NOTE: Energy Offset: $parms->{ENG_OFFSET}\n" if ($parms->{MODE} eq "ENG" and exists($parms->{ENG_OFFSET}));
	open OUTFILE, "> ${outprefix}.compare.dat" or die "ERROR: Cannot write to ${outprefix}.compare.dat: $!\n";
	print OUTFILE "\n================================\n";
	print OUTFILE "STRUCTURE INFORMATION Comparison\n";
	print OUTFILE "================================\n";
	if($parms->{MODE} ne "ATOMQ" and ($parms->{MODE} ne "QMOMENT" or $parms->{QMOMENT_FIT} == 1)) {
		printf OUTFILE "%-40s %10s %20s %20s %20s %20s %20s\n","Structure","ID#","Weight","Initial","Final","Target","||ERROR(%)||";
	} elsif ($parms->{MODE} eq "QMOMENT") {
		printf OUTFILE "%-40s %10s %20s %36s %36s %36s %20s\n","Structure","ID#","Weight","Initial","Final","Target","||ERROR(%)||";
	} else {
		printf OUTFILE "%-40s %10s %20s %20s %20s %20s %20s\n","Structure","ATOM#","","Initial","Final","Target","||ERROR(%)||";
	}
					
	for $i (1 .. 156) { print OUTFILE "-"; }
	printf OUTFILE "\n";
	for $i (sort {$DATA->{$a}{DISP} <=> $DATA->{$b}{DISP} } keys %{ $DATA }) {
		$fname = basename($DATA->{$i}{FILE});
		$fname =~ s/\.\w+$//;
		if ($parms->{MODE} eq "ENG") {
			$error = 100*($DATA->{$i}{CURR}*$DATA->{$i}{WEIGHT});
			$error = 100*($DATA->{$i}{CURR}-($DATA->{$i}{ENG}-$parms->{ENG_OFFSET}))*$DATA->{$i}{WEIGHT}/($DATA->{$i}{ENG}-$parms->{ENG_OFFSET}) 
				if ($DATA->{$i}{ENG} != 0);
			printf OUTFILE "%-40s %10.3f %20.5G %20.5G %20.5G %20.5G %20.5f\n", $fname, $DATA->{$i}{DISP}, $DATA->{$i}{WEIGHT},
															$DATA->{$i}{INIT}, $DATA->{$i}{CURR}, $DATA->{$i}{ENG}-$parms->{ENG_OFFSET}, $error;
		} elsif ($parms->{MODE} eq "QMOMENT") {
			if($parms->{QMOMENT_FIT} == 1) {
				$error = 100*($DATA->{$i}{'DIPOLE'}-$DATA->{$i}{QMOMENT}{'DIPOLE'})**2;
				printf OUTFILE "%-40s %10.3f %20.5G %20s %20.5G  %20.5G %20.5f\n", $fname, $DATA->{$i}{DISP}, $DATA->{$i}{WEIGHT},
															$DATA->{$i}{INIT}, $DATA->{$i}{QMOMENT}{'DIPOLE'}, $DATA->{$i}{'DIPOLE'}, $error;
			} elsif($parms->{QMOMENT_FIT} == 2) {
				$error = 100*(  ($DATA->{$i}{'Qxx'}-$DATA->{$i}{QMOMENT}{'Qxx'})**2 +
								($DATA->{$i}{'Qxy'}-$DATA->{$i}{QMOMENT}{'Qxy'})**2 +
								($DATA->{$i}{'Qxz'}-$DATA->{$i}{QMOMENT}{'Qxz'})**2 +
								($DATA->{$i}{'Qyy'}-$DATA->{$i}{QMOMENT}{'Qyy'})**2 +
								($DATA->{$i}{'Qyz'}-$DATA->{$i}{QMOMENT}{'Qyz'})**2 +
								($DATA->{$i}{'Qzz'}-$DATA->{$i}{QMOMENT}{'Qzz'})**2);
				printf OUTFILE "%-40s %10.3f %20.5G %20s   %5.3G/%-5.3G/%5.3G/%-5.3G/%5.3G/%-5.3G %5.3G/%-5.3G/%5.3G/%-5.3G/%5.3G/%-5.3G %5.3G/%-5.3G/%5.3G/%-5.3G/%5.3G/%-5.3G %20.5f\n", 
								$fname, $DATA->{$i}{DISP}, $DATA->{$i}{WEIGHT}, 
								$DATA->{$i}{INIT}, $DATA->{$i}{QMOMENT}{'Qxx'}, $DATA->{$i}{QMOMENT}{'Qyy'}, 
															$DATA->{$i}{'Qxx'}, $DATA->{$i}{'Qyy'}, $error;
			}
		} elsif ($parms->{MODE} eq "IP-EA") {
		} elsif ($parms->{MODE} eq "ATOMQ") {
			$tmp = $fname;
			for $q (@{ $DATA->{$i}{aList} }) {
				$fftype = $DATA->{$i}{FFTYPE}[$q];
				$error = 100*(($DATA->{$i}{CURR}[$q] - $DATA->{$i}{REFQ}[$q])/$DATA->{$i}{REFQ}[$q])**2;
				printf OUTFILE "%-40s %10d %20s %20.5G %20.5G %20.5G %20.5f\n", $fname, ($q+1), $fftype, 
				$DATA->{$i}{INITQ}[$q], $DATA->{$i}{CURR}[$q], $DATA->{$i}{REFQ}[$q], $error;
			}
		}
	}
	close OUTFILE;
	&printFile2Screen("${outprefix}.compare.dat");

	open OUTFILE, "> ${outprefix}.compare.plt" or die "ERROR: Cannot write to ${outprefix}.compare.plt: $!\n";
	print OUTFILE "pl '${outprefix}.compare.dat' u 2:4 w lp pt 6 ps 1.5 t 'inital'," .
				  " '' u 2:5 w lp pt 8 ps 1.5 t 'optimized', '' u 2:6 w lp pt 4 ps 1.5 t 'target', 0 w l dt 2 not\n";
	close OUTFILE;			  
	print "\n";
	die "\n" if(exists($parms->{CHECK_ONLY}) and $parms->{CHECK_ONLY} = 1);

}

sub calcResiduals {
	my ($i, $cmds, $tot, $changeList); 

	$count++;

	@{ $PARMS->{OPT}{pCurr} } = @_;
	&setCoupledVals($PARMS) if (exists($PARMS->{OPT}{couple}));
	$cmds = updateCoeffs($ffData, \@{ $PARMS->{OPT}{pCurr} });
	&updateQEqParms($PARMS, $ffData, \@{ $PARMS->{OPT}{pCurr} }) if(exists($PARMS->{QEQ_OPT}));
	&updateSWparms($PARMS, $ffData, \@{ $PARMS->{OPT}{pCurr} }) if (exists($PARMS->{SW}));
	$PARMS->{ENG_OFFSET} = $_[0] if (exists($PARMS->{VARY}{ENG_OFFSET}));
	push @{ $cmds }, $PARMS->{SW} if(exists($PARMS->{SW}));
	push @{ $cmds }, "run 10";

	$tot = calcError($PARMS, $DATA, $cmds);	
	
	$changeList = compareParms($oldParms, $PARMS);
	printf $DATFILE "%-20d %20.5G\n",$count, $tot;
	for $i (0 .. $#{ $changeList }) {
		if($i == 0) {
			printf "%-10d %12.5G %-50s\n", $count, $tot, $changeList->[$i];
		} else {
			printf "%-10s %12s %-50s\n", "", "", $changeList->[$i];
		}
	}
	@{ $oldParms } =@{ $PARMS->{OPT}{pCurr} };
	return $tot;
}

sub setCoupledVals {
	my ($parms) = $_[0];
	my ($i);

	for $i (@{ $parms->{OPT}{couple}}) {
		eval($i);
	}
	print "";
}
sub calcError {
	my ($parms, $data, $cmds) = @_;
	my ($tot, $j, $eng, $emin, $i, $q, $mode, $tmp, $isNB);

	$isNB = 0;
	$isNB = 1 if (exists($parms->{VDW_OPT}) or exists($parms->{QEQ_OPT}) or exists($parms->{"3BODY_OPT"}));

	$tot = 0;
	$mode = $parms->{MODE};
	for $j (keys %{ $parms }) {
		#first get energies and compare
		next if ($j !~ /(vdw|qeq|bond|angle|torsion|inversion|3body)_OPT/i);
		$j = uc $1;
		$eng = 0;
		undef $emin;
		for $i (keys %{ $data }) {
			$eng = getLMPeng($data->{$i}{FULL}, $cmds);
			if($mode eq "ENG") {
				if($j =~ /VDW|QEQ|3BODY/ and exists($parms->{MOVABLE})) {
					$eng -= getLMPeng($data->{$i}{FRAG1}, $cmds);
					$eng -= getLMPeng($data->{$i}{FRAG2}, $cmds);
				}
				$emin = $eng if (!defined($emin) or $eng<$emin);
				$data->{$i}{CURR} = $eng;
			} elsif ($mode eq "QMOMENT") {
				if ($parms->{QMOMENT_FIT} == 1) { 
					#dipole
					$data->{$i}{QMOMENT}{"DIPOLE"}=$data->{$i}{FULL}{LMP}->extract_variable("qt",0);
					$tot += ($data->{$i}{"DIPOLE"}-$data->{$i}{QMOMENT}{"DIPOLE"})**2;
				} else { 
					#quadrupole
					$data->{$i}{QMOMENT}{"Qxx"}=$data->{$i}{FULL}{LMP}->extract_variable("Qxx",0);
					$data->{$i}{QMOMENT}{"Qxy"}=$data->{$i}{FULL}{LMP}->extract_variable("Qxy",0);
					$data->{$i}{QMOMENT}{"Qxz"}=$data->{$i}{FULL}{LMP}->extract_variable("Qxz",0);
					$data->{$i}{QMOMENT}{"Qyy"}=$data->{$i}{FULL}{LMP}->extract_variable("Qyy",0);
					$data->{$i}{QMOMENT}{"Qyz"}=$data->{$i}{FULL}{LMP}->extract_variable("Qyz",0);
					$data->{$i}{QMOMENT}{"Qzz"}=$data->{$i}{FULL}{LMP}->extract_variable("Qzz",0);
					$tmp = (($data->{$i}{"Qxx"}-$data->{$i}{QMOMENT}{"Qxx"})**2 +
							($data->{$i}{"Qxy"}-$data->{$i}{QMOMENT}{"Qxy"})**2 +
							($data->{$i}{"Qxz"}-$data->{$i}{QMOMENT}{"Qxz"})**2 +
							($data->{$i}{"Qyy"}-$data->{$i}{QMOMENT}{"Qyy"})**2 +
							($data->{$i}{"Qyz"}-$data->{$i}{QMOMENT}{"Qyz"})**2 +
							($data->{$i}{"Qzz"}-$data->{$i}{QMOMENT}{"Qzz"})**2);
					$tot += $tmp;
				}
			} elsif ($mode eq "PRESSURE") {
				if ($parms->{PRESSURE_FIT} == 1) { #total pressure
					$data->{$i}{PRESSURE}{TOT}=$data->{$i}{FULL}{LMP}->get_thermo("press");;
					$tot += ($data->{$i}{PRESSURE}-$data->{$i}{PRESSURE}{TOT})**2;
				} elsif($parms->{PRESSURE_FIT} == 2) { #diagonal stress tensor
					for $tmp ("Pxx", "Pyy", "Pzz") {
						$data->{$i}{PRESSURE}{$tmp}=$data->{$i}{FULL}{LMP}->get_thermo(lcfirst $tmp);
						$tot += ($data->{$i}{$tmp}-$data->{$i}{PRESSURE}{$tmp})**2;
					}
				} elsif($parms->{PRESSURE_FIT} == 3) { #full stress tensor
					for $tmp ("Pxx", "Pyy", "Pzz", "Pxz", "Pyz", "Pxy") {
						$data->{$i}{PRESSURE}{$tmp}=$data->{$i}{FULL}{LMP}->get_thermo(lcfirst $tmp);
						$tot += ($data->{$i}{$tmp}-$data->{$i}{PRESSURE}{$tmp})**2;
					}
				}
			} elsif ($mode eq "IP-EA") {
				if(exists($parms->{IP_FIT})) {
					$tmp = $eng - getLMPeng($data->{$i}{LMPIP}, $cmds);
					$tot += ($data->{$i}{IP}-$tmp)**2;
				}
				if(exists($parms->{EA_FIT})) {
					$tmp = $eng - getLMPeng($data->{$i}{LMPEA}, $cmds);
					$tot += ($data->{$i}{EA}-$tmp)**2;
				}
			} elsif ($mode eq "ATOMQ") {
				$data->{$i}{CURR} = getLMPAtomVar($data->{$i}{FULL},$cmds,"q");
				#pop @{ $data->{$i}{CURR} };
				for $q (@{ $data->{$i}{aList} }) {
					$tot += ($data->{$i}{REFQ}[$q] - $data->{$i}{CURR}[$q])**2;
				}
			}
		}
		for $i (keys %{ $data }) {
			if($mode eq "ENG") {
				if(!exists($parms->{VDW_OPT}) and ! $isNB) {
					$data->{$i}{CURR} -= $emin;
				}
				$tot += $data->{$i}{WEIGHT}*($data->{$i}{ENG} - $data->{$i}{CURR} - $parms->{ENG_OFFSET})**2;
			}
		}
	}
	return $tot;
}

sub getLMPeng {
	my ($lmp, $cmds) = @_;
	my ($i, $val, $ret);

	for $i (@{ $cmds }) {
		$ret = $lmp->{LMP}->command($i);
	}
	$val = $lmp->{LMP}->get_thermo("etotal");
	return $val;
}

sub getLMPAtomVar_old {
	my ($lmp, $cmds, $var) = @_;
	my ($i, $atomVars, $ret);

	for $i (@{ $cmds }) {
		$ret = $lmp->{LMP}->command($i);
	}
	$atomVars = $lmp->{LMP}->extract_atom("$var");
	return $atomVars;
}

sub getLMPAtomVar {
	my ($lmp, $cmds, $var) = @_;
	my ($i, $atomVars, $ret, $tmp);

	for $i (@{ $cmds }) {
		$ret = $lmp->{LMP}->command($i);
	}
	$lmp->{LMP}->command("run 10");
	@{ $tmp } = `gawk -v f="$var" -f ${Bin}/getLastLammpsTrjAtomVals.awk $lmp->{LMPTRJ}`;
	for $i (@{ $tmp }) {
		$i =~ s/\n$//;
		push @{ $atomVars }, $i;
	}
	#$atomVars = $lmp->{LMP}->extract_atom("q");
	return $atomVars;
}

sub calcResiduals_old {
	my ($execStr, $i, $pair_coeff, $inp, $valid, $tot, $operation);

	$count++;
	$valid = $tot = 0;

	if(exists($PARMS->{VDW_OPT})) {
		$execStr = "$lammpsBinary -var datafile ";
		$pair_coeff = updatePairCoeffs($PARMS, \@_);
	}
	if(exists($PARMS->{QEQ_OPT})) {
		&updateQEqParms($PARMS, \@_);
	}

	for $i (keys %{ $DATA }) {
		$inp = "_lmptmp/in.tmp";
		`cp $DATA->{$i}{FULL}{LMPINP} $inp`;
		open DAT, ">> $inp" or die "ERROR: Cannot write to $inp: $!\n";
		print DAT $pair_coeff;
		close DAT;
		$execStr = "$lammpsBinary -in $inp |";
		open LMP, "${execStr} |" or die "Cannot execute '${execStr}': $!\n";
		$valid = 0;
		while (<LMP>) {
			chomp;
			if ($_ =~ /interact\s+=\s+(\-?\d+\.\d+)/i) {
				$DATA->{$i}{INIT} = $1 if(!exists($DATA->{$i}{INIT}));
				$DATA->{$i}{CURR} = $1;
				$tot += sprintf("%.5f", ($DATA->{$i}{WEIGHT} * ($1 - $DATA->{$i}{ENG})**2));
				$valid = 1;
			} elsif ($_ =~ /forces\s+=\w+(\d+\.\d+)/) {
				$DATA->{$i}{INIT} = $1 if(!exists($DATA->{$i}{INIT}));
				$DATA->{$i}{CURR} = $1;
				$tot += sprintf("%.5f", ($DATA->{$i}{WEIGHT} * $1));
				$valid = 1;
			}
		}
		close LMP;
		die "ERROR: No valid data found file executing '$execStr'\n" if (! $valid);
	}
	$operation = compareParms($oldParms, \@_, $PARMS);
	print $DATFILE "$count $tot\n";
	printf "%-10d %12.5G $operation\n",$count,$tot;
	@{ $oldParms } = @_;
	return $tot;
}

sub updateQEqParms {
	my ($parms, $newVals, $cmds) = @_;
	my ($i, $curr);

	for $i (0 .. $#{ $parms->{OPT}{map} }) {
		next if (! defined($newVals->[$i]));
		next if ($parms->{OPT}{map}[$i]{otype} ne "qeq");
		$curr = $parms->{OPT}{map}[$i]{name}{1};
		$curr->{VALS}[$parms->{OPT}{map}[$i]{pOrder}] = $newVals->[$i];
	}
	&writeQEqFile($parms->{QEq},"_lmptmp/tmp.qeq.dat");
	push @{ $cmds }, "unfix charge";
	push @{ $cmds }, "fix charge all qeq/rg 1 8 1.0e-6 100 _lmptmp/tmp.qeq.dat";
}

sub writeQEqFile {
	my ($parms, $outfile) = @_;
	my ($tmp, $i);

	for $i (keys %{ $parms }) {
		$tmp->{$parms->{$i}{typeID}} = $parms->{$i};
	}

	open QEqOutFile, "> $outfile" or die "ERROR: Cannot create $outfile: $!\n";
	for $i (sort numerically keys %{ $tmp }) {
		print QEqOutFile "$i @{ $tmp->{$i}->{1}{VALS} }\n";
	}
	close QEqOutFile;
}

sub updateSWparms {
	my ($parms, $ff, $newParms) = @_;
	my ($i, $j, $idx, $curr, $tmp);

	for $i (keys %{ $ff->{"3BODYS"} }) {
		next if (!exists($ff->{"3BODYS"}{$i}{sw}));
		$curr = $ff->{"3BODYS"}{$i}{sw};
		next if(!exists($curr->{USED}));
		for $j (0 .. $#{ $curr->{VALS} }) {
			if(defined($curr->{VALSIDXPTR}[$j])) {
				$idx = $curr->{VALSIDXPTR}[$j];
				$curr->{VALS}[$j] = $newParms->[$idx];
				$curr->{VALS}[$j] = 0
					if($curr->{VALS}[$j] < 0);
			}
		}
	}
	open SWFILE, "> $parms->{SW_file}" or die "ERROR: Cannot create $parms->{SW_file}: $!\n";
	for $i (keys %{ $ff->{"3BODYS"} }) {
		next if (!exists($ff->{"3BODYS"}{$i}{sw}));
		$tmp = ();
		$curr = $ff->{"3BODYS"}{$i}{sw};
		&flattenHash($curr->{KEY}, \@{ $tmp });
		printf SWFILE "@{ $tmp } ";
		for $j (0 .. $#{ $curr->{VALS} }) {
			printf SWFILE "%s ", $curr->{VALS}[$j];
		}
		printf SWFILE "\n";
	}
	close SWFILE;
}

sub updateCoeffs {
	my ($ff, $newParms) = @_;
	my ($i, $j, $k, $l, $curr, $hybridName, $cmds, $cmdStr, $val, $idx);

	for $i (keys %{ $ff->{VDWS} }) {
		for $j (keys %{ $ff->{VDWS}{$i} }) {
			for $k (keys %{ $ff->{VDWS}{$i}{$j} }) {
				$curr = $ff->{VDWS}{$i}{$j}{$k};
				next if(!exists($curr->{USED}));
				$hybridName = "";
				$hybridName = $k if (exists($hybridTypes->{"VDW"}));
				$cmdStr = "pair_coeff $j $i $hybridName ";
				$cmdStr = "pair_coeff $i $j $hybridName " if($i<$j);
				#first we update the curr array with the newParms array
				for $l (0 .. $#{ $curr->{VALS} }) {
					if(defined($curr->{VALSIDXPTR}[$l])) {
						$idx = $curr->{VALSIDXPTR}[$l];
						$curr->{VALS}[$l] = $newParms->[$idx];
					}
				}
				#here we update the independent variable if we are optimizing dependent ones
				&updateIndpVars($curr) if (exists($curr->{DEPENDENCE}));
				#finally we populate the cmdStr
				for $l (0 .. $#{ $curr->{VALS} }) {
					$val = $curr->{VALS}[$l];
					$cmdStr .= "$val " 
						if(defined($val) and (!exists($curr->{DEPENDENCE}) or ! defined($curr->{DEPENDENCE}[$l])));
				}
				push @{ $cmds }, $cmdStr;
			}
		}
	}

	for $i ("BOND", "ANGLE", "DIHEDRAL", "TORSION") {
		next if (! exists($ff->{"${i}S"}));
		for $j (keys %{ $ff->{"${i}S"} }) {
			for $k (keys %{ $ff->{"${i}S"}{$j} }) {
				$curr = $ff->{"${i}S"}{$j}{$k};
				$hybridName = "";
				$hybridName = $k if (exists($hybridTypes->{$i}));
				$cmdStr = lc $i . "_coeff $j $hybridName ";
				for $l (0 .. $#{ $curr->{VALS} }) {
					$val = $curr->{VALS}[$l];
					if(defined($curr->{VALSIDXPTR}[$l])) {
						$idx = $curr->{VALSIDXPTR}[$l];
						$val = $newParms->[$idx];
					}
					$cmdStr .= "$val ";
				}
				push @{ $cmds }, $cmdStr;
			}
		}
	}

	return $cmds;
}

sub updateIndpVars {
	my ($curr, $newParm) = @_;
	my ($i);	

	for $i (0 .. $#{ $curr->{DEPENDENCE} }) {
		next if (! defined($curr->{DEPENDENCE}[$i]));
		eval($curr->{DEPENDENCE}[$i]);
	}
}

sub updatePairCoeffs {
	my ($ffData, $newParms) = @_;
	my ($lmpCmds, $i, $j, @vdws, @tmp, $curr); 
	my ($pair, $vdwType, $type1, $type2, $hybridStr);
	
	for $i (0 .. $#{ $ffData->{OPT}{map} }) {
		next if (! defined($newParms->[$i]));
		next if ($ffData->{OPT}{map}[$i]{otype} ne "vdw");
		$curr = $ffData->{OPT}{map}[$i]{name}{1}; # j points to $PARMS{VDW}{type1}{type2}
		$curr->{VALS}[$ffData->{OPT}{map}[$i]{pOrder}] = $newParms->[$i]; # update j val to reflect new val from simplex
		$vdwType = uc $curr->{TYPE};
		($type1, $type2) = split /\-/,$ffData->{OPT}{map}[$i]{label};
		$type1 = $ffData->{TYPE_MAP}{$type1};
		$type2 = $ffData->{TYPE_MAP}{$type2};
		@vdws = @{ $curr->{VALS} };
		if ($vdwType =~ /^MORSE|STRETCH_MORSE$/) {
			($vdws[0], $vdws[1]) = ($vdws[1], $vdws[0]);
			($vdws[1], $vdws[2]) = ($vdws[2], $vdws[1]);
			$vdws[1] /= ($vdws[2] * 2); # changed from div to multi 07/28/2007
			$vdws[3] /= ($vdws[2] * 2) if ($vdwType =~ /^STRETCH_MORSE/); #aplha2 for stretch morse
		} elsif ($vdwType =~ /^LJ_6_12$/) {
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			$vdws[1] = $vdws[1] / (2**(1/6));
			if ($#vdws > 1) {
				($vdws[2],$vdws[3]) = ($vdws[3],$vdws[2]);
				$vdws[3] = $vdws[3] /(2**(1/6));
			}
		} elsif ($vdwType =~ /^EXPO_6$/) {
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			@tmp = @vdws;
			@vdws = ();
			$vdws[0] = $tmp[0] * (6/($tmp[2]-6)) * exp($tmp[2]);
			$vdws[1] = $tmp[1]/$tmp[2];
			$vdws[2] = $tmp[1]**6 * $tmp[0] * ($tmp[2]/($tmp[2] - 6));
		} elsif ($vdwType =~ /^BUCKINGHAM/) {
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			$vdws[1] = 1/$vdws[1];
		}
		$pair->{$type1}{$type2}{TYPE} = $curr->{TYPE};
		@{ $pair->{$type1}{$type2}{VALS} } = @vdws;
	}

	$hybridStr = "";
	for $type1 (keys %{ $pair }) {
		for $type2 (keys %{ $pair->{$type1} }) {
			#next if ($type1 > $type2);
			#$pair->{$type1}{$type2}{VALS} = mix($type1, $type2, $pair) if (! exists($pair->{$type1}{$type2}));
			$hybridStr = getPairName($pair->{$type1}{$type2}{TYPE}) if ($isHybrid);
			push @{ $lmpCmds }, "pair_coeff $type1 $type2 ${hybridStr} @{ $pair->{$type1}{$type2}{VALS} }" if($type1<$type2);
			push @{ $lmpCmds }, "pair_coeff $type2 $type1 ${hybridStr} @{ $pair->{$type1}{$type2}{VALS} }" if($type1>$type2);
		}
	}
	return ($lmpCmds);
}

sub updatePairCoeffs_old {
	my ($ffData, $newParms) = @_;
	my ($outStr, $i, $j, @vdws, @tmp, $curr); 
	my ($pair, $vdwType, $type1, $type2, $hybridStr);
	
	for $i (0 .. $#{ $ffData->{OPT}{map} }) {
		next if (! defined($newParms->[$i]));
		next if ($ffData->{OPT}{map}[$i]{otype} ne "vdw");
		$curr = $ffData->{OPT}{map}[$i]{name}{1}; # j points to $PARMS{VDW}{type1}{type2}
		$curr->{VALS}[$ffData->{OPT}{map}[$i]{pOrder}] = $newParms->[$i]; # update j val to reflect new val from simplex
		$vdwType = uc $curr->{TYPE};
		($type1, $type2) = split /\-/,$ffData->{OPT}{map}[$i]{label};
		$type1 = $ffData->{TYPE_MAP}{$type1};
		$type2 = $ffData->{TYPE_MAP}{$type2};
		@vdws = @{ $curr->{VALS} };
		if ($vdwType =~ /^MORSE|STRETCH_MORSE$/) {
			($vdws[0], $vdws[1]) = ($vdws[1], $vdws[0]);
			($vdws[1], $vdws[2]) = ($vdws[2], $vdws[1]);
			$vdws[1] /= ($vdws[2] * 2); # changed from div to multi 07/28/2007
			$vdws[3] /= ($vdws[2] * 2) if ($vdwType =~ /^STRETCH_MORSE/); #aplha2 for stretch morse
		} elsif ($vdwType =~ /^LJ_6_12$/) {
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			$vdws[1] = $vdws[1] / (2**(1/6));
			if ($#vdws > 1) {
				($vdws[2],$vdws[3]) = ($vdws[3],$vdws[2]);
				$vdws[3] = $vdws[3] /(2**(1/6));
			}
		} elsif ($vdwType =~ /^EXPO_6$/) {
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			@tmp = @vdws;
			@vdws = ();
			$vdws[0] = $tmp[0] * (6/($tmp[2]-6)) * exp($tmp[2]);
			$vdws[1] = $tmp[1]/$tmp[2];
			$vdws[2] = $tmp[1]**6 * $tmp[0] * ($tmp[2]/($tmp[2] - 6));
		} elsif ($vdwType =~ /^BUCKINGHAM/) {
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			$vdws[1] = 1/$vdws[1];
		}
		$pair->{$type1}{$type2}{TYPE} = $curr->{TYPE};
		@{ $pair->{$type1}{$type2}{VALS} } = @vdws;
	}

	$hybridStr = "";
	for $type1 (keys %{ $pair }) {
		for $type2 (keys %{ $pair->{$type1} }) {
			#next if ($type1 > $type2);
			#$pair->{$type1}{$type2}{VALS} = mix($type1, $type2, $pair) if (! exists($pair->{$type1}{$type2}));
			$hybridStr = getPairName($pair->{$type1}{$type2}{TYPE}) if ($isHybrid);
			$outStr .= sprintf("%-12s%8d%8d	 ${hybridStr} @{ $pair->{$type1}{$type2}{VALS} }\n", "pair_coeff", $type1, $type2) if ($type1<=$type2);
			$outStr .= sprintf("%-12s%8d%8d	 ${hybridStr} @{ $pair->{$type1}{$type2}{VALS} }\n", "pair_coeff", $type2, $type1) if ($type1>$type2);
		}
	}
	$outStr .= "\ncompute		 interact group1 group/group group2\n";
	#$outStr .= "compute		 forceX all reduce sum fx\n";
	#$outStr .= "compute		 forceY all reduce sum fy\n";
	#$outStr .= "compute		 forceZ all reduce sum fz\n";
	#$outStr .= "variable		forces equal c_forceX^2+c_forceY^2+c_forceZ^2\n";
	$outStr .= "thermo_style	custom step cpu pe c_interact evdwl\n";
	$outStr .= "thermo_modify   line multi\n";
	$outStr .= "run			 0\n";

	return ($outStr);
}

sub getPairName {
	my ($name) = $_[0];
	my ($return_str);

	$return_str = "lj/charmm/coul/long/opt";
	$return_str = "lj/gromacs" if ($name eq 'LJ_6_12' and $isPQEQ);
	$return_str = "buck/coul/long" if ($name eq "BUCKINGHAM");
	$return_str = "morse/coul/long/opt" if ($name eq "VDW_MORSE");

	return $return_str;
}

sub compareParms {
	my ($pOld, $parms) = @_;
	my ($i, $pNew, $rlist, $pKey);

	$pNew = $parms->{OPT}{pCurr};
	$pKey = $parms->{OPT}{key};
	for $i (0 .. $#{ $pNew }) {
		next if (exists($parms->{OPT}{pSkip}) and exists($parms->{OPT}{pSkip}{$i})); 
		push @{ $rlist }, $pKey->[$i] . sprintf("%10.5G -> %10.5G",$pOld->[$i],$pNew->[$i]) if ($pOld->[$i] != $pNew->[$i]);
	}

	return $rlist;
}

sub mix  {
	my ($t1, $t2, $pairData) = @_;
	my ($aa, $bb, @vals);

	$aa = $pairData->{$t1}{$t1};
	$bb = $pairData->{$t2}{$t2};
	return () if ($aa->{TYPE} ne $bb->{TYPE});
	
	if (uc($aa->{TYPE}) eq "VDW_MORSE") {
		$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
		$vals[1] = 0.5*($aa->{VALS}[1]+$bb->{VALS}[1]);
		$vals[2] = 0.5*($aa->{VALS}[2]+$bb->{VALS}[2]);
	} elsif (uc($aa->{TYPE}) eq "LJ_6_12") {
		$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
		$vals[1] = 0.5*($aa->{VALS}[1]+$bb->{VALS}[1]);
		if ($#{ $aa->{VALS} } > 1) {
			$vals[2] = sqrt($aa->{VALS}[2]*$bb->{VALS}[2]);
			$vals[3] = 0.5*($aa->{VALS}[3]+$bb->{VALS}[3]);
		}
	} elsif (uc($aa->{TYPE}) =~ /EXPO_6|BUCKINGHAM/) {
		$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
		$vals[1] = 0.5*($aa->{VALS}[1]+$bb->{VALS}[1]);
		$vals[2] = sqrt($aa->{VALS}[2]*$bb->{VALS}[2]);
	}
	return \@vals;
}


sub assignBGFS {
	my ($data, $bgfs, $dataFileSstr, $nprefix, $nsuffix) = @_;
	my ($i, $j, $k, $isValid, $count, @tmp, $sstr); 
	my ($engVals, $offset, $matchArrayFields, $engStr, $testStr);

	if(scalar(keys %{ $data }) == 1 and scalar(@{ $bgfs }) == 1) {
		$data->{0}{FILE} = abs_path($bgfs->[0]);
		$count = 1;
		print "read in $count datavalues and assigned corresponding bgffiles...";
		return;
	}

	$count = 0;
	$dataFileSstr = basename($dataFileSstr);
	$dataFileSstr =~ s/\.\w+$//;
	#$dataFileSstr =~ s/\.\w+$//;
	if(defined($nprefix) or defined($nsuffix)) {
		$matchArrayFields = 1;
		@{ $engVals } = split /\./, $dataFileSstr;
		$engStr = "";
		$offset = 0;
		if(defined($nprefix)) {
			for $i (0 .. $nprefix-1) {
				$engStr = "${engStr}" . $engVals->[$i] . ".";
			}
			$offset = $nprefix;
		}
		$engStr .= "dist_here(\.0+A){0,1}\.";
		if(defined($nsuffix)) {
			for $i ($offset .. ($offset + $nsuffix - 1)) {
				$engStr = "${engStr}" . $engVals->[$i] . ".";
			}
		}
	}

	for $i (reverse sort numerically keys %{ $data }) {
		$isValid = 0;
		for $j (0 .. $#{ $bgfs }) {
			if($matchArrayFields) {
				$testStr = $engStr;
				$testStr =~ s/dist_here/${i}/;
				if($bgfs->[$j] =~ /$testStr/) {
					$data->{$i}{FILE} = abs_path($bgfs->[$j]);
					splice @{ $bgfs }, $j, 1;
					$isValid = 1;
					last;
				}
			} else {
				if ($bgfs->[$j] =~ /\.${i}A?\./) {
					$sstr = basename($bgfs->[$j]);
					$sstr =~ s/\.\w+$//;
					$sstr =~ s/\d*${i}A?.*//;
					$sstr =~ s/\.$//;
					if($dataFileSstr =~ /$sstr/i) {
						$data->{$i}{FILE} = abs_path($bgfs->[$j]);
						splice @{ $bgfs }, $j, 1;
						$isValid = 1;
						last;
					}
				}
			}
		}
		die "ERROR: No corresponding bgf file found while search for $i!\n"
			if (! $isValid);
		$count++;
	}
	print "read in $count datavalues and assigned corresponding bgffiles...";
}

sub getObjectiveFunctionData {
	my ($parms, $eFile, $factor, $emin) = @_;
	my (%EFUNC, $rec, $i, $elast, $rlast, $eng, $dist, $isNB);

	$isNB = 0;
	$isNB = 1 if (exists($parms->{VDW_OPT}) or exists($parms->{QEQ_OPT}) or exists($parms->{"3BODY_OPT"}));
	open ENGFILE, $eFile or die "ERROR: Cannot open energy file $eFile: $!\n";
	$parms->{EWEIGHT_FACT} = 1 if ( ! exists($PARMS->{EWEIGHT_FACT}));
	while (<ENGFILE>) {
		chomp;
		if ($_ =~ /^\s*(\-?\d+\.?\d*)\s+(\-?\d*\.\d+E?\-?\+?\d*)\s+(\d*\.?\d*E?\-?\+?\d*)/i) {
			$eng = $factor*$2;
			$dist = $1;
			$dist = $1 + $parms->{DIST_OFFSET} if ($parms->{DIST_OFFSET});
			$rec = (
					{
						"DISP"   => $dist,
						"ENG"	=> $eng,
					}
					);
			$rec->{WEIGHT} = $3 if ($3 ne "" and !exists($parms->{IGNORE_FILE_WEIGHTS}));		
			$$emin = $eng if (!exists($parms->{EMIN_ADJUSTED}) and (!defined($$emin) or $$emin > $eng));
			$EFUNC{$dist} = $rec;
			if (!defined($rlast) or $1>$rlast) {
				$rlast = $dist;
				$elast = $eng;
			}
		} elsif ($_ =~ /^\s*(\-?\d+\.?\d*)\s+(\-?\d*\.\d+E?\+?\-?\d*)/) {
			$eng = $factor*$2;
			$dist = $1;
			$dist = $1+$parms->{DIST_OFFSET} if ($parms->{DIST_OFFSET});
			$rec = (
					{
						"DISP"   => $dist,
						"ENG"	=> $eng,
					}
					);
			$EFUNC{$dist} = $rec;
			$$emin = $eng if (!exists($parms->{EMIN_ADJUSTED}) and (!defined($$emin) or $$emin > $eng));
			if (!defined($rlast) or $1>$rlast) {
				$rlast = $dist;
				$elast = $eng;
			}
		}
	}
	close ENGFILE;
	die "ERROR: Energy file $eFile does not contain any valid information!\n"
		if (! %EFUNC);

	if(!exists($parms->{EMIN_ADJUSTED}) and defined($parms->{ENG_OFFSET}) and $isNB) {
		#$$emin -= $parms->{ENG_OFFSET};
		$parms->{EMIN_ADJUSTED} = 1;
	}
	foreach $i (keys %EFUNC) {
		if(defined($parms->{ENG_OFFSET})) {
			#$EFUNC{$i}{ENG} -= $parms->{ENG_OFFSET};
		} elsif (defined($parms->{ESHIFT}) and $parms->{ESHIFT} == 1) {
			$EFUNC{$i}{ENG} -= $elast;
		}
		if($isNB) {
			#$EFUNC{$i}{WEIGHT}=1;
			next if(exists($EFUNC{$i}{WEIGHT}));
			$EFUNC{$i}{WEIGHT} = exp(-($EFUNC{$i}{ENG}-$$emin)/.592126/$parms->{EWEIGHT_FACT});
		} else {
			if(exists($EFUNC{$i}{WEIGHT})) { 
				$EFUNC{$i}{ENG} -= $$emin;
				next;
			}
			$EFUNC{$i}{WEIGHT} = exp(-($EFUNC{$i}{ENG}-$$emin)/.592126);
			#$EFUNC{$i}{WEIGHT} = 1/($EFUNC{$i}{ENG}-$$emin) if($EFUNC{$i}{ENG} != $$emin);
			$EFUNC{$i}{ENG} -= $$emin;
		}
	}
	return \%EFUNC;
}

sub getBounds {
	my ($ff, $parms) = @_;
	my ($headerStr, $count, $factor, $fMap, $type, $header, $pPtr, $fPtr);

	$headerStr = "vdw|bond|angle|torsion|inversion|3body";
	$count = 0;
	$factor = 0.05;
	$fMap = \%{ $parms->{OPT} };

	if (exists($parms->{VARY}{ENG_OFFSET})) {
		#$parms->{ENG_OFFSET} = 1;
		&getEngOffsetPlist($parms,\%{ $fMap });
	}
	for $type (keys %{ $parms }) {
		next if ($type !~ /($headerStr)_OPT/i or !exists($ff->{uc $1 . "S"}));
		$header = uc $1;
		$pPtr = $parms->{$type};
		$fPtr = $ff->{"${header}S"};
		if($header eq "VDW") {
			&getVDWPlist($pPtr, $fPtr, $header, \%{ $fMap });
		} else {
			&getValencePlist($pPtr, $fPtr, $header, \%{ $fMap });
		}
	}
	if (exists($parms->{QEQ_OPT}) and exists($parms->{QEQ_PARM_FILE})) {
		&parseQEqFile($parms);
		&getQEqPlist($parms, \%{ $fMap });
	}
	&coupleParms($parms, $fMap) if(exists($parms->{COUPLE}));

	die "ERROR: No valid info found while creating objective function!\n" if (! $fMap);
}

sub coupleParms {
	my ($parms, $fMap) = @_;
	my ($i, $j, $k, $v1, $v2, $mode, $n);

	$n = $#{ $fMap->{init} };
	for $i (keys %{ $parms->{COUPLE}}) {
		if($i =~ /^(diff|scale)/i) {
			$mode = lc $1;
			for $j (keys %{ $parms->{COUPLE}{$i}}) {
				$v1 = $j - 1;
				next if ($v1 < 0 or $v1 > $n);
				for $k (keys %{ $parms->{COUPLE}{$i}{$j}}) {
					$v2 = $k - 1;
					next if ($v2 < 0 or $v2 > $n);
					if ($mode eq "diff") {
						push @{ $fMap->{couple} }, "\$parms->{OPT}{pCurr}[$v2]=\$parms->{OPT}{init}[$v2]+\$parms->{OPT}{pCurr}[$v1]-\$parms->{OPT}{init}[$v1]";
					} elsif ($mode eq "scale") {
						push @{ $fMap->{couple} }, "\$parms->{OPT}{pCurr}[$v2]=\$parms->{OPT}{init}[$v2]*\$parms->{OPT}{pCurr}[$v1]/\$parms->{OPT}{init}[$v1]";
					}
					$fMap->{pGuess}[$v2] = $fMap->{pScale}[$v2] = 0;
					$fMap->{pSkip}{$v2} = 1;
				}
			}
		}
	}
}

sub getEngOffsetPlist {
	my ($parms, $plist) = @_;

	$parms->{ENG_OFFSET} = 1 if (! exists($parms->{ENG_OFFSET}) or ! $parms->{ENG_OFFSET});
	push @{ $plist->{pCurr} }, $parms->{ENG_OFFSET};
	push @{ $plist->{pScale} }, $parms->{ENG_OFFSET};
	push @{ $plist->{init} }, $parms->{ENG_OFFSET};
	push @{ $plist->{pGuess} }, $parms->{ENG_OFFSET};
	push @{ $plist->{key} }, sprintf("%-53s","ENG_OFFSET");
}

sub getQEqPlist {
	my ($parms, $plist) = @_;
	my ($validQEq, $i, $j, $pOrder, $currParm, $incr, $used, $factor);

	$validQEq = 0;
	for $i (keys %{ $parms->{QEQ_OPT} }) {
		next if (!exists($parms->{QEq}{$i}));
		for $j (keys %{$parms->{QEQ_OPT}{$i} }) {
			$pOrder = getQEqParmOrder($j);
			next if (exists($used->{qeq}{$i}{$j}));
			next if (!defined($pOrder));
			$validQEq = 1;
			$currParm = $parms->{QEq}{$i};
			$incr = $factor * $currParm->{1}{VALS}[$pOrder];
			if($currParm->{1}{VALS}[$pOrder] == 0) {
				$incr = 0.1;
				$incr = 1.0 if ($i eq "chi");
				$incr = 10.0 if ($i eq "gamma");
			}
			$incr = 0.1 if ($incr == 0);
			push @{ $plist->{pGuess} }, $currParm->{1}{VALS}[$pOrder];
			push @{ $plist->{pScale} }, $incr;
			#push @{ $plist->{pGuess} }, $fPtr->{$i}{$j}{VALS}[$indx];
			#push @{ $plist->{pCurr} },  $fPtr->{$i}{$j}{VALS}[$indx];
			#push @{ $plist->{pScale} }, $fPtr->{$i}{$j}{VALS}[$indx]*$factor;
			#push @{ $plist->{init} },   $fPtr->{$i}{$j}{VALS}[$indx];
			$used->{qeq}{$i}{$j} = 1;
		}
	}
	if ($validQEq) {
		`sed -i '/fix.*qeq/d' _lmptmp/in.tmp.qeq`;
		`echo 'fix charge all qeq/rg 1 8 1.0e-6 100 _lmptmp/tmp.qeq.dat' >> _lmptmp/in.tmp.qeq`;
	} else {
		delete $parms->{QEQ_OPT};
	}
}

sub getValencePlist {
	my ($pPtr, $fPtr, $type, $plist) = @_;
	my ($i, $j, $k, $curr, $indx, $nKeys); 
	my ($factor, $lmpName, $c2Name, $tmp);

	$factor = 0.15;
	for $i (keys %{ $fPtr }) {
		for $j (keys %{ $fPtr->{$i} }) {
			#search pPtr ofr matching key
			$curr = ();
			&findVal($pPtr, $fPtr->{$i}{$j}{KEY}, $j, \%{ $curr });
			if (! keys %{ $curr } and ($type !~ /3BODY/ or $j !~ /sw/)) {
				$curr = ();
				$nKeys = getrevKeys($fPtr->{$i}{$j}{KEY});
				&findVal($pPtr, $nKeys, $j, \%{ $curr });
			}
			next if (! keys %{ $curr });
			($lmpName, $c2Name) = getTypeName($type,$j);
			for $k (keys %{ $curr->{VAL} }) {
				next if(!exists($typeMap->{$type}{$lmpName}{PARMS}{lc $k}));
				$indx = $typeMap->{$type}{$lmpName}{PARMS}{lc $k};
				push @{ $plist->{pGuess} }, $fPtr->{$i}{$j}{VALS}[$indx];
				push @{ $plist->{pCurr} },  $fPtr->{$i}{$j}{VALS}[$indx];
				push @{ $plist->{pScale} }, $fPtr->{$i}{$j}{VALS}[$indx]*$factor;
				push @{ $plist->{init} },   $fPtr->{$i}{$j}{VALS}[$indx];
				$tmp = ();
				&flattenHash($fPtr->{$i}{$j}{KEY}, \@{ $tmp });
				push @{ $plist->{key} }, sprintf("%-20s %10s %10s @{ $tmp }",$type,$j,$k);
				$fPtr->{$i}{$j}{USED} = 1;
				$fPtr->{$i}{$j}{VALSIDXPTR}[$indx] = $#{ $plist->{pCurr} };
			}
		}
	}
}

sub flattenHash {
	my ($kStr, $ret) = @_;
	my ($i);

	for $i (keys %{ $kStr }) {
		push @{ $ret }, $i;
		&flattenHash($kStr->{$i}, \@{ $ret });
	}
}

sub getrevKeys {
	my ($k) = $_[0];
	my ($tmp, $rKeys, $curr, $i);

	$tmp = ();
	&reverseKeys($k, \@{ $tmp });
	$rKeys = ();
	$curr = \%{ $rKeys };
	for $i (@{ $tmp }) {
		$curr->{$i} = ();
		$curr = \%{ $curr->{$i} };
	}

	return $rKeys;
}

sub reverseKeys {
	my ($kList, $holder) = @_;
	my ($i);

	for $i (keys %{ $kList }) {
		unshift @{ $holder }, $i;
		&reverseKeys($kList->{$i}, \@{ $holder });
	}
}

sub findVal {
	my ($pPtr, $kPtr, $sstr, $curr) = @_;
	my ($i, $tmp, $kStr);

	for $i (keys %{ $pPtr }) {
		if($i eq $sstr) {
			$curr->{VAL} = $pPtr->{$i};
			last;
		}
		last if (keys %{ $curr });
		@{ $tmp } = keys %{ $kPtr };
		$kStr = $tmp->[0];
		if($i eq $kStr) {
			&findVal($pPtr->{$i}, $kPtr->{$i}, $sstr, \%{ $curr });
		}
	}
}

sub getVDWPlist {
	my ($pPtr, $fPtr, $type, $plist) = @_;
	my ($i, $j, $k, $l, $iType, $jType, $indx, $curr);
	my ($lmpName, $c2Name, $valid);

	$valid = 0;
	for $i (keys %{ $pPtr }) {
		$iType = getTypeIndex($i);
		next if (!defined($iType) or ! exists($fPtr->{$iType}));
		for $j (keys %{ $pPtr->{$i} }) {
			$jType = getTypeIndex($j);
			next if (! defined($jType));
			for $k (keys %{ $pPtr->{$i}{$j} }) {
				if(exists($fPtr->{$iType}) and exists($fPtr->{$iType}{$jType}) and exists($fPtr->{$iType}{$jType}{$k})) {
					$curr = \%{ $fPtr->{$iType}{$jType}{$k} };
				} elsif(exists($fPtr->{$jType}) and exists($fPtr->{$jType}{$iType}) and exists($fPtr->{$jType}{$iType}{$k})) {
					$curr = \%{ $fPtr->{$jType}{$iType}{$k} };
				} else {
					next;
				}
				($lmpName, $c2Name) = getTypeName("VDW",$k);
				for $l (keys %{ $pPtr->{$i}{$j}{$k} }) {
					next if(!exists($typeMap->{$type}{$lmpName}{PARMS}{lc $l}));
					$indx = setplistVal($plist,$pPtr->{$i}{$j}{$k},$curr,$typeMap->{$type}{$lmpName},$l);
					if(exists($hybridTypes->{VDW})) {
						push @{ $plist->{key} }, sprintf("%-20s %10s %10s %10s",$k, $i, $j, $l);
					} else {
						push @{ $plist->{key} }, sprintf("%10s %10s %10s",$i, $j, $l);
					}
					$curr->{USED} = 1;
					$valid++;
				}
			}
		}
	}

	die "ERROR: Cannot find any relevant search file in VDW_OPT. Check parm file!\n"
		if (! $valid);

}

sub setplistVal {
	my ($plist, $pPtr, $curr, $tMap, $tKey) = @_;
	my ($indx, $currVal, $i, $factor, $cVal);

	$indx = $tMap->{PARMS}{lc $tKey};
	$factor = 0.5;
	#first check for dependence
	if(exists($tMap->{DEPENDENCE}) and exists($tMap->{DEPENDENCE}{$tKey})) {
		#now check that the dependent variable and the independent variable are not both specified
		for $i (grep {!/^calc/} keys %{ $tMap->{DEPENDENCE}{$tKey} }) {
			die "ERROR: Both the dependent variable '$tKey' and independent variable '$i' specified for optimization!\n"
				if(exists($pPtr->{$i}));
		}
		#calculate the value and set the flag
		$curr->{VALS}[$indx] = eval($tMap->{DEPENDENCE}{$tKey}{val});
		$curr->{DEPENDENCE}[$indx] = $tMap->{DEPENDENCE}{$tKey}{calc};
	}
	$cVal = $curr->{VALS}[$indx];
	push @{ $plist->{pGuess} }, $cVal;
	push @{ $plist->{pCurr} }, $cVal;
	push @{ $plist->{pScale} }, $cVal*$factor;
	push @{ $plist->{init} }, $cVal;

	$curr->{VALSIDXPTR}[$indx] = $#{ $plist->{pCurr} };
	if(exists($tMap->{EQUILV}) && exists($tMap->{EQUILV}{lc $tKey})) { #equivalent parameters
		$indx = $tMap->{PARMS}{ $tMap->{EQUILV}{lc $tKey} };
		$curr->{VALSIDXPTR}[$indx] = $#{ $plist->{pCurr} };
	}
	return $indx;
}

sub getTypeIndex {
	my ($sStr) = $_[0];
	my ($i, $typeIndx, $found, $key);

	$found = 0;
	for $i (keys %{ $ffData->{ATOMTYPES} }) {
		@{ $key } = keys %{ $ffData->{ATOMTYPES}{$i}{KEY} };
		if ($key->[0] eq $sStr) {
			$found = 1;
			$typeIndx = $i;
			last;
		}
	}
	return $typeIndx;
}

sub getQEqParmOrder {
	my ($key) = $_[0];

	return 0 if ($key =~ /chi/i);
	return 1 if ($key =~ /gamma/i);
	return 3 if ($key =~ /zeta/i);
	return;
}

sub parseQEqFile {
	my ($parms) = $_[0];
	my ($qeqfile, $typeID);

	$qeqfile = $parms->{QEQ_PARM_FILE};
	die "ERROR: Cannot access $qeqfile: $!\n"
		if(! -e $qeqfile or ! -r $qeqfile or ! -T $qeqfile);
	open QEQFILE, $qeqfile or die "ERROR: Cannot open $qeqfile: $!\n";
	while(<QEQFILE>) {
		chomp;
		if ($_ =~ /^(\S+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
			next if (!exists($parms->{TYPE_MAP}{$1}));
			$typeID = $parms->{TYPE_MAP}{$1};
			$parms->{QEq}{$1}{name} = $1;
			$parms->{QEq}{$1}{typeID} = $typeID;
			@{ $parms->{QEq}{$1}{1}{VALS} } = ($2,$3,$4,$5,$6);
		}
	}
	close QEQFILE;

	die "ERROR: No valid data found while reading $qeqfile!\n"
		if (!exists($parms->{QEq}));

}

sub getVDwParmOrder {
	my ($vdwParm) = $_[0];

	return 0 if ($vdwParm =~ /r0/i); 
	return 1 if ($vdwParm =~ /d0/i); 
	return 2 if ($vdwParm =~ /alpha/i);
	return undef

}

sub cleanup {
	my ($data) = $_[0];
	my ($i);

	for $i (keys %{ $data }) {
		$data->{$i}{FULL}{LMP}->close();
		$data->{$i}{FRAG1}{LMP}->close() if(exists($data->{$i}{FRAG1}));
		$data->{$i}{FRAG2}{LMP}->close() if(exists($data->{$i}{FRAG1}));
	}
}

sub printFile2Screen {
	open INFILE, $_[0] or die "ERROR: Cannot access $_[0]:$!\n";
	while(<INFILE>) {
		print $_;
	}
	close INFILE;
}

sub getC2Header {
	my ($sstr) = $_[0];

	if($sstr eq "BONDS") {
		return ("BOND_STRETCH", 2);
	} elsif($sstr eq "ANGLES") {
		return ("ANGLE_BEND", 3);
	} elsif($sstr eq "TORSIONS") {
		return ("TORSIONS", 4);
	} elsif($sstr eq "INVERSIONS") {
		return ("INVERSIONS",4);
	} elsif($sstr eq "3BODY") {
		return ("3BODY", 3);
	}
}

sub getElement {
	my ($fftype) = $_[0];
	my ($eleNum);
	
	$fftype =~ s/_.*//;
	$fftype =~ s/\d+.*//;
	$eleNum = FindElement($ELEMENTS, $fftype);
	die "ERROR: Cannot locate element for fftype $fftype!\n"
		if (! defined($eleNum));
	return $ELEMENTS->{$eleNum}{SYMBOL};
}

sub calcIP_EA {
	my ($curr, $ipIndx, $eaIndx) = @_;
	my ($ip, $ea);

	#first get the ip
	if(!exists($curr->{VALS}[$ipIndx]) or ! defined($curr->{VALS}[$ipIndx])) {
		$ip = $curr->{VALS}[0]+$curr->{VALS}[1]/2;
	} else {
		$ip = $curr->{VALS}[$ipIndx];
	}

	#now get the ea
	if(!exists($curr->{VALS}[$eaIndx]) or ! defined($curr->{VALS}[$eaIndx])) {
		$ea = $curr->{VALS}[0]-$curr->{VALS}[1]/2;
	} else {
		$ea = $curr->{VALS}[$eaIndx];
	}
	
	#chi
	$curr->{VALS}[0] = ($ip + $ea)/2;
	#eta
	$curr->{VALS}[1] = $ip - $ea;
}

sub init {
	my (%OPTS);
	my ($tmp, $findCmd, $ffList, $singleBGF, $flist, $rec);

	getopt('psc',\%OPTS);
	($paramFile, $saveName, $tmp) = ($OPTS{p}, $OPTS{s}, $OPTS{c});
	die "usage: $0 -p param_file -s (save_name) -c (check_only=0)\n"
		if (! exists($OPTS{p}));
	print "Initializing...";
	FileTester($OPTS{p});
	if(defined($tmp)) {
		$tmp = 1 if ($tmp =~ /yes|1/i);
		$tmp = 0 if ($tmp !~ /yes|1/i);
		$PARMS->{CHECK_ONLY} = $tmp;
	}
	$saveName = "ffopt.dat" if (! defined($saveName));
	$lammpsBinary = "$Bin/../codes/bin/lmp_mpi";
	$oldParms = ();
	$isHybrid = 0;
	print "Done\n";
	$ELEMENTS = &LoadElements;

	$typeMap = {
		'VDW' => {
			'hbond/dreiding/morse' => {
				'c2Name' => "MORSE_COSP",
				'c2INDX' => [0,1,2],
				'c2FACT' => [1,1,1],
				'c2FACU' => [1,1,1],
				'PARMS'  => {
					'd0'    => 2,
					'alpha' => 3,
					'r0'    => 4,
				},
			},
			'yukawa'     => {
				'c2Name' => "YUKAWA",
				'PARMS'  => {
						'a' => 0,
						'kappa' => 1,
						'cutoff' => 2,
				},
				'c2INDX' => [0, 1, 2],
				'c2FACT' => [1, 1, 1],
				'c2FACU' => [1, 1, 1],
			},
			'rexpon/unb'     => {
				'c2Name' => "REXPON_UNB",
				'PARMS'  => {
						'd0' => 0,
						'r0' => 1,
						'l' => 2,
				},
				'c2INDX' => [0, 1, 2],
				'c2FACT' => [1, 1, 1],
				'c2FACU' => [1, 1, 1],
			},
			'lj/gromacs' => {
				'c2Name' => 'LJ_6_12',
				'PARMS'  => {
						'd0' => 0,
						'r0' => 1,
				},
				'c2INDX' => [1, 0],
				'c2FACT' => [2**(1/6), 1],
				'c2FACU' => [2**(1/6), 1],
			},
			'buck'       => {
				'c2Name' => 'BUCKINGHAM',
				'PARMS'  => {
						'a'  => 0,
						'r'  => 1,
						'c'  => 2,
				},
				'c2INDX' => [0, 1, 2],
				'c2FACT' => [1, 'my $tdx=$typeMap->{VDW}{$lmpName}{c2INDX}[2]; 1/($parms->{OPT}{final}[$tdx]**2)', 1],
				'c2FACU' => [1, 'my $tmp=$ff->{VDWS}{$iType}{$jType}{$k}{VALS}[2]; 1/($tmp**2)', 1],
			},
			'morse/opt'     => {
				'c2Name' => 'VDW_MORSE',
				'PARMS'  => {
						'd0'    => 0,
						'alpha' => 1,
						'r0'    => 2,
				},
				'c2INDX' => [2, 0, 1],
				'c2FACT' => [1, 1, 'my $n = $typeMap->{VDW}{$lmpName}{c2INDX}[0]; my $tdx = $ff->{VDWS}{$iType}{$jType}{$k}{VALSIDXPTR}[$n]; 2*$parms->{OPT}{final}[$tdx]'],
				'c2FACU' => [1, 1, 'my $n = $ff->{VDWS}{$iType}{$jType}{$k}{VALS}[0]; my $tdx = $ff->{VDWS}{$iType}{$jType}{$k}{VALSIDXPTR}[$n]; 2*$parms->{OPT}{final}[$tdx]'],
			},
			'lj/charmm' => {
				'c2Name' => 'LJ_6_12',
				'PARMS'  => {
						'd0' => 0,
						'r0' => 1,
				},
				'c2INDX' => [1, 0],
				'c2FACT' => [2**(1/6), 1],
				'c2FACU' => [2**(1/6), 1],
			},
			'coul/pqeqgauss' => {
				'c2Name' => 'PQEQGAUSS',
				'PARMS'  => {
					'chi'   => 0,
					'eta'   => 1,
					'rc'    => 2,
					'p'     => 3,
					'zcore' => 4,
					'rs'    => 5,
					'ks'    => 6,
					'ks1'   => 7,
					'ip'    => 8, #note that the IP is a dependent variable, used to determine chi and eta
					'ea'    => 9, #also a dependent variable
				},
				'EQUILV' => {
					'rc' => 'rs',
					'rs' => 'rc',
				},
				'DEPENDENCE' => {
					'ip'  => {
						'chi'  => 1,
						'eta'  => 1,
						'val'  => '$curr->{VALS}[0]+$curr->{VALS}[1]/2',
						'calc' => '&calcIP_EA($curr,8,9)',
					},
					'ea' => {
						'chi'  => 1,
						'eta'  => 1,
						'val'  => '$curr->{VALS}[0]-$curr->{VALS}[1]/2',
						'calc' => '&calcIP_EA($curr,8,9)',
					}
				},
				'c2INDX' => [3, 0, 1, 4, 2, 5, 6, 7, -1, -1], #note that the index for dependent variables is -1
				'c2FACT' => [1, 1, 1, 1, 1, 1, 1, 1, -1, -1],
				'c2FACU' => [1, 1, 1, 1, 1, 1, 1, 1, -1, -1],
			},
		},
		'BOND' => {
			'harmonic' => {
				'c2Name' => 'HARMONIC',
				'PARMS'  => {
						'k'  => 0,
						'r0' => 1,
				},
				'c2INDX' => [0, 1],
				'c2FACT' => [2, 1],
				'c2FACU' => [2, 1],
			},
			'morse'    => {
				'c2Name' => 'MORSE',
				'PARMS'  => {
						'k'     => 0,
						'alpha' => 1,
						'r0'    => 2,
				},
				'c2INDX' => [1, 2, 0],
				'c2FACT' => ['my $n = $typeMap->{$i}{$k}{c2INDX}[2]; my $tdx = $ff->{"${i}S"}{$j}{$k}{VALSIDXPTR}[$n]; $parms->{OPT}{final}[$idx]*2*$parms->{OPT}{final}[$tdx]',1,1],
				'c2FACU' => ['$ff->{BONDS}{$iType}{$jType}{$k}{VALS}[2]*2*$ff->{BONDS}{$iType}{$jType}{$k}{VALS}[0]',1,1],
			},
		},
		'ANGLE' => {
			'cosine/squared' => {
				'c2Name' => 'COS_HARMON',
				'PARMS' => {
						'k'      => 0,
						'theta0' => 1,
				},
				'c2INDX' => [0, 1],
				'c2FACT' => [1, 1],
				'c2FACU' => [1, 1],
			},
			'harmonic' => {
				'c2Name' => 'THETA_HARM',
				'PARMS' => {
						'k'      => 0,
						'theta0' => 1,
				},
				'c2INDX' => [0, 1],
				'c2FACT' => [1, 1],
				'c2FACU' => [1, 1],
			},
			'quartic' => {
				'c2Name' => 'QUARTIC',
				'PARMS' => {
						'theta0' => 0,
						'k1'     => 1,
						'k2'     => 2,
						'k3'     => 3,
				},
				'c2INDX' => [0, 1, 2, 3],
				'c2FACT' => [1, 1, 1, 1],
				'c2FACU' => [1, 1, 1, 1],
			},
		},
		'3BODY' => {
			'sw' => {
				'c2Name' => 'SW',
				'c2Head' => "ANGLE_BEND",
				'PARMS'  => {
					     'epsilon'    => 0,
						 'sigma'      => 1,
						 'acut'       => 2,
						 'lambda'     => 3,
						 'gamma'      => 4,
						 'costheta0'  => 5,
						 'a'          => 6,
						 'b'          => 7,
						 'p'          => 8,
						 'q'          => 9,
						 'tol'        => 10,
				 },
				 'c2INDX' => [0,1,2,3,4,5,6,7,8,9,10],
				 'c2FACT' => [0.5,1,1,1,1,1,1,1,1,1,1],
				 'c2FACU' => [0.5,1,1,1,1,1,1,1,1,1,1],
			},
		},
	};
}

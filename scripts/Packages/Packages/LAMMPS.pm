package LAMMPS;

require Exporter;
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use General qw(PrintProgress GetTime);
use File::Basename qw(basename);
use Math::Trig;

our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);
our $polarizable;
my ($rigidMols, $rigidType, $rigidOpt, $rigidNVE, $flexMols);
use strict;

my $scripts_dir = "$Bin"; #change here as necessary
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(ParseLammpsTrajectoryFile ReadDataFile ParseLAMMPSLogFile CreateLAMMPSTrj 
				ParseLAMMPSTrj CreateInputFile GetLammpsByteOffset GetLammpsTrjType 
				ConvertLammpsBox ParseLAMMPSTrjSave);
$VERSION = "1.00";

sub determineIfScaled {
	my ($DATA, $OPTS, $tmp2) = @_;
	my ($i, $j);
	my (@dim) = ("XCOORD", "YCOORD", "ZCOORD");

	$OPTS->{ISIFT} = 0;
	$OPTS->{scaled} = $OPTS->{imaged} = 1;
	for $i (keys %{ $DATA->{ATOMS} }) {
		$OPTS->{ISIFT} = 1 if ($DATA->{ATOMS}{$i}{SYZ});
		for $j (0 .. $#dim) {
			$OPTS->{scaled} = 0 if ($DATA->{ATOMS}{$i}{$dim[$j]} > 2 or $DATA->{ATOMS}{$i}{$dim[$j]} < -2);
			if (! $OPTS->{scaled}) {
				$OPTS->{imaged} = 0 if ($DATA->{ATOMS}{$i}{$dim[$j]} > $DATA->{"BOX BOUNDS"}[$j]{hi} or
										$DATA->{ATOMS}{$i}{$dim[$j]} < $DATA->{"BOX BOUNDS"}[$j]{lo});
				last;
			}
		}
	}
	print "reading ";
	if (! $OPTS->{scaled}) {
		print "unscaled ";
	} else {
		print "scaled ";
	}
	if ($OPTS->{imaged}) {
		print "imaged ";
	} else {
		print "unimaged ";
	}
	print "coordinates...";
}

sub GetLammpsTrjType {
	my ($SELECT, $trjFile, $field, $OPTS) = @_;
	my ($i, $tmp1);

	print "Determining LAMMPS trajectory type...";
	for $i (keys %{ $SELECT }) {
		$tmp1->{$i} = $SELECT->{$i};
		last;
	}
	ParseLAMMPSTrj($OPTS, $trjFile, $tmp1, $field, \&determineIfScaled, undef, undef);
	print "Done\n";
}

sub ReadDataFile {
	my ($inFile) = $_[0];
	my (%DATA, $var, @val, $i);

	$var = "";
	open DATAFILE, $inFile or die "ERROR: Cannot open LAMMPS data file $inFile: $!\n";
	while (<DATAFILE>) {
		chomp;
		if ($_ =~ /^\s*(\d+) (atoms|bonds|angles|dihedrals|impropers)\s*$/) {
			$DATA{TOTAL}{uc($2)} = $1;
			$DATA{TOTAL}{count}++;
		} elsif ($_ =~ /^\s*(\d+) (atom|bond|angle|dihedral|improper) types\s*$/) {
			$DATA{TYPE}{uc($2)} = $1;
			$DATA{TYPE}{count}++;
		} elsif ($_ =~ /^\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(x|y|z)lo (x|y|z)hi\s*$/) {
			$DATA{BOX}{uc($3)}{lo} = $1;
			$DATA{BOX}{uc($3)}{hi} = $2;
			$DATA{BOX}{count}++;
		} elsif ($_ =~ /^Masses\s*$/) {
			$var = "MASS";
		} elsif ($_ =~ /^(\w+) Coeffs/) {
			$var = uc($1) . " COEFFS";
		} elsif ($_ =~ /^(Atoms|Bonds|Angles|Dihedrals|Impropers|Velocities)\s*$/) {
			$var = uc($1);
		} elsif (defined($var) and $_ =~ /^\s*(\d+)\s+(.+)$/) {
			@val = split /\s+/, $2;
			for $i (0 .. $#val) {
				$DATA{$var}{$1}{$i} = $val[$i];
			}
		} elsif ($_ =~ /^\S+/) {
			undef($var);
		}
	}

	for $var ("TYPE", "TOTAL", "BOX", "MASS", "ATOMS") {
		die "ERROR: Relevant header (" . uc($var) . ") not found when reading file!\n"
			if (! exists($DATA{$var})); 
	}

	die "ERROR: Invalid format in data file\n"
		if ($DATA{TYPE}{count} < 5 || $DATA{TOTAL}{count} < 5 || $DATA{BOX}{count} < 3);

	return \%DATA;
			
}

sub ParseLammpsTrajectoryFile {
	my ($fileName, $dType, $selection) = @_;
	my (%LINE, $timeStep, $totAtoms, %DATA, $atomC, $scaled);
	my (@dataVals, $patern, $BYTEINFO, $i);

	$timeStep = $atomC = 0;
	$scaled = 1;
	$patern = '^(\d+)\s+(\d+\s+\-?\d+\.?\d*\s+\-?\d+\.?\d*\s+\-?\d+\.?\d*)\s+';
	#$patern .= '(\-?\d+\s+\-?\d+\s+\-?\d+)?\s+';
	$patern .= '(\-?\d+\.?\d*e?\-?\d*)?\s$';

	#getByteOffset($selection, $fileName);
	open TRAJFILE, $fileName or die "ERROR: Cannot open trajectory file $fileName: $!\n";
	while (<TRAJFILE>) {
		chomp;
		if ($_ =~ /^\s*ITEM: TIMESTEP/) {
			setZero(\%LINE);
			$LINE{"tstep"} = 1;
			if ($timeStep > 0) {
				delete $DATA{$timeStep} if (! exists($DATA{$timeStep}{"ATOMS"}) and lc($dType) ne "density");
			}
			$timeStep = $atomC = 0;
		} elsif ($_ =~ /^ITEM: NUMBER OF ATOMS/ and $timeStep > 0) {
			setZero(\%LINE);
			$LINE{"num_atoms"} = 1;
		} elsif ($_ =~ /^ITEM: BOX BOUNDS/ and $timeStep > 0) {
			setZero(\%LINE);
			$LINE{"box"} = 1;
				$LINE{triclinic} = 0;
				$LINE{triclinic} = 1 if ($_ =~ /xy xz yz/);
		} elsif ($_ =~ /^ITEM: ATOMS/ and $timeStep > 0) {
			setZero(\%LINE);
			$LINE{"atoms"} = 1;
		} elsif ($_ =~ /^(\d+)$/ and $LINE{"tstep"}) {
			setZero(\%LINE);
			$timeStep = $1;
			#print "GOT TIMESTEP $timeStep\n";
		} elsif ($_ =~ /^(\d+)$/ and $LINE{"num_atoms"} and $timeStep > 0) {
			setZero(\%LINE);
			$totAtoms = $1;
			$DATA{$timeStep}{"TOTAL_ATOMS"} = $totAtoms;
		} elsif ($_ =~ /^(\-?\d+\.?\d*)\s+(\-?\d+\.?\d*)/ and $LINE{"box"} and $timeStep > 0) {
			setZero(\%LINE);
			$LINE{"box"} = 1;
			push @{ $DATA{$timeStep}{"BOX"} }, (
												{
													"lo" => $1,
													"hi" => $2,
												}
												);
		} elsif ($LINE{"atoms"} and $timeStep > 0 and $_ =~ /^(\d+)\s(.+)/) {
			$atomC = $1;
			next if (defined($selection) and $selection ne "*" and ! exists($selection->{$atomC}));
			@dataVals = split /\s/, $2;
			if ($dType =~ /bgf|com|avg|rms/) {
				die "ERROR: Atom data not found for atom $atomC\n" if ($#dataVals < 3);
				$DATA{$timeStep}{"ATOMS"}{$atomC}{"TYPEID"} = shift @dataVals;
				$DATA{$timeStep}{"ATOMS"}{$atomC}{"XCOORD"} = $dataVals[0];
				$DATA{$timeStep}{"ATOMS"}{$atomC}{"YCOORD"} = $dataVals[1];
				$DATA{$timeStep}{"ATOMS"}{$atomC}{"ZCOORD"} = $dataVals[2];
				if ($scaled) {
					for $i (0 .. 2) {
						if ($dataVals[$i] > 2) {
							$scaled = 0;
							last;
						}
					}
				}
				if ($#dataVals > 4) {
					$DATA{$timeStep}{"ATOMS"}{$atomC}{"XINDEX"} = $dataVals[3];
					$DATA{$timeStep}{"ATOMS"}{$atomC}{"YINDEX"} = $dataVals[4];
					$DATA{$timeStep}{"ATOMS"}{$atomC}{"ZINDEX"} = $dataVals[5];
				}
			}

			if ($dType =~ /vel/) {
				die "ERROR: Veleocity data not found for atom $atomC\n" if ($#dataVals < 2);
				$DATA{$timeStep}{"ATOMS"}{$atomC}{"XCOORD"} = $dataVals[0];
				$DATA{$timeStep}{"ATOMS"}{$atomC}{"YCOORD"} = $dataVals[1];
				$DATA{$timeStep}{"ATOMS"}{$atomC}{"ZCOORD"} = $dataVals[2];
			}

			if ($dType =~ /eng/) {
				$DATA{$timeStep}{"ATOMS"}{$atomC}{"ENERGY"} = pop @dataVals;
				$scaled = 0;
			}
		}
	}
	close TRAJFILE or die "ERROR: Cannot close file $fileName: $!\n";

	die "ERROR: $fileName does not contain any valid informaton\n"
		if (!%DATA);
	
	return (\%DATA, $scaled);
}

sub setZero {
	my ($LINE, $exception) = @_;
	
	for (keys %{ $LINE }) {
		next
			if ($exception and $exception =~ /\s+$_/);
		$LINE->{$_} = 0;
	}
}

sub numerically {
	($a<=>$b);
}

sub ParseLAMMPSLogFile {
	my ($inFile, $selection, $saveFunc, $OUTDATA) = @_;
	my ($tStep, %DATA, $rem, $count);

	$tStep = -1;
	$count = 0;
	open LOGFILE, $inFile or die "ERROR: Cannot read from LAMMPS log file $inFile: $!\n";
	while (<LOGFILE>) {
		chomp;
		if ($_ =~ /^---------------- Step\s+(\d+)\s----- CPU =\s+(\d+\.\d+)/) {
			if (keys %DATA) {
				$saveFunc->(\%DATA, $tStep, $OUTDATA);
				%DATA = ();
				$tStep = -1;
			}
			next if (defined($selection) and ! exists($selection->{$1}));
			$tStep = $1;
			$DATA{CPU} = $2;
			$count++;
		} elsif ($tStep > -1 && $_ =~ /^(\S+)\s+\=\s+(\-?\d+\.?\d*E?\-?\d*)(.+)/) {
			$DATA{lc($1)} = $2;
			$rem = $3;
			while ($rem =~ /(\S+)\s+\=\s+(\-?\d+\.?\d*E?\-?\d*)/g) {
				$DATA{lc($1)} = $2;
			}
		} else {
			if (keys %DATA) {
				$saveFunc->(\%DATA, $tStep, $OUTDATA);
				%DATA = ();
			}
			$tStep = -1;
		}
	}

	$saveFunc->(\%DATA, $tStep, $OUTDATA) if (keys %DATA);
	close LOGFILE;
	
	die "ERROR: LAMMPS log file $inFile does not contain any valid information!\n"
		if (! $count);
}

sub CreateLAMMPSTrj {
	my ($ATOMS, $HEADERS, $OUTFILE) = @_;
	my ($i, $atomC, $atom, @dim, $index, $j, $bLen);

	@dim = ("TYPE", 
			"XCOORD", "YCOORD", "ZCOORD", 
			"XINDEX", "YINDEX", "ZINDEX",
			"XVEL", "YVEL", "ZVEL");

	for $i ("TIMESTEP", "NUMBER OF ATOMS") {
		print $OUTFILE "ITEM: $i\n";
		for $j (@{ $HEADERS->{$i} }) {
			print $OUTFILE "$j\n";
		}
	}
	# BOX
	print $OUTFILE "ITEM: BOX BOUNDS\n";
	for $i (@{ $HEADERS->{"BOX BOUNDS"} }) {
		print $OUTFILE "$i->{lo} $i->{hi}\n";
	}
	print $OUTFILE "ITEM: ATOMS id xu yu zu\n";
	for $atomC (sort numerically keys %{ $ATOMS }) {
		$atom = \%{ $ATOMS->{$atomC} };
		print $OUTFILE "$atomC ";
		for $i (@dim) {
			next if (! exists($atom->{$i}));
			print $OUTFILE "$atom->{$i} ";
		}
		print $OUTFILE "\n";
	}
}	

sub getFormat {
	my ($type) = $_[0];
	my (@fields);

	$type = lc($type);
	if ($type eq "vel") {
		@fields = ("XVEL", "YVEL", "ZVEL");
	} elsif ($type eq "atom") {
		@fields = ("TYPE", "XCOORD", "YCOORD", "ZCOORD", "XINDEX", "YINDEX", "ZINDEX", "CHARGE");
	} elsif ($type eq "ift") {
		#@fields = ("TYPE", "XCOORD", "YCOORD", "ZCOORD", "XINDEX", "YINDEX", 
						#"ZINDEX", "SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ");
		@fields = ("TYPE", "XCOORD", "YCOORD", "ZCOORD", "SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ");
	} elsif ($type eq "energy") {
		@fields = ("TYPES", "ENERGY");
	}

	return \@fields;
}

sub getTotBytes {
	my ($inFile) = $_[0];
	my ($count) = 0;
	if (open(INFILE, $inFile)) {
		while (<INFILE>) {
			chomp;
			if ($_ =~ /^totalbytes: (\d+)/) {
				$count = $1;
				last;
			}
		}
		close INFILE;
	}
	return $count;
}

sub GetLammpsByteOffset {
	my ($selection, $trjFile, $junk) = @_;
	my ($counter, $isValid, $hasSelections, $count, $tot, $saveName, $offset, $last, $bSelect, $shouldWrite);
	my ($grepCmd) = "grep -b 'ITEM: TIMESTEP' $trjFile";

	$shouldWrite = 1;
	$saveName = basename($trjFile);
	$saveName = ".byte.offset.${saveName}";
	print "Computing byte offset for LAMMPS trajectory $trjFile...";
	for $counter (keys %{ $selection }) {
		$bSelect->{$counter}{start} = -1;
	}
	$counter = $isValid = $hasSelections = $count = $tot = 0;
	$tot = (-s $trjFile);
	$hasSelections = 1 if (keys %{ $selection });
	if (! -e $saveName) {
		open GREPCMD, "$grepCmd |" || die "ERROR: Cannot execute $grepCmd: $!\n";
	} else {
		$count = getTotBytes($saveName);
		if (($count != $tot) or ! open(GREPCMD, $saveName)) {
			open GREPCMD, "$grepCmd |" || die "ERROR: Cannot execute $grepCmd: $!\n";
		} else {
			print "reading file...";
			$shouldWrite = 0;
		}
	}
	while (<GREPCMD>) {
		chomp;
		if ($_ =~ /^(\d+)/) {
			$counter++;
			$offset->{$counter} = $1;
			if (! $hasSelections || exists($selection->{$counter}) || exists($selection->{all})) {
				$bSelect->{$counter}{start} = $1;
				$isValid = 1;
			}
			$last = $1;
		}
	}
	close GREPCMD;

	if (exists($selection->{"-1"})) {
		$bSelect->{-1}{start} = $last;
		$isValid = 1;
	}
	if ($shouldWrite and open(BYTEOFFSET, "> $saveName")) {
		print BYTEOFFSET "totalbytes: $tot\n";
		for $count (sort { $a<=>$b} keys %{ $offset }) {
			print BYTEOFFSET "$offset->{$count}\n";
		}
		close BYTEOFFSET;
	}
	$tot = $counter;
	$count = 0;
	for $counter (keys %{ $selection }) {
		if ($bSelect->{$counter}{start} == -1) {
			delete $bSelect->{$counter};
			delete $selection->{$counter};
		} else {
			if (exists($offset->{ ($counter + 1) })) {
				$bSelect->{$counter}{end} = $offset->{ ($counter + 1) } -1;
			} else {
				$bSelect->{$counter}{end} = $last;
			}
			$count++;
		}
	}

	die "ERROR: No valid frame found in trajectory!\n" if (! $isValid);

	$count = $tot if (! $hasSelections or exists($selection->{all}));
	print "using $count of $tot snapshots...Done\n";
	%{ $selection } = %{ $bSelect };
}

sub getNewFormat {
	my ($fmtStr, $oldFmt) = @_;
	my (@fields, $curr);
	while ($fmtStr =~ /(\w+)/g) {
		$curr = uc $1;
		if ($curr =~ /type/i) {
			push @fields, "TYPE";
		} elsif ($curr =~ /^(x|y|z)(u|s)?/i) {
			push @fields, "${1}COORD";
		} elsif ($curr =~ /^i(x|y|z)/i) {
			push @fields, "${1}INDEX";
		} elsif ($curr =~ /^v(x|y|z)/i) {
			push @fields, "${1}VEL";
		} elsif ($curr =~ /^(s)(x|y|z)(x|y|z)/i) {
			push @fields, "${1}${2}${3}";
		} elsif ($curr =~/q/i) {
			push @fields, "CHARGE";
		} elsif ($curr !~ /^id/i) {
			push @fields, $curr;
		}
	}
	if (! @fields) {
		@fields = @{ $oldFmt };
	}
	return \@fields;
}

sub ParseLAMMPSTrj {
	my ($LOGDATA, $lammpsFile, $SELECT, $trjType, $doAnal, $printStr, $fileHandle) = @_;
	my ($inStr, $tStepData, $fields, $counter, $tot, $field, $filesize, $frame, $totF);
	my ($atomC, $currPos, $start, $strLen, $tStep, %header, $rec, $coords, $i, @vals);
	my ($isTri);

	print "${printStr}Calculating time remaining\r" if (defined($printStr));
	$strLen = length("Calculating time remaining");
	$fields =  getFormat($trjType);
	$totF = scalar(@{ $fields });
	%header = ("TIMESTEP"=>1,"NUMBER OF ATOMS"=>1,"BOX BOUNDS"=>1,"ATOMS"=>1);
	#getByteOffset($SELECT, $lammpsFile);

	open TRAJFILE, $lammpsFile or die "ERROR: Cannot open trajectory file $lammpsFile: $!\n";
	$tot = scalar keys %{ $SELECT };
	$start = time();
	$currPos = $isTri = 0;

	for $i (sort numerically keys %{ $SELECT }) {
		next if ( ! seek(TRAJFILE, $SELECT->{$i}{start}, 0));
		$currPos++;
		while (<TRAJFILE>) {
			chomp;
			$inStr = $_;
			study;
			$isTri = 1 if ($inStr =~ /BOX BOUNDS/ && $inStr =~ /xy xz yz/);
			if ($inStr =~ /^ITEM:\s+(TIMESTEP|NUMBER OF ATOMS|BOX BOUNDS|ATOMS)(.*)$/io) {
				if (exists($header{$1})) {
					$field = $1;
					if ($field eq "TIMESTEP" && keys %{ $tStepData }) {
						$tStepData->{FRAME} = $currPos;
						&fixBox($tStepData->{"BOX BOUNDS"}) if ($isTri);
						$doAnal->($tStepData, $LOGDATA, $fileHandle);
						$strLen = PrintProgress($currPos, $tot, $start, $printStr);
						undef($tStepData);
						last;
					}elsif ($field eq "ATOMS" and defined($2)) {
						$fields = getNewFormat($2, $fields);
						$totF = scalar(@{ $fields });
					}
				} else {
					undef($field);
				}
			} elsif (defined($field) && $inStr =~ /^\s*(\-?\d+\.?\d*e?[\-|\+]?\d*)\s*(\-?\d*.*)$/o) {
				$atomC = $1;
				$coords = $2;
				if ($field eq "BOX BOUNDS") {
					$rec = (
							{
								"lo" => $1,
								"hi" => $2,
							}
							);
					push @{ $tStepData->{$field} }, $rec;
				} elsif ($field =~ /^ATOMS/o) {
					@vals = split /\s+/, $coords;
					$counter = 0;
					while ($counter < $totF) {
						#$vals[$counter] *= 10E8 if ($fields->[$counter] eq "CHARGE");
						$tStepData->{ATOMS}{$atomC}{ $fields->[$counter] } = $vals[$counter];
						$counter++;
					}
					#while ($coords =~ /\s(\-?\d+\.?\d*e?[\-|\+]?\d*)/g && ($counter <= $#{ $fields })) {
						#$tStepData->{ATOMS}{$atomC}{ $fields->[$counter] } = $1;
						#$counter++;
					#}
				} else {
					$tStepData->{$field}[0] = $1;
				}
			}
		}
	}
	close TRAJFILE or die "ERROR: Cannot close $lammpsFile: $!\n";
	if (keys %{ $tStepData }) {
		&fixBox($tStepData->{"BOX BOUNDS"}) if ($isTri);
		$doAnal->($tStepData, $LOGDATA, $fileHandle);
	}
	$i = GetTime(time() - $start);
	printf "$printStr%-${strLen}s\n", "${i}s elapsed..Done" if (defined($printStr));
}

sub fixBox {
	my ($box) = $_[0];
	my ($xlo,$xhi,$xy,$ylo,$yhi,$xz,$zlo,$zhi,$yz);
	my ($lx,$ly,$lz,$a,$b,$c,$alpha,$beta,$gamma);

	$xlo = $box->[0]{lo};
	($xhi,$xy) = split / /,$box->[0]{hi};
	$lx = $xhi - $xlo;
	$ylo = $box->[1]{lo};
	($yhi,$xz) = split / /,$box->[1]{hi};
	$ly = $yhi - $ylo;
	$zlo = $box->[2]{lo};
	($zhi,$yz) = split / /,$box->[2]{hi};
	$lz = $zhi - $zlo;

	$a = $lx;
	$b = sqrt($ly*$ly+$xy*$xy);
	$c = sqrt($lz*$lz+$xz*$xz+$yz*$yz);
	$alpha = rad2deg(acos(($xy*$xz+$ly*$yz)/$b/$c));
	$beta = rad2deg(acos($xz/$c));
	$gamma = rad2deg(acos($xy/$b));

	$box->[0]{lo} = $box->[1]{lo} = $box->[2]{lo} = 0;
	$box->[0]{hi} = $a; $box->[1]{hi} = $b; $box->[2]{hi}=$c;
	$box->[3]{alpha} = $alpha; $box->[3]{gamma} = $gamma; $box->[3]{beta} = $beta;
}

sub ParseLAMMPSTrj_new {
	my ($LOGDATA, $lammpsFile, $SELECT, $trjType, $doAnal, $printStr, $fileHandle) = @_;
	my ($strLen, $fields, $totF, $tot, $start, $currPos, $headerPtn, $aC, $field, $totL);
	my ($i, $j, $k, $tStepData, $bLen, $buf, $counter, @tmp, @vals, $atomC, $aT, $vArrayTot);

	print "${printStr}Calculating time remaining\r" if (defined($printStr));
	$strLen = length("Calculating time remaining");
	$fields =  getFormat($trjType);
	$totF = scalar(@{ $fields });
	open TRAJFILE, $lammpsFile or die "ERROR: Cannot open trajectory file $lammpsFile: $!\n";
	$tot = scalar keys %{ $SELECT };
	$start = time();
	$currPos = $vArrayTot = 0;
	$headerPtn = qr/^ITEM: TIMESTEP\n(\d+)\nITEM: NUMBER OF ATOMS\n(\d+)\nITEM: BOX BOUNDS\n([^\n]+)\n([^\n]+)\n([^\n]+)\nITEM: ATOMS(?:[^\n]*)\n/;
	$field = qr/(\S+)/;
	for $i (sort numerically keys %{ $SELECT }) {
		$tStepData = ();
		$bLen = $SELECT->{$i}{end} - $SELECT->{$i}{start};
		next if ( ! sysseek(TRAJFILE, $SELECT->{$i}{start}, 0));		
		next if ( ! sysread(TRAJFILE, $buf, $bLen,0));
		next if ($buf !~ /$headerPtn/g);
		$currPos++;
		$tStepData->{TIMESTEP}[0] = $1;
		$tStepData->{"NUMBER OF ATOMS"}[0] = $2;
		$aT = $2;
		$counter = 0;
		for $i ($3, $4, $5) {
			 @tmp = split /\s+/, $i;
			$tStepData->{"BOX BOUNDS"}[$counter]{lo} = $tmp[0];
			$tStepData->{"BOX BOUNDS"}[$counter]{hi} = $tmp[1];
			$counter++;
		}
		$counter = -1;
		$aC = 1;
		if (! $vArrayTot) {
			$buf =~ /(.+)\n/g;
			@tmp = split /\s+/, $1;
			$totL = scalar(@tmp); # total number of fields per line
			@vals = (@tmp, split /\s+/, $');
			$vArrayTot = scalar(@vals) - $totL;
			$k = $totF;
			$k = $totL if ($totL < $k);
		} else {
			@vals = split /\s+/, $';
		}

		$j = 0;
		while ($j <= $vArrayTot) {
			$aC++;
			#last if ($aC > $aT);
			$atomC = $vals[$j];
			$counter = 0;
			while ($counter < $k) {
				$tStepData->{ATOMS}{$atomC}{ $fields->[$counter] } = $vals[($j + $counter + 1)];
				$counter++;
			}
			$j += $totL;
		}
		$tStepData->{FRAME} = $currPos;
		$doAnal->($tStepData, $LOGDATA, $fileHandle);
		$strLen = PrintProgress($currPos, $tot, $start, $printStr);
	}
	close TRAJFILE or die "ERROR: Cannot close $lammpsFile: $!\n";

	$i = GetTime(time() - $start);
	printf "$printStr%-${strLen}s\n", "${i}s elapsed..Done" if (defined($printStr));
}

sub CreateInputFile {
	my ($parms) = $_[0];
	my ($header, $middle, $tail, $fileName, $i, $tmp, $line); 
	my ($qeqDump, $electrodeQStr, $electrodeStr, $cmapStr, $pStr);

	$polarizable = 0; 
	$rigidMols = "none";
	$flexMols = "all";
	$cmapStr = "";
	$cmapStr = "thermo_style    custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong f_CMAP press vol\nthermo_modify   line multi\n" 
		if(exists($parms->{PARMS}{CMAP}));
	$polarizable = 1 if(exists($parms->{PARMS}{POLARIZATION}) and keys %{ $parms->{PARMS}{POLARIZATION} });
	($electrodeQStr, $electrodeStr) = setElectrodeOpts($parms);
	$qeqDump = "";
	if(exists($parms->{QEq}) and exists($parms->{QEq}{Opt}) and $parms->{QEq}{Opt} > 0) {
		$qeqDump = <<"DATAdump";
dump                 1 all custom 1 \${sname}.lammps id type xu yu zu q d2_rs[1] d2_rs[2] d2_rs[3]
dump_modify          1 first yes sort id
DATAdump

	}
	#polarizable option
	$pStr = "";
	if($polarizable and $parms->{PARMS}{POLARIZATION}{type}>0) {
		$pStr = <<"DATA";
thermo 10
print                          .
print ==========================================
print "SHELL MINIMIZATION"
print ==========================================
print    				        .
dump                 1 all custom 10 \${sname}.shell-min.lammps id type xu yu zu
fix                  shells_only_move ATOMS spring/self 500000000
min_style            cg
minimize             1.0e-4 1.0e-4 500 5000
unfix                shells_only_move
undump               1
DATA
	}
	$qeqDump =~ s/ d2_rs.*// if ($polarizable); 
	$fileName = "in." . $parms->{PARMS}{SUFFIX};
	$header = createInputHeader($parms, $parms->{PARMS}{FFTYPE}, $parms->{PARMS}{NUM_FILES});
	$middle = createMiddle($parms, $parms->{PARMS}{SOLUTE}, $parms->{PARMS}{SUFFIX});
	$tail = readTail($parms->{PARMS}{SOLUTE}, $parms);
	$tail =~ s/qeqElectrodeOpts/$electrodeStr/;
	$tail =~ s/qElectrodeOpts/$electrodeQStr/;
	open LAMMPSINPUT, "> $fileName" || die "ERROR: Cannot create $fileName:$!\n";
	@{ $tmp } = split /\n/,"$header\n$middle\n${cmapStr}${qeqDump}\n${pStr}\n$tail";
	for $i (@{ $tmp }) {
		@{ $line } = split /\s+/, $i;
		printf LAMMPSINPUT "%-20s %s\n",shift @{ $line }, "@{ $line }";
	}
	close LAMMPSINPUT;
	open SINGLEPOINT, "> ${fileName}_singlepoint" || die "ERROR: Cannot create ${fileName}_singlepoint:$!\n";
	@{ $tmp } = split /\n/,"$header\n$middle\n${cmapStr}${qeqDump}\n${pStr}\nrun 0";
	for $i (@{ $tmp }) {
		@{ $line } = split /\s+/, $i;
		printf SINGLEPOINT "%-20s %s\n",shift @{ $line }, "@{ $line }";
	}
	close SINGLEPOINT;
}

sub createMiddle {
	my ($parms, $soluteAtms, $suffix) = @_;
	my ($middle, $numAtms, $i, $j, $totAtms, $tmp);
	
	$parms->{PARMS}{OFF_DIAG} = "$parms->{PARMS}{PAIR_COEFFS}\n$parms->{PARMS}{OFF_DIAG}" 
		if (exists($parms->{PARMS}{PAIR_COEFFS}));
	$middle = "";
	if(exists($parms->{PARMS}{ALL_TINKER}) and $parms->{PARMS}{ALL_TINKER} == 1) {
		$middle .= <<DATA;
# per-atom properties required by AMOEBA or HIPPO
fix             amtype all property/atom i_amtype ghost yes
fix             extra all property/atom i_amgroup d_redID d_pval ghost yes
fix             extra2 all property/atom i_polaxe d2_xyzaxis 3

DATA

		if(exists($parms->{PI_TORSIONS}) and $parms->{PI_TORSIONS}{counter}> 0) {
			$middle .= "fix             pit all amoeba/pitorsion\nfix_modify      pit energy yes";
		}
		if(exists($parms->{BI_TORSIONS}) and $parms->{BI_TORSIONS}{counter}> 0) {
			$middle .= "fix             bit all amoeba/bitorsion $parms->{PARMS}{BI_TORSIONS}{FILE}\nfix_modify      bit energy yes\n";
		}
	}
	$middle .= "fix            CMAP all cmap $parms->{PARMS}{CMAP}{FILE}\n" if(exists($parms->{PARMS}{CMAP}));	
	$middle .= "read_data	   data.${suffix}";
	$middle .= " extra/special/per/atom 8" if ($polarizable and $parms->{PARMS}{POLARIZATION}{type}> 0);
	$middle .= " fix Isotopes NULL Isotopes" if(exists($parms->{PARMS}{OPTIONS}{ISOTOPES}));
	$middle .= " fix shells NULL \"Shell Atoms\"" if($parms->{QEq}{Opt} == 2 and ! exists($parms->{PARMS}{POLARIZATION}{type}));
	$middle .= " fix amtype NULL \"Tinker Types\"" if(exists($parms->{PARMS}{ALL_TINKER}) and $parms->{PARMS}{ALL_TINKER} == 1);
	$middle .= "&\nfix pit \"pitorsion types\" \"PiTorsion Coeffs\" &\nfix pit pitorsions PiTorsions"
				if(exists($parms->{PI_TORSIONS}) and $parms->{PI_TORSIONS}{counter}> 0);
	$middle .= "&\nfix bit bitorsions BiTorsions\n"
				if(exists($parms->{BI_TORSIONS}) and $parms->{BI_TORSIONS}{counter}> 0);			
	$middle .= "\ninclude $parms->{PARMS}{include_file}\n" if (exists($parms->{PARMS}{include_file}));
	$middle .= $parms->{PARMS}{kspace_str} if($parms->{PARMS}{has_kspace} and $parms->{PARMS}{isTriclinic} == 1);

	if(exists($parms->{PARMS}{CMAP})) {
		$middle .= " fix CMAP crossterm CMAP" ;
		$middle .= "\nfix_modify     CMAP energy yes";
	}
	$middle .= "\n\n$parms->{PARMS}{OFF_DIAG}";
	$middle .= "pair_modify	 mix $parms->{PARMS}{mix_rule}\n" if (exists($parms->{PARMS}{mix_rule}));
	if($parms->{QEq}{Opt} == 2) {
		$middle .= "neighbor		2.0 bin\n";
	} else {
		$middle .= "neighbor		2.0 multi\n";
	}
	if (! $polarizable) {
		$middle .= "neigh_modify	every 2 delay 4 check yes\n";
	} else {
		$middle .= "neigh_modify	every 2 delay 4 check yes one 10000\n";
	}
	$middle .= "thermo_style	multi\nthermo_modify		line multi format float %14.6f flush yes\n";
	$parms->{HYBRID_VALENCE} = "" if (! exists($parms->{HYBRID_VALENCE}));
	$middle .= "$parms->{HYBRID_VALENCE}";
	$middle .= "variable		input string in.${suffix}\n";
	$middle .= "variable		sname string $suffix\n";
	$numAtms = scalar keys %{ $soluteAtms };
	$totAtms = $parms->{PARMS}{NUM_ATOMS};
	if ($numAtms > 0 and $numAtms < $totAtms) {
		$middle .= <<DATA;
variable		sAtoms index $numAtms
group		   solute id <> 1 \${sAtoms}
group		   solvent subtract all solute
DATA
	}
	$middle .= "\n$parms->{PARMS}{POLARIZATION}{input}\n" if ($polarizable);
	$middle .= "\n$parms->{MESODNA}" if ($parms->{MESODNA});
	if($parms->{PARMS}{USE_HBOND}) {
		$middle .= "\ncompute   hb all pair hbond/dreiding/lj\n" if ($parms->{PARMS}{USE_HBOND} == 1);
		$middle .= "\ncompute   hb all pair hbond/dreiding/morse\n" if ($parms->{PARMS}{USE_HBOND} == 2);
		$middle .= <<DATA;
variable	nHB equal c_hb[1]
variable	eHB equal c_hb[2]
thermo_style		 custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong v_eHB v_nHB press vol
thermo_modify   line multi flush yes
DATA
	}
	if(exists($parms->{PARMS}{REAX})) {
		$middle .= "\nfix		charge all qeq/reax 1 0.0 10.0 1.0e-6 reaxff\n";
	} elsif($parms->{QEq}{Opt} == 1) {
		$middle .= "\nfix		charge all qeq/rg 1 0.0 10.0 1.0e-06 ${suffix}.param.qeq\n";
	} elsif($parms->{QEq}{Opt} == 2 and ! $polarizable) {
		$middle .= "\nfix		pqeq all pqeq method 2 nevery 1 charge $parms->{QEq}{sys_charge} tolerance 1.0e-6 damp 1.0\n";
	} elsif($parms->{QEq}{Opt} == 2 and $polarizable) {
		#drude induced dipole with QEq
		$middle .= "fix	pqeq CORES pqeq method 0 nevery 1 charge $parms->{QEq}{sys_charge} tolerance 1.0e-6 damp 1.0 drude yes\n";
	}
	if(exists($parms->{PARMS}{OPTIONS}{ISOTOPES})) {
		#$middle .= $parms->{PARMS}{OPTIONS}{ISOTOPES}{OUTSTR};
		print "";
	}
	$parms->{PARMS}{OPTIONS}{MD_TIMESTEP} /=2 if ($polarizable);
	$middle .= "\n$parms->{PARMS}{TYPE_CHARGES}\ntimestep $parms->{PARMS}{OPTIONS}{MD_TIMESTEP}\n";
	if(exists($parms->{PARMS}{RIGID}) and keys %{ $parms->{PARMS}{RIGID}{MOLS}}) {
		$rigidType = "molecule";
		$rigidType = "single" if (scalar(keys %{ $parms->{PARMS}{MOLS} } == 1));
		#here we write rigid molecules id
		if(exists($parms->{PARMS}{RIGID}{MOLS}{all})) {
			$rigidMols = "all";
			$rigidMols = "ATOMS" if($polarizable and $parms->{PARMS}{POLARIZATION}{type}>0);	
			$flexMols = "";	
			#search to see any of the rigid bodies are larger than 10 atoms (treshold for "small")
			$rigidOpt = "";
			undef $tmp;
			$tmp = 0;
			map { $tmp++ if(scalar(keys %{ $parms->{PARMS}{MOLS}{$_} }) > 9) } keys %{ $parms->{PARMS}{MOLS} }; 
			$rigidOpt = "/small" if (! $tmp);
		} else {
			#search to see any of the rigid bodies are larger than 10 atoms (treshold for "small")
			$rigidOpt = "";
			undef $tmp;
			$tmp = 0;
			map { $tmp++ if(scalar(keys %{ $parms->{PARMS}{RIGID}{MOLS}{$_} }) > 9) } keys %{ $parms->{PARMS}{RIGID}{MOLS} }; 
			$rigidOpt = "/small" if (! $tmp);
			#now see if a larger number of atoms are rigid or flexible
			$tmp = scalar(keys %{ $parms->{PARMS}{RIGID}{ATOMS} });
			$parms->{PARMS}{majorRigid} = 1;
			$parms->{PARMS}{majorRigid} = 0 if($tmp < $parms->{PARMS}{NUM_ATOMS}/2);
			undef $tmp;

			@{ $tmp } = sort numerically keys %{ $parms->{PARMS}{RIGID}{MOLS} };
			if($polarizable and $parms->{PARMS}{POLARIZATION}{type}>0) { 
				if($parms->{PARMS}{POLARIZATION}{rigidOpt} ne "ATOMS") {
					$middle .= "\n" . getRange($tmp, "molecule", "rtmp") . 
								"\ngroup RIGID subtract rtmp SHELLS" . 
								"\ngroup FLEX subtract ATOMS RIGID\n";
				} else {
					$middle .= "\n" . getRange($tmp, "molecule", "RIGID") . "\ngroup FLEX subtract ATOMS RIGID\n";
				}
			} else {
				$middle .= "\n" . getRange($tmp, "molecule", "RIGID") . "\ngroup FLEX subtract all RIGID\n";
			}
			$rigidMols = "RIGID";
			$flexMols = "FLEX";				
		}
		#				$middle .= <<DATA;
		#neigh_modify exclude molecule/intra $rigidMols
		#delete_bonds $rigidMols multi
		#run 10			
		#DATA
	}
	return $middle;
}

sub getRange {
	my ($molList, $itype, $vName) = @_;
	my ($sorted_list, $i, $j, $rStr, $range, $fStr);

	@{ $sorted_list } = sort numerically @{ $molList };
	return "group $vName $itype $sorted_list->[0]" if ($#{ $sorted_list } == 0);
	$i = 0;
	$j = 1;
	while (($i+$j) <= $#{ $sorted_list }) {
		if (($sorted_list->[$i+$j] - $sorted_list->[$i]) <= $j) {
			$j++;
			next;
		}
		if (($sorted_list->[$i]-$sorted_list->[$i+$j-1]) == 0) {
			push @{ $range }, "$sorted_list->[$i] ";
		} else {
			push @{ $range }, "<> $sorted_list->[$i] $sorted_list->[$i+$j-1] ";
		}
		$i += $j;
		$j = 1;
	}
	if (($i-$j) == 1) {
		push @{ $range }, "$sorted_list->[$i]";
	} else {
		push @{ $range }, "<> $sorted_list->[$i] $sorted_list->[$i+$j-1]";
	}
	foreach $i (0 .. $#{ $range }) {
		$rStr .= "group tmp${i} $itype $range->[$i]\n";
		$fStr .= " tmp${i}";
	}
	if($i>0) {
		$rStr .= "group $vName union $fStr\n";
	} else {
		$rStr .= "group $vName union tmp0\n";
	}
	return $rStr;

}

sub readTail {
	my ($soluteAtms, $parms) = @_;
	my ($hasSolute, $outStr, $inStr, $inputFile, $shakeOpts, $solAtms, $numAtms, $fepStr);
	my ($tmpVals, $unfixShellid, $sstr, $unfixList, $pStr, $polarBondUpdate);

	$hasSolute = 1;
	if (keys %{ $soluteAtms }) {
		$solAtms = scalar(keys %{ $soluteAtms });
		$numAtms = $parms->{PARMS}{NUM_ATOMS};
		$hasSolute = 0 if ($numAtms == $solAtms);
	} else {
		$hasSolute = 0;
	}

	if(exists($parms->{PARMS}{FEP})) {
		if (! exists($parms->{PARMS}{FEP}{var_list})) {
			$parms->{PARMS}{FEP}{var_list} = "scaling";
		}
		$fepStr = <<"DATA";

#FEP
$parms->{PARMS}{FEP}{var_str}

fix                  ADAPT all adapt/fep 100000 &
$parms->{PARMS}{FEP}{fix_adapt_str} $parms->{PARMS}{FEP}{atom_fix_str}                     after yes

compute              FEP all fep \${rtemp} &					 
$parms->{PARMS}{FEP}{compute_fep_str} $parms->{PARMS}{FEP}{atom_compute_str}                      volume no

fix                  PRINT all print 100000 "FEP variables\\n & 
                                            lambda = \${lambda}\\n & 
$parms->{PARMS}{FEP}{var_list}"
DATA
	}

	foreach $inputFile (@{ $parms->{PARMS}{INPUTLOC} }) {
		$polarBondUpdate = 0;
		$shakeOpts = "";
		if($parms->{PARMS}{OPTIONS}{SHAKE}) {
			if (exists($parms->{PARMS}{SHAKE_MASS}) and $parms->{PARMS}{SHAKE_MASS} ne "") {
				$shakeOpts .= " m $parms->{PARMS}{SHAKE_MASS}";
			} elsif (exists($parms->{PARMS}{SHAKE_BOND}) and $parms->{PARMS}{SHAKE_BOND} ne "") {
				$shakeOpts .= " b $parms->{PARMS}{SHAKE_BOND}";
			}
			if ($shakeOpts and exists($parms->{PARMS}{SHAKE_ANGLE})) {
				$shakeOpts .= " $parms->{PARMS}{SHAKE_ANGLE}";
			}
		}
		$rigidNVE = 0;		
		open TAIL, $inputFile || die "ERROR: Cannot open $inputFile: $!\n";
		while (<TAIL>) {
			chomp;
			$inStr = $_;
			if (! $hasSolute) {
				next if $inStr =~ /solute|restraint/;
				$inStr =~ s/solvent/all/g;
			}
			#polarization options
			if($polarizable) {
			}
			#rigid options
			if($rigidMols ne "none") {
				$inStr = "" if ($inStr =~ /shakeOpts/i and $flexMols eq "");
				$inStr = "" if ($inStr =~ /(restraint|minimize)/i);
				$inStr =~ s/nreset\s+\d+//i;
				$inStr =~ s/(t|p)loop\s+\d+//gi;
				if($inStr =~ /^\s*fix\s+(\S+)\s+(\S+)\s+(n.t|nph|nve)\s*(.*)$/i) {
					if ($3 eq "nve") {
						$inStr = "fix $1 $rigidMols rigid/${3}${rigidOpt} $rigidType";
						$rigidNVE = 1;
						$inStr .= "\nfix ${1}_flex $flexMols $3 $4" if ($flexMols ne "");
					} elsif ($3 eq "npt" and $flexMols ne "") {
						$tmpVals = ();
						push @{ $tmpVals }, ($1,$2,$3,$4);
						if($parms->{PARMS}{majorRigid} == 1) {
							#most of the atoms are rigid, so invoke npt for rigid...
							$inStr = "fix $1 $rigidMols rigid/$tmpVals->[2]${rigidOpt} $rigidType $tmpVals->[3] dilate all";
							# and nvt for flexible
							$inStr .= "\nfix $tmpVals->[0]_flex $flexMols nvt";
							$tmpVals->[3] =~ s/^(.*) (iso|aniso|x|y|z)\S?\s+\d+.*//;
							$inStr .= " $1";
						} else {
							$inStr = "fix $tmpVals->[0]_flex $flexMols $tmpVals->[2] $tmpVals->[3]";
							$tmpVals->[2] =~ s/^(.*) (iso|aniso|x|y|z)\S?\s+\d+.*//;
							$inStr .= "\nfix $tmpVals->[0] $rigidMols rigid/nvt${rigidOpt} $rigidType $1";
						}
					} else {
						$inStr = "fix $1 $rigidMols rigid/${3}${rigidOpt} $rigidType $4";
						$inStr .= "\nfix ${1}_flex $flexMols $3 $4" if ($flexMols ne "");
					}
					$unfixList->{$1} = 1;
				} elsif ($inStr =~ /^\s+fix\s+\S+\s+\S+\s+langevin (.*)$/ and $rigidNVE == 1) {
					$outStr =~ s/\R$//;
					$inStr = "langevin $1";
					$rigidNVE = 0;
				} elsif ($inStr =~ /^\s+fix\s+\S+\s+\S+\s+langevin (.*)$/) {
					$inStr = "";
				} elsif ($inStr =~ /^unfix\s+(\S+)/ and $flexMols ne "" and exists($unfixList->{$1})) {
					$inStr .= "\nunfix ${1}_flex";
					delete $unfixList->{$1};
				}
			}
			if ($inStr =~ /shakeOpts/) {
				$inStr =~ s/ all / $parms->{PARMS}{OPTIONS}{SHAKE_WHAT} / if($hasSolute && exists($parms->{PARMS}{OPTIONS}{SHAKE_WHAT}));
				$inStr =~ s/ all / $flexMols / if($rigidMols ne "none" and $flexMols == "");
				if ($shakeOpts ne "") {
					$inStr =~ s/shakeOpts/$shakeOpts/;
				 } else {
					$inStr = "";
				 }
			} elsif ($inStr =~ /shakeH/ and $shakeOpts eq "") {
				$inStr = "";
			}
			$inStr = "thermo_style		custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong v_E_hbond v_n_hbond press vol" 
				if ($inStr =~ /^thermo_style/ and $parms->{PARMS}{USE_HBOND});
			$inStr =~ s/all npt.*$/all npt temp 300.0 300.0 100.0 x 1.0 1.0 2000.0 y 1.0 1.0 2000.0 couple xy/ 
				if ($parms->{PARMS}{OPTIONS}{PERIODICITY} == 2);
			$inStr =~ s/wall_str/fix zconfine all wall\/harmonic zlo EDGE 10.0 1.0 2.5 zhi EDGE 10.0 1.0 2.5 units box/
				if ($parms->{PARMS}{OPTIONS}{PERIODICITY} == 2);
			$inStr =~ s/wall_str// if($parms->{PARMS}{OPTIONS}{PERIODICITY} != 2);
			$inStr = "" if ($inStr =~ /unfix\s+shakeH/ && $shakeOpts eq "");
			$inStr .= " q" if($inStr =~ /dump.*custom/ and $parms->{QEq}{Opt}>0);
			#polarization options
			if($polarizable) {
				$inStr =~ s/^velocity\s*all\s*create.*/velocity ATOMS create \${rtemp} \${seed} dist gaussian/;
				if($inStr =~ /^\s*fix\s+(\S+).* (rigid.n.t\S+|temp.berendsen|temp.rescale|temp.csvr|n.t|langevin) /i) {
					$inStr .= "\nfix_modify $1 temp TATOMS";
				} elsif ($inStr =~ /^\s*fix\s+(\S+).* (rigid.npt\S+|npt|nph|press.berendsen) /i) {
					$inStr .= "\nfix_modify $1 press thermo_press_lmp";
				} elsif (($parms->{PARMS}{POLARIZATION}{type}==0) and $inStr =~ /^\s*velocity.*create/) {
					$inStr .= " bias yes temp TATOMS";
				} elsif ($inStr =~ / rsx/) {
					$inStr =~ s/ rsx.*//;
				}
				if ($parms->{PARMS}{POLARIZATION}{type}>0) {
					#drude particles
					if ($inStr =~ /^(.*) (langevin) (.*)/) {
						$inStr = "$1 langevin/drude $2 1.00 20 13977 zero yes";
					}
					#$inStr =~ s/ all / ATOMS /i if ($inStr =~ /^\s*fix.*all\s+shake/);
					$inStr .= " v_Tshell v_Tcore" if ($inStr =~ /^\s*thermo_style\s+custom/i);
					$sstr = "all|ATOMS";
					$sstr = $parms->{PARMS}{POLARIZATION}{rigidOpt}
						if(exists($parms->{PARMS}{POLARIZATION}{rigidOpt}));
					if($inStr =~ /fix\s+(\S+)\s+($sstr)\s+(rigid.npt|rigid.nvt|npt|nvt)/) {
						$tmpVals = ();
						push @{ $tmpVals }, ($1,$2, $3);
						$inStr =~ s/ all / ATOMS /;
						$inStr = "fix DIRECT all drude/transform/direct\n${inStr}";
						$inStr .= "\nfix $tmpVals->[0]_shells SHELLS nvt temp \${drude_T} \${drude_T} 20";
						$inStr .= "\nfix INVERSE all drude/transform/inverse";
						$unfixShellid = $tmpVals->[0];
					} elsif($inStr =~ /fix\s+(\S+)\s+($sstr)\s+(rigid|nve)/) {
						$tmpVals = ();
						push @{ $tmpVals }, ($1,$2, $3);
						$inStr =~ s/all/ATOMS/;
						$inStr .= "\nfix $tmpVals->[0]_shells SHELLS nve";
						$unfixShellid = $tmpVals->[0];
					} elsif(defined($unfixShellid) and $inStr =~ /unfix\s+${unfixShellid}/) {
						$inStr .= "\nunfix ${unfixShellid}_shells\n";
						$inStr .= "unfix DIRECT\nunfix INVERSE";
						undef $unfixShellid;
					}
					if($inStr =~ /velocity ATOMS/) {
						$polarBondUpdate = 1;
						$inStr .= updateCoreShellBondK($parms->{BONDS}, 2);
					}elsif($polarBondUpdate and $inStr =~ /^\s*run/) {
						$polarBondUpdate = 0;
						$inStr .= updateCoreShellBondK($parms->{BONDS}, 1) . $inStr;
					}
				}
				$inStr =~ s/ temp\s+1.0\s+\$\{rtemp\}/ temp \$\{rtemp\} \$\{rtemp\}/g;
			}
			$outStr .= "$inStr\n";
			$outStr =~ s/fep_str/$fepStr/;
		}

		close TAIL;
	}
	return $outStr;
}

sub updateCoreShellBondK {
	my ($bondParms, $sf) = @_;
	my ($i, $j, $k, $csData, $fc, $hybrid, $outStr);
	
	for $i (keys %{ $bondParms }) {
		next if ($i =~ /TYPE|counter/);
		for $j (values %{ $bondParms->{$i} }) {
			for $k (values %{ $j }) {
				next if (! exists($k->{CORE_SHELL}));
				$fc = $k->{VALS}[0] * $sf;
				$csData->{ $k->{Lammps}{name} }{ $k->{INDEX} } = "$fc $k->{VALS}[1]";
			}
		}
	}

	$hybrid = 0;
	$hybrid = 1 if(scalar(keys %{ $csData }) > 1);
	for $i (keys %{ $csData }) {
		for $j (sort numerically keys %{ $csData->{$i} }) {
			$outStr .= "bond_coeff ";
			$outStr .= "$i " if ($hybrid);
			$outStr .= "$j $csData->{$i}{$j}\n";
		}
	}

	return "\n$outStr";
}

sub createInputHeader {
	my ($parms, $ffType, $numFiles) = @_;
	my ($header, $i, $parm, $j, $tmp, $hasCharmmTorsions, $hasCharmmLJ);

	#figure out if the torsions are of type charmm
	@{ $tmp } = keys %{ $parms->{TORSIONS}{TYPE} }; 
	$hasCharmmTorsions = $hasCharmmLJ = 0;
	$hasCharmmTorsions = 1 if ("@{ $tmp }" =~ /charmm/);

	$parms->{PARMS}{dielectric} += 0;
	$header = "units		   real\n";
#	if($parms->{QEq}{Opt} == 2 and ! exists($parms->{PARMS}{POLARIZATION}{type})) {
#		$header .= "atom_style	  pqeq\n";
#	} elsif (exists($parms->{PARMS}{ALL_TINKER}) and $parms->{PARMS}{ALL_TINKER} == 1) {
	if (exists($parms->{PARMS}{ALL_TINKER}) and $parms->{PARMS}{ALL_TINKER} == 1) {
		$header .= "atom_style	  amoeba\n";
	} else {
		$header .= "atom_style	  full\n";
	}
	if($parms->{PARMS}{OPTIONS}{PERIODICITY} == 0) { 
		$header .= "boundary		f f f\n";
	}elsif ($parms->{PARMS}{OPTIONS}{PERIODICITY} == 2) {
		$header .= "boundary		p p f\n";
	}else {
		$header .= "boundary		p p p\n";
	}
	if(exists($parms->{PARMS}{ALL_TINKER}) and $parms->{PARMS}{ALL_TINKER} == 1) {
		$header .= "atom_modify     sort 0 0.0\n";
	}
	$header .= "dielectric	  $parms->{PARMS}{dielectric}\n";
	#$header .= "newton		  off\n" if ($parms->{PARMS}{HAS_HBONDS});

	#if ($parms->{PARMS}{isTriclinic}) {
		#$header .= "special_bonds   lj/coul 0 0 0.5\n\n";
		#$header .= "pair_style	  lj/coul long long 10.0";
#	} elsif ($ffType == 1) { #amber
	if (exists($parms->{PARMS}{ALL_TINKER}) and $parms->{PARMS}{ALL_TINKER} == 1) {
		$header .= "special_bonds   lj/coul 0.5 0.5 0.5 one/five yes\n\n";
		$header .= "pair_style amoeba";
	} elsif ($ffType == 1 and ! $parms->{PARMS}{is_tip4p}) {
		$header .= "special_bonds   amber\n\n";
		$header .= "pair_style	  lj/charmm/coul/long/opt 9.0 10.0";
	} elsif ($ffType == 2 && ! $parms->{PARMS}{is_tip4p} && scalar(keys %{ $parms->{VDW}{TYPE} }) == 1) { #charmm
		$header .= "special_bonds   charmm\n\n";
		$hasCharmmLJ = 1;
		if($parms->{PARMS}{OPTIONS}{PERIODICITY} == 0) {
			$header .= "pair_style	  lj/charmmfsw/coul/charmmfsh 10 12";
		} else {
			$header .= "pair_style	  lj/charmmfsw/coul/long 10 12";
		}
#	} elsif ($numFiles == 1 && $ffType == 4) { #dreiding
#		$header .= "special_bonds   0.0 0.0 0.5\n\n";
#		$header .= "pair_style	  hybrid/overlay dreiding/hb 5.0 90 lj/cut/coul/long " . 
#					"$parms->{PARMS}{cut_coul} $parms->{PARMS}{cut_vdw}";
		#} elsif ($ffType == 3) { #MESODNA
		#$header .= "special_bonds   amber\n\n";
		#$header .= "pair_style	  hybrid yukawa 3.0 6.0  dreiding/hb 6.0 120.0 morse/opt 12.0";
	} elsif ($ffType == 5) { #reax
		#$header .= "\npair_style	  reax 10.0 1.0e-5";
		$header .= "\npair_style	  reaxff NULL";
	} elsif ($ffType == 6  and ! $parms->{PARMS}{is_tip4p}) { #opls
		$header .= "special_bonds   lj/coul 0.0 0.0 0.5\n\n";
		$header .= "pair_style	  lj/charmm/coul/long/opt 9.0 10.0";
		$header .= " coul 1.0 1.0 1.0 ";
	} else {
		if ($parms->{QEq}{Opt} == 2) {
			$tmp = " coul ";
			$header .= "special_bonds lj ";
			for $i (2 .. 4) {
				$parms->{PARMS}{"scale_vdw_1${i}"} = 0.00000001
					if($parms->{PARMS}{"scale_vdw_1${i}"} == 0);
				$parms->{PARMS}{"scale_cou_1${i}"} = 0.00000001
					if($parms->{PARMS}{"scale_cou_1${i}"} == 0);
				$header .= sprintf("%s ", $parms->{PARMS}{"scale_vdw_1${i}"});
				$tmp .= sprintf("%s ", $parms->{PARMS}{"scale_cou_1${i}"});
			}
			$header .= $tmp;
		} else {	
			if ($ffType == 2 && ! $parms->{PARMS}{is_tip4p}) {
				$header .= "special_bonds   charmm";
				$hasCharmmLJ = 1;
			} else {
				$header .= "special_bonds   lj";
				$header .= "/coul" if ($parms->{PARMS}{same_scale});
				$header .= " ";
				for $i (2 .. 4) {
					$header .= sprintf("%2.1f ", $parms->{PARMS}{"scale_cou_1${i}"});
				}
				if (! $parms->{PARMS}{same_scale}) {
					$header .= "coul ";
					for $i (2 .. 4) {
						$header .= sprintf("%2.1f ", $parms->{PARMS}{"scale_vdw_1${i}"});
					}
				}
			}
		}
		$header .= "\n\npair_style	  ";
		if (scalar keys %{ $parms->{VDW}{TYPE} } > 1) {
			if ($parms->{PARMS}{HYBRID_VDW} == 2) {
				$header .= "hybrid/overlay ";
			} else {
				$header .= "hybrid ";
			}
		}elsif (scalar keys %{ $parms->{VDW}{TYPE } } == 0) {
			$header .= "none";
		}
		for $j (sort { $a cmp $b } keys %{ $parms->{VDW}{TYPE} }) {
			if($hasCharmmTorsions or $hasCharmmLJ) {
				$parms->{VDW}{TYPE}{$j}{LMPNAME} =~ s/lj\/charmm\/coul\/long.*/lj\/charmmfsw\/coul\/long/;
				$parms->{VDW}{TYPE}{$j}{LMPNAME} =~ s/lj\/charmm\/coul\/charmm.*/lj\/charmmfsw\/coul\/charmmfsh/;
				if(exists($parms->{PARMS}{OFF_DIAG})) {
					$parms->{PARMS}{OFF_DIAG} =~ s/lj\/charmm\/coul\/long\/opt/lj\/charmmfsw\/coul\/long/g;
					$parms->{PARMS}{OFF_DIAG} =~ s/lj\/charmm\/coul\/long/lj\/charmmfsw\/coul\/long/g;
					$parms->{PARMS}{OFF_DIAG} =~ s/lj\/charmm\/coul\/charmm/lj\/charmmfsw\/coul\/charmmfsh/g;
				}
			}
			$header .= "$parms->{VDW}{TYPE}{$j}{LMPNAME} $parms->{VDW}{TYPE}{$j}{OPTS} ";
		}
	}
	
	for $i ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS") {
		$parm = lc(substr($i, 0, -1));
		$parm = "dihedral" if ($parm eq "torsion");
		$parm = "improper" if ($parm eq "inversion");
		$header .= sprintf("\n%-16s", $parm . "_style");
		if (scalar keys %{ $parms->{$i}{TYPE} } > 1) {
			$header .= "hybrid ";
		} elsif (scalar keys %{ $parms->{$i}{TYPE } } == 0) {
			$header .= "none";
		}
		for $j (sort { $a cmp $b } keys %{ $parms->{$i}{TYPE} }) {
			if ($i eq "TORSIONS" and $j =~ /charmm/) {
				$header  .= "charmmfsw ";
			} else {
				$header .= "$j ";
			}
		}
	}
	$parms->{PARMS}{has_kspace} = 1;
	if ($parms->{PARMS}{isTriclinic} and $parms->{PARMS}{OPTIONS}{PERIODICITY} and $header =~ /coul.*long/ and ($parms->{PARMS}{HAS_CHARGE} or $parms->{PARMS}{OPTIONS}{VDW_EWALD}) ) {
		$parms->{PARMS}{kspace_str} = sprintf("\nkspace_style	pppm %-8.4g\n",$parms->{PARMS}{coul_accuracy});
	} elsif ($header =~ /tip4p/ and $parms->{PARMS}{OPTIONS}{PERIODICITY} and $parms->{PARMS}{HAS_CHARGE}) {
		$parms->{PARMS}{kspace_str} = sprintf("\nkspace_style	pppm/tip4p %-8.4g\n",$parms->{PARMS}{coul_accuracy});
	} elsif ($header =~ /coul.*long/ and $parms->{PARMS}{OPTIONS}{PERIODICITY} and $parms->{PARMS}{HAS_CHARGE} and ! $parms->{PARMS}{OPTIONS}{VDW_EWALD}) {
		$parms->{PARMS}{kspace_str} = sprintf("\nkspace_style	pppm %-8.4g\n",$parms->{PARMS}{coul_accuracy});
	} elsif (($header =~ /coul.*long/ and $parms->{PARMS}{OPTIONS}{PERIODICITY} and $parms->{PARMS}{HAS_CHARGE}) or $parms->{PARMS}{OPTIONS}{VDW_EWALD}) {
		$parms->{PARMS}{kspace_str} = sprintf("\nkspace_style	pppm %-8.4g\n",$parms->{PARMS}{coul_accuracy});
	} else {
		$parms->{PARMS}{kspace_str} = "\nkspace_style	none\n";
		$parms->{PARMS}{has_kspace} = 0;
	}
	$header .= $parms->{PARMS}{kspace_str};
	$header .= "kspace_modify		slab 2.0\n" if($parms->{PARMS}{OPTIONS}{PERIODICITY} == 2 and $header =~ /coul\/long/);
	$header .= "fix             Isotopes all property/atom rmass ghost yes\n" if(exists($parms->{PARMS}{OPTIONS}{ISOTOPES}));
	$header .= "fix                  shells all property/atom d2_rs 3\n" if($parms->{QEq}{Opt} == 2 and ! exists($parms->{PARMS}{POLARIZATION}{type}));
	return $header;
}

sub ConvertLammpsBox {
	my ($box) = $_[0];
	my ($generalBox, $i, @dim, %BOX);

	@dim = ("X", "Y", "Z");

	for $i (0 .. 2) {
		$BOX{$dim[$i]}{lo} = $BOX{$dim[$i] . "COORD"}{lo} = $box->[$i]{lo};
		$BOX{$dim[$i]}{hi} = $BOX{$dim[$i] . "COORD"}{hi} = $box->[$i]{hi};
		$BOX{$dim[$i]}{len} = $BOX{$dim[$i] . "COORD"}{len} = $box->[$i]{hi} - $box->[$i]{lo};
		$BOX{$dim[$i]}{CENTER} = $BOX{$dim[$i] . "COORD"}{CENTER} = $BOX{$dim[$i]}{len}/2 - $box->[$i]{lo};
	}
	$BOX{ALPHA} = $BOX{BETA} = $BOX{GAMMA} = $BOX{X}{angle} = $BOX{Y}{angle} = $BOX{Z}{angle} = 90;
	$BOX{ALPHA} = $BOX{X}{angle} = $box->[3]{alpha} if ($#{ $box } > 2);
	$BOX{BETA} = $BOX{Y}{angle} = $box->[3]{beta} if ($#{ $box } > 2);
	$BOX{GAMMA} = $BOX{Z}{angle} = $box->[3]{gamma} if ($#{ $box } > 2);
	return \%BOX;
}

sub setElectrodeOpts {
	my ($parms) = $_[0];
	my ($i, $electrodeQStr, $electrodeStr, $chiStr, $electrodeMolStr, $tmp);

	#first set the moleculeids for top and botton electrodes
	@{ $tmp->{1} } = keys %{ $parms->{PARMS}{ELECTRODE}{top}{MOLID} };
	@{ $tmp->{2} } = keys %{ $parms->{PARMS}{ELECTRODE}{bot}{MOLID} };
	$electrodeQStr = <<"DATAQ";
group  topElectrode molecule @{ $tmp->{1} }
group  botElectrode molecule @{ $tmp->{2} }
group  electrodes union topElectrode botElectrode
group  electrolyte subtract all electrodes 
DATAQ

	if(exists($parms->{PARMS}{ELECTRODE}) and $parms->{PARMS}{ELECTRODE}{type} eq "qeq") {
		$electrodeQStr .= "\nvariable        volt equal 4.0 #NOTE: Change here as necessary";
		$chiStr = "";
		for $i (sort numerically keys %{ $parms->{PARMS}{ELECTRODE}{CHI} }) {
			$chiStr .= <<"DATAChi";
group          g${i} type $i			
set            group g${i} d_loc $parms->{PARMS}{ELECTRODE}{CHI}{$i}
DATAChi

		}
		$electrodeStr = <<"DATAQEq";
group                electrodes union topElectrode botElectrode
print                .
print ================================================
print "Constant Potential Simulation (\${volt}V)"
print ================================================
print                .
# Initialize the property
fix            echemdid_lap all property/atom d_loc d_lap d_Wi
compute        echemdid_lap all property/atom d_loc d_lap d_Wi
# Initialize electronegativity in each electrodes
variable       topChi equal $parms->{PARMS}{ELECTRODE}{top}+(\${volt}/2.0)
variable       botChi equal $parms->{PARMS}{ELECTRODE}{bot}-(\${volt}/2.0)

${chiStr}
set            group all d_lap 0.0
fix            biasedElectrodes electrodes echemdid 1 k 6.0 cut 4.0 norm 0.737106 nelec 10 boundary topElectrode \${topChi} botElectrode \${botChi}
fix            pqeq all pqeq method 2 cycle 1 nevery 1 charge 0 tolerance 1.0e-6 damp 1.0
group    movable union electrodes electrolyte
DATAQEq

	} elsif (exists($parms->{PARMS}{ELECTRODE}) and $parms->{PARMS}{ELECTRODE}{type} eq "conp") { #fix CONP
		#this fix reqqires a unique moleculeid for each electrode, so create here
		$electrodeMolStr =  "variable electrode1MolID equal " . ($parms->{PARMS}{ELECTRODE}{max_res} + 1) . "\n";
		$electrodeMolStr .= "variable electrode2MolID equal " . ($parms->{PARMS}{ELECTRODE}{max_res} + 2) . "\n";
		$electrodeMolStr .= "#set left electrode molid\n";
		for $i (keys %{ $parms->{PARMS}{ELECTRODE}{top}{MOLID} }) {
			$electrodeMolStr .= "set mol $i mol \${electrode1MolID}\n";
		}
		$electrodeMolStr .= "#set right electrode molid\n";
		for $i (keys %{ $parms->{PARMS}{ELECTRODE}{bot}{MOLID} }) {
			$electrodeMolStr .= "set mol $i mol \${electrode2MolID}\n";
		}

		$electrodeQStr .= <<"DATAConp";

${electrodeMolStr}
variable        mu equal 1.979; #**CHANGE HERE**
variable        volt equal 1.0; #**CHARGE HERE**

DATAConp

		$electrodeStr = <<"DATAConp";
print                .
print ================================================
print "Constant Potential (CONP) Simulation (\${volt}V)"
print ================================================
print                .
fix      biasedElectrodes all conp 1 \${mu} \${electrode1MolID} \${electrode2MolID} \${volt} -\${volt} cg 
group    movable subtract all electrodes
DATAConp

	} elsif (exists($parms->{PARMS}{ELECTRODE})) { #fixed Q
		$electrodeQStr .= <<"DATAQ";
variable             electrodeAtomQ index 0.001 #NOTE: Change here as necessary
set                  group topElectrode charge \${electrodeAtomQ}        
set                  group botElectrode charge -\${electrodeAtomQ}
DATAQ

		$electrodeStr = <<"DATAfixedQ";
print                          .
print ==================================================================================
print "NVT production dynamics with \${electrodeAtomQ}e- charges on each electrode atom"
print ==================================================================================
print                       .
#slowly ramp the charge of the electrode atoms from 0 to \${electrodeQ} over 10ps
variable qstepP equal ramp(0.0,\${electrodeAtomQ})
variable qstepN equal ramp(0.0,-\${electrodeAtomQ})
fix topElectrodeQchange topElectrode adapt 10 atom charge v_qstepP
fix botElectrodeQchange botElectrode adapt 10 atom charge v_qstepN
dump    1 all custom 100 \${sname}.\${rtemp}K.electrodeQchage.lammps id type xu yu zu vx vy vz q
fix                  2 all nvt temp \${rtemp} \${rtemp} 100.0
run                  10000
undump 1		
group    movable union electrodes electrolyte
DATAfixedQ

	}
	return ($electrodeQStr, $electrodeStr);
}
1;

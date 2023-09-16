package GROMACS;

require Exporter;
use General qw(PrintProgress FindElementByMass LoadConverter LoadElements);
use ManipAtoms qw(CreateBondsByDistance);
use CERIUS2 qw(getLammpsOpts findDuplicate);
use FileFormats qw(GetXYZFileInfo);
use constant q2e => 18.2223;
use Math::Trig qw(pi acos);

use strict;

our (@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(ParseGroFF UpdateGroFFData parseGroFile parseGroTopFile);
$VERSION = "1.00";

sub numerically { ($a<=>$b); }

sub UpdateGroFFData {
	my ($parms,$elements,$atoms) = @_;
	my ($i, $j, $curr, $ff, $c, $rec, $fftype, @vals, $prev);

	for $i (keys %{ $parms->{ATOMTYPES} }) {
		$j = $parms->{ATOMTYPES}{$i}{FFTYPE};
		next if($j =~ /OW|HW/);
		$ff->{atoms}{$j}{VALS}[0]{name} = $parms->{ATOMTYPES}{$i}{LABEL};
		$ff->{atoms}{$j}{VALS}[0]{mass} = $parms->{ATOMTYPES}{$i}{MASS};
		$ff->{atoms}{$j}{VALS}[0]{r} = $parms->{ATOMTYPES}{$i}{sigma}*10*2**(1/6)/2;
		$ff->{atoms}{$j}{VALS}[0]{e} = $parms->{ATOMTYPES}{$i}{epsilon}/4.184;
		$ff->{atoms}{$j}{VALS}[0]{hybrid} = 0;
		$ff->{atoms}{$j}{VALS}[0]{charge} = $parms->{ATOMTYPES}{$i}{CHARGE};
		$ff->{atoms}{$j}{VALS}[0]{element} = General::FindElementByMass($parms->{ATOMTYPES}{$i}{MASS}, $elements);
	}
	$c = 1; #4-body counter
	for $i (keys %{ $parms }) {
		for $j (sort numerically keys %{ $parms->{$i} }) {
			next if (!exists($parms->{$i}{$j}{ATOMS}));
			$curr = \%{ $ff->{lc $i} };
			if($i ne "INVERSIONS") {
				for (@{ $parms->{$i}{$j}{ATOMS} }) {
					$fftype = $_;
					$prev = $curr;
					$curr = \%{ $curr->{$fftype} };
				}
			} else {
				for (@{ $parms->{$i}{$j}{ATOMS} }) {
					$fftype = $atoms->{$_}{FFTYPE};
					$prev = $curr;
					$curr = \%{ $curr->{$fftype} };
				}
			}
			if($i eq "BONDS") {
				$curr->{VALS}[0]{r0} = $parms->{$i}{$j}{VALS}[1]*10;
				$curr->{VALS}[0]{kb} = $parms->{$i}{$j}{VALS}[0]/418.2/2;
			} elsif($i eq "ANGLES") {
				$curr->{VALS}[0]{t0} = $parms->{$i}{$j}{VALS}[1]*pi/180;
				$curr->{VALS}[0]{kt} = $parms->{$i}{$j}{VALS}[0]/4.182/2;
			} elsif($i eq "TORSIONS") {
				@vals = ();
				$vals[0] = -2*$parms->{$i}{$j}{VALS}[1]-3*$parms->{$i}{$j}{VALS}[3]/2;
				$vals[1] = -$parms->{$i}{$j}{VALS}[2] + $parms->{$i}{$j}{VALS}[4];
				$vals[2] = -$parms->{$i}{$j}{VALS}[3]/2;
				$vals[3] = -$parms->{$i}{$j}{VALS}[4]/4;
				for (0 .. $#vals) {
					next if ($vals[$_] < 0.0001);
					$rec = ();
					$rec->{n} = $_ + 1;
					$rec->{kp} = $vals[$_]/4.184/4;
					$rec->{p0} = pi;
					push @{ $curr->{VALS} }, $rec;
				}
				if(!exists($curr->{VALS})) {
					delete $prev->{$fftype};
					next;
				}
				$curr->{counter} = $c;
				$ff->{torsionOrders}{$c} = "0 1 2 3";
				$c++;
			} elsif($i eq "INVERSIONS") {
				$curr->{VALS}[0]{p0} = $parms->{$i}{$j}{VALS}[0];
				$curr->{VALS}[0]{kp} = $parms->{$i}{$j}{VALS}[1]/4.184/2;
				$curr->{VALS}[0]{n} = $parms->{$i}{$j}{VALS}[2];
				$curr->{counter} = $c;
				$ff->{inversionOrders}{$c} = "2 0 1 3";
				$c++;
			}
		}
	}

	return $ff;
}

sub getTopTypeOpts {
	my ($header) = $_[0];
	my ($fieldStr, $num);

	if ($header eq "PARMS") {
		$fieldStr = '{VDWTYPE} {mix_rule} {gen_pairs} {scale_vdw} {scale_cou}';
		$num = 5;
	}elsif($header eq "ATOMTYPES") {
		$fieldStr = '{TYPE} {LABEL} {MASS} {CHARGE} {UNK} {sigma} {epsilon}';
		$num = 7;
	}elsif($header eq "BONDS") {
		$fieldStr = '{ATOMS}[0] {ATOMS}[1] {TYPE} {VALS}[1] {VALS}[0]';
		$num = 5;
	}elsif($header eq "ANGLES") {
		$fieldStr = '{ATOMS}[0] {ATOMS}[1] {ATOMS}[2] {TYPE} {VALS}[1] {VALS}[0]';
		$num = 6;
	}elsif($header eq "TORSIONS") {
		$fieldStr = '{ATOMS}[0] {ATOMS}[1] {ATOMS}[2] {ATOMS}[3] {TYPE} {VALS}[0] {VALS}[1] {VALS}[2] {VALS}[3] {VALS}[4] {VALS}[5]';
		$num = 11;
	}elsif($header eq "INVERSIONS") {
		$fieldStr = '{ATOMS}[0] {ATOMS}[1] {ATOMS}[2] {ATOMS}[3] {TYPE} {VALS}[0] {VALS}[1] {VALS}[2]';
		$num = 8;
	}elsif($header eq "ATOMS") {
		$fieldStr = '{INDEX} {TYPE} {RESNUM} {RESNAME} {ATMNAME} {UNK} {CHARGE}';
		$num = 7;
	}elsif($header eq "BONDLIST") {
		$fieldStr = '{ATOMS}[0] {ATOMS}[1] {ORDER}';
		$num = 3;
	} elsif($header eq "SETTLES") {
		$fieldStr = '{ATOMS}[0] {ATOMS}[1] {R1} {R2}';
		$num = 4;
	} elsif($header eq "MOLECULE") {
		$fieldStr = '{MOLNAME} {NREXCL}';
		$num = 2;
	}

	return($num, $fieldStr);
}

sub getATPTypeOpts {
	my ($header) = $_[0];
	my ($fieldStr, $num);

	if ($header eq "PARMS") {
		$fieldStr = '{VDWTYPE} {mix_rule} {gen_pairs} {scale_vdw} {scale_cou}';
		$num = 5;
	}elsif($header eq "ATOMTYPES") {
		$fieldStr = '{FFTYPE} {EQUIVALENCE} {ELEMENT} {MASS} {CHARGE} {UNK} {sigma} {epsilon}';
		$num = 7;
	}elsif($header eq "BONDS") {
		$fieldStr = '{ATOMS}[0] {ATOMS}[1] {FUNCT} {VALS}[1] {VALS}[0]';
		$num = 5;
	}elsif($header eq "ANGLES") {
		$fieldStr = '{ATOMS}[0] {ATOMS}[1] {ATOMS}[2] {FUNCT} {VALS}[1] {VALS}[0]';
		$num = 6;
	}elsif($header eq "TORSIONS") {
		$fieldStr = '{ATOMS}[0] {ATOMS}[1] {ATOMS}[2] {ATOMS}[3] {FUNCT} {VALS}[0] {VALS}[1] {VALS}[2] {VALS}[3] {VALS}[4] {VALS}[5]';
		$num = 11;
	}elsif($header eq "INVERSIONS") {
		$fieldStr = '{ATOMS}[0] {ATOMS}[1] {ATOMS}[2] {ATOMS}[3] {VALS}[0] {VALS}[1] {VALS}[2]';
		$num = 7;
	}elsif($header eq "ATOMS") {
		$fieldStr = '{INDEX} {TYPE} {RESNUM} {RESNAME} {ATMNAME} {UNK} {CHARGE}';
		$num = 7;
	}elsif($header eq "BONDLIST") {
	$fieldStr = '{ATOMS}[0] {ATOMS}[1] {ORDER}';
		$num = 3;
	}

	return($num, $fieldStr);
}

sub getItpTypeOpts {
	my ($header) = $_[0];
	my ($fieldStr, $num);

	if($header eq "ATOMS") {
		$fieldStr = '{INDEX} {TYPE} {RESNUM} {RESNAME} {ATMNAME} {UNK} {CHARGE} {MASS}';
		$num = 8;
	} elsif ($header eq "BONDS") {
		$fieldStr = '{I} {J} {FUNCT} {C0} {C1}';
		$num = 5;
	} elsif ($header eq "ANGLES") {
		$fieldStr = '{I} {J} {K} {FUNCT} {THETA} {K0}';
		$num = 6;
	} elsif ($header eq "DIHEDRALS") {
		$fieldStr = '{I} {J} {K} {L} {FUNCT} {ANGLE} {K0} {MULT}';
		$num = 8;
	}

	return($num, $fieldStr);
}

sub ParseGroFF {
	my ($inFile, $mode, $isITP, $ref) = @_;
	my ($topData, $header, $fields, $num, @vals, $counter, $molName);
	my ($getOpts, $tmp, $parms, $atoms, $bonds, $iOffset, $i, $j, $n);

	$header = "";
	$getOpts = \&getTopTypeOpts if (! defined($isITP) or $isITP==0);
	$getOpts = \&getItpTypeOpts if ($isITP==1);
	$getOpts = \&getATPTypeOpts if ($isITP==2);
	$inFile = $inFile->{FF} if ($mode == 1);
	open GROMACSTOPFILE, $inFile or die "ERROR: Cannot open GROMACS topology file $inFile: $!\n";
	while(<GROMACSTOPFILE>) {
		chomp;
		if ($header eq "DIHEDRALS" && $_ =~ /(improper dihedrals|; improper_)/) {
			$header = "INVERSIONS";
			$fields = '{I} {J} {K} {L} {FUNCT} {ANGLE} {K0}';
			$num = 7;
			$counter = 1;
		}
		#$_ =~ s/\;.*$//;
		$_ =~ s/\t/ /g;
		if ($_ =~ /^\[ (\w+) \]/) {
			next if !($header = parseHeading($1, $isITP));
			($num, $fields) = $getOpts->($header);
			$counter = 1;
			$counter = scalar(keys %{ $topData->{$header} }) + 1 if (exists($topData->{$header}));
		} elsif ($header =~ /TORSION/ and $_ =~ /define improper_([^;]*);*.*$/)  {
			$tmp = $1;
			$tmp =~ s/[Y|Z]/X/g;
			$tmp =~ s/_+/ /g;
			@vals= split /\s+/,$tmp;
			($num, $fields) = $getOpts->("INVERSIONS");
			$counter = 1;
			$counter = scalar(keys %{ $topData->{INVERSIONS} }) + 1 if (exists($topData->{INVERSIONS}));
			while ($fields =~ /(\S+)/g) {
				die "ERROR: header: INVERSIONS counter: $counter field: $1 vals: @vals error: $!\n" 
					if (! eval('$topData->{INVERSIONS}{' . $counter . '}' . $1 . ' = shift @vals'));
			}
			$counter++;
		} elsif ($header and $_ =~ /^\s*([^;]*);*.*$/) {
			@vals = split /\s+/, $1;
			if($header =~ /INVERSION/) { #offset the bond/torsion/improper indices to match numbering atom numbering
				$n = 3;
				$n = 1 if ($header =~ /BONDLIST/);
				for $j (0 .. $n) {
					$vals[$j] += $iOffset;
				}
			}
			next if (scalar(@vals) < $num); # or ($header eq "ATOMS" and $vals[3] eq "SOL"));
			next if ($header eq "INVERSIONS" and ! $isITP and $_ !~ /improper/);
			next if ($header eq "ATOMTYPES" and $_ =~ / 0 /);
			$molName = $vals[0] if ($header eq "MOLECULE");
			while ($fields =~ /(\S+)/g) {
				die "ERROR: header: $header counter: $counter field: $1 vals: @vals error: $!\n" 
					if (! eval('$topData->{' . $header . '}{' . $counter . '}' . $1 . ' = shift @vals'));
			}
			$iOffset = $counter - $topData->{ATOMS}{$counter}{INDEX} if($header eq "ATOMS");
			$topData->{$header}{$counter}{MOLNAME} = $molName if($header eq "SETTLES");
			$counter++;
		}
	}
	close GROMACSTOPFILE;
	die "ERROR: $inFile is not valid!\n" if (! $topData);

	$parms = genGroFF($topData, $ref->{FF});
	return $parms if ($mode == 1);

	($atoms, $bonds) = genGroStruct($topData, $ref->{atoms}, $ref->{box}, $parms);
	return ($atoms, $bonds);
}

sub genGroStruct {
	my ($topData, $atoms, $box, $FF) = @_;
	my ($i, $j, $parm, $typeMap, $bonds, $atomSel);

	for $i (values %{ $topData->{ATOMS} }) {
		$typeMap->{$i->{ATMNAME}}{aData} = $i;
		for $j (values %{ $topData->{ATOMTYPES}}) {
			if($i->{TYPE} eq $j->{TYPE}) {
				$typeMap->{$i->{ATMNAME}}{pData} = $j;
				last;
			}
		}
	}
	for $i (sort numerically keys %{ $atoms }) {
		undef $parm;
		die "ERROR: Cannot locate type data for $atoms->{$i}{ATMNAME} in topology file... Cannot continue...\n"
			if (! exists($typeMap->{ $atoms->{$i}{ATMNAME} }));
		$parm = $typeMap->{ $atoms->{$i}{ATMNAME} };
		$atoms->{$i}{CHARGE} = $parm->{aData}{CHARGE};
		$atoms->{$i}{MASS} = $parm->{pData}{MASS};
		$atoms->{$i}{FFTYPE} = $parm->{pData}{LABEL};
		$atomSel->{$i} = 1;
	}	

	#generate bondlist
	for $i (values %{ $topData->{BONDLIST}}) {
		push @{ $bonds->{ $i->{ATOMS}[0] } }, $i->{ATOMS}[1];
		push @{ $bonds->{ $i->{ATOMS}[1] } }, $i->{ATOMS}[0];
		delete $atomSel->{ $i->{ATOMS}[0] };
		delete $atomSel->{ $i->{ATOMS}[1] };
	}
	if(scalar(keys %{ $atomSel }) > 0) {
		print "Generating bonding info for " . scalar(keys %{ $atomSel }) . "/" . scalar(keys %{ $atoms }) . " atoms...";
		&ManipAtoms::CreateBondsByDistance(\%{ $atoms }, \%{ $bonds }, $box, $atomSel, $FF);
	}
	return ($atoms, $bonds);
}
sub genGroFF {
	my ($topData, $parms) = @_;
	my ($i, $j, $k, $l, $v, $vType, $rec, $tmp, $akey, $aTypes);
	my ($counter, $curr,  $CNV, $elements, $eleNum, $vMap, $tMap); 

	$CNV = &General::LoadConverter();
	$elements = &General::LoadElements();
	#set some defaults
	$parms->{PARMS}{dielectric} = 1;
	$parms->{PARMS}{cut_coul} = 10;
	$parms->{PARMS}{cut_vdw} = 10;
	$parms->{PARMS}{coul_accuracy} = 1e-06;
	#map atom numbers to atom fftypes
	for $i (keys %{ $topData->{ATOMS} }) {
		for $j (values %{ $topData->{ATOMTYPES}}) {
			if($topData->{ATOMS}{$i}{TYPE} eq $j->{TYPE}) {
				$tMap->{$i} = $j->{LABEL};
				last;
			}
		}
	}

	$rec = (
			{
				"TYPE"		 => "LJ_6_12",
				"Lammps"	 => getLammpsOpts("LJ_6_12","vdw", $CNV, 0 , $parms->{PARMS}),
				"USED"		 => 0,
				"IT"		 => "vdw",
				"IGNORE"     => 0,
				"COUL"       => "0010.5",
				"VDW"        => "0010.5",
				"USE_CHARGE" => 0,
			}
	);	
	$rec->{Lammps}{opts} = "9 10";
	for $l (sort numerically keys %{ $topData->{ATOMTYPES} }) {
		#atom types
		$i = $topData->{ATOMTYPES}{$l};
		$eleNum = General::FindElementByMass($i->{MASS}, $elements);
		$k = $i->{LABEL};
		$curr = \%{ $parms->{ATOMTYPES}{$k} };
		for $j (keys %{ $i }) {
			$curr->{$j} = $i->{$j};
		}
		$curr->{ELEMENT} = $elements->{$eleNum};
		$akey = "$i->{LABEL} $i->{LABEL} ";
		for $j (values %{ $topData->{ATOMS} }) {
			if ($j->{TYPE} eq $i->{TYPE}) {
				$curr->{ATOM} = $j->{ATMNAME};
				last;
			}
		}
		$curr->{USE_CHARGE} = 0; $curr->{TYPEID} = $l;
		#vdw
		$counter = findDuplicate($parms->{VDW}{$k}{$k}, "LJ_6_12");
		$counter = 1 if (! $counter);
		$curr = \%{ $parms->{VDW}{$k}{$k}{$counter} };
		for $j (keys %{ $rec }) {
			$curr->{$j} = $rec->{$j};
		}
		$curr->{ATOM} = $curr->{KEY} = $akey;
		@{ $curr->{VALS} } = ($i->{epsilon}/4.184, $i->{sigma}*10);
	}
	#settles to make O-H bond and H-O-H angles
	&createH2Oopts($topData) if(exists($topData->{SETTLES}));

	#bonds/angles/torsions/imporpers
	delete $rec->{IT}; delete $rec->{TYPE}; delete $rec->{Lammps};
	$rec->{NUM} = 1;
	$vMap->{BONDS}{T}      = "HARMONIC";   $vMap->{BONDS}{S}      = [836.8,0.1];
	$vMap->{ANGLES}{T}     = "THETA_HARM"; $vMap->{ANGLES}{S}     = [8.368,1];
	$vMap->{TORSIONS}{T}   = "MULTIHAR";   $vMap->{TORSIONS}{S}   = [4.184,-4.184,4.184,-4.184,4.184];
	$vMap->{INVERSIONS}{T} = "IT_JIKL";    $vMap->{INVERSIONS}{S} = [1,8.368,1];
	for $k (keys %{ $topData }) {
		next if ($k !~ /BONDS|ANGLES|TORSIONS|INVERSIONS/);
		$v = $vMap->{$k};
		$rec->{TYPE} = $v->{T};
		$vType = lc(substr($k,0,-1));
		$vType = "dihedral" if ($k eq "TORSIONS");
		$rec->{Lammps} = getLammpsOpts($v->{T}, $vType, $CNV);
		for $l (values %{ $topData->{$k} }) {
			$curr = \%{ $parms->{$k} };
			@{ $aTypes } = @{ $l->{ATOMS} };
			if ($k eq "INVERSIONS") {
				@{ $aTypes } = map { $tMap->{$_} } @{ $l->{ATOMS} };
				($aTypes->[2], $aTypes->[0]) = ($aTypes->[0], $aTypes->[2]);
			}
			for $j (@{ $aTypes }) {
				$curr = \%{ $curr->{$j} };
			}
			$counter = findDuplicate($curr, $v->{T});
			$counter = 1 if (! $counter);
			$curr = \%{ $curr->{$counter} };
			for $j (keys %{ $rec }) {
				$curr->{$j} = $rec->{$j};
			}
			$curr->{KEY} = $curr->{ATOM} = "@{ $aTypes } ";
			for $j (0 .. $#{ $v->{S} }) {
				$curr->{VALS}[$j] = $l->{VALS}[$j]/$v->{S}[$j];
			}
			if ($k eq "INVERSIONS") {
				($curr->{VALS}[0], $curr->{VALS}[1]) = ($curr->{VALS}[1], $curr->{VALS}[0]);
				if($curr->{VALS}[1] == 0) {
					$curr->{VALS}[1] = 1;
				} else {
					$curr->{VALS}[1] = -1;
				}
			}
		}
	}
	return $parms;
}

sub createH2Oopts {
	my ($topData) = $_[0];
	my ($i, $ff1, $ff2, $oh_bond, $hh_bond, $hoh_angle, $counter);

	#get fftypes
	for $i (values %{ $topData->{ATOMS} }) {
		$ff1 = $i->{TYPE} if (($i->{RESNAME} eq $topData->{SETTLES}{1}{MOLNAME}) and $i->{ATMNAME} =~ /^O/);
		$ff2 = $i->{TYPE} if (($i->{RESNAME} eq $topData->{SETTLES}{1}{MOLNAME}) and $i->{ATMNAME} =~ /^H/);
	}
	die "ERROR: Cannot figure out oxygen ($ff1) or hydrogen ($ff2) fftypes!"
		if (! defined($ff1) or ! defined($ff2));
	for $i (values %{ $topData->{ATOMTYPES}}) {
		$ff1 = $i->{LABEL} if($i->{TYPE} eq $ff1);
		$ff2 = $i->{LABEL} if($i->{TYPE} eq $ff2);
	}

	#bond
	$oh_bond = $topData->{SETTLES}{1}{R1};
	$counter = scalar(keys %{ $topData->{BONDS} }) + 1;
	$topData->{BONDS}{$counter} = (
									{
										ATOMS => [$ff1,$ff2],
										TYPE  => 1,
										VALS  => [4184,$oh_bond],
									},
	);

	#angle
	$hh_bond = $topData->{SETTLES}{1}{R2};
	$hoh_angle = acos(($oh_bond**2-$hh_bond**2/2)/$oh_bond**2)*180/pi;
	$counter = scalar(keys %{ $topData->{ANGLES}}) + 1;
	$topData->{ANGLES}{$counter} = (
									{
										ATOMS => [$ff2,$ff1,$ff2],
										TYPE  => 1,
										VALS  => [418.4,$hoh_angle],
									},
	);
	
}
sub parseHeading {
	my ($inStr, $isITP) = @_;
	my ($header);

	$header = "";
	if($inStr =~/(defaults|atomtypes|bondtypes|angletypes|dihedraltypes|atoms|bonds|angles|dihedrals|settles|moleculetype)/) {
		if($1 eq "atomtypes") { $header = "ATOMTYPES"; }
		elsif($1 eq "dihedraltypes") { $header = "TORSIONS"; }
		elsif($1 eq "defaults") { $header = "PARMS"; }
		elsif($1 eq "bonds" and ! $isITP) { $header = "BONDLIST"; }
		elsif($1 eq "dihedrals" and ! $isITP) { $header = "INVERSIONS"; }
		else { $header = uc $1; $header =~ s/type//i; }
	}
	return $header;
}

sub parseGroFile {
	my ($parms, $inFile, $ffTypeOpt) = @_;
	my ($lc, $atoms, $i, $id, $atom1, $atom2, $bonds, $bStr, $box);

	$ffTypeOpt = 1 if (! defined($ffTypeOpt));
	open GROFILE, $inFile or die "ERROR: Cannot access GRO file $inFile: $!\n";
	$lc = 0;
	while (<GROFILE>) {
		chomp;
		$lc++;
		next if ($lc < 3);
		if ($_ =~ /^\s*(\d+)(\S+)\s+(\S+)\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
			$atoms->{$4} = (
							{
								RESNUM  => $1,
								RESNAME => $2,
								ATMNAME => $3,
								INDEX   => $4,
								XCOORD  => $5*10,
								YCOORD  => $6*10,
								ZCOORD  => $7*10,
								LABEL   => "ATOM",
							}
						);
		} else {
			$bStr = $_;
		}
	}
	close GROFILE or die "ERROR: Cannot close $inFile: $!\n";
	die "ERROR: No valid ATOM data read from $inFile"
		if (! defined($atoms));

	$box = createBox($atoms, $bStr) if (defined($bStr));
	return ($atoms, undef, $box) if (! $ffTypeOpt);


	#add fftype
	for $i (keys %{ $atoms } ) {
		$id = findTypeID($atoms->{$i}{ATMNAME}, $parms);
		die "ERROR: Cannot find fftype for atom # $i ($atoms->{$i}{ATMNAME})!\n"
			if (!defined($id));
		$atoms->{$i}{FFTYPE} = $parms->{ATOMTYPES}{$id}{FFTYPE};
		$atoms->{$i}{CHARGE} = $parms->{ATOMTYPES}{$id}{CHARGE};
	}

	for $i (keys %{ $parms->{BONDLIST} }) {
		($atom1, $atom2) = @{ $parms->{BONDLIST}{$i}{ATOMS} };
		push @{ $bonds->{$atom1} }, $atom2;
		push @{ $bonds->{$atom2} }, $atom1;
	}

	return ($atoms, $bonds, $box);
}

sub createBox {
	my ($atoms, $bStr) = @_;
	my ($box);

	if ($bStr =~ /^\s*(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s*$/) {
		$box->{X}{len} = $box->{XCOORD}{len} = $box->{X}{hi} = $box->{XCOORD}{hi} = $1*10;
		$box->{Y}{len} = $box->{YCOORD}{len} = $box->{Y}{hi} = $box->{YCOORD}{hi} = $2*10;
		$box->{Z}{len} = $box->{ZCOORD}{len} = $box->{Z}{hi} = $box->{ZCOORD}{hi} = $3*10;
		$box->{X}{lo} = $box->{Y}{lo} = $box->{Z}{lo} = $box->{XCOORD}{lo} = $box->{YCOORD}{lo} = $box->{ZCOORD}{lo} = 0;
		$box->{ALPHA} = $box->{BETA} = $box->{GAMMA} = $box->{X}{angle} = $box->{Y}{angle} = $box->{Z}{angle} = 90;
	}
	return $box;
}

sub findTypeID {
	my ($atmName, $parms) = @_;
	my ($ffid, $i);

	for $i (values %{ $parms->{ATOMS} }) {
		if ($i->{ATMNAME} eq $atmName) {
			$ffid = findFFType($i->{TYPE},$parms->{ATOMTYPES});
			if (defined($ffid)) {
				$parms->{ATOMTYPES}{$ffid}{CHARGE} = $i->{CHARGE};
				last;
			}
		}
	}

	return $ffid;
}

sub findFFType {
	my ($sStr, $types) = @_;
	my ($i, $ffid);

	for $i (keys %{ $types }) {
		if ($types->{$i}{LABEL} eq $sStr) {
			$ffid = $i;
			last;
		}
	}
	return $ffid;
}

sub parseGroTopFile {
	my ($topFile,$atoms,$bonds,$writeBonds) = @_;
	my ($i, $mode, $noAtoms, $indxFromName, $j, $valid);

	$noAtoms = scalar(keys %{ $atoms });
	$writeBonds = 1 if (! defined($writeBonds));
	%{ $bonds } = () if (defined($bonds) and $writeBonds);
	$valid=1;
	$mode = 0;
	open TOPFILE, $topFile or die "ERROR: Cannot open $topFile: $!\n";
	while(<TOPFILE>) {
		chomp;
		if ($_ =~ /^\{\s+atomType\s+mass/i) {
			$mode = 1; #masses
			$i = 1;
		}elsif ($_ =~ /^\{\s+atomName\s+atomType\s+charge/i) {
			$mode = 2; #atom types and charge
			$i = 1;
		}elsif ($_ =~ /^\{\s+Bonds:/i && $writeBonds) {
			$mode = 3; #bonds
		}elsif ($mode==1 and $_ =~ /^MASS\s+(\S+)\s+(\d+\.\d+)/) {
			die "ERROR: Atom # $i does not exists in atoms array!\n" if(! $noAtoms && ! exists($atoms->{$i}));
			$atoms->{$i}{ATMNAME}  = $1; $atoms->{$i}{MASS} = $2;
			$indxFromName->{$1} = $i; $i++;
		}elsif ($mode==2 and $_ =~ /^ATOM\s+(\S+)\s+TYPE=\s+(\S+)\s+CHARGE=\s+(\-?\d+\.\d+)/) {
			$atoms->{$i}{ATMNAME} = $1; $atoms->{$i}{FFTYPE} = $2; $atoms->{$i}{CHARGE} = $3;
			$indxFromName->{$1} = $i; $i++;
		}elsif($mode==3 and $_ =~ /^BOND\s+(\S+)\s+(\S+)/) {
			$i = $indxFromName->{$1}; $j = $indxFromName->{$2};
			push @{ $bonds->{$i} }, $j; push @{ $bonds->{$j} }, $i;
		}
	}
	close TOPFILE;
}
1;

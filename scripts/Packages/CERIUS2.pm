package CERIUS2;

require Exporter;

use FindBin qw($Bin);
use lib "$FindBin::Bin";

use strict;
use Math::Trig qw(pi);
use Storable qw(dclone);
use General qw(LoadElements LoadConverter);

use constant PI => atan2(1,1) * 4;

our(@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(parseCerius2FF saveCeriusFF parseMPSimFF parseReaxFF GenerateUFFParms getLammpsOpts findDuplicate);
$VERSION = "1.00";

sub getLammpsOpts {
	my ($parmName, $parmType, $SHASH, $istip4p, $parms) = @_;
	my (%RET, $same_cut);
	my ($searchStr) = lc($parmType) . "_" . lc($parmName);

	$istip4p = 0 if (! defined($istip4p));
	$searchStr .= "_tip4p" if ($parmType eq "vdw" and $istip4p and $parmName eq "LJ_6_12");
	
	$same_cut = 1;
	#$same_cut = 1 if ($parms->{cut_coul} == $parms->{cut_vdw});
	if (exists($SHASH->{$searchStr})) {
		%RET = %{ $SHASH->{$searchStr} };
	} else {
		$RET{name} = lc($parmName);
		$RET{opts} = "";
		$RET{missing} = 1;
	}

	if($parmType eq "vdw") {
		if($RET{name} =~ /charmm/) {
			if (! $same_cut) {
				$RET{opts} = ($parms->{cut_vdw}-1) . " $parms->{cut_vdw} " . ($parms->{cut_coul}-1) . " $parms->{cut_coul}";
			} else {
				$RET{opts} = ($parms->{cut_coul}-1) . " $parms->{cut_coul}";
			}
		} elsif ($RET{name} =~ /lj|gromacs/ && $RET{name} !~ /dreiding/) {
			$RET{opts} = " $parms->{cut_vdw}  $parms->{cut_coul}";
			$RET{opts} = $parms->{cut_coul} if ($same_cut);
		}
	}
	return \%RET;
}

sub findDuplicate {
	my ($parm, $currType) = @_;
	my ($i, $typeCounter);
	
	$typeCounter = 0;
	$currType = lc($currType);
	for $i (keys %{ $parm }) {
		if (lc($parm->{$i}{TYPE}) eq $currType) {
			$typeCounter = $i;
			last;
		}
	}

	return $typeCounter;
}

sub parseReaxFF {
	my ($ff_file, $elements, $alter, $oldFF, $rexponOpt) = @_;
	my (%PARMS, $type_counter, $currParm);

	$type_counter = 0;
	if (defined($oldFF)) {
		%PARMS = %{ $oldFF };
	}

	open FORCEFIELD, $ff_file->{FF} or die "Cannot open $ff_file->{FF}: $!\n";
	while (<FORCEFIELD>) {
		chomp;
		if ($_ =~ /^\s*([A-Za-z][A-Za-z]?)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) { # atom line
			$type_counter++;
			$PARMS{"ATOMTYPES"}{$1} = (
									   {
										"TYPEID"		=>	$type_counter,
										"ATOM"			=>	$1,
										"MASS"			=>	$4,
										"CHARGE"		=>	undef,
										"NUMBONDS"		=>	$3,
										"LONEPAIRS"		=>	undef,
										"OTHER"			=>	undef,
										"LABEL"			=>	$1,
										"USED"			=>	0,
										"USE_CHARGE"	=>	0,
									   }
									);
#		   print "ATOMTYPES: $1: $type_counter\n";
				$currParm = \%{ $PARMS{"VDW"}{$1}{$1} };
				$currParm->{1} =(
									{
										"TYPE"		=>	"REAX",
										"Lammps"	=>	{
														"name"	=>	"reax",
														"opts"	=>	"",
														},
										"KEY"		=>	"$1 $1 ",
										"VALS"		=>	[0],
										"ATOM"		=>	"$1 $1 ",
										"USED"		=>	0,
										"IT"		=>	"vdw",
									}
								);
#				print "VDW: $1: $type_id\n"
				if($rexponOpt) {
						$currParm->{1}{Lammps} = (
													{
														"name" => "rexpon",
														"opts" => "NULL lgvdw yes coulomb_off yes checkqeq no safezone 2.0 mincap 200",
													}
												);
				}

		}
	}
	close FORCEFIELD;
	die "ERROR: Force Field file $ff_file->{FF} is invalid!\n"
		if (! %PARMS);
	return (\%PARMS);
}

sub findElementNum {
	my ($symbol, $elements) = @_;
	my ($i, $elementNum);

	$symbol = lc $symbol;
	$elementNum = 0;
	for $i (keys %{ $elements }) {
		if(lc($elements->{$i}{SYMBOL}) eq $symbol) {
			$elementNum = $i;
			last;
		}
	}

	return $elementNum;
}

sub getHbondList {
	my ($parms, $donor, $acceptor, $hydrogen, $elements) = @_;
	my ($hbondList, $i, $j, $k, $iNum, $jNum, $dList); 
	my ($aList, $hList, $rec, $stored);

	$hbondList = ();

	if($donor ne "X" && $acceptor ne "X" && $hydrogen ne "X") {
		$hbondList = [{ "donor" => $donor, "acceptor" => $acceptor, "hydrogen" => $hydrogen, }];
		return $hbondList;
	}

	push @{ $dList }, $donor if($donor ne "X");
	push @{ $hList }, $hydrogen if($hydrogen ne "X");
	push @{ $aList }, $acceptor if($acceptor ne "X");

	#first search for donors and hydrogens
	for $i (keys %{ $parms->{BONDS} }) {
		next if($i =~ /^\s*H_\s*$/);
		$iNum = findElementNum($parms->{ATOMTYPES}{$i}{ATOM}, $elements);
		for $j (keys %{ $parms->{BONDS}{$i} }) {
			next if($j =~ /^\s*H_\s*$/);
			$jNum = findElementNum($parms->{ATOMTYPES}{$j}{ATOM}, $elements);
			next if((! $iNum && ! $jNum) || ($iNum != 1 && $jNum != 1));
			next if($donor ne "X" && ($i ne $donor && $j ne $donor));
			if($iNum == 1) {
				push @{ $hList }, $i if($hydrogen eq "X");
				push @{ $dList }, $j if($donor eq "X" && $jNum =~ /^(7|8|9)$/);
			} else {
				push @{ $hList }, $j if($hydrogen eq "X");
				push @{ $dList }, $i if($donor eq "X" && $iNum =~ /^(7|8|9)$/);
			}
		}
	}

	#now search for acceptors
	if($acceptor eq "X") {
		for $i (keys %{ $parms->{ATOMTYPES} }) {
			$iNum = findElementNum($parms->{ATOMTYPES}{$i}{ATOM}, $elements);
			next if ($iNum !~ /^(6|7|8|9|16|17)$/); # not C,N,O,F,S,Cl
			next if ($parms->{ATOMTYPES}{$i}{LONEPAIRS} == 0); # must have lonepairs
			push @{ $aList }, $i;
		}
	}

	for $i (@{ $dList }) {
		for $j (@{ $aList }) {
		   #next if ($i eq $j);
			for $k (@{ $hList }) {
				next if ($k !~ /H___|H_F/);
				next if (exists($stored->{$i}{$j}{$k}));
				$stored->{$i}{$j}{$k} = 1;
				$rec = ({ "donor" => $i, "acceptor" => $j, "hydrogen" => $k, });
				push @{ $hbondList }, $rec;
		   }
		}
	}

	return $hbondList;
}

sub searchForHBparm {
	my ($currParm, $hatom) = @_;
	my ($index, $i, $count);

	$index = $count = 0;
	for $i (keys %{ $currParm }) {
		$count = $i;
		next if (! exists($currParm->{$i}{HATOM}));
		if($currParm->{$i}{HATOM} eq $hatom) {
			$index = $i;
			last;
		}
	}

	$index = $count + 1 if(! $index );

	return $index;
}

sub getVals {
	my ($parmStr) = $_[0];
	my (@ret);

	$parmStr =~ s/#.*$//;
	@ret = split /\s+/,$parmStr;

	return @ret;
}

sub parseCerius2FF {
	my ($ff_file, $alter, $oldFF) = @_;
	my (%PARMS, $which_var, $in_data, $type_counter, @vdws, $hb_counter);
	my (@bonds, @angles, $tmp1, @inversions, $atom4, $currParm, $calcExp6, $istip4p);
	my (@torsions, $parmType, $counter, $use_hb, $use_charge, $ty_tmp);
	my ($bond_counter, $angle_counter, $torsion_counter, @tmp, $inversion_counter, $bool, $pref, $special_vdw, $ignore_ex);
	my ($atom1, $atom2, $atom3, $parmDat, $i, $dihdr_scale, $CNV, $vdw_counter, $crossType, $rec, $ignore, $special_coul);
	my ($donors, $acceptors, $hydrogens, $elements, $hbtmp, $scale_flag, $curr_atype, $parmParent, $parmCounter, $tmpParm);

	$CNV = &General::LoadConverter();
	$elements = &General::LoadElements();
	$bond_counter = $angle_counter = $torsion_counter = $counter = 0;
	$use_hb = $hb_counter = $crossType = 0;
	$calcExp6 = 1;
	$istip4p = $which_var = $ignore_ex = 0;
	#$istip4p = 1 if ($ff_file =~ /tip4/);
	if (defined($oldFF)) {
		%PARMS = %{ $oldFF };
	} else {
		$PARMS{PARMS} = (
							{
								SW                                   => 0,
								cut_vdw                              => 14.0,
								cut_coul                             => 15.0,
								coul_accuracy                        => 0.00001,
								hbond_distance                       => 2.5,
								hbond_angle                          => 90,
								scale_torsions                       => 0,
								scale_torsions_by_n_defined_torsions => 0,
								single_inversion                     => 0,
								dielectric                           => 1,
								EXO_CYCLIC                           => 1.0,
								QEq                                  => 0,
								USE_HBOND                            => 0,
								HAS_HBONDS                           => 0,
							}
		);
	}
	$istip4p =1 if (exists($PARMS{"PARMS"}{"tip4_om_dist"}));
	$type_counter = 0;
	$type_counter = ($PARMS{"PARMS"}{"type_counter"} + 1) if (exists($PARMS{"PARMS"}{"type_counter"}));
	$use_charge = 0;

	open FORCEFIELD, $ff_file->{FF} or die "Cannot open force field file $ff_file->{FF}: $!\n";
	while (<FORCEFIELD>) {
		chomp;
		$in_data = $_;
		$in_data =~ s/#.*$// if ($which_var != 3 and $which_var != 8);
		$PARMS{PARMS}{SW} = 1 if ($in_data =~ / SW /);
		if ($in_data =~ /^END/) {
			undef $pref if (defined($pref));
			$which_var = 0;
		} elsif ($in_data =~ /^PREFERENCES/) {
			$pref = 1;
		} elsif ($in_data =~ /^\s+MD_TIMESTEP\s+(\d+\.?\d*)/) {
			$PARMS{PARMS}{OPTIONS}{MD_TIMESTEP} = $1 
			if (! exists($PARMS{PARMS}{OPTIONS}{MD_TIMESTEP}) or $1 < $PARMS{PARMS}{OPTIONS}{MD_TIMESTEP});
		} elsif ($in_data =~ /^\s+EXOCYCLIC_TORSIONS_SCALE_FACTOR\s+(\d+\.?\d*)/) {
			$PARMS{"PARMS"}{EXO_CYCLIC} = $1;
		} elsif ($in_data =~ /^\s+USE_CURRENT_EXP6_VALS\s+T/) {
			$calcExp6 = 0;
		} elsif ($in_data =~ /^\s*COU_DIELETRIC_CONSTANT\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"dielectric"} = $1;
		} elsif ($in_data =~ /^\s+IGNORE_COUL_VDW_EXCLUSIONS\s+(T)/i) {
			$ignore_ex = 1;
		} elsif ($in_data =~ /^\s*SINGLE_INVERSION\s+T/) {
			$PARMS{"PARMS"}{"single_inversion"} = 1;
		 } elsif ($in_data =~ /^\s*SCALE_TORSIONS_ABOUT_COMMON_BOND\s+T/) {
			$PARMS{"PARMS"}{"scale_torsions"} = 1;
		 } elsif ($in_data =~ /^\s*SCALE_BY_N_DEFINED_TORSIONS\s+T/) {
			$PARMS{"PARMS"}{"scale_torsions_by_n_defined_torsions"} = 1;
		} elsif ($in_data =~ /^\s*COU_SPLINE_OFF\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"cut_coul"} = $1;
		} elsif ($in_data =~ /^\s*VDW_SPLINE_OFF\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"cut_vdw"} = $1;
		} elsif ($in_data =~ /^\s*VDW_COMBINATION_RULE\s+(\w+)/) {
			$PARMS{"PARMS"}{"mix_rule"} = lc($1);
		} elsif ($in_data =~ /^\s*EWALD_SUM_COU_ACCURACY\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"coul_accuracy"} = $1;
		} elsif ($in_data =~ /^\s*(COU|VDW)_EXCLUDE_(\d)\-(\d)\s+(\w)/) {
			if (lc($4) eq "t") {
				$bool = 0;
			} else {
				$bool = 1;
			}
			$PARMS{"PARMS"}{"scale_" . lc($1) . "_" . $2 . $3} = $bool;
			$special_coul .= $bool if(lc($1) eq "cou");
			$special_vdw .= $bool if (lc($1) eq "vdw");
		 } elsif ($in_data =~ /^\s*COU_1-4_SCALE_FACTOR\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"scale_cou"} = $1;
			$PARMS{"PARMS"}{"scale_cou_14"} = $1 if ($PARMS{"PARMS"}{"scale_cou_14"});
			$special_coul .= $1;
		 } elsif ($in_data =~ /^\s*VDW_1-4_SCALE_FACTOR\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"scale_vdw"} = $1;
			$PARMS{"PARMS"}{"scale_vdw_14"} = $1 if ($PARMS{"PARMS"}{"scale_vdw_14"});
			if($PARMS{"PARMS"}{"scale_cou"} ne $PARMS{"PARMS"}{"scale_vdw"}) { 
				$PARMS{"PARMS"}{"same_scale"} = 0;
			} else {
				$PARMS{"PARMS"}{"same_scale"} = 1;
			}
			$special_vdw .= $1;
		} elsif ($in_data =~ /^ HYDROGEN_BONDS\s+T/) {
			$use_hb = 1;
			$PARMS{"PARMS"}{"hbond_distance"} = 4.5;
			$PARMS{"PARMS"}{"hbond_angle"} = 60;
			$PARMS{PARMS}{USE_HBOND} = 1;
		} elsif ($use_hb and $in_data =~ /^ H-BOND_LIST_DISTANCE_OFF\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"hbond_distance"} = $1;
		} elsif ($use_hb and $in_data =~ /^ H-BOND_LIST_ANGLE_OFF\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"hbond_angle"} = $1;
		} elsif ($in_data =~ /^ ASSIGN_CHARGE\s+T/) {
			$use_charge  = 1;
		} elsif ($in_data =~ /^ TIP4_OM_DIST\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"tip4_om_dist"} = $1;
		} elsif ($in_data =~ /^ATOMTYPES/) {
			$which_var = 1;
		} elsif ($in_data =~ /HYDROGEN_BONDS/) {
				$which_var = 2;			
		} elsif ($in_data =~ /^DIAGONAL_VDW/) {
				$which_var = 3;
		} elsif ($in_data =~ /^BOND_STRETCH/) {
			$which_var = 4;
		} elsif ($in_data =~ /^ANGLE_BEND/) {
			$which_var = 5;
			$crossType = "";
		} elsif ($in_data =~ /^TORSIONS/) {
			$which_var = 6;
			$crossType = "";
		} elsif ($in_data =~ /^INVERSIONS/) {
			$which_var = 7;
			$crossType = "";
		} elsif ($in_data =~ /^STRETCH_STRETCH/) {
			$which_var = 5;
			$crossType = "BondBond";
		} elsif ($in_data =~ /^STRETCH_BEND_STRETCH/) {
			$which_var = 5;
			$crossType = "BondAngle";
		} elsif ($in_data =~ /^OFF_DIAGONAL_VDW/) {
			$which_var = 8;
		} elsif ($in_data =~ /^UREY_BRADLEY/) {
			$which_var = 9;
			$crossType = "";
		} elsif ($in_data =~ /^SEPARATED_STRETCH_STRETCH/) {
			$which_var = 6;
			$crossType = "13BondBond";
		} elsif ($in_data =~ /^TORSION_BEND_BEND/) {
			$which_var = 6;
			$crossType = "AngleAngle";
		} elsif ($in_data =~ /^BEND_BEND/) {
			$which_var = 7;
			$crossType = "AngleAngle";
		} elsif ($in_data =~ /^GENERATOR/) {
			$which_var = 11;
		} elsif ($in_data =~ /^QEq/) {
			$which_var = 10;
		} elsif ($in_data =~ /^ATOM_TYPING_RULES/) {
			$which_var = 202;
		} elsif ($in_data =~ /^EQUIVALENCE/) {
			$which_var = 303;
		} elsif ($which_var == 303 && $in_data =~ /^\s*(\S+)\s+(.*)$/) {
			$atom1 = $1;
			$PARMS{ATOMTYPES}{$atom1}{USED} = 0;
			$PARMS{ATOMTYPES}{$atom1}{TYPEID} = 0;
			@tmp = split /\s+/, $2;
			for $i (@tmp) {
				$PARMS{EQUIVALENCE}{$i}{$atom1} = 1;
			}
		} elsif ($which_var == 202 && $in_data =~ /^\s+(.*)$/) { #atomtyping
			@tmp = split /\s+/, $1;
			next if ($#tmp < 4);
			$parmCounter = $tmp[4];
			if ($#tmp == 5) {
				$curr_atype = $tmp[0];
				$rec = (
						{
							"ELE"     => $tmp[1],
							"BOND"    => undef,
							"HYBRID"  => $tmp[2],
							"RING"    => $tmp[3],
							"FFTYPE"  => $curr_atype,
							"NBOND"   => $parmCounter,
							"NBCOUNT" => 0,
							"PARENT"  => undef,
						}
					);
				push @{ $PARMS{TYPING}{$tmp[1]} }, $rec;
				if($parmCounter > 0) {
					$currParm = \%{ $PARMS{TYPING}{$tmp[1]}[$#{ $PARMS{TYPING}{$tmp[1]} }]};
				}
			} else {
				$parmCounter = $tmp[3];
				$tmp[0] = "any" if ($tmp[0] =~ /\*\*/);
				$rec = (
						{
							"ELE" => $tmp[0],
							"BOND" => undef,
							"HYBRID" => $tmp[1],
							"RING" => $tmp[2],
							"PARENT" => \%{ $currParm },
							"FFTYPE" => $curr_atype,
							"NBOND" => $parmCounter,
							"NBCOUNT" => 0,
						}
					);
				push @{ $currParm->{BOND} }, $rec;
				$currParm->{NBCOUNT}++;
				if ($parmCounter > 0) {
					$tmpParm = $currParm;
					$currParm = \%{ $tmpParm->{BOND}[$#{ $tmpParm->{BOND} } ] };
				}
				while(defined($currParm->{PARENT}) and $currParm->{NBCOUNT} == $currParm->{NBOND}) {
					$tmpParm = $currParm->{PARENT};
					$currParm = $tmpParm;
				}
			}
		} elsif ($in_data =~ /^\s*(\S+)\s+(.+)$/ and ($which_var == 11)) { #UFF Generator
			@tmp = split /\s+/, $2;
			next if ($#tmp < 9);
			$PARMS{PARMS}{UFF}{$1} = (
										{
										"radius"		=> $tmp[0],
										"angle"			=> $tmp[1],
										"zstar"			=> $tmp[2],
										"zeta"			=> $tmp[3],
										"uenerg"		=> $tmp[4],
										"uang"			=> $tmp[5],
										"prd"			=> $tmp[6],
										"cis"			=> $tmp[7],
										"torbar"		=> $tmp[8],
										"elecneg"		=> $tmp[9],
										}
									);
		} elsif ($in_data =~ /^\s*(\S+)\s+(\d+\.\d+\e?\-?\d*.+)/ and ($which_var == 10)) { #QEq/PQEq?
			@tmp = split /\s+/, $2;
			next if ($#tmp < 2);
			$PARMS{PARMS}{QEq} = 1;
			$PARMS{QEq}{$1}{ZETA} = $PARMS{QEq}{$1}{CHRG} = 0;
			$PARMS{QEq}{$1} = (
								{
									"CHI"   => $tmp[0],
									"ETA"   => $tmp[1],
									"GAMMA" => $tmp[2],
									"USED"  => 0,
								}
							);
			if ($#tmp > 2) {
				$PARMS{QEq}{$1}{ZETA} = $tmp[3];
				$PARMS{QEq}{$1}{CHRG} = $tmp[4] if ($#tmp > 4);
			}
		} elsif ($in_data =~ /^\s*(\S+)\s+(\w+)\s+(\d+\.\d+e?[+-]?\d*)\s+(\-?\d+\.\d+)\s+(\d+)\s+(\d+)\s+(\d+)/i and ($which_var == 1)) { #atom types
			$type_counter += 1;
			$PARMS{"ATOMTYPES"}{$1} = (
										{
											"TYPEID"		=> $type_counter,
											"ATOM"			=> $2,
											"MASS"			=> $3,
											"CHARGE"		=> $4,
											"NUMBONDS"		=> $5,
											"LONEPAIRS"		=> $7,
											"OTHER"			=> $6,
											"LABEL"			=> $1,
											"USED"			=> 0,
											"USE_CHARGE"	=> $use_charge,
										}
									);
#			print "ATOMTYPES: $1: $type_counter\n";
		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+.*)/i and $which_var == 2) { #hbonds
			$donors = $hydrogens = $acceptors = "";
			$donors = $1 if(exists($PARMS{"ATOMTYPES"}{$1}) || $1 eq "X");
			$hbtmp = getLammpsOpts($3,"hbond", $CNV);
			if(! exists($PARMS{ATOMTYPES}{$3}) && $3 ne "X" && ! exists($hbtmp->{missing})) { #old style hbonds, only donor and acceptor
				$acceptors = $2 if(exists($PARMS{ATOMTYPES}{$2}) || $2 eq "X");
				$hydrogens = "X";
				$parmType = $3;
				@vdws = getVals("$4 $5");
			} else { #new hbond: donor-hydrogen ... aceptor
				$acceptors = $3 if(exists($PARMS{ATOMTYPES}{$3}) || $3 eq "X");
				$hydrogens = $2 if(exists($PARMS{ATOMTYPES}{$2}) || $2 eq "X");
				$parmType = $4;
				@vdws = getVals($5);
			}
			next if (! $donors || ! $acceptors || ! $hydrogens);
			$hbtmp = getHbondList(\%PARMS, $donors, $hydrogens, $acceptors, $elements);
			next if (! $hbtmp);
			push(@vdws, "4") if($parmType eq "LJ_12_10"); # append the power for the cosine - 4 - to the end of the data
			for $i (@{ $hbtmp }) {
				$atom1 = $i->{donor};
				$atom2 = $i->{acceptor};
				if (! exists($PARMS{VDW}{$atom1}{$atom2})) {
					$vdw_counter = 1;
				} else {
					$vdw_counter = searchForHBparm($PARMS{VDW}{$atom1}{$atom2}, $i->{hydrogen});
				}
				$currParm = \%{ $PARMS{"VDW"}{$atom1}{$atom2} };
				$currParm->{$vdw_counter} = (
												{
													"TYPE"		=> $parmType,
													"Lammps"	=> getLammpsOpts($parmType,"hbond", $CNV),
													"KEY"		=> "$atom1 $atom2 ",
													"VALS"		=> [$i->{hydrogen},@vdws],
													"ATOM"		=> "$atom1 $atom2 ",
													"USED"		=> 0,
													"IT"		=> "hbond",
													"HATOM"		=> $i->{hydrogen},
													"ACCEPTOR"	=> $i->{acceptor},
													"DONOR"		=> $i->{donor},
													"IGNORE"    => 0,
												}
											);
			}

		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s*(\S*)\s*(.*)/i and ($which_var == 3 or $which_var == 8)) { #vdw
			$atom1 = $1;
			$special_coul = $special_vdw = undef if ($ignore_ex);
			if ($which_var == 3) {
				$atom2 = $1;
				$parmType = $2;
				$parmDat = "$3 $4 $5";
			} else {
				$atom2 = $2;
				$parmType = $3;
				$parmDat = "$4 $5";
			}
			if ($atom1 lt $atom2) {
				($atom1, $atom2) = ($atom2, $atom1);
			}
			$ignore = 0;
			$ignore = 1 if (uc($parmDat) =~ /IGNORE/);
			$ignore = 2 if (uc($parmDat) =~ /DELETE/);
			$ignore = 3 if (uc($parmDat) =~ /REPLACE/);
			if (uc($parmType) !~ /UFF_GEN/ and (defined($4) or $ignore>0)) {
				next if ((! validType(\%PARMS, $atom1) or ! validType(\%PARMS, $atom2)) and
						($atom1 !~ /_shell/ and $atom2 !~ /shell/));
				if($ignore > 0) {
					@vdws = ()
				}elsif ($parmDat =~ /(.*)\# 1\-4 scaling: (.*)/) {
					@vdws = getVals("$1 $2");
				} else {
					@vdws = getVals($parmDat);
				}
				@vdws = GetAbs(\@vdws) if ($parmType !~ /EAM/);

				if ($ignore == 0 and (! defined($alter) or $alter != 0)) {
					if (uc($parmType) eq "VDW_MORSE") {
						($vdws[0], $vdws[1]) = ($vdws[1], $vdws[0]);
						($vdws[1], $vdws[2]) = ($vdws[2], $vdws[1]);
						$vdws[1] /= ($vdws[2] * 2); # changed from div to multi 07/28/2007
						#$vdws[0] *= 2;
						#$vdws[1] /= ($vdws[2]/2);
					} elsif (uc($parmType) eq "LJ_6_12") {
						($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
						$vdws[1] = $vdws[1] / (2**(1/6));
						if ($#vdws > 1) {
							($vdws[2],$vdws[3]) = ($vdws[3],$vdws[2]);
							$vdws[3] = $vdws[3] /(2**(1/6));
						}
#						$vdws[1] = $vdws[1]/2;
					} elsif (uc($parmType) eq "EXPO_6" and $calcExp6) {
						($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
						@tmp = @vdws;
						@vdws = ();
						#$vdws[0] = 6 * $tmp[1] * exp($tmp[2]) / ($tmp[2] - 6);
						#$vdws[1] = $tmp[0]/$tmp[2];
						#$vdws[2] = $tmp[1] * $tmp[2] * $tmp[0]**6/($tmp[2] - 6);
						$vdws[0] = $tmp[0] * (6/($tmp[2]-6)) * exp($tmp[2]);
						$vdws[1] = $tmp[1]/$tmp[2];
						$vdws[2] = $tmp[1]**6 * $tmp[0] * ($tmp[2]/($tmp[2] - 6));
					} elsif (uc($parmType) eq "BUCKINGHAM") {
						#($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
						$vdws[1] = 1/$vdws[1];
					} elsif (uc($parmType) eq "LJ_6_9") {
						($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
						$vdws[1] = $vdws[1] / 1.14471438354;
					} elsif (uc($parmType) eq "BORN") {
						splice(@vdws,2,0,0);
						$vdws[$#vdws] *= -1;
					}
				}
				if ($atom1 eq "*" or $atom2 eq "*") {
					&updateAllVDW($PARMS{VDW},$atom1,$atom2,$parmType,\@vdws,$ignore, $special_coul, $CNV, $istip4p, \%PARMS);
					next;
				}
				if (! exists($PARMS{VDW}{$atom1}{$atom2})) {
					$vdw_counter = 1;
				} else {
					$vdw_counter = findDuplicate($PARMS{VDW}{$atom1}{$atom2}, $parmType);
					$vdw_counter = scalar(keys %{ $PARMS{VDW}{$atom1}{$atom2} }) + 1 if (! $vdw_counter);
				}
				$tmp1 = "vdw";
				$currParm = \%{ $PARMS{"VDW"}{$atom1}{$atom2} };
				@vdws = @{ $currParm->{$vdw_counter}{VALS} } 
					if ($ignore == 1 and exists($currParm->{$vdw_counter}{VALS}));
				$ignore = 1 if ($ignore == 2);
				&createVDWParm($currParm, $vdw_counter, $atom1, $atom2, $parmType, \@vdws, $ignore, $special_coul, $CNV, $istip4p, \%PARMS);	
				if ($parmType eq "YUKAWA") {
					$currParm->{$vdw_counter}{Lammps}{opts} = "$vdws[1] $vdws[2]";
				}
			}
		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\w+)\s+(.+)/ and ($which_var == 4)) { #bonds
			($atom1, $atom2, $parmType, $parmDat) = ($1,$2,$3,$4);
			($atom1, $atom2) = ($atom2, $atom1) 
				if (exists($PARMS{BONDS}{$atom2}) and 
					exists($PARMS{BONDS}{$atom2}{$atom1}));
			$ignore = 0;
			$ignore = 1 if (uc($parmDat) =~ /IGNORE/);
			if ($parmType !~ /UFF_GEN/i) {
				next if (! validType(\%PARMS, $atom1) or ! validType(\%PARMS, $atom2));
				if(!$ignore) {
					@bonds = getVals($parmDat);
					next if (! @bonds);
				} else {
					@bonds = ($parmDat);
				}
				if (! $ignore and (! defined($alter) or $alter != 0)) {
					$bonds[0] = $bonds[0] / 2; #fix for the 1/2 factor in harmonic eqns between cerius and lammps
					if ($parmType =~ /^MORSE/i) {
						$bonds[0] *= 2;
						($bonds[0], $bonds[2]) = ($bonds[2], $bonds[0]); # kb is now bonds[0]
						($bonds[1], $bonds[2]) = ($bonds[2], $bonds[1]); # r0 is now bonds[2]
						$bonds[1] = sqrt($bonds[1]/(2 * $bonds[0]));
						pop @bonds if ($#bonds > 2);
					} elsif ($parmType =~ /^CLASS2/i) {
						$bonds[0]*=2;
					}
				}
				$bond_counter = 1;
				if (exists($PARMS{BONDS}{$atom1}{$atom2})) {
					$bond_counter = findDuplicate($PARMS{BONDS}{$atom1}{$atom2},$parmType);
					$bond_counter = scalar(keys %{ $PARMS{BONDS}{$atom1}{$atom2} }) + 1 if (! $bond_counter);
				}
				$currParm = \%{ $PARMS{"BONDS"}{$atom1}{$atom2} };
				$currParm->{$bond_counter} = (
												{
													"INDEX"		=> $bond_counter,
													"TYPE"		=> $parmType,
													"Lammps"	=> getLammpsOpts($parmType,"bond"),
													"VALS"		=> [@bonds],
													"USED"		=> 0,
													"KEY"		=> "$atom1 $atom2 ",
													"IGNORE"    => $ignore,
												}
											);

			}
		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)/ and ($which_var == 5 or $which_var == 9)) { #angles/urey-bradley
			($atom1, $atom2, $atom3, $parmType, $parmDat) = ($1,$2,$3,$4,$5);
			($atom1, $atom2, $atom3) = ($atom3, $atom2, $atom1) 
				if (exists($PARMS{ANGLES}{$atom3}) and 
					exists($PARMS{ATOMS}{$atom3}{$atom2}) and 
					exists($PARMS{ATOMS}{$atom3}{$atom2}{$atom1}));
			$ignore = 0;
			$ignore = 1 if (uc($parmDat) =~ /IGNORE/);
			if ($parmDat !~ /UFF_GEN/i) {
				next if (! validType(\%PARMS, $atom1) or ! validType(\%PARMS, $atom2) or ! validType(\%PARMS, $atom3));
				if (! $ignore) {
					@angles = getVals($parmDat);
					next if (! @angles);
				} else {
					@angles = ($parmDat);
				}
				if (! $ignore and $crossType eq "" and (! defined($alter) or $alter != 0)) {
					$angles[0] = $angles[0]/2; #same fix as bonds
				}
				if ($parmDat =~ /^COS_HARMON/i) {
					if (sin($angles[1] * pi/180) > 0.001) {
						$angles[0] /= sin($angles[1] * pi/180)**2;
					}
				} elsif ($parmDat =~ /^R-COSINE/i) { #bond angle cosine cross term
					$angles[3] = -$angles[3]/sin($angles[2] * pi/180);
					$angles[4] = -$angles[4]/sin($angles[2] * pi/180);
				} elsif ($parmDat =~ /^SW/i) {
					$angles[0] *= 2;
				} elsif ($parmDat =~ /^QUARTIC/i) {
					$angles[0] *= 2;
				}
				$angle_counter = 1;
				if($which_var == 9) { #urey-bradley
					next if (! exists($PARMS{ANGLES}{$atom1}) or
							 ! exists($PARMS{ANGLES}{$atom1}{$atom2}) or
							 ! exists($PARMS{ANGLES}{$atom1}{$atom2}{$atom3}));
					for $i (values %{ $PARMS{ANGLES}{$atom1}{$atom2}{$atom3} }) {
						next if ($i->{Lammps}{name} !~ /harmonic/);
						push @{ $i->{VALS } }, @angles;
					}
				} else {
					if (exists($PARMS{ANGLES}{$atom1}{$atom2}{$atom3})) {
						$angle_counter = findDuplicate($PARMS{ANGLES}{$atom1}{$atom2}{$atom3},$parmType);
						$angle_counter = scalar(keys %{ $PARMS{ANGLES}{$atom1}{$atom2}{$atom3} }) + 1 if (! $angle_counter);
					}
					$currParm = \%{ $PARMS{"ANGLES"}{$atom1}{$atom2}{$atom3} };
					$currParm->{$angle_counter} =  (
													{
														"INDEX"		=> $angle_counter,
														"TYPE"		=> $parmType,
														"Lammps"	=> getLammpsOpts($parmType,"angle", $CNV),
														"VALS"		=> [@angles],
														"USED"		=> 0,
														"KEY"		=> "$atom1 $atom2 $atom3 ",
														"IGNORE"	=> $ignore,
														"CTYPE"		=> $crossType,
													}
												   );
				}
			}
		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)/ and ($which_var == 6)) { #torsions
			($atom1, $atom2, $atom3, $atom4, $parmType, $parmDat) = ($1,$2,$3,$4,$5,$6);
			($atom1, $atom2, $atom3, $atom4) = ($atom4, $atom3, $atom2, $atom1) 
				if (exists($PARMS{TORSIONS}{$atom4}) and
					exists($PARMS{TORSIONS}{$atom4}{$atom3}) and 
					exists($PARMS{TORSIONS}{$atom4}{$atom3}{$atom2}) and 
					exists($PARMS{TORSIONS}{$atom4}{$atom3}{$atom2}{$atom1}));
			$ignore = 0;
			$ignore = 1 if (uc($parmDat) =~ /IGNORE/);
			$currParm = ();
			if ($parmType !~ /^UFF_GEN/i) {
				if (! validType(\%PARMS, $atom1) or 
						 ! validType(\%PARMS, $atom2) or 
						 ! validType(\%PARMS, $atom3) or 
						 ! validType(\%PARMS, $atom4)) {
							 undef $atom1;
							 next;
						 }
				@torsions = getVals($parmDat);
				if (! $ignore and $crossType eq "" and (! defined($alter) or $alter != 0)) {
					if ($parmType =~ /^DIHEDRAL/i) {
						if ($torsions[2] == 1) {
							$torsions[2] = 180;
						} else {
							$torsions[2] = 0;
						}
					}
					if($parmType !~ /^MULTIHAR|OPLS/i) {
						$torsions[0] = $torsions[0] / 2; #1/2 fix again
					}
				}
				$torsions[0] *= 2 if ($parmType =~ /^FOURIER/i);
				$torsion_counter = 1;
				$torsion_counter = scalar(keys %{ $PARMS{TORSIONS}{$atom1}{$atom2}{$atom3}{$atom4} }) + 1
					if (exists($PARMS{TORSIONS}{$atom1}{$atom2}{$atom3}{$atom4}));
				$scale_flag = $PARMS{PARMS}{scale_torsions};
				$currParm = \%{ $PARMS{"TORSIONS"}{$atom1}{$atom2}{$atom3}{$atom4} };
				$currParm->{$torsion_counter} =  (
													{
														"TYPE"		=> $parmType,
														"Lammps"	=> getLammpsOpts($parmType,"dihedral", $CNV),
														"INDEX"		=> $torsion_counter,
														"VALS"		=> [@torsions],
														"USED"		=> 0,
														"KEY"		=> "$atom1 $atom2 $atom3 $atom4 ",
														"NUM"		=> 1,
														"PER"		=> ($#torsions + 1),
														"1_4scale"	=> $PARMS{"PARMS"}{"scale_vdw_14"},
														"CTYPE"		=> $crossType,
														"do_scale"	=> $scale_flag,
														"COUL"      => $special_coul,
														"VDW"       => $special_vdw,
														"IGNORE"	=> $ignore,
														"FFTYPEID"  => $ff_file->{FFTYPEID},
													}
												);
			}
		} elsif ($in_data =~ /^\s+(\-?\d+\.\d+)\s+(.+)/ and $which_var == 6) { #multiple torsions
			next if (! defined($atom1) or ! validType(\%PARMS, $atom1)); #prevents recording of spurious lines
			@torsions = ();
			@torsions = getVals($2);
#			GetAbs(\@torsions);
			if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
				if (lc($parmType) eq "dihedral") {
					if ($torsions[1] == 1) {
						$torsions[1] = 180;
					} else {
						$torsions[1] = 0;
					}
				}

				$tmp1 = $1 / 2; #1/2 fix again
			} else {
				$tmp1 = $1;
			}
			$currParm = \%{ $PARMS{"TORSIONS"}{$atom1}{$atom2}{$atom3}{$atom4} };
			push @{ $currParm->{$torsion_counter}{"VALS"} }, $tmp1;
			push @{ $currParm->{$torsion_counter}{"VALS"} }, @torsions;
			$currParm->{$torsion_counter}{NUM}++;
			$currParm->{$torsion_counter}{BLANK} = "$atom1 $atom2 $atom3 $atom4";
			#print "FOUND MULTIPLE TORSIONS FOR $atom1 $atom2 $atom3 $atom4 # $currParm->{$torsion_counter}{NUM}\n";
		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)/ and ($which_var == 7)) { #inversions
			($atom1, $atom2, $atom3, $atom4, $parmType, $parmDat) = ($1,$2,$3,$4,$5,$6);
			$ignore = 0;
			$ignore = 1 if (uc($parmDat) =~ /IGNORE/);
			if ($parmType !~ /^UFF_GEN/i) {
				next if (! validType(\%PARMS, $atom1) or 
					     ! validType(\%PARMS, $atom2) or 
						 ! validType(\%PARMS, $atom3) or 
						 ! validType(\%PARMS, $atom4));
				$inversion_counter++;
				@inversions = getVals($parmDat);
				next if (! @inversions);
#				GetAbs(\@torsions);
 
				if (! $ignore and $crossType eq "" and (! defined($alter) or $alter != 0)) {
					if (lc($parmType) eq "it_jikl") {
						if ($inversions[1] == 180) {
							$inversions[1] = -1;
						} else {
							$inversions[1] = 1;
						}
					}
					$inversions[0] = $inversions[0] / 2 if ($parmType !~ /UMBRELLA|DISTANCE/); #1/2 fix again
				}
				$currParm = \%{ $PARMS{"INVERSIONS"}{$atom1}{$atom2}{$atom3}{$atom4} };
				$currParm->{$inversion_counter} = (
														{
															"TYPE"		=> $parmType,
															"Lammps"	=> getLammpsOpts($parmType,"inversion", $CNV),
															"INDEX"		=> $inversion_counter,
															"VALS"		=> [@inversions],
															"USED"		=> 0,
															"KEY"		=> "$atom1 $atom2 $atom3 $atom4 ",
															"CTYPE"		=> $crossType,
															"IGNORE"	=> $ignore,
														}
												);
			}
		}			
		
	}
	close FORCEFIELD;

	die "ERROR: Force Field file $ff_file->{FF} is invalid!\n"
		if (! %PARMS);
	
	if ($PARMS{"PARMS"}{"cut_vdw"} > $PARMS{"PARMS"}{"cut_coul"}) {
		$PARMS{"PARMS"}{"cut_vdw"} = $PARMS{"PARMS"}{"cut_coul"};
	}
	$PARMS{"PARMS"}{"type_counter"} = $type_counter;

	return (\%PARMS);
}

sub validType {
	my ($parms, $sType) = @_;
	my ($rVal);

	$rVal = 0;
	$rVal = 1 if(exists($parms->{ATOMTYPES}{$sType}) or (lc($sType) eq "x") or $sType eq "*" or
		        (exists($parms->{EQUIVALENCE}) and exists($parms->{EQUIVALENCE}{$sType})));
	return $rVal;
}

sub updateAllVDW {
	my ($vdws,$atom1,$atom2,$parmType,$vals,$ignore, $special_coul, $CNV, $istip4p, $PARMS) = @_;
	my ($i, $j, $curr, $count, $tot);

	($atom1,$atom2) = ($atom2, $atom1) if ($atom1 eq "*");
	for $i (keys %{ $vdws }) {
		if (exists($vdws->{$i}{$atom1})) {
			$curr = \%{ $vdws->{$i}{$atom1} };
		} else {
			$curr = \%{ $vdws->{$atom1}{$i} };
		}
		if (! keys %{ $curr }) {
			&createVDWParm($curr, 1, $atom1, $atom2, $parmType, $vals, $ignore, $special_coul, $CNV, $istip4p, $PARMS);
		} else {
			$count = 0;
			for $j (keys %{ $curr }) {
				if ($curr->{$j}{TYPE} eq $parmType) {
					$count = $j;
					last;
				}
			}
			if (! $count) {
				$tot = scalar(keys %{ $curr });
				&createVDWParm($curr, ($tot+1), $atom1, $atom2, $parmType, $vals, $ignore, $special_coul, $CNV, $istip4p, $PARMS);
			} else {
				if($ignore == 1 || $ignore == 3) {
					$curr->{$count}{IGNORE} = $ignore;
				} elsif ($ignore == 2) {
					delete $curr->{$count};
				} elsif (@{ $vals }) {
					@{ $curr->{$count}{VALS} } = @{ $vals };
				}
			}
		}
	}
}

sub createVDWParm {
	my ($curr, $count, $atom1, $atom2, $parmType, $vals, $ignore, $special_coul, $CNV, $istip4p, $PARMS) = @_;
	$curr->{$count} = (
					{
						"TYPE"		=> $parmType,
						"KEY"		=> "$atom1 $atom2 ",
						"VALS"		=> [@{$vals}],
						"ATOM"		=> "$atom1 $atom2 ",
						"USED"		=> 0,
						"IGNORE"	=> $ignore,
						"IT"		=> "vdw",
						"COUL"      => $special_coul,
						"Lammps"    => getLammpsOpts($parmType,"vdw", $CNV, $istip4p, $PARMS->{PARMS}),
					}
	);
	$curr->{$count}{VALS} = [] if ($ignore == 3 and ($atom1 ne $atom2)); #this takes care of the * creation placeholders
}

sub updateAllVDWold {
	my ($vdws,$atom1,$atom2,$parmType,$vals,$ignore) = @_;
	my ($i, $j, $curr);

	($atom1,$atom2) = ($atom2, $atom1) if ($atom1 eq "*");
	for $i (keys %{ $vdws }) {
		if (exists($vdws->{$i}{$atom1})) {
			$curr = \%{ $vdws->{$i}{$atom1} };
		} elsif (exists($vdws->{$atom1}) and exists($vdws->{$atom1}{$i})) {
			$curr = \%{ $vdws->{$atom1}{$i} };
		} else {
			next;
		}
		for $j (keys %{ $curr } ) {
			next if ($curr->{$j}{TYPE} ne $parmType);
			if($ignore == 1) {
				$curr->{$j}{IGNORE} = 1;
			} elsif ($ignore == 2) {
				delete $curr->{$j};
			} elsif (@{ $vals }) {
				@{ $curr->{$j}{VALS} } = @{ $vals };
			}
		}
	}
}

sub saveCeriusFF {
	my ($ffData, $save_name, $ELEMENTS) = @_;
	my ($c, $shft_angle, $torsions, $counter, $element, $vdws, $i, $chrg, $torsionType);
	my ($atom1, $atom2, $atom3, $atom4, @ATMS, $hk, $inversions, @TMP, $equiv, $equivStr);

	$ELEMENTS = &General::LoadElements();
	$vdws = "";
	open OUTFILE, "> $save_name" or die "Cannot create file $save_name: $!\n";
	print OUTFILE "ATOMTYPES\n";
	for $atom1 (keys %{ $ffData->{"atoms"} }) {
		next
			if ($atom1 eq "?");
		$c = \%{ $ffData->{"atoms"}{$atom1}{"VALS"}[0] };
		next if (! exists($c->{r}));
		$element = $ELEMENTS->{ $c->{"element"} }{"SYMBOL"};
		$chrg = 0;
		$chrg = $c->{charge} if (exists($c->{charge}));
		printf OUTFILE " %-11s%-3s%12.5f%8.4f%4s%4s%4s",
		$atom1,$element,$c->{"mass"},$chrg,$c->{"hybrid"},0,0;
		$vdws .= sprintf(" %-11s LJ_6_12%14.4f%8.4f",
		$atom1, ($c->{"r"}  * 2), $c->{"e"});
		$c->{r} = 0.00001 if ($c->{r} == 0);
		#if (! exists($c->{"1_4"})) {
			#$c->{"1_4"}{"r"} = $c->{"r"};
			#$c->{"1_4"}{"e"} = $c->{"e"};
		#}
		if (exists($c->{"1_4"})) {
			$vdws .= sprintf(" # 1-4 scaling: %8.5f%8.5f\n",($c->{"1_4"}{"r"}*2),abs($c->{"1_4"}{"e"}));
		} else {
			$vdws .= "\n";
		}
		if (exists($c->{"name"})) {
			printf OUTFILE " # $c->{name}";
		}
		printf OUTFILE "\n";
		$equiv->{$c->{equivalence}}{$atom1} = 1
			if(exists($c->{equivalence}));
	}
	if(defined($equiv)) {
		for $i (keys %{ $equiv }) {
			@TMP = keys %{ $equiv->{$i} };
			$equivStr .= sprintf(" %-11s @TMP\n",$i);
		}
		print OUTFILE "END\nEQUIVALENCE\n$equivStr";
	}
	print OUTFILE "END\n#\nDIAGONAL_VDW\n$vdws";
	print OUTFILE "END\n#\nATOM_TYPING_RULES\nEND\n\#\nOFF_DIAGONAL_VDW\nEND\n#\nBOND_STRETCH\n";
	for $atom1 (keys %{ $ffData->{"bonds"} }) {
		for $atom2 (keys %{ $ffData->{"bonds"}{$atom1} }) {
			$c = \%{ $ffData->{"bonds"}{$atom1}{$atom2}{"VALS"}[0] };
			printf OUTFILE " %-9s%-11s HARMONIC%13.4f%10.4f\n",
			$atom1, $atom2, ($c->{"kb"}  * 2), $c->{"r0"};
		}
	}
	print OUTFILE "END\n#\nANGLE_BEND\n";
	for $atom1 (keys %{ $ffData->{"angles"} }) {
		for $atom2 (keys %{ $ffData->{"angles"}{$atom1} }) {
			for $atom3 (keys %{ $ffData->{"angles"}{$atom1}{$atom2} }) {
				$c = \%{ $ffData->{"angles"}{$atom1}{$atom2}{$atom3}{"VALS"}[0] };
				$c->{"t0"} = ($c->{"t0"} * 180/pi);
				printf OUTFILE " %-9s%-9s%-11s THETA_HARM%10.4f%10.4f\n",
				$atom1, $atom2, $atom3, ($c->{"kt"} * 2), $c->{"t0"};
			}
		}
	}
	print OUTFILE "END\n#\nUREY_BRADLEY\n";
	for $atom1 (keys %{ $ffData->{"urey_bradley"} }) {
		for $atom2 (keys %{ $ffData->{"urey_bradley"}{$atom1} }) {
			for $atom3 (keys %{ $ffData->{"urey_bradley"}{$atom1}{$atom2} }) {
				$c = \%{ $ffData->{"urey_bradley"}{$atom1}{$atom2}{$atom3}{"VALS"}[0] };
				printf OUTFILE " %-9s%-9s%-11s HARMONIC%13.4f%10.4f\n",
				$atom1, $atom2, $atom3, ($c->{"ku"} * 2), $c->{"su"};
			}
		}
	}
	print OUTFILE "END\n#\nTORSIONS\n";
	for $atom1 (keys %{ $ffData->{"torsions"} }) {
		for $atom2 (keys %{ $ffData->{"torsions"}{$atom1} }) {
			for $atom3 (keys %{ $ffData->{"torsions"}{$atom1}{$atom2} }) {
				for $atom4 (keys %{ $ffData->{"torsions"}{$atom1}{$atom2}{$atom3} }) {
					@TMP = @ATMS = ();
					$torsions = \@{ $ffData->{"torsions"}{$atom1}{$atom2}{$atom3}{$atom4}{"VALS"} };
					@TMP = ( $atom1, $atom2, $atom3, $atom4 );
					$i = $ffData->{"torsions"}{$atom1}{$atom2}{$atom3}{$atom4}{"counter"};
					while ($ffData->{"torsionOrders"}{$i} =~ /(\d)/g) {
						push @ATMS, $TMP[$1];
					}
					undef $torsionType;
					$torsionType = $ffData->{"torsions"}{$atom1}{$atom2}{$atom3}{$atom4}{TYPE}
						if(exists($ffData->{"torsions"}{$atom1}{$atom2}{$atom3}{$atom4}{TYPE}));
					if(!defined($torsionType) or $torsionType ne "MHARMONIC") {	
						printf OUTFILE " %-9s%-9s%-9s%-11s SHFT_DIHDR", 
						$ATMS[0], $ATMS[1], $ATMS[2], $ATMS[3];
						for $counter (0 .. $#{ $torsions }) {
							$c = \%{ $torsions->[$counter] };
							$c->{"p0"} = (($c->{"p0"} * 180/pi));# % 360);
							if ($counter == 0) {
								printf OUTFILE "%11.6f%10.4f%10.4f\n", ($c->{"kp"} * 2), $c->{"n"}, $c->{"p0"};
							} else {
								printf OUTFILE "%50s%11.6f%10.4f%10.4f\n", "", ($c->{"kp"} * 2), $c->{"n"}, $c->{"p0"};
							}
						}
					} else {
						printf OUTFILE " %-9s%-9s%-9s%-11s MULTIHAR  ", 
						$ATMS[0], $ATMS[1], $ATMS[2], $ATMS[3];
						for $counter(0 .. 4) {
							printf OUTFILE "%11.6f", $torsions->[0][$counter];
						}
						printf OUTFILE "\n";
					}
				}
			}
		}
	}
	
	print OUTFILE "END\n#\nINVERSIONS\n";
	for $atom1 (keys %{ $ffData->{"inversions"} }) {
		for $atom2 (keys %{ $ffData->{"inversions"}{$atom1} }) {
			for $atom3 (keys %{ $ffData->{"inversions"}{$atom1}{$atom2} }) {
				for $atom4 (keys %{ $ffData->{"inversions"}{$atom1}{$atom2}{$atom3} }) {
					@TMP = @ATMS = ();
					$torsions = \@{ $ffData->{"inversions"}{$atom1}{$atom2}{$atom3}{$atom4}{"VALS"} };
					@TMP = ( $atom1, $atom2, $atom3, $atom4 );
					$i = $ffData->{"inversions"}{$atom1}{$atom2}{$atom3}{$atom4}{"counter"};
					while ($ffData->{"inversionOrders"}{$i} =~ /(\d)/g) {
						push @ATMS, $TMP[$1];
					}
					$c = \%{ $torsions->[0] };
					if (! defined($c->{type}) || $c->{type} eq "IT_JIKL") {
						printf OUTFILE " %-9s%-9s%-9s%-11s IT_JIKL%14.4f%10.4f%10.4f\n",
							 $ATMS[0], $ATMS[1], $ATMS[2], $ATMS[3],
							 ($c->{"kp"} * 2), $c->{"p0"},$c->{"n"};
					} else {
						printf OUTFILE " %-9s%-9s%-9s%-11s IT_IJKL%14.4f%10.4f\n",
							 $ATMS[0], $ATMS[1], $ATMS[2], $ATMS[3],
							 ($c->{"kp"} * 2), $c->{"p0"};
					}
				}  
			}
		}
	}
	print OUTFILE "END\n\#\nCOULOMBIC\n X		X		   CONST-EPS\nEND\n";
	close OUTFILE;
}

sub GetAbs {
	my ($inVals) = @_;
	my ($counter, @ret);

	for $counter (@{ $inVals }) { 
		if ($counter =~ /(.+)e(.+)/i) {
			push @ret, $1 * 10 ** $2;
		} elsif ($counter =~ /\d+\.?\d*/) {
			push @ret, $counter;
		}
	}
	return (@ret);
}


sub GenerateUFFParms {
	my ($atoms,$bonds,$parms) = @_;
	my ($i,$j,$k,$l,@aList,$uff,$bo,$rbo,$ren,$rij,$rjk);
	my ($s1,$s2,$index,$n,$ang,$ang_rad,$fc,$curr,$inv);
	my ($rik,$beta,$kij,$Kijk,$tor,$n1,$n2,$n3,$nbond);

	for $i (keys %{ $parms->{PARMS}{UFF} }) {
		delete $parms->{PARMS}{UFF}{$i} if (! exists($parms->{ATOMTYPES}{$i}) or (! $parms->{ATOMTYPES}{$i}{USED}));
	}
	@aList = keys %{ $parms->{PARMS}{UFF} };
	$uff = $parms->{PARMS}{UFF};
	$bo = getBondOrders($atoms, $bonds);
	for $i (@aList) {
		for $j (@aList) {
			next if (! exists($bo->{$i}{$j}));
			$ang = ();
			$rbo = -0.1332*($uff->{$i}{radius} + $uff->{$j}{radius})*log($bo->{$i}{$j}); #see eqn 3 of uff paper
			$ren = $uff->{$i}{radius}*$uff->{$j}{radius}*(sqrt($uff->{$i}{elecneg})-sqrt($uff->{$j}{elecneg}))**2 /
				($uff->{$i}{elecneg}*$uff->{$i}{radius} + $uff->{$j}{elecneg}*$uff->{$j}{radius}); #eqn 4
			$rij = $uff->{$i}{radius} + $uff->{$j}{radius} + $rbo + $ren; #eqn 2
			$kij = 664.12*$uff->{$i}{zstar}*$uff->{$j}{zstar}/($rij**3); #eqn 6
			$ang_rad = $uff->{$j}{angle}*PI/180;
			$n = ();
			$ang->{lammps} = "cosine/periodic";
			$ang->{type} = "COS_PERIOD";
			if ($uff->{$j}{angle} == 180) { $ang->{k_factor} = 1/4;  $ang->{n} = 2; $ang->{b} =  1; } #linear = 180 }
			elsif ($uff->{$j}{angle} == 120) { $ang->{k_factor} = 1/9;  $ang->{n} = 3; $ang->{b} = -1; } #trigonal-planar
			elsif ($uff->{$j}{angle} ==  90) { $ang->{k_factor} = 1/16; $ang->{n} = 4; $ang->{b} =  1; } #square-planar
			elsif ($uff->{$j}{angle} ==  60) { $ang->{k_factor} = 1/36; $ang->{n} = 6; $ang->{b} =  1; } #octahedral
			else {
				$ang->{k_factor} = 1/(2*sin($ang_rad)**2); #using identity (cos(theta)-cos(theta0)**2/(sin(2*theta0)) = C0 + C1*cos(theta)+C2*cos(2*theta))
				$ang->{b} = $uff->{$j}{angle};
				$ang->{lammps} = "cosine/squared";
				$ang->{type} = "COS_HARMON";
			}
			for $k (@aList) {
				next if (! exists($bo->{$j}{$k}));
				$rbo = -0.1332*($uff->{$j}{radius} + $uff->{$k}{radius})*log($bo->{$j}{$k});
				$ren = $uff->{$j}{radius}*$uff->{$k}{radius}*(sqrt($uff->{$j}{elecneg})-sqrt($uff->{$k}{elecneg}))**2 /
				($uff->{$j}{elecneg}*$uff->{$j}{radius} + $uff->{$k}{elecneg}*$uff->{$k}{radius}); 
				$rjk = $uff->{$j}{radius} + $uff->{$k}{radius} + $rbo + $ren;
				$rik = sqrt($rij**2 + $rjk**2 - 2*$rij*$rjk*cos($ang_rad)); 
				for $l (@aList) {
					next if (! exists($bo->{$k}{$l}));
					next if (exists($parms->{TORSIONS}{$l}) and 
						exists($parms->{TORSIONS}{$l}{$k}) and
						exists($parms->{TORSIONS}{$l}{$k}{$j}) and
						exists($parms->{TORSIONS}{$l}{$k}{$j}{$i}));
					#generate torsions
					$curr = \%{ $parms->{TORSIONS}{$i}{$j}{$k}{$l} };
					$index = 1;
					$index = scalar(keys %{ $curr }) + 1 if(keys %{ $curr });
					$tor = ();
					$n1 = $parms->{ATOMTYPES}{$j}{NUMBONDS}; #hybridization
					$n2 = $parms->{ATOMTYPES}{$k}{NUMBONDS};
					if ($n1 == $n2 and $n2 == 3) { #case 1: X -- sp3 -- sp3 -- X
						$tor->{v0} = sqrt($uff->{$j}{torbar}*$uff->{$k}{torbar}); $tor->{n} = 3; $tor->{d} = 1; #eqn 16
					} elsif ($n1 == $n2 and $n2 == 2) { #case 2: X -- sp2 --- sp2 -- X
						$tor->{v0} = 5*sqrt($uff->{$j}{torbar}*$uff->{$k}{torbar})*(1+4.18*log($bo->{$j}{$k}));
						$tor->{n} = 2; $tor->{d} = -1; #eqn 17
					} elsif (($n1 == 2 and $n2 == 3) or ($n1 == 3 and $n2 == 2)) {
						$tor->{v0} = 2; $tor->{n} = 3; $tor->{d} = -1;
					}
					#special cases for torsions involving group 6 elements
					$n1 = isgroup6($parms, $j);
					$n2 = isgroup6($parms, $k);
					if($n1 or $n2) {
						$tor->{v0} = 5*sqrt($uff->{$j}{torbar}*$uff->{$k}{torbar})*(1+4.18*log($bo->{$j}{$k}));
						$tor->{n} = 2; $tor->{d} = 1;
						if ($n1 == $n2) {
							$n1 = $n2 = 6.8; #set torbar = 6.8 for all group 6 elements
							$n1 = 2 if ($parms->{ATOMTYPES}{$j}{ELENUM} == 6); #special case for oxygen
							$n2 = 2 if ($parms->{ATOMTYPES}{$k}{ELENUM} == 6); #special case for oxygen
							$tor->{v0} = sqrt($n1*$n2); $tor->{n} = 2; $tor->{d} = 1;
						}
					}
					next if (!$tor->{v0});
					$curr->{$index} = (
										{
										"INDEX"			=> $index,
										"TYPE"			=> "HARMONIC",
										"VALS"			=> [$tor->{v0}/2,$tor->{d},$tor->{n}],
										"USED"			=> 0,
										"KEY"			=> "$i $j $k $l ",
										"PER"			=> 4,
										"1_4scale"		=> $parms->{PARMS}{scale_vdw_14},
										"CTYPE"			=> "",
										"do_scale"		=> $parms->{PARMS}{scale_torsions},
										"COUL"          => "0011",
										"VDW"           => "0011",
										}
									);
					$curr->{$index}{Lammps}{name} = "harmonic"; $curr->{$index}{Lammps}{opts} = "";
					$curr->{$index}{NUM}++;
				}
				# generate bond angles
				next if (exists($parms->{ANGLES}{$k}) and
					exists($parms->{ANGLES}{$k}{$j}) and
					exists($parms->{ANGLES}{$k}{$j}{$i}));
				$beta = 664.12/$rij/$rjk;
				$fc = $beta*$uff->{$i}{zstar}*$uff->{$k}{zstar}*$rij*$rjk*(3*$rij*$rjk*
					((1-cos($ang_rad)**2)-$rik*$rik*cos($ang_rad)))/$rik**5; #eqn 13
				$curr = \%{ $parms->{ANGLES}{$i}{$j}{$k} };
				$index = 1;
				$index = scalar(keys %{ $curr }) + 1 if (keys %{ $curr });
				$Kijk = $ang->{k_factor} * $fc/2;
				@{ $ang->{vals} } = ($Kijk,$ang->{b});
				push @{ $ang->{vals} }, $ang->{n} if ($ang->{lammps} eq "cosine/periodic");
				$curr->{$index} = (
									{
									"INDEX"		=> $index,
									"TYPE"		=> $ang->{type},
									"VALS"		=> [@{ $ang->{vals} }],
									"USED"		=> 0,
									"KEY"		=> "$i $j $k ",
									"CTYPE"		=> "",
									}
								);
				$curr->{$index}{Lammps}{name} = $ang->{lammps};
				$curr->{$index}{Lammps}{opts} = "";
			}
			# generate bonds
			next if (exists($parms->{BONDS}{$j}) and 
				exists($parms->{BONDS}{$j}{$i}));
			$curr = \%{ $parms->{BONDS}{$i}{$j} };
			$index = 1;
			$index = scalar(keys %{ $curr }) + 1 if(keys %{ $curr });
			$curr->{$index} = (
								{
								"INDEX"	=> $index,
								"TYPE"	=> "HARMONIC",
								"VALS"	=> [$kij/2,$rij],
								"USED"	=> 0,
								"KEY"	=> "$i $j ",
								}
							);
			$curr->{$index}{Lammps}{name} = "harmonic"; $curr->{$index}{Lammps}{opts} = "";
			print "";
		}
		$nbond = 0;
		for (keys %{ $bo->{$i} }) {
			#$nbond += ceil($bo->{$i}{$_});
		}
		next if ($uff->{$i}{uenerg}==0 or $nbond < 3);
		# generate inversions
		$inv = genInvList($bo->{$i});
		for $j (keys %{ $inv }) {
			for $k (keys %{ $inv->{$j} }) {
				for $l (keys %{ $inv->{$j}{$k} }) {
					$curr = \%{ $parms->{INVERSIONS}{$i}{$j}{$k}{$l} };
					$index = 1;
					$index = scalar(keys %{ $curr }) + 1 if (keys %{ $curr });
					$curr->{$index} = (
										{
										"INDEX"	=> $index,
										"TYPE"	=> "UMBRELLA",
										"VALS"	=> [$uff->{$i}{uenerg},$uff->{$i}{uang}],
										"USED"	=> 0,
										"KEY"	=> "$i $j $k $l",
										"CTYPE"			=> "",
										}
									);
					$curr->{$index}{Lammps}{name} = "umbrella"; $curr->{$index}{Lammps}{opts} = "";
				}
			}
		}
	}
}

sub getBondOrders {
	my ($atoms, $bonds) = @_;
	my ($bondorders, $i, $j, $k, $typei, $typej, $order);

	for $i (keys %{ $atoms }) {
		$typei = $atoms->{$i}{FFTYPE};
		next if (!defined($bonds->{$i}) or ! @{ $bonds->{$i} });
		for $j (0 .. $#{ $bonds->{$i} }) {
			$k = $bonds->{$i}[$j];
			$typej = $atoms->{$k}{FFTYPE};
			next if (exists($bondorders->{$typei}{$typej}) and defined($bondorders->{$typei}{$typej}));
			$order = 1;
			$order = $atoms->{$i}{ORDER}[$j] if(exists($atoms->{$i}{ORDER}) and $#{ $atoms->{$i}{ORDER} } <= $j);
			$order = 1.41 if (($atoms->{$i}{FFTYPE} eq "C_R" and $atoms->{$k}{FFTYPE} eq "N_R") or
				($atoms->{$i}{FFTYPE} eq "N_R" and $atoms->{$k}{FFTYPE} eq "C_R"));
			$order = 1.5 if($atoms->{$i}{FFTYPE} eq "C_R" and $atoms->{$k}{FFTYPE} eq "C_R");
			$order = 1 if (!defined($order));
			$bondorders->{$typei}{$typej} = $bondorders->{$typej}{$typei} = $order;
		}
	}

	return $bondorders;
}

sub isgroup6 {
	my ($parms, $atomtype) = @_;
	my ($curr);

	$curr = $parms->{ATOMTYPES}{$atomtype}{ELENUM};
	return 0 if (! defined($curr));
	return 1 if ($curr == 8 or $curr == 16 or $curr == 34 or $curr == 52 or $curr == 84 or $curr == 116);
	return 0;
}

1;

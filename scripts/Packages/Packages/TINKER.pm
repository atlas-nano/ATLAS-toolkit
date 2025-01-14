package TINKER;

require Exporter;
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use Math::Trig qw(pi);
use Storable qw(dclone); 
use General qw(FileTester FindElementByMass LoadConverter LoadElements);
use CERIUS2 qw(getLammpsOpts findDuplicate);

use strict;

our (@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(parseTinkerFF);
$VERSION = "1.00";

sub numerically { ($a<=>$b); }

sub parseTinkerFF {
	my ($ff_file, $alter, $oldFF) = @_;
	my ($validFF, $CNV, $elements, %PARMS, $headers, $in_data, @tmp, $TYPE_MAP, $rec);
	my ($i, $j, $k, $l, $m, $n, $c, $curr, @vals, $pflag, $ubflag, $r2d, $tortors, $multipole);
	my ($lmp1, $lmp2, $lmp3, $lmp4, $lmp5, $lmp6, $nxny, $eleNum, $type_counter, $MPOLE, $descrp);

	$headers = "forcefield|bond-(cub|quart)ic|angle-(cub|quart|pent|sext)ic|opbendtype|opbend-(cub|quart|pent|sext)ic|" . 
				"torsionunit|vdwtype|radiusrule|radiustype|radiussize|epsilonrule|dielectric|polarization|vdw-1[2-5]-scale|" .
				"mpole-1[2-5]-scale|polar-1[2-5]-scale|polar-14-intra|direct-1[1-4]-scale|mutual-1[1-4]-scale";
	$r2d = 180/pi;
	$CNV = General::LoadConverter();
	$elements = General::LoadElements();
	$tortors = $multipole = 0;
	$MPOLE = (
				{
					"1" => {
						"n"  => 3,
						"l" => "DIPOLE"
					},
					"2" => {
						"n" => 1,
						"l" => "QUADRUPOLE",
					},
					"3" => {
						"n" => 2,
						"l" => "QUADRUPOLE",
					},
					"4" => {
						"n" => 3,
						"l" => "QUADRUPOLE",
					},
				}
	);
	if(defined($oldFF)) {
		%PARMS = %{ $oldFF };
	} else {
		$PARMS{PARMS} = (
							{
								'EXO_CYCLIC'							=> 1.00000,
								'HAS_HBONDS'							=> 0,
								'QEq'									=> 0,
								'SW'									=> 0,
								'USE_HBOND'								=> 0,
								'coul_accuracy'							=> 0.00001,
								'cut_coul'								=> 10.00000,
								'cut_vdw'								=> 10.00000,
								'dielectric'							=> 1.00000000000000000,
								'hbond_angle'							=> 90,
								'hbond_distance'						=> 2.5,
								'mix_rule'								=> 'geometric',
								'same_scale'							=> 1,
								'scale_cou'								=> 0.50000,
								'scale_cou_12'							=> 0.50000,
								'scale_cou_13'							=> 0.50000,
								'scale_cou_14'							=> 0.50000,
								'scale_torsions'						=> 0,
								'scale_torsions_by_n_defined_torsions'	=> 0,
								'scale_vdw'								=> 0.50000,
								'scale_vdw_12'							=> 0.50000,
								'scale_vdw_13'							=> 0.50000,
								'scale_vdw_14'							=> 0.50000,
								'single_inversion'						=> 1,
							}
						);
	}

	open FORCEFIELD, $ff_file->{FF} or die "Cannot open force field file $ff_file->{FF}: $!\n";
	while (<FORCEFIELD>) {
		chomp;
		$in_data = $_;
		if(@tmp = $in_data =~ /^\s*($headers)\s+(\S+)/) {
			$PARMS{TINKER}{HEADERS}{$tmp[0]} = $tmp[$#tmp];
			$tortors = $multipole = 0;
		} elsif($in_data =~ /^\s*atom /) { #atom types
			$tortors = $multipole = 0;
			$in_data =~ s/^\s*//;
			$in_data =~ /.*("[^*]*").*/;
			$descrp = $1;
			$descrp =~ s/\"//g;
			$descrp =~ s/^\s+//;
			$descrp =~ s/\s+$//;
			@tmp = split /\s+"[^"]+"\s+|\s+/, $in_data;
			next if ($#tmp < 4);
			$eleNum = General::FindElementByMass($tmp[5], $elements);
			$rec = (
						{
							"TYPEID"        => $tmp[1],
							"ELEMENT"       => $eleNum,
							"ATOM"          => $elements->{$eleNum}{SYMBOL},
							"MASS"          => $tmp[5],
							"USE_CHARGE"    => 1,
							"USED"          => 0,
							"LABEL"         => $tmp[3],
							"CHARGE"        => 0,
							"TINKER_ID"     => $tmp[1],
							"TINKER_PRM"    => $ff_file->{FF},
							"AMOEBA_FLAG"   => 0,
							"TINKER_DESCRP" => $descrp,
							"TINKER_NBOND"  => $tmp[6],
							}
			);
			$in_data =~ /"([^"]*)"/;
			$rec->{OTHER} = $1;
			$PARMS{ATOMTYPES}{$tmp[3]} = $rec;
			$TYPE_MAP->{$tmp[1]} = $tmp[3];
		} elsif ($in_data =~ /^\s*(vdw|vdwpr) /) { #vdw
			$tortors = $multipole = 0;
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/,$in_data;
			$i = getAtmType($tmp[1],$TYPE_MAP);
			if($in_data =~ /^\s*vdwpr/) {
				$j = getAtmType($tmp[2],$TYPE_MAP);
				shift @tmp;
			} else {
				$j = $i;
			}
			@vals = ($tmp[3],$tmp[2]);
			push @vals, $tmp[4] if ($#tmp > 3);
			$PARMS{VDW}{$i}{$j}{1} = (
										{
											"TYPE"    => "AMOEBA",
											"Lammps"  => getLammpsOpts("AMOEBA","vdw", $CNV,0,$PARMS{PARMS}),
											"KEY"     => "$i $j ",
											"VALS"    => [@vals],
											"ATOM"    => "$i $j ",
											"USED"    => 0,
											"IGNORE"  => 0,
											"IT"      => "vdw",
										}
									);
		} elsif($in_data =~ /^\s*bond /) {
			$tortors = $multipole = 0;
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/,$in_data;
			$i = getAtmType($tmp[1],$TYPE_MAP); $j = getAtmType($tmp[2],$TYPE_MAP);
			if(exists($PARMS{BONDS}{$j}) and exists($PARMS{BONDS}{$j}{$i})) {
				$curr=\%{$PARMS{BONDS}{$j}{$i}};
			} else {
				$curr=\%{$PARMS{BONDS}{$i}{$j}};
			}
			$type_counter = findDuplicate($curr,"CLASS2");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			$lmp1 = $tmp[4];
			$lmp2 = $tmp[3];
			&checkHeader(\%PARMS, "bond-cubic",   $ff_file->{FF}); 
			&checkHeader(\%PARMS, "bond-quartic", $ff_file->{FF});
			$lmp3 = $PARMS{TINKER}{HEADERS}{"bond-cubic"}  *$tmp[3];
			$lmp4 = $PARMS{TINKER}{HEADERS}{"bond-quartic"}*$tmp[3];
			$curr->{$type_counter} = (
										{
											"INDEX"   => $type_counter,
											"TYPE"    => "CLASS2",
											"Lammps"  => getLammpsOpts("CLASS2","bond", $CNV),
											"VALS"    => [$lmp1, $lmp2, $lmp3, $lmp4],
											"USED"    => 0,
											"KEY"     => "$i $j ",
											"IGNORE"  => 0,
										}
									);
		} elsif($in_data =~ /^\s*(angle|anglep)/) {
			$tortors = $multipole = 0;
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/, $in_data;
			$i = getAtmType($tmp[1], $TYPE_MAP); $j = getAtmType($tmp[2],$TYPE_MAP); $k = getAtmType($tmp[3],$TYPE_MAP);
			if(exists($PARMS{ANGLES}{$k}) and exists($PARMS{ANGLES}{$k}{$j}) and exists($PARMS{ANGLES}{$k}{$j}{$i})) {
				$curr=\%{ $PARMS{ANGLES}{$k}{$j}{$i} };
			} else {
				$curr=\%{ $PARMS{ANGLES}{$i}{$j}{$k} };
			}
			$validFF = 1;
			$pflag = $ubflag = 0;
			$pflag = 1 if ($in_data =~ /anglep/);
			$type_counter = findDuplicate($curr,"AMOEBA");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			$curr->{$type_counter} = (
										{
											"INDEX"   => $type_counter,
											"TYPE"    => "AMOEBA",
											"Lammps"  => getLammpsOpts("AMOEBA","angle", $CNV),
											"USED"    => 0,
											"KEY"     => "$i $j $k ",
											"IGNORE"  => 0,
											"CTYPE"   => "",
											"TINKER"  => 1,
										}
			);
			$lmp1 = $tmp[5];
			$lmp2 = $tmp[4];
			&checkHeader(\%PARMS, "angle-cubic",   $ff_file->{FF});
			&checkHeader(\%PARMS, "angle-quartic", $ff_file->{FF});
			&checkHeader(\%PARMS, "angle-pentic",  $ff_file->{FF});
			&checkHeader(\%PARMS, "angle-sextic",  $ff_file->{FF});    
			$lmp3 = $PARMS{TINKER}{HEADERS}{"angle-cubic"}   * $tmp[4] * $r2d;
			$lmp4 = $PARMS{TINKER}{HEADERS}{"angle-quartic"} * $tmp[4] * $r2d * $r2d;
			$lmp5 = $PARMS{TINKER}{HEADERS}{"angle-pentic"}  * $tmp[4] * $r2d * $r2d * $r2d;
			$lmp6 = $PARMS{TINKER}{HEADERS}{"angle-sextic"}  * $tmp[4] * $r2d * $r2d * $r2d * $r2d;

			#check to see if urey_bradley is defined for current angle
			$ubflag = 1 if(exists($curr->{$type_counter}{UREY_BRADLEY}));

			$curr->{$type_counter}{VALS} = [ $pflag, $ubflag, $lmp1, $lmp2, $lmp3, $lmp4, $lmp5, $lmp6 ];
			if ($#tmp > 5) { #bond angle
				push @{ $curr->{$type_counter}{VALS} }, [ $pflag, $ubflag, $tmp[6], $lmp2, $lmp3, $lmp4, $lmp5, $lmp6 ];
			} 
			if ($#tmp > 6 ) { #urey bradley 
				push @{ $curr->{$type_counter}{VALS} }, [ $pflag, $ubflag, $tmp[7], $lmp2, $lmp3, $lmp4, $lmp5, $lmp6 ];
			}
		} elsif ($in_data =~ /^\s*strbnd /) {
			$tortors = $multipole = 0;
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/, $in_data;
			$i = getAtmType($tmp[1], $TYPE_MAP); $j = getAtmType($tmp[2],$TYPE_MAP); $k = getAtmType($tmp[3],$TYPE_MAP);
			if(exists($PARMS{ANGLES}{$k}) and exists($PARMS{ANGLES}{$k}{$j}) and exists($PARMS{ANGLES}{$k}{$j}{$i})) {
				$curr=\%{ $PARMS{ANGLES}{$k}{$j}{$i} };
			} else {
				$curr=\%{ $PARMS{ANGLES}{$i}{$j}{$k} };
			}

			$type_counter = findDuplicate($curr,"AMOEBA");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			$curr->{$type_counter}{BOND_ANGLE}{TYPE} = "CLASS2/P6";
			$curr->{$type_counter}{BOND_ANGLE}{TYPE} = getLammpsOpts("CLASS2/P6","angle", $CNV);
			$lmp1 = $tmp[5];
			$lmp2 = $tmp[4];

			#i - j bond
			if(exists($PARMS{BONDS}{$i}) and exists($PARMS{BONDS}{$i}{$j})) {
				$type_counter = findDuplicate($PARMS{BONDS}{$i}{$j}, "CLASS2");
				$lmp3 = $PARMS{BONDS}{$i}{$j}{$type_counter}{VALS}[0];
			} elsif (exists($PARMS{BONDS}{$j}) and exists($PARMS{BONDS}{$j}{$i})) {
				$type_counter = findDuplicate($PARMS{BONDS}{$j}{$i}, "CLASS2");
				$lmp3 = $PARMS{BONDS}{$j}{$i}{$type_counter}{VALS}[0];
			} else {
				die "ERROR: Cannot find bond parameters for $i - $i\n";
			}
			#j - l bond
			if(exists($PARMS{BONDS}{$j}) and exists($PARMS{BONDS}{$j}{$k})) {
				$type_counter = findDuplicate($PARMS{BONDS}{$j}{$k}, "CLASS2");
				$lmp4 = $PARMS{BONDS}{$j}{$k}{$type_counter}{VALS}[0];
			} elsif (exists($PARMS{BONDS}{$k}) and exists($PARMS{BONDS}{$k}{$j})) {
				$type_counter = findDuplicate($PARMS{BONDS}{$k}{$j}, "CLASS2");
				$lmp4 = $PARMS{BONDS}{$k}{$j}{$type_counter}{VALS}[0];
			} else {
				die "ERROR: Cannot find bond parameters for $j - $k\n";
			}

			@{ $curr->{$type_counter}{BOND_ANGLE}{VALS} } = ( $lmp1, $lmp2, $lmp3, $lmp4 );
		} elsif ($in_data =~ /^\s*ureybrad /) {
			$tortors = $multipole = 0;
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/, $in_data;
			$i = getAtmType($tmp[1], $TYPE_MAP); $j = getAtmType($tmp[2],$TYPE_MAP); $k = getAtmType($tmp[3],$TYPE_MAP);
			if(exists($PARMS{ANGLES}{$k}) and exists($PARMS{ANGLES}{$k}{$j}) and exists($PARMS{ANGLES}{$k}{$j}{$i})) {
				$curr=\%{ $PARMS{ANGLES}{$k}{$j}{$i} };
			} else {
				$curr=\%{ $PARMS{ANGLES}{$i}{$j}{$k} };
			}
			$type_counter = findDuplicate($curr,"AMOEBA");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			$curr->{$type_counter}{UREY_BRADLEY}{TYPE} = "UREY_BRADLEY";
			$curr->{$type_counter}{UREY_BRADLEY}{Lammps} = getLammpsOpts("UREY_BRADLEY","angle", $CNV);
			@{ $curr->{$type_counter}{UREY_BRADLEY}{VALS} } = ($tmp[4], $tmp[5]);
			#now find the angle equivalent and update the ubflag indicator
			if(exists($curr->{$type_counter}{VALS})) {
				$i = 1;
				while($i<$#{ $curr->{$type_counter}{VALS} }) {
					$curr->{$type_counter}{VALS}[$i] = 1;
					$i += 8;
				}
			}
		} elsif($in_data =~ /^\s*torsion /) { #dihedral (proper) torsion
			$tortors = $multipole = 0;
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/, $in_data;
			die "ERROR: torsion line is invalid! \"$in_data\"\n" if($#tmp < 7);
			die "ERROR: Expected torsion i j k l [eng angle period]n. Got \"$in_data\"\n" if (($#tmp - 4) %3 > 0);
			$i = getAtmType($tmp[1],$TYPE_MAP); $j = getAtmType($tmp[2],$TYPE_MAP); 
			$k = getAtmType($tmp[3],$TYPE_MAP); $l = getAtmType($tmp[4],$TYPE_MAP);
			if(exists($PARMS{TORSIONS}{$l}) and exists($PARMS{TORSIONS}{$l}{$k}) and exists($PARMS{TORSIONS}{$l}{$k}{$j}) and exists($PARMS{TORSIONS}{$l}{$k}{$j}{$i}) and exists($PARMS{TORSIONS}{$l}{$k}{$j}{$i})) {
				$curr=\%{ $PARMS{TORSIONS}{$l}{$k}{$j}{$i} };
			} else {
				$curr=\%{ $PARMS{TORSIONS}{$i}{$j}{$k}{$l} };
			}
			$type_counter = findDuplicate($curr,"FOURIER");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			$curr->{$type_counter} = (
										{
											"INDEX"   => $type_counter,
											"TYPE"    => "FOURIER",
											"Lammps"  => getLammpsOpts("FOURIER","dihedral",$CNV),
											"USED"    => 0,
											"KEY"     => "$i $j $k $l ",
											"IGNORE"  => 0,
											"CTYPE"   => "",
											"NUM"	  => 1,
											"COUL"    => "11111",
											"VDW"     => "11111",
										}
									);
			$curr->{$type_counter}{VALS} = [($#tmp - 4)/3]; #ntorsion
			$m = 5;
			while ($m<$#tmp) {
				$lmp1 = $tmp[$m];
				$lmp2 = $tmp[$m+2];
				$lmp3 = $tmp[$m+1];
				push @{ $curr->{$type_counter}{VALS} }, ( $lmp1, $lmp2, $lmp3 );
				$m +=3;
			}
		} elsif($in_data =~ /^\s*opbend /) { #improper torsion
			$tortors = $multipole = 0;
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/,$in_data;
			$i = getAtmType($tmp[1],$TYPE_MAP); $j = getAtmType($tmp[2],$TYPE_MAP); 
			$k = getAtmType($tmp[3],$TYPE_MAP); $l = getAtmType($tmp[4],$TYPE_MAP);
			if(exists($PARMS{INVERSIONS}{$j})) {
				$curr = \%{ $PARMS{INVERSIONS}{$j} };
				if(exists($curr->{$i}) and exists($curr->{$i}{$k}) and exists($curr->{$i}{$k}{$l})) {
					$curr = \%{ $curr->{$i}{$k}{$l} };
				}elsif(exists($curr->{$i}) and exists($curr->{$i}{$l}) and exists($curr->{$i}{$l}{$k})) {
					$curr = \%{ $curr->{$i}{$l}{$k} };
				}elsif(exists($curr->{$k}) and exists($curr->{$k}{$i}) and exists($curr->{$k}{$i}{$l})) {
					$curr = \%{ $curr->{$k}{$i}{$l} };
				}elsif(exists($curr->{$k}) and exists($curr->{$k}{$l}) and exists($curr->{$k}{$l}{$i})) {
					$curr = \%{ $curr->{$k}{$l}{$i} };
				}elsif(exists($curr->{$l}) and exists($curr->{$l}{$i}) and exists($curr->{$l}{$i}{$k})) {
					$curr = \%{ $curr->{$l}{$i}{$k} }
				}elsif(exists($curr->{$l}) and exists($curr->{$l}{$k}) and exists($curr->{$l}{$k}{$i})) {
					$curr = \%{ $curr->{$l}{$i}{$k} }
				}else {
					$curr = \%{ $curr->{$i}{$k}{$l} };
				}
			} else {
				$curr=\%{ $PARMS{INVERSIONS}{$j}{$i}{$k}{$l} };
			}
			$type_counter = findDuplicate($curr,"AMOEBA");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			$curr->{$type_counter} = (
										{
											"INDEX"   => $type_counter,
											"TYPE"    => "AMOEBA",
											"Lammps"  => getLammpsOpts("AMOEBA","inversion", $CNV),
											"VALS"    => [$tmp[5]],
											"USED"    => 0,
											"KEY"     => "$j $i $k $l ",
											"IGNORE"  => 0,
											"CTYPE"   => "",
										}
									);
		} elsif($in_data =~ /^\s*pitors /) {
			$tortors = $multipole = 0;
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/,$in_data;
			$i = getAtmType($tmp[1],$TYPE_MAP); $j = getAtmType($tmp[2],$TYPE_MAP);
			if(exists($PARMS{PI_TORSIONS}{$j}) and exists($PARMS{PI_TORSIONS}{$j}{$i})) {
				$curr=\%{$PARMS{PI_TORSIONS}{$j}{$i}};
			} else {
				$curr=\%{$PARMS{PI_TORSIONS}{$i}{$j}};
			}
			$type_counter = findDuplicate($curr,"PITORSION");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			$curr->{$type_counter} = (
										{
											"INDEX"   => $type_counter,
											"TYPE"    => "PITORSION",
											"Lammps"  => getLammpsOpts("PITORSION","fix", $CNV),
											"VALS"    => [$tmp[3]],
											"USED"    => 0,
											"KEY"     => "$i $j ",
											"IGNORE"  => 0,
										}
									);
		} elsif($in_data =~ /^\s*tortors /) {
			$tortors = 1; $multipole = 0;
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/, $in_data;
			die "ERROR: tortors line is invalid! \"$in_data\"\n" if($#tmp < 7);
			shift @tmp;
			$i = getAtmType(shift @tmp,$TYPE_MAP); $j = getAtmType(shift @tmp,$TYPE_MAP); $k = getAtmType(shift @tmp,$TYPE_MAP); 
			$l = getAtmType(shift @tmp,$TYPE_MAP); $m = getAtmType(shift @tmp,$TYPE_MAP);
			if(exists($PARMS{BI_TORSIONS}{$m}) and exists($PARMS{BI_TORSIONS}{$m}{$l}) and exists($PARMS{BI_TORSIONS}{$m}{$l}{$k}) 
			and exists($PARMS{BI_TORSIONS}{$m}{$l}{$k}{$j}) and exists($PARMS{BI_TORSIONS}{$m}{$l}{$k}{$j}{$i}) and exists($PARMS{BI_TORSIONS}{$m}{$l}{$k}{$j}{$i})) {
				$curr=\%{ $PARMS{BI_TORSIONS}{$m}{$l}{$k}{$j}{$i} };
			} else {
				$curr=\%{ $PARMS{BI_TORSIONS}{$i}{$j}{$k}{$l}{$m} };
			}
			$type_counter = findDuplicate($curr,"BITORSION");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			$curr->{$type_counter} = (
										{
											"INDEX"   => $type_counter,
											"TYPE"    => "BITORSION",
											"Lammps"  => getLammpsOpts("BITORSION","fix", $CNV),
											"USED"    => 0,
											"KEY"     => "$i $j $k $l $m ",
											"IGNORE"  => 0,
											"NX"      => $tmp[0],
											"NY"      => $tmp[1],
										}
									);
			$nxny = $tmp[0]*$tmp[1];
			$n = 0;
		} elsif ($tortors == 1 and $_ =~ /^\s*(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
			@{ $curr->{$type_counter}{VALS}[$n++] } = ([$1,$2,$3]);
			die "ERROR: Read more nx,ny data for tortor ($n) than expected ($nxny)!\n" if ($n > $nxny);
		} elsif ($tortors == 1 and $_ =~ /^\s*$/) {
			die "ERROR: Read less nx,ny data for tortor ($n) than expected ($nxny)!\n" if ($n < $nxny);
		} elsif($in_data =~ /^\s*multipole /) {
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/, $in_data;
			die "ERROR: multipole line is invalid!\n" if ($#tmp < 2);
			#check previous entry for validity
			if($multipole > 0 and defined($curr)) {
				die "ERROR: Invalid multipole line for $curr->{LABEL}!\n"
					if(!exists($curr->{MULTIPOLE}{DIPOLE}) or ! exists($curr->{MULTIPOLE}{QUADRUPOLE}));
			}
			$tortors = 0; $multipole = 1;
			shift @tmp;
			$i = getAtmType(shift @tmp, $TYPE_MAP);
			$curr = \%{ $PARMS{ATOMTYPES}{ $i} };
			$curr->{MULTIPOLE}{MONOPOLE} = pop @tmp; 
			$curr->{MULTIPOLE}{AXIS} = ();
			for $j (@tmp) {
				if ($j == 0) {
					push @{ $curr->{MULTIPOLE}{AXIS} }, 1;
					push @{ $curr->{MULTIPOLE}{ATMTYPES } }, "X";
				} elsif($j < 0) {
					push @{ $curr->{MULTIPOLE}{AXIS} }, -1;
					push @{ $curr->{MULTIPOLE}{ATMTYPES} }, getAtmType(-$j,$TYPE_MAP);
				} else {
					push @{ $curr->{MULTIPOLE}{AXIS} }, 1;
					push @{ $curr->{MULTIPOLE}{ATMTYPES} }, getAtmType($j,$TYPE_MAP);
				}
			}
		} elsif($multipole > 0 and $multipole < 5 and defined($curr)) { #multipole
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/, $in_data;
			die "ERROR: multipole line for $curr->{LABEL} is invalid for specifiying dipole (line " . ($multipole+1) . "). " . 
				"Expected " . $MPOLE->{$multipole}{n} . " floats. got \"$in_data\"!\n" if(($#tmp + 1) != $MPOLE->{$multipole}{n});
			for $i (@tmp) {
				die "ERROR: Expected xx.xxx for each multipole value. Got \"$i\"!\n" 
					if ($i !~ /^\-?\d+\.\d+/);
			}	
			push @{ $curr->{MULTIPOLE}{ $MPOLE->{$multipole}{l} } }, @tmp;
			$multipole++;	
		} elsif($in_data =~ /^\s*polarize /) {
			$tortors = $multipole = 0;
			$in_data =~ s/^\s*//;
			@tmp = split /\s+/, $in_data;
			die "ERROR: polarize line is invalid!\n" if ($#tmp < 3);
			shift @tmp;
			$i = getAtmType(shift @tmp, $TYPE_MAP);
			$PARMS{ATOMTYPES}{$i}{AMOEBA_FLAG} = 1;
			$lmp1 = shift @tmp;
			$lmp2 = shift @tmp;
			$curr = \%{ $PARMS{ATOMTYPES}{ $i} };
			$curr->{POLARIZE} = (
									{
										"p1"    => $lmp1,
										"p2"    => $lmp2,
										"ATOMS" => [map { getAtmType($_, $TYPE_MAP) } @tmp],
									}
			);
		}
	}
	close FORCEFIELD;
	return \%PARMS;	
}

sub checkHeader {
	my ($parms, $i, $ff) = @_;

	die "ERROR: Cannot find keyword $i in Tinker forcefield $ff!\n" 
		if (! exists($parms->{TINKER}{HEADERS}{$i}));
}

sub getAtmType {
	my ($type_id, $type_list) = @_;

	if(exists($type_list->{$type_id})) {
		return $type_list->{$type_id};
	} elsif ($type_id < 1) {
		return "X";
	} else {
		die "ERROR: Cannot find corresponding atom type for type index: $type_id\n";
	}
}

1;

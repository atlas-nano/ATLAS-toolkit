#!/usr/bin/perl -w
use strict;
use constant PI => atan2(1,1) * 4;

sub Initialize();
sub ReadParmFile();
sub AssignParms(@);
sub CalcCylDimensions();
sub HexagonalPacking();
sub nint(@);
sub PlaceSolute();
sub PlaceFluid();
sub IsOverLap(@);
sub AddAtm(@);
sub WriteData();
sub Createfile(@);
sub WriteLammpsHeader(@);
sub DetermineXtreme(@);
sub CalcDebye();
sub GetNumCoIons();
sub OrderRes(@);
sub CorrectKeyName(@);
sub CreateLammpsIn();

die "usage: $0 parameterfile [sim_details]\n"
    if (! @ARGV);

my ($parmfile, $sim_data) = @ARGV;
my (%PARMS, $ANG_cnv, %CYLINDER, %DATA, %HEXAGON, $atm_counter, $num_overlap);

$num_overlap = 0;
Initialize();
ReadParmFile();

# Convert charge to Angstroms

$PARMS{"CYLINDER"}{"charge"} *= $ANG_cnv;
CalcCylDimensions();

$atm_counter = 1;
PlaceSolute()
    if ($PARMS{"SOLUTE"}{"incS"});

PlaceFluid();
WriteData();

CreateLammpsIn()
    if ($sim_data);

sub Initialize() {
    $ANG_cnv = 1/(1.6022E-19 * 1.0E20);    # convertion to Angstoms

    if ($sim_data) {
	if (! -e $sim_data or ! -T $sim_data) {
	    $sim_data = 0;
	}
    } else {
	$sim_data = 0;
    }

    die "Cannot locate $parmfile: $!\n"
	if (! -e $parmfile);
    die "Error: $parmfile is not readable\n"
	if (! -r $parmfile or ! -T $parmfile);

    %PARMS = (
	      "SOLVENT"      => {
				 "itypef"       => -9999.99,
				 "qf"           => -9999.99,
				 "rf"           => -9999.99,
				 "mf"           => -9999.99,
				 },
	      "IONS"          => {
				  "itypefion"    => -9999.99,
				  "qfion"        => -9999.99,
				  "rfion"        => -9999.99,
				  "mfion"        => -9999.99,
			      },
	      "COIONS"        =>  {
		                  "incCoIon"     => 0,
		                  "itypecoion"   => -9999.99,
				  "qcoion"       => -9999.99,
				  "rcoion"       => -9999.99,
				  "mcoion"       => -9999.99,
				  "ion_str"      => -9999.99,
			         },
	      "SOLUTE"       => {
		                  "incS"        => 0,
				  "itypeSol"    => -9999.99,
				  "qSol"        => -9999.99,
				  "rSol"        => -9999.99,
				  "buffer"      => -9999.99,
				  "mSol"        => -9999.99,
				 },
	      "CYLINDER"     => {
				 "rhof_targ"   => -9999.99,
				 "sigma"       => -9999.99,
				 "rc"          => -9999.99,
				 "xc"          => -9999.99,
				 "qw"          => -9999.99,
			         },
	      "SIMULATION"   => {
		                 "cut_off"     => -9999.99,
				 "well_depth"  => -9999.99,
				 "alpha_coeff" => -9999.99,
			         },
	      );
}

sub ReadParmFile() {
    my ($header, $pkeys, $invalid);

    open INFILE, $parmfile or die "Cannot read $parmfile: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^\s*(\w+)\s*=\s*(\-?\d+\.?\d*)/) {
	    if (! AssignParms($1,$2)) {
		print "Invalid key: $1\n";
	    }
	}
    }
    close INFILE;
    
    $invalid = 0;
    for $header (keys %PARMS) {
	if ($header ne "SOLUTE" or ($PARMS{"SOLUTE"}{"incS"} == 1)) {
	    for $pkeys (keys %{ $PARMS{$header} }) {
		if ($PARMS{$header}{$pkeys} <= -9999) {
		    die "Error in parameter file. $pkeys is either missing or invalid\n";
		    $invalid = 1;
		    last;
		}
	    }
	}
	if ($invalid) {
	    last;
	}
    }
}

sub AssignParms(@) {
    my ($keynm, $keyval) = @_;
    my ($header, $pkeys, $isvalid);

    $isvalid = 0;
    for $header (keys %PARMS) {
	for $pkeys (keys %{ $PARMS{$header} }) {
	    if ($pkeys eq $keynm) {
		delete $PARMS{$header}{$pkeys};
		$pkeys = CorrectKeyName($pkeys);
		$PARMS{$header}{$pkeys} = $keyval;
		$isvalid = 1;
		last;
	    }
	}
	if ($isvalid) {
	    last;
	}
    }

    return $isvalid;
}

sub CalcCylDimensions() {
    my ($s_area_calc, $num_ions, $s_area_real);
    my ($pos_ions, $add_ions, $adj_vol);

    $adj_vol = 0; # The adjusted volume to compensate for vol of the solute and coions
    $s_area_calc = 2 * PI * $PARMS{"CYLINDER"}{"radii"} * $PARMS{"CYLINDER"}{"xc"}; # 2piR*H
    $num_ions = -1 * nint($s_area_calc * $PARMS{"CYLINDER"}{"charge"}/ $PARMS{"IONS"}{"charge"});
    $s_area_real = -1 * $PARMS{"IONS"}{"charge"} * $num_ions/ $PARMS{"CYLINDER"}{"charge"};
   
    $add_ions = 0;
    $CYLINDER{"COIONS"} = 0;

    die "Solvent and Ions have same charge!!\n"
	if ($num_ions < 0);
    if ($PARMS{"COIONS"}{"incCoIon"}) {
	die "Ions and Counter Ions cannot have the same charge!\n"
	    if ($PARMS{"COIONS"}{"charge"}/$PARMS{"IONS"}{"charge"} > 0);
    }

#   Adjust the number of ions to include the fact that there is a charged solute
    if ($PARMS{"SOLUTE"}{"incS"}) {
	$num_ions -= $PARMS{"SOLUTE"}{"charge"};
	$adj_vol = (2 * PI * $PARMS{"SOLUTE"}{"radii"} * $PARMS{"CYLINDER"}{"xc"})/3;
	$adj_vol += $num_ions * (2 * PI * $PARMS{"IONS"}{"radii"} * $PARMS{"CYLINDER"}{"xc"})/3;
    }


    $CYLINDER{"HEIGHT"} = $s_area_real/(2 * PI * $PARMS{"CYLINDER"}{"radii"});
    $CYLINDER{"RADIUS"} = $PARMS{"CYLINDER"}{"radii"} - ($PARMS{"CYLINDER"}{"sigma"}/2);
    $CYLINDER{"VOL"} = PI * $CYLINDER{"RADIUS"}**2 * $CYLINDER{"HEIGHT"};

    if ($PARMS{"COIONS"}{"incCoIon"}) {
	$PARMS{"COIONS"}{"num"} = GetNumCoIons();
	$add_ions = $PARMS{"COIONS"}{"num"} * abs( $PARMS{"COIONS"}{"charge"} / $PARMS{"IONS"}{"charge"});
	$adj_vol += $PARMS{"COIONS"}{"num"} * (2 * PI * $PARMS{"COIONS"}{"radii"} * $PARMS{"CYLINDER"}{"xc"})/3;
	$adj_vol += $add_ions * (2 * PI * $PARMS{"IONS"}{"radii"} * $PARMS{"CYLINDER"}{"xc"})/3;
    }

    $CYLINDER{"NSOLV"} = nint(($CYLINDER{"VOL"} - $adj_vol) * $PARMS{"CYLINDER"}{"rhof_targ"});
    $CYLINDER{"FLUID_DENSITY"} = $CYLINDER{"NSOLV"}/($CYLINDER{"VOL"} - $adj_vol);
    $CYLINDER{"SAREA"} =  $s_area_real;


    $CYLINDER{"NIONS"} = $num_ions;

#   CYLINDRICAL DIMENSIONS
    $CYLINDER{"DIMENSIONS"} = (
			      {
				  "XHi" => -99.0,
				  "XLo" => 99.0,
				  "YHi" => $PARMS{"CYLINDER"}{"radii"} * 2,
				  "YLo" => -1 * $PARMS{"CYLINDER"}{"radii"} * 2,
				  "ZHi" => $PARMS{"CYLINDER"}{"radii"} * 2,
				  "ZLo" => -1 * $PARMS{"CYLINDER"}{"radii"} * 2,
			      }
			      );

    print "--==STATS===--\n";
    print "Effective Cylindrical radius is " . $CYLINDER{"RADIUS"} . "\n";
    print "Effective Cylindrical length is " . $CYLINDER{"HEIGHT"} . "\n";
    print "Effective Cylindrical Surface Area is " . $CYLINDER{"SAREA"} . "\n";
    print "Effective Cylindrical Volume is " . $CYLINDER{"VOL"} . "\n";
    print "Wall charge density in C/m^2 is " . $PARMS{"CYLINDER"}{"charge"} / $ANG_cnv . "\n";
    print "Wall charge density in Angstroms is " . $PARMS{"CYLINDER"}{"charge"} . "\n";
    print "Number of ions to neutralize wall is " . $CYLINDER{"NIONS"} . "\n";

    if ($PARMS{"COIONS"}{"incCoIon"}) {
	$pos_ions = $CYLINDER{"NIONS"} + $add_ions;
	$CYLINDER{"NIONS"} = $pos_ions;
	$CYLINDER{"COIONS"} = $PARMS{"COIONS"}{"num"};
	print "Total Number of positive ions is " . $pos_ions . "\n";
	print "Total Number of co-ions is " . $PARMS{"COIONS"}{"num"} . "\n";
    }

    CalcDebye();
    HexagonalPacking();
}

sub HexagonalPacking() {
    my ($hex_rank, $atm_counter, $hex_spacing, $hex_num, $hex_layers);
    my ($h_one, $h_zero, $h_half, $h_rt32);

#    find the largest hexagon spacing, as well as the 
    $hex_rank = $atm_counter = $hex_spacing = $hex_num = $hex_layers = 0;

    while ($atm_counter < $CYLINDER{"NSOLV"}) {
	$hex_rank++;
	$hex_spacing = $CYLINDER{"RADIUS"}/$hex_rank;
	$hex_num = nint($CYLINDER{"HEIGHT"}/$hex_spacing);
	$hex_layers = $CYLINDER{"HEIGHT"}/$hex_num;
	$atm_counter = (1 + 6 * $hex_rank * ($hex_rank + 1)/2) * $hex_num;
	printf "%8d%8d%8d%8d\n", $hex_rank, $hex_num, $atm_counter, $CYLINDER{"NSOLV"};
    }

    print "Number of fluid atoms is " . $CYLINDER{"NSOLV"} . "\n";
    print "Maximum number of fluid atoms is " . $atm_counter . "\n";
    print "Rank of hexagon is " . $hex_rank . "\n";
    print "Number of hexagon layers is " . $hex_num . "\n";
    print "Spacing of hexagon sites is " . $hex_spacing . "\n";
    print "Spacing of hexagon layers is " . $hex_layers . "\n";
    
    $h_one = $hex_spacing;
    $h_zero = 0.0;
    $h_half = $hex_spacing/2;
    $h_rt32 = $hex_spacing * sqrt(3.0)/2;

    %HEXAGON = (
		"RANK"          => $hex_rank,
		"NUM_LAYERS"    => $hex_num,
		"LAYER_SPACING" => $hex_layers,
		"RING_SPACING"  => $hex_spacing,
		"YCOORDS"       => [$h_half, $h_one, $h_half, -1*$h_half, -1*$h_one,-1*$h_half],
		"ZCOORDS"       => [$h_rt32, $h_zero, -1*$h_rt32,-1*$h_rt32, $h_zero, $h_rt32], 
		    
		);
}

sub nint(@) {
    return sprintf("%.0f", $_[0]);
}

sub PlaceSolute() {

    my ($rz, $ry, $rx);


# Deposit the solute at the top center of the box.

    die "Error: solvent radius exceeds cylinder radius\n"
	if ( ($PARMS{"SOLUTE"}{"radii"} + $PARMS{"SOLUTE"}{"buffer"}) > ($CYLINDER{"RADIUS"}) );
    
    die "Error: solvent radius exceeds cylinder length\n"
	if ( ($PARMS{"SOLUTE"}{"radii"} + $PARMS{"SOLUTE"}{"buffer"}) > ($CYLINDER{"HEIGHT"}/2) );
    
    $ry = $rz = 0.0;
    $rx = $PARMS{"SOLUTE"}{"radii"} + $PARMS{"SOLVENT"}{"radii"};
    DetermineXtreme($rx);
    
    print "Placed Solute at $rx $ry $rz\n";
    $DATA{"SOLUTE"}[0] = (
			 {
			     "XCOORD" => $rx,
			     "YCOORD" => $ry,
			     "ZCOORD" => $rz,
			     "CHARGE" => $PARMS{"SOLUTE"}{"charge"},
			     "MASS"   => $PARMS{"SOLUTE"}{"mass"},
			     "RADII"  => $PARMS{"SOLUTE"}{"radii"},
			 }
		       );
    $atm_counter++;
}

sub PlaceFluid() {
    my ($ion_counter, $hex_counter, $counter, $layer_counter, $rank_counter, $solv_particles);
    my ($rx, $ry, $rz) = (0,0,0);

    my ($total_atms);

    $total_atms = $CYLINDER{"NSOLV"}  + $CYLINDER{"NIONS"} + $CYLINDER{"COIONS"};
    $total_atms++
	if ($PARMS{"SOLUTE"}{"incS"});
#   Place the first set of ions/solvent as a single line along the middle of the tube
#   After will place rings of ions/fluids around center line in hexagonal packing arrangement

    open OVERLAP, " > ./overlap.dat" or die "Cannot write to regular file overlap.dat: $!\n";

    $ion_counter = $solv_particles = 0;

    for $layer_counter (1 .. $HEXAGON{"NUM_LAYERS"}) {
	if ($ion_counter < $CYLINDER{"NSOLV"}) {
	    if ($ion_counter >= ($CYLINDER{"NIONS"} + $CYLINDER{"COIONS"})) {     # Ran out of ions so place solvent
		$solv_particles += AddAtm(2,$rx,$ry,$rz,$atm_counter);
	    }elsif ($ion_counter >= $CYLINDER{"NIONS"}) { # Ran out of positive ions so place negative
		$ion_counter += AddAtm(1, $rx, $ry, $rz, $atm_counter);
	    } else {                                     # Place an ion
		$ion_counter += AddAtm(0,$rx,$ry,$rz,$atm_counter);
	    }
	    $rx += $HEXAGON{"LAYER_SPACING"};
	    DetermineXtreme($rx);
	}
    }

    for $rank_counter (1 .. $HEXAGON{"RANK"}) {
	$ry += $HEXAGON{"YCOORDS"}[4];
	$rz += $HEXAGON{"ZCOORDS"}[4];
	for $hex_counter (0 .. 5) {
	    for $counter (1 .. $rank_counter) {
		$ry += $HEXAGON{"YCOORDS"}[$hex_counter];
		$rz += $HEXAGON{"ZCOORDS"}[$hex_counter];
		$rx = 0.0;

		for $layer_counter (1 .. $HEXAGON{"NUM_LAYERS"}) {
		    if ($atm_counter <= $total_atms) {
			if ($ion_counter >= ($CYLINDER{"NIONS"} + $CYLINDER{"COIONS"})) {     # Ran out of ions so place solvent
			    $solv_particles += AddAtm(2,$rx,$ry,$rz,$atm_counter);
			}elsif ($ion_counter >= $CYLINDER{"NIONS"}) { # Ran out of positive ions so place negative
			    $ion_counter += AddAtm(1, $rx, $ry, $rz, $atm_counter);
			} else {                                     # Place a positive ion
			    $ion_counter += AddAtm(0,$rx,$ry,$rz,$atm_counter);
			}
			$rx += $HEXAGON{"LAYER_SPACING"};
			DetermineXtreme($rx);
		    } else {
			last;
			last;
			last;
		    }
		}
	    }
	}
    }

    close OVERLAP;
    $atm_counter = 0;
    print "\n--===RESULTS===---\n";
    if ($PARMS{"SOLUTE"}{"incS"}) {
	print "Placed 1 Solute Molecule\n";
	$atm_counter = 1;
    } 

    print "Placed $ion_counter Ions\n";
    print "Placed $solv_particles Solvent Molecules\n";
    print "Target fluid density is " . $PARMS{"CYLINDER"}{"rhof_targ"} . "\n";
    print "Actual fluid density is " . ($CYLINDER{"FLUID_DENSITY"}  * $solv_particles/$CYLINDER{"NSOLV"}) . "\n";
    
    if ($num_overlap > 0) {
	print "Recorded " . ($num_overlap + 1) . " overlaps with the solute\n";
	print "Check file overlap.dat for more information\n";
    }

    $atm_counter += ($ion_counter + $solv_particles);
}

sub AddAtm(@) {
    my ($is_solvent, $xc, $yc, $zc) = @_;
    my ($return_val, $fluid_nm, $charge, $mass, $radii);

    if ($is_solvent == 2) {
	$fluid_nm = "SOLVENT";
    } elsif ($is_solvent == 1) {
	$fluid_nm = "COIONS";
    } else {
	$fluid_nm = "IONS";
    }


    $charge = $PARMS{$fluid_nm}{"charge"};
    $mass = $PARMS{$fluid_nm}{"mass"};
    $radii = $PARMS{$fluid_nm}{"radii"};

    $return_val = 1;

    $return_val = IsOverlap($radii, $xc, $yc, $zc, $is_solvent)
	if ($PARMS{"SOLUTE"}{"incS"});

    if ($return_val) {
	$atm_counter++;
	push @{ $DATA{$fluid_nm} }, (
			    {
				"XCOORD" => $xc,
				"YCOORD" => $yc,
				"ZCOORD" => $zc,
				"CHARGE" => $charge,
				"MASS"   => $mass,
				"RADII"  => $radii,
			    }
			    );
    }

    return $return_val;
}

sub IsOverlap(@) {
    my ($radii, $rx, $ry, $rz, $isSolv) = @_;
    my ($returnval, $dist);


    $dist = sqrt( ($DATA{"SOLUTE"}[0]{"XCOORD"} - $rx)**2 + ($DATA{"SOLUTE"}[0]{"YCOORD"} - $ry)**2 +
	($DATA{"SOLUTE"}[0]{"ZCOORD"} - $rz)**2 );

    $radii += $PARMS{"SOLUTE"}{"radii"};
    
    if ( ($dist - $radii) > 0) {
	$returnval = 1;
    } else {
	$num_overlap++;
	if ($isSolv == 2) {
	    print OVERLAP "Overlap of solvent ";
	} else {
	    print OVERLAP "Overlap of ion ";
	}
	print OVERLAP "with Solute at $rx $ry $rz\n";
	$returnval = 0;
    }

    return $returnval;
}

sub WriteData() {
    my ($counter, $fmt1, $string_1, $fmt2, $string_2, $res_type, $atm_name);
    my ($type_counter, $temp, $arry_temp, $itype, $out_file_1, $out_file_2);
    my ($type_holder);

    $counter = 1;
    $fmt1 = "%8d%8d%8d%10.4f%15.6f%15.6f%15.6f\n";
    $fmt2 = "%-6s%5d  %4s %3s  %4d    %8.3f%8.3f%8.3f %5.1f %5.1f\n";
    
    $type_holder = ["SOLUTE", "COIONS", "IONS", "SOLVENT"];
    
    $type_holder = OrderRes($type_holder);
    $string_1 = WriteLammpsHeader($type_holder);

    for $type_counter (@{ $type_holder }) {
	if ($type_counter eq "SOLVENT") {
	    $res_type = "WAT";
	    $atm_name = "Wat";
	} elsif ($type_counter eq "COIONS") {
	    $res_type = "CL-";
	    $atm_name = "Cl-";
	} elsif ($type_counter eq "IONS") {
	    $res_type = "NA+";
	    $atm_name = "Na+";
	} else {
	    $res_type = "SVT";
	    $atm_name = "Svt";
	}		
	
	$itype = $PARMS{$type_counter}{"itype"};
	$temp = \%{ $DATA{$type_counter} };
	for $arry_temp ( @{ $temp } ) {
	    $string_1 .= sprintf($fmt1, $counter, 0, $itype, $arry_temp->{"CHARGE"},
				 $arry_temp->{"XCOORD"}, $arry_temp->{"YCOORD"}, $arry_temp->{"ZCOORD"});
	    $string_2 .= sprintf($fmt2, "ATOM", $counter, $atm_name, $res_type,  $itype,
				 $arry_temp->{"XCOORD"}, $arry_temp->{"YCOORD"}, $arry_temp->{"ZCOORD"},
				 $arry_temp->{"CHARGE"}, $arry_temp->{"RADII"});
	    
	    $counter++;
	}
	
    }
    
    $out_file_1 = "cyl_builder_fluid.lmp";
    Createfile($out_file_1, $string_1);
    
    $out_file_2 = "cyl_builder_fluid.pqr";
    Createfile($out_file_2, $string_2);
}

sub Createfile(@) {
    my ($file_nm, $file_text) = @_;

    print "Writing to $file_nm....";
    open OUTFILE, "> $file_nm" or die "Cannot write to $file_nm: $!\n";
    print OUTFILE $file_text;
    close OUTFILE;
    print "SUCESS\n";
}

sub DetermineXtreme(@) {
    my ($xVal) = $_[0];

    $CYLINDER{"DIMENSIONS"}{"XHi"} = $xVal
	if ($xVal > $CYLINDER{"DIMENSIONS"}{"XHi"});

    $CYLINDER{"DIMENSIONS"}{"XLo"} = $xVal
	if ($xVal < $CYLINDER{"DIMENSIONS"}{"XLo"});

}

sub WriteLammpsHeader(@) {
    my ($out_string, $num_types, $counter);
    my ($type_list) = $_[0];

    $out_string = "LAMMPS Description\n\n";
    $out_string .= $atm_counter . " atoms\n";
    $out_string .= "0 bonds\n0 angles\n0 dihedrals\n0 impropers\n\n";

    $num_types = $#{ $type_list }  + 1;
    $out_string .= "$num_types atom types\n";

    $out_string .= "\n" . sprintf("%8.2f", 0);
    $out_string .= sprintf("%8.2f", $CYLINDER{"DIMENSIONS"}{"XHi"}) . " xlo xhi";
    $out_string .= "\n" . sprintf("%8.2f", $CYLINDER{"DIMENSIONS"}{"YLo"});
    $out_string .= sprintf("%8.2f", $CYLINDER{"DIMENSIONS"}{"YHi"}) . " ylo yhi";
    $out_string .= "\n" . sprintf("%8.2f", $CYLINDER{"DIMENSIONS"}{"ZLo"});
    $out_string .= sprintf("%8.2f", $CYLINDER{"DIMENSIONS"}{"ZHi"}) . " zlo zhi\n\n";

    $out_string .= "Masses\n\n";
    for $counter (@{ $type_list }) {
	$out_string .= sprintf("%-5d%8.2f\n",$PARMS{$counter}{"itype"},$PARMS{$counter}{"mass"});
    }

    $out_string .= "\nAtoms\n\n";

    return $out_string;
    

}

sub CalcDebye() {

#   This subroutine will calculate the Debye shielding distance for this system and determine
#   whether this Debye length is greater than the radius of the solute

    my ($eta, $gas_const, $temp, $faraday, $conc, $mol_liter);
    my ($debye, $is_valid);

    $is_valid = 1;

    $mol_liter = (1E27)/(6.023E23); # convert from particles/Angstroms^3 to moles/liter    
    $eta = 78.3; # dielectric constant of water at 298.15 Kelvin
    $eta *= 8.854E-12; # eta = eta * permitivity of vacuum
    $gas_const = 8.3144; # J mol-1 K-1
    $temp = 298.15; # Kelvin   
    $faraday = 9.64846E4; # C mol-1

    if ( ($PARMS{"SOLUTE"}{"charge"}/$PARMS{"IONS"}{"charge"} ) < 0) {
    # Solute is oppositely charged to +'tive Ions
	$faraday *= abs($PARMS{"IONS"}{"charge"}); 
	$conc = abs($CYLINDER{"NIONS"}/$CYLINDER{"VOL"}) * $mol_liter;
	print "Concentration (mol/liter) of positive ions is " . $conc . "\n";
	$conc *= ($PARMS{"IONS"}{"charge"}**2);
    } elsif ($PARMS{"COIONS"}{"incCoIon"}) {
	$faraday *= abs($PARMS{"COIONS"}{"charge"}); # C mol-1
	$conc = abs($CYLINDER{"COIONS"}/$CYLINDER{"VOL"}) * $mol_liter;
	print "Concentration (mol/liter) of coions is " . $conc . "\n";
	$conc *= ($PARMS{"COIONS"}{"charge"}**2);
    } else {
	print "No ions are found of opposite charge to the solute: Infinite Debye Length\n";
	$is_valid = 0;
    }

    if ($is_valid) {
	$conc *= 1E3; # Convert Molarity to moles/m^3
	$debye = sqrt( ($eta * $gas_const * $temp)/(2 * $faraday**2 * $conc));
	print "DeBye Length of system is " . sprintf("%12.8G", $debye) . "\n";
    }
					  
}

sub GetNumCoIons() {

    my ($returnval);

    $returnval = ($PARMS{"COIONS"}{"ion_str"} * $CYLINDER{"VOL"} * 6.023E23)/ 1E27;
    return nint($returnval);

}

sub OrderRes(@) {
    my ($in_res) = $_[0];
    my ($counter, @result);

    for $counter (@{ $in_res  }) {
	if ($DATA{$counter}) {
	    $result[($PARMS{$counter}{"itype"} - 1)] = $counter;
	}
    }

    return \@result;

}

sub CorrectKeyName(@) {

    my ($pkey) = $_[0];

    if ($pkey =~ /^itype/) {
        $pkey = "itype";
    } elsif ($pkey =~ /^q/) {
	$pkey = "charge";
    } elsif ($pkey =~ /^r/ and $pkey !~ /^rh/) {
	$pkey = "radii";
    } elsif ($pkey =~ /^m/) {
	$pkey = "mass";
    }

    return $pkey;
}

sub CreateLammpsIn() {
    my ($type_holder, $counter, $counter1, $out_string, $alpha, $radii, $OUT_DATA, $out_file);
    my ($res1, $res2);

    print "Creating Lammps input files....";

    $type_holder = ["SOLUTE", "COIONS", "IONS", "SOLVENT"];
    
    $type_holder = OrderRes($type_holder);

    for $counter (0 .. $#{ $type_holder }) {
	for $counter1 ($counter .. $#{ $type_holder }) {
	    $radii = $PARMS{$type_holder->[$counter]}{"radii"} + $PARMS{$type_holder->[$counter1]}{"radii"};
	    $res1 = $PARMS{$type_holder->[$counter]}{"itype"};
	    $res2 = $PARMS{$type_holder->[$counter1]}{"itype"};
	    
	    $alpha = $PARMS{"SIMULATION"}{"alpha_coeff"}/$radii;
	    $out_string .= "nonbond coeff";
	    $out_string .= sprintf("%4d%8d%13.4f%6.1f%8.1f%13.7f\n", $res1, $res2,
				   $PARMS{"SIMULATION"}{"well_depth"}, $radii, $PARMS{"SIMULATION"}{"cut_off"}, 
				   $alpha);
	}
    }

    $OUT_DATA = GetTemplateData($sim_data);
    for $counter (@{ $OUT_DATA }) {
	$counter->{"OUTDATA"} = $counter->{"HEADER"} . "\n$out_string\n" . $counter->{"FOOTER"};
	$out_file = $counter->{"FILENAME"};
	open OUTDATA, "> $out_file" or die "Cannot write to $out_file: $!\n";
	print OUTDATA $counter->{"OUTDATA"};
	close OUTDATA;
    }
    print "Done\n";
}

sub GetTemplateData(@) {
    my ($in_file) = $_[0];
    my (@DATA, $which_section, $header, $footer, $which_file);
    my ($fixstyle, $assignfix, $counter);

    $which_section = $which_file = 0;
    open INFILE, $in_file or die "Cannot open $in_file: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^HEADER/) {
	    $which_section  = 1;
	} elsif ($which_section == 1) {
	    if ($_ =~ /^FILEDATA/) {
		$which_section  = 2;
	    } else {
		if ($_ =~ /^\[fix style\]\s+(.+)/) {
		    ($fixstyle, $assignfix) = GenerateField($1);
		}elsif ($_ =~ /^\[(.+)\]\s+(.+)$/) {
		    $header .= sprintf("%-30s", $1) . $2 . "\n\n";
		}
	    }
	} elsif ($which_section == 2) {
	    if ($_ =~ /^FILE (\d+)/) {
		$which_file = $1;
	    } elsif (defined($which_file) and $which_file > 0) {
		$DATA[($which_file - 1)]{"FILENAME"} = "input_lammps_sol_0" . $which_file;
		if ($_ =~  /^\[read data\]\s+(\d+)/) {
		    if ($1) {
			$DATA[($which_file - 1)]{"HEADER"} = sprintf("%-30s", "read data") .
			    "cyl_builder_fluid.lmp\n\n";
			$DATA[($which_file - 1)]{"HEADER"} .= sprintf("%-30s", "create temp") .
			    "gaussian 300.0 688\n";
		    } else {
			$DATA[($which_file - 1)]{"HEADER"} = sprintf("%-30s", "read restart") .
			    "restart_lammps_0" . ($which_file - 1) . "\n";
;
		    }
		    $DATA[($which_file - 1)]{"HEADER"} .= sprintf("%-30s", "dielectric") . "80.0\n";
		} elsif ($_ =~ /^\[fixinfo\]\s+(\d+)/) {
		    if ($1 > 0) {
			$DATA[($which_file - 1)]{"HEADER"} .= $fixstyle;
			if ($1 == 1) {
			    $DATA[($which_file - 1)]{"HEADER"} .= $assignfix;
			}
		    }
		}elsif ($_ =~ /^\[thermo flag\]\s+(\d+)/) {
		    $DATA[($which_file - 1)]{"FOOTER"} .= sprintf("%-30s", "thermo flag") . $1 . "\n\n";
		}elsif ($_ =~ /^\[dump_(\w+)\]\s+(\d+)/) {
		    $DATA[($which_file - 1)]{"FOOTER"} .= sprintf("%-30s", "dump $1") .
			$2 . "   ldump_" . $1 . "_0" . $which_file . "\n\n";
		} elsif ($_ =~ /^\[restart\]\s+(\d+)/) {
		    $DATA[($which_file - 1)]{"FOOTER"} .= sprintf("%-30s", "restart") . 
			$1 . "  2   restart_lammps_0" . $which_file . "   restart_lammps_0" . $which_file . "a\n\n";
		}elsif ($_ =~ /^\[(.+)\]\s+(.+)$/) {
		    $DATA[($which_file - 1)]{"FOOTER"} .= sprintf("%-30s", $1) . $2 . "\n\n";
		}
	    }
	}
    }
    close INFILE;

    for $counter (0 .. $#DATA) {
	$DATA[$counter]{"HEADER"} = $header . $DATA[$counter]{"HEADER"};
    }

    return \@DATA;

}

sub GenerateField(@) {
    my ($fieldstr) = $_[0];
    my ($addfield, $addstyle, $counter, $charge, $restype, $type_holder);
    my (@field_data) = split /\s+/, $fieldstr;

    $addstyle  = "\n" . sprintf("%-30s", "fix style 1   addforce");
    for $counter (@field_data) {
	$addstyle .= sprintf("%10.3G", $counter);
    }
    $addstyle .= "\n";

    $addstyle  .= sprintf("%-30s", "fix style 2  addforce");
    for $counter (@field_data) {
	$addstyle .= sprintf("%-10.3G", (-1 * $counter));
    }
    $addstyle .="\n\n";

    $type_holder = ["SOLUTE", "COIONS", "IONS", "SOLVENT"];

    $type_holder = OrderRes($type_holder);

    for $counter (0 .. $#{ $type_holder }) {
	$charge = $PARMS{$type_holder->[$counter]}{"charge"};
	$restype = $PARMS{$type_holder->[$counter]}{"itype"};
	if($charge > 0) {
	    $addfield .= sprintf("%-30s", "assign fix 1") . "type $restype\n";
	} elsif ($charge < 0) {
	    $addfield .= sprintf("%-30s", "assign fix 2") . "type $restype\n";
	}
    }

    return ($addstyle, $addfield);
}

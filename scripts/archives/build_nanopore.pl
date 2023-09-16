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
sub WriteLammpsHeader();
sub DetermineXtreme(@);

die "usage: $0 parameterfile\n"
    if (! @ARGV);

my ($parmfile) = $ARGV[0];
my (%PARMS, $ANG_cnv, %CYLINDER, %DATA, %HEXAGON, $atm_counter);

Initialize();
ReadParmFile();

# Convert charge to Angstroms

$PARMS{"CYLINDER"}{"qw"} *= $ANG_cnv;
CalcCylDimensions();

$atm_counter = 1;
PlaceSolute()
    if ($PARMS{"SOLUTE"}{"incS"});

PlaceFluid();
WriteData();
 
sub Initialize() {
    $ANG_cnv = 1/(1.6022E-19 * 1.0E20);    # convertion to Angstoms

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

    $s_area_calc = 2 * PI * $PARMS{"CYLINDER"}{"rc"} * $PARMS{"CYLINDER"}{"xc"}; # 2piR*H
    $num_ions = -1 * nint($s_area_calc * $PARMS{"CYLINDER"}{"qw"}/ $PARMS{"IONS"}{"qfion"});
    

    die "Solvent and Ions have same charge!!\n"
	if ($num_ions < 0);

#   Adjust the number of ions to include the fact that there is a charged solute
    $num_ions -= abs($PARMS{"SOLUTE"}{"qSol"})
	if ($PARMS{"SOLUTE"}{"incS"});

    $s_area_real = -1 * $PARMS{"IONS"}{"qfion"} * $num_ions/ $PARMS{"CYLINDER"}{"qw"};
    $CYLINDER{"HEIGHT"} = $s_area_real/(2 * PI * $PARMS{"CYLINDER"}{"rc"});
    $CYLINDER{"RADIUS"} = $PARMS{"CYLINDER"}{"rc"} - ($PARMS{"CYLINDER"}{"sigma"}/2);
    $CYLINDER{"VOL"} = PI * $CYLINDER{"RADIUS"}**2 * $CYLINDER{"HEIGHT"};
    $CYLINDER{"NSOLV"} = nint($CYLINDER{"VOL"} * $PARMS{"CYLINDER"}{"rhof_targ"});
    $CYLINDER{"FLUID_DENSITY"} = $CYLINDER{"NSOLV"}/$CYLINDER{"VOL"};
    $CYLINDER{"SAREA"} =  $s_area_real;
    $CYLINDER{"NIONS"} = $num_ions;

#   CYLINDRICAL DIMENSIONS
    $CYLINDER{"DIMENSIONS"} = (
			      {
				  "XHi" => -99.0,
				  "XLo" => 99.0,
				  "YHi" => $PARMS{"CYLINDER"}{"rc"} * 2,
				  "YLo" => -1 * $PARMS{"CYLINDER"}{"rc"} * 2,
				  "ZHi" => $PARMS{"CYLINDER"}{"rc"} * 2,
				  "ZLo" => -1 * $PARMS{"CYLINDER"}{"rc"} * 2,
			      }
			      );

    print "--==STATS===--\n";
    print "Cylinder radius is " . $PARMS{"CYLINDER"}{"rc"} . "\n";
    print "Cylinder length is " . $CYLINDER{"HEIGHT"} . "\n";
    print "Wall charge density in C/m^2 is " . $PARMS{"CYLINDER"}{"qw"} / $ANG_cnv . "\n";
    print "Wall charge density in Angstroms is " . $PARMS{"CYLINDER"}{"qw"} . "\n";
    print "Number of fluid ions is " . $CYLINDER{"NIONS"} . "\n";
    print "Cylinder Surface Area is " . $CYLINDER{"SAREA"} . "\n";
    print "Cylinder Volume is " . $CYLINDER{"VOL"} . "\n";

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

    print "Target fluid density is " . $PARMS{"CYLINDER"}{"rhof_targ"} . "\n";
    print "Actual fluid density is " . $CYLINDER{"FLUID_DENSITY"} . "\n";
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

    die "Error: Solute and counter ions cannot be oppositely charged!"
	if ($PARMS{"SOLUTE"}{"qSol"}/$PARMS{"IONS"}{"qfion"} < 0);

# Deposit the solute at the top center of the box.

    die "Error: solvent radius exceeds cylinder radius\n"
	if ( ($PARMS{"SOLUTE"}{"rSol"} + $PARMS{"SOLUTE"}{"buffer"}) > ($CYLINDER{"RADIUS"}) );
    
    die "Error: solvent radius exceeds cylinder length\n"
	if ( ($PARMS{"SOLUTE"}{"rSol"} + $PARMS{"SOLUTE"}{"buffer"}) > ($CYLINDER{"HEIGHT"}/2) );
    
    $ry = $rz = 0.0;
    $rx = $PARMS{"SOLUTE"}{"buffer"} + $PARMS{"SOLVENT"}{"rf"};
    DetermineXtreme($rx);
    
    print "Placed Solute at $rx $ry $rz\n";
    $DATA{"SOLUTE"} = (
			 {
			     "XCOORD" => $rx,
			     "YCOORD" => $ry,
			     "ZCOORD" => $rz,
			     "CHARGE" => $PARMS{"SOLUTE"}{"qSol"},
			     "MASS"   => $PARMS{"SOLUTE"}{"mSol"},
			     "RADII"  => $PARMS{"SOLUTE"}{"rSol"},
			 }
		       );
    $atm_counter++;
}

sub PlaceFluid() {
    my ($ion_counter, $hex_counter, $counter, $layer_counter, $rank_counter, $solv_particles);
    my ($rx, $ry, $rz) = (0,0,0);

#   Place the first set of ions/solvent as a single line along the middle of the tube
#   After will place rings of ions/fluids around center line in hexagonal packing arrangement

    $ion_counter = $solv_particles = 0;

    for $layer_counter (1 .. $HEXAGON{"NUM_LAYERS"}) {
	if ($ion_counter < $CYLINDER{"NSOLV"}) {
	    if ($ion_counter > $CYLINDER{"NIONS"}) {     # Ran out of ions so place solvent
		$solv_particles += AddAtm(1,$rx,$ry,$rz,$atm_counter);
	    } else {                                     # Place an ion
		$ion_counter += AddAtm(0,$rx,$ry,$rz,$atm_counter);
	    }
	    $rx += $HEXAGON{"LAYER_SPACING"};
	    DetermineXtreme($rx);
	    $atm_counter++;
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
		    if ($atm_counter <= $CYLINDER{"NSOLV"}) {
			if ($ion_counter > $CYLINDER{"NIONS"}) {     # Ran out of ions so place solvent
			    $solv_particles += AddAtm(1,$rx,$ry,$rz,$atm_counter);
			} else {                                     # Place an ion
			    $ion_counter += AddAtm(0,$rx,$ry,$rz,$atm_counter);
			}
			$rx += $HEXAGON{"LAYER_SPACING"};
			DetermineXtreme($rx);
			$atm_counter++;
		    }
		}
	    }
	}
    }

    $atm_counter = 0;
    print "\n--===RESULTS===---\n";
    if ($PARMS{"SOLUTE"}{"incS"}) {
	print "Placed 1 Solute Molecule\n";
	$atm_counter = 1;
    } 

    print "Placed $ion_counter Ions\n";
    print "Placed $solv_particles Solvent Molecules\n";
    
    $atm_counter += ($ion_counter + $solv_particles);
}

sub AddAtm(@) {
    my ($is_solvent, $xc, $yc, $zc) = @_;
    my ($return_val, $fluid_nm, $charge, $mass, $radii);

    if ($is_solvent) {
	$fluid_nm = "SOLVENT";
	$charge = "qf";
	$mass = "mf";
	$radii = "rf";
    } else {
	$fluid_nm = "IONS";
	$charge = "qfion";
	$mass = "mfion";
	$radii = "rfion";
    }

    $return_val = 1;

    $return_val = IsOverlap($is_solvent, $xc, $yc, $zc)
	if ($PARMS{"SOLUTE"}{"incS"});

    if ($return_val) {
	push @{ $DATA{$fluid_nm} }, (
			    {
				"XCOORD" => $xc,
				"YCOORD" => $yc,
				"ZCOORD" => $zc,
				"CHARGE" => $PARMS{$fluid_nm}{$charge},
				"MASS"   => $PARMS{$fluid_nm}{$mass},
				"RADII"  => $PARMS{$fluid_nm}{$radii},
			    }
			    );
    }

    return $return_val;
}

sub IsOverlap(@) {
    my ($isSolv, $rx, $ry, $rz) = @_;
    my ($radii, $returnval, $dist);

    if ($isSolv) {
	$radii = $PARMS{"SOLVENT"}{"rf"};
    }else {
	$radii = $PARMS{"IONS"}{"rfion"};
    }

    $dist = sqrt( ($DATA{"SOLUTE"}{"XCOORD"} - $rx)**2 + ($DATA{"SOLUTE"}{"YCOORD"} - $ry)**2 +
	($DATA{"SOLUTE"}{"ZCOORD"} - $rz)**2 );

    $radii += $PARMS{"SOLUTE"}{"rSol"};
    
    if ( ($dist - $radii) > 0) {
	$returnval = 1;
    } else {
	if ($isSolv) {
	    print "Overlap of solvent ";
	} else {
	    print "Overlap of ion ";
	}
	print "with Solvent at $rx $ry $rz\n";
	$returnval = 0;
    }

    return $returnval;
}

sub WriteData() {
    my ($counter, $fmt1, $string_1, $fmt2, $string_2, $res_type, $atm_name);
    my ($type_counter, $temp, $arry_temp, $itype, $out_file_1, $out_file_2);

    $counter = 1;
    $fmt1 = "%8d%8d%8d%10.4f%15.6f%15.6f%15.6f\n";
    $fmt2 = "%-6s%5d  %4s %3s  %4d    %8.3f%8.3f%8.3f %5.1f %5.1f\n";

    $string_1 = WriteLammpsHeader();
    if ($PARMS{"SOLUTE"}{"incS"}) {
	$string_1 .= sprintf($fmt1, $counter, 0, $PARMS{"SOLUTE"}{"itypeSol"}, $DATA{"SOLUTE"}{"CHARGE"},
			    $DATA{"SOLUTE"}{"XCOORD"}, $DATA{"SOLUTE"}{"YCOORD"}, $DATA{"SOLUTE"}{"ZCOORD"});
	$string_2 = sprintf($fmt2, "ATOM", $counter, "SO+", "SO+", $PARMS{"SOLUTE"}{"itypeSol"},
			    $DATA{"SOLUTE"}{"XCOORD"}, $DATA{"SOLUTE"}{"YCOORD"}, $DATA{"SOLUTE"}{"ZCOORD"},
			    $DATA{"SOLUTE"}{"CHARGE"}, $DATA{"SOLUTE"}{"RADII"});
 	$counter++;
   }

    for $type_counter ("IONS", "SOLVENT") {
	if ($type_counter ne "SOLUTE") {
	    if ($type_counter eq "SOLVENT") {
		$itype = $PARMS{"SOLVENT"}{"itypef"}; 
		$res_type = "WAT";
		$atm_name = "WAT";
	    } else {
		$itype = $PARMS{"IONS"}{"itypefion"};
		$res_type = "Na+";
		$atm_name = "Na+";
	    }

	    $temp = \%{ $DATA{$type_counter} };
	    for $arry_temp ( @{ $temp } ) {
		$string_1 .= sprintf($fmt1, $counter, 0, $itype, $arry_temp->{"CHARGE"},
				     $arry_temp->{"XCOORD"}, $arry_temp->{"YCOORD"}, $arry_temp->{"ZCOORD"});
		$string_2 .= sprintf($fmt2, "ATOM", $counter, $atm_name, $res_type, $itype,
				    $arry_temp->{"XCOORD"}, $arry_temp->{"YCOORD"}, $arry_temp->{"ZCOORD"},
				    $arry_temp->{"CHARGE"}, $arry_temp->{"RADII"});
		
		$counter++;
	    }
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

sub WriteLammpsHeader() {
    my ($out_string);

    $out_string = "LAMMPS Description\n\n";
    $out_string .= $atm_counter . " atoms\n";
    $out_string .= "0 bonds\n0 angles\n0 dihedrals\n0 impropers\n\n";

    if ($PARMS{"SOLUTE"}{"incS"}) {
	$out_string .= "3 atom types\n";
    } else {
	$out_string .= "2 atom types\n";
    }

    $out_string .= "\n" . sprintf("%8.2f", 0);
    $out_string .= sprintf("%8.2f", $CYLINDER{"DIMENSIONS"}{"XHi"}) . " xlo xhi";
    $out_string .= "\n" . sprintf("%8.2f", $CYLINDER{"DIMENSIONS"}{"YLo"});
    $out_string .= sprintf("%8.2f", $CYLINDER{"DIMENSIONS"}{"YHi"}) . " ylo yhi";
    $out_string .= "\n" . sprintf("%8.2f", $CYLINDER{"DIMENSIONS"}{"ZLo"});
    $out_string .= sprintf("%8.2f", $CYLINDER{"DIMENSIONS"}{"ZHi"}) . " zlo zhi\n\n";

    $out_string .= "Masses\n\n";
    $out_string .= sprintf("%-5d%8.2f\n",$PARMS{"IONS"}{"itypefion"},$PARMS{"IONS"}{"mfion"});
    $out_string .= sprintf("%-5d%8.2f\n",$PARMS{"SOLVENT"}{"itypef"},$PARMS{"SOLVENT"}{"mf"});
    $out_string .= sprintf("%-5d%8.2f\n\n",$PARMS{"SOLUTE"}{"itypeSol"},$PARMS{"SOLUTE"}{"mSol"})
	if ($PARMS{"SOLUTE"}{"incS"});

    $out_string .= "Atoms\n\n";

    return $out_string;
    

}

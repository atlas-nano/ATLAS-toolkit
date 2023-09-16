#!/usr/bin/perl -w

use strict;
use constant PI => 3.14159265358979;
# This program will perform apbs for the specified molecule and report the
# electrostatics


sub VerifyFile(@);
sub GetDimensions(@);
sub CreateInput(@);
sub RunApbs(@);
sub do_round(@);
sub SoluteDiElectric();
sub GetLargest(@);
sub CalcGP(@);
sub GetParms(@);
sub ExtractParms(@);
sub GenerateInputs();
sub Initialize();
sub WriteVals(@);

my ($pqr_file, $parm_file, $cubic) = @ARGV;
my (%ATOMINFO, @DATA, $counter);
my ($Grid, $Grid_Pts, $in_file_name);
my ($was_sucess, $duration, $elec, $index);
my ($fine_grid, $mol_volume, $input_template, $calc_epsilon);
my ($curr_variant, $new_variant, $my_space);

die "usage: $0 pqr_file parm_file [cubic]\n"
    if (! @ARGV or $#ARGV < 1);

$cubic = 0
    if (! $cubic);

my (%PARMS); 

GetParms($parm_file);

$index = 1;


for $counter (@DATA) {
    print "\nCalculation $index of ", ($#DATA +1), "\n";
    ($counter->{"Grid"}, $counter->{"Grid_Pts"}, $my_space) = GetDimensions($pqr_file, $counter);
    $in_file_name = CreateInput($counter, $counter->{"variant"});
    ($was_sucess, $duration, $elec) = RunApbs($counter, $in_file_name);
    
    if ($index ==  1) {
	WriteVals($counter,0,0,1, 0);
	$curr_variant = $counter->{"variant"};
    } 
    $new_variant = $counter->{"variant"};
    if (! $new_variant or ! $curr_variant or $new_variant ne $curr_variant or $index == 1) {
	$curr_variant = $new_variant;
	$new_variant = 1;
    } else {
	$new_variant = 0;
    }

    WriteVals($counter, $duration, $elec, 0, $new_variant);
 	
    if (! $was_sucess) {
	system "rm -fr core*";
    }
   
    print "Done Calculation $index\n";
    $index++;
}

print "Check the comparison.dat file for a comparison\n";

sub VerifyFile(@) {
    die "Cannot locate $_[0]: $!\n"
	if (! -e $_[0]);
    die "$_[0] is not a text file!\n"
	if (! -T $_[0]);
    die "$_[0] is not readable!\n"
	if (! -r $_[0]);
    
}

sub GetDimensions(@) {
    my ($in_file, $curr_input) = @_;
    my (@BBOX, @MDIMS, $flen_factor, $clen_factor, @F_MESH, @C_MESH, @my_spacing);
    my ($GRID, $GRID_PTS, $counter, $spacing, @TN, $holder, $largest_dim, $nlev, $val1);

    $flen_factor = $curr_input->{"grid_factor"};  # The number of amgstrongs to add to each box coordinate for the fine mesh
    $clen_factor = 1.0; # Factor to multiply the box dimensions by for the coarse length
    $spacing = $curr_input->{"grid_spacing"};     # Grid spacing

    @BBOX = (-999, 999, -999, 999, -999, 999);
    $counter = 0;
    open PQRFILE, $in_file or die "Cannot open $in_file: $!\n";
    while (<PQRFILE>) {
	chomp;
	if ($_ =~ /^ATOM\s+\d+\s+(\w+)\S*\s+\w+\s+(\d+)\s+(.+)/) {
	    $counter++;
	    @TN = split /\s+/, $3;	    
	    if (! $ATOMINFO{$1}) {
		$ATOMINFO{$1} = (
				     {
					 "ATOMS"  => $1 . $2,
					 "NUM_ATOMS" => 1,
					 "RADII"     => $TN[4],
				     }
				     );
	    } else {
		$ATOMINFO{$1}{"NUM_ATOMS"}++;
		$ATOMINFO{$1}{"ATOMS"} .= " $1" . $2;
	    }

	    $BBOX[0] = $TN[0]
		if ($TN[0] > $BBOX[0]);
	    $BBOX[1] = $TN[0]
		if ($TN[0] < $BBOX[1]);
            $BBOX[2] = $TN[1]
                if ($TN[1] > $BBOX[2]);
            $BBOX[3] = $TN[1]
                if ($TN[1] < $BBOX[3]);
            $BBOX[4] = $TN[2]
                if ($TN[2] > $BBOX[4]);
            $BBOX[5] = $TN[2]
                if ($TN[2] < $BBOX[5]);
	}
    }
    close PQRFILE;

    print "Found $counter atoms\n";
    @MDIMS = (($BBOX[0] - $BBOX[1]), ($BBOX[2] - $BBOX[3]), ($BBOX[4] - $BBOX[5]));
    
    $largest_dim = GetLargest(\@MDIMS);
    for $counter (0 .. 2) {
	$C_MESH[$counter] = $clen_factor + $MDIMS[$counter];
	if (! $cubic) {
	    $F_MESH[$counter] = $flen_factor * $MDIMS[$counter];
        } else {
	    $F_MESH[$counter] = $flen_factor * $largest_dim;
        }
 	
        if ($fine_grid) {   
            $GRID .= sprintf("%-11.2f", $F_MESH[$counter]);
	    ($val1, $my_spacing[$counter]) = CalcGP($F_MESH[$counter], $spacing, 4);
	} else {
	    $GRID .= sprintf("-%11.2f", $C_MESH[$counter]);
	    ($val1, $my_spacing[$counter]) = CalcGP($C_MESH[$counter], $spacing, 4);
	}
	$GRID_PTS .= $val1;
    }

    return ($GRID, $GRID_PTS, \@my_spacing);
}

sub CreateInput(@) {

    my ($curr_input, $which_variant) = @_;
    my ($in_string, $epsilon, $water_perm, $fle_nm);

    $fle_nm = $curr_input->{$which_variant} . "_" . $which_variant;
    if ($calc_epsilon) {
	$epsilon = SoluteDiElectric();
    } else {
	$epsilon = $curr_input->{"solute_episolon"};
    }

    $water_perm = $curr_input->{"solvent_epsilon"};

#    $input_template = "./apbs.in";
    open INFILE, "$input_template" or die "Cannot open $input_template: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /mol pqr/) {
	    $in_string .= "    mol pqr $pqr_file\n";
	}elsif ($_ =~ /glen/) {
	    $in_string .= "    glen " . $curr_input->{"Grid"} . "\n";
	} elsif ($_ =~ /(\s+dime)/) {
	    $in_string .= "$1 " . $curr_input->{"Grid_Pts"} . "\n";
	} elsif ($_ =~ /pdie/) {
	    $in_string .= "    pdie " . sprintf("%5.2f", $epsilon) . "\n";
	} elsif ($_ =~ /_sdie/) {
	    $in_string .= "    sdie " . sprintf("%5.2f", $water_perm) . "\n";
	} elsif ($_ =~ /(\s+nlev)/) {
	    $in_string .= "$1 4\n";
	} else {
	    $in_string .= $_ . "\n";
	}
    }
    close INFILE;

    open OUTFILE, "> ./$fle_nm" . ".in" or die "Cannot create apbs input file $fle_nm " . ".in : $!\n";
    print OUTFILE $in_string;
    close OUTFILE;

    return $fle_nm;
}

sub RunApbs(@) {
    my ($curr_in, $in_file) = @_;
    my ($out_text, $is_valid, $start, $end, $duration, $elec, $counter, $out_file);

    $out_file = $in_file . ".out";
    $in_file .= ".in";

    $is_valid = 0;

    print "\n--==INPUT PARAMETERS==--\n\n";
    for $counter (keys %{ $curr_in }) {
	printf "%-15s%15s\n", $counter, $curr_in->{$counter};
    }
    
    print "Actual Grid Spacing: ";
    for $counter (@{ $my_space } ) {
	printf "%12.3f x", $counter;
    }
    print "\n";

    delete $curr_in->{"Grid_Pts"};
    delete $curr_in->{"Grid"};

    $elec = 0;
    print "\n--==RESULTS==--\n\n";
    $start = time();
    print "APBS started: ", scalar localtime($start), "\n";

    open APBS, "/ul/tpascal/apbs-0.3.1/bin/i686-intel-linux/apbs ./$in_file |" or 
	die "Error while executing apbs: $!\n";
    while (<APBS>) {
	$out_text .= $_;
	chomp;
	if ($_ =~ /Global net energy\s+.\s+(\-?\d+\.\d+.+)\s+kJ/) {
	    $is_valid = 1;
	    $elec = $1/4.184;
	    print "Electrostatic energy: " . $elec . " kcal/mol\n";
	}
	    
    }
    close APBS;
    $end = time();
    print "APBS ended: ", scalar localtime($end), "\n";

    open OUTFILE, "> $out_file" or die "Cannot create $out_file: $!\n";
    print OUTFILE $out_text;
    close OUTFILE;

    if (! $is_valid) {
	print "Error while running apbs. Check apbs.out file for more details\n";
    }else {
	print "Check apbs.out file for details of calculation\n";
    }
    
    $duration = $end - $start;
    if ($duration > 1) {
	print "Apbs Run time: ",$duration, " sec\n";
    }

    return ($is_valid, $duration, $elec);
}

sub do_round(@) {
    my ($c1, $c2, $c3) = @_;
    my ($return_val, $grid_factor);

    $grid_factor = 1;
    if ($c1 > $c2) {
	if ($c1 > $c3) {
	   $return_val = $c1;
        }else {
	   $return_val = $c3;
	} 
    }elsif ($c2 > $c3) {
	   $return_val = $c2;
    } else {
	   $return_val = $c3;
    }

   $return_val = sprintf("%-11.2f", $return_val);
   return "$return_val $return_val $return_val";
}

sub SoluteDiElectric() {
# Nag B R Appl. Phys. Lett. 65 1938, 1994
# Xue, Betzler and Hesse J. Phys.: Condens. Matter 12, 3113, 2000
# Marquez R and Rincon C Phys. Status Solidi b 191, 115, 1995
# mol_eps = 1+(s/(eps0-gamma*s))
#
#       s = SUM N(i)*alpha(i)
# where
#
#       N(i) = number of atoms of species i per unit volume (1/A^3)
#
# An estimate of the molar volume Vm may be obtained from the empirical expression
#
#        Vm = 1.1001*Mv
#
# where Vm is the molar volume Vm = cc/mol and Mv is the molecular volume in A^3
# Thus, the effective volume per molecule is  V = No Vm, where No is avogadros number
# In A^3 this expression becomes  V = 6.022 x 10^23 * 1.1001 * Mv * 1.0 10^-24
#                                 V = 0.6022 x 1.1001 x Mv = 0.66248 Mv  A^3
#
# In calculating the Molecular Volume:
# V = sum (4/3 pi * radius(i)^3)

    my ($sum_alpha) = 0; # alpha is the molecular dielectric
    my ($lorentz, $epsO) = (1, 1); #lorentz factor = 1 since v = 0, epso is permitivity of vacuum
    my ($epsilon) = 0;
    my ($curr_atm, %ALPHA_INFO, $dat_file);


    $dat_file = "/ul/tpascal/scripts/apbs/ubm.dat";
    open INFILE, $dat_file or die "Cannot locate $dat_file: $!\n";
    while (<INFILE>) {
	if ($_ =~ /^([a-zA-Z]+)_?\s+\d+\.\d+\s+(\d+\.\d+)/) {
	    $ALPHA_INFO{$1} = $2;
	}
    }
    close INFILE;

    if (! $mol_volume) {
	$epsilon = 1;
	return $epsilon;
    }

    for $curr_atm (keys %ATOMINFO) {
	if ($ALPHA_INFO{$curr_atm}) {
	    $sum_alpha += $ALPHA_INFO{$curr_atm} * $ATOMINFO{$curr_atm}{"NUM_ATOMS"};
	} else {
	    print "BAD ATOM: $curr_atm";
	}
    }

    $sum_alpha = $sum_alpha / (0.66248 * $mol_volume);
    $epsilon = 1 + ($sum_alpha/($epsO - $lorentz * $sum_alpha));
    return $epsilon;

}

sub GetLargest(@) {

    my ($bdim) = $_[0];
    my ($largest_dim, $counter);

    $largest_dim = 0;
    for $counter (@{ $bdim }) {
	if ($counter > $largest_dim) {
	    $largest_dim = $counter;
	}
    }

    return $largest_dim;

}

sub CalcGP(@) {

    my ($largest_dim, $spacing, $nlev) = @_;
    my ($needed_pts, $counter, $real_spacing);

    $needed_pts = int(($largest_dim)/$spacing) + 1;

    $counter = 1;
    
    if ((($needed_pts -1) % 32) > 0) {
	while ((32 * $counter) < $needed_pts) {
	    $counter++;
	}
	
	$counter = (32 * $counter) + 1;
    } else {
	$counter = $needed_pts;
    }
    
    $real_spacing = $largest_dim / $counter;
    print "Pts needed for spacing: $needed_pts, Pts used by APBS: $counter\n";
    print "Spacing Gotten: " . $real_spacing . "\n";
    return (sprintf("%-12d", $counter), $real_spacing);
    

}

sub GetParms(@) {
    my ($in_file) = $_[0];
    my ($result);
    die "Cannot locate parameter file $in_file: $!\n"
	if (! -e $in_file);

    open INFILE, $in_file or die "Error opening $in_file: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^(\w+)\s+(.+)$/) {
	    $result = ExtractParms($2);
	    $PARMS{$1} = $result
		if ($result);
	} elsif ($_ =~ /^\s+(\w+)\s+(.+)$/) {
	    $PARMS{$1} = $2;
	}
	    
    }
    close INFILE;

    Initialize();
    GenerateInputs();

}

sub ExtractParms(@) {
    my ($invals) = $_[0];
    my (@headers) = ("average", "hi", "low", "increment");
    my ($counter, %rec, $key);

    $counter = 0;

    while ($invals =~ /(\d+\.?\d*)\s*/g) {
	$key = $headers[$counter];
	$rec{$key} = $1;
	$counter++;
    }

    if (! $counter) {
	return 0;
    } else {
	return \%rec;
    }
}

sub GenerateInputs() {
    my ($key1, $key2, $counter);
    my ($hi, $lo, $increment, $avg, $rec);

    for $key1 (keys %PARMS) {
	$hi = $PARMS{$key1}{"hi"};
	$lo = $PARMS{$key1}{"low"};
	$avg = $PARMS{$key1}{"average"};
	$increment = $PARMS{$key1}{"increment"};

	
	next
	    if ($increment == "9999");
	while ($lo <= $hi) {
	    $rec = {};
	    $rec->{$key1} = $lo;
	    $rec->{"variant"} = $key1;
	    for $key2 (keys %PARMS) {
		if ($key1 ne $key2) {
		    $rec->{$key2} = $PARMS{$key2}{"average"};
		}
	    }
	    push @DATA, $rec;
	    $lo += $increment;
	}
    }
}

sub Initialize() {

    if ($PARMS{"calc_epislon"} ne "") {
	$calc_epsilon = $PARMS{"calc_epislon"};
	delete $PARMS{"calc_epislon"};
    } else {
	die "Error: Parameter calc_epislon is missing from the parameter file\n";
    }

    if ($PARMS{"mol_volume"} and $PARMS{"mol_volume"} ne "na") {
	$mol_volume = $PARMS{"mol_volume"};
    } else {
	$mol_volume = "";
    }
    delete $PARMS{"mol_volume"};
    
    
    if ($PARMS{"fine_grid"} ne "") {
	$fine_grid = $PARMS{"fine_grid"};
	delete $PARMS{"fine_grid"};
    } else {
	die "Error: Parameter fine_grid is missing from the parameter file\n";
    }
    
    
    if ($PARMS{"input_template"} and $PARMS{"input_template"} ne "na") {
	$input_template = $PARMS{"input_template"};
    } else {
	$input_template = "/ul/tpascal/scripts/apbs/apbs.in";
    }
    delete $PARMS{"input_template"};
    
    die "use_fine_grid expected 1 or 0"
	if ($fine_grid !~ /^[1|0]$/);
    
    
    VerifyFile($pqr_file);
    
    
    VerifyFile($input_template);

}

sub WriteVals(@) {
    my ($curr_file, $length, $elec, $header, $variant) = @_;
    my ($out_string);

    

    if ($header) {
	open OUTDATA, "> comparison.dat" or die "Cannot write to comparison.dat: $!\n";
    } else {
	open OUTDATA, ">> comparison.dat" or die "Cannot write to comparison.dat: $!\n";
    }
    
    if ($variant) {
	printf OUTDATA "%60s%20s%40s\n", " ", $curr_file->{"variant"}, " ";
    }
    
    delete $curr_file->{"variant"};
    
    for $counter (keys %{ $curr_file }) {
	if ($header) {
	    $out_string .= sprintf("%20s", $counter);
	} else {
	    $out_string .= sprintf("%20.6f", $curr_file->{$counter});
	}
    }
    
    if ($header) {
	$out_string .= sprintf("%20s%20s%12s%12s%12s\n", "Elec", "Duration", "X Spacing", "Y Spacing", "Z Spacing");
    } else {
	if ($elec == 0) {
	    $out_string .= sprintf("%20s%20s", "Mem Error", "Mem Error");

	} else {
            $out_string .= sprintf("%20.6f%20.6G", $elec, $length);
	}
	for $counter (@{ $my_space }) {
	    $out_string .= sprintf("%12.6f", $counter);
	}
	$out_string .= "\n";
    }
    
    print OUTDATA $out_string;
    close OUTDATA;
}

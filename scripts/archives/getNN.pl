#!/usr/bin/perl -w
#    This script will march along the specified double helices and collect the
#    NN triples according to 
#       E(BP) = Internal(BPx) + 0.5*NonValence(BPx:Bpx-1 + BPx:BPx+1)
#    It assumes that the energy of a particular BP is only due to it's nearest
#    neighbours. It then collects a system of linear equations and then uses 
#    Mathematica to solve the corresponding linear system for the variables
BEGIN {
    push @INC, "/ul/tpascal/.libs/Math-ematica-1.108";
    push (@INC, "/ul/tpascal/scripts/");
}

use strict;
use File::Basename qw(basename);
use Packages::General;
use Packages::FileFormats;
use Packages::HelixLayout;
use Packages::GetParms;
use Math::ematica qw(:PACKET :TYPE :FUNC);

sub Numerically { ($a<=>$b); }

my ($anal_order, $save_name);
my ($FILES, $bgfFile, $datFile); $strand_length, $out_name, $energyfile, $error_norm);
my ($PDB_Info, $Helix_Info, $Energy_Info, $EQNS, $Vars, $result, $f_c);

$FILES = &init;

$f_c = 0;
for $i (@{ $file_list }) {
    print "Creating Linear System from BGF file $i->{BGFFILE}...";
    ($ATOMS, $BONDS) = GetBGFFileInfo($i->{BGFFILE}, 0);
	$strand_length = 0;
	print "Processing $pdbfile....";
	$out_name .= ".anal";
	$energyfile = $pdbfile;
	$energyfile =~ s/\.pdb/\.dat/g;
	($strand_length, $PDB_Info) = ParsePDBFile($pdbfile);
	next
	    if (! $strand_length);
	($Helix_Info) = GetHelixLayout($strand_length);
	($Energy_Info) = GetMPSimEnergy($energyfile, ($strand_length * 2));
	CreateLinearSystem(\%{ $EQNS }, $Helix_Info, $Energy_Info, $anal_order, $f_c, \%{ $Vars });
	print "Done\n";
	$f_c++;
    }

    die "ERROR: No Valid data found\n"
	if (! defined($EQNS));
    
    ($result, $error_norm) = SolveLinearSystem(\%{ $EQNS }, \%{ $Vars });
    if (! $save_name) {
	$save_name = "results.txt";
    }
    WriteOutput($save_name, $result, $error_norm, $Vars);
}

sub init {
    my (@bgffiles, $curr_file, $dat_file, @valid_files, %OPTS, $loc, $work_dir);

    $usage = "usage: $0 -d directory|bgf file -l level -s [save name]\n" .
	"Options:\n\t-d directory|bgf file: location of bgf file or directory of files\n" . 
	"\t-l level: 1-5. Level 1: No nearest neighbor correction. Level 2: Nearest neighbor correction\n" .
	"\t\tLevel 3: Next Nearest neighbor correction etc.\n" .
	"\t-s [save name]: (Optional) file name to store data\n";

    getopt('dls',\%OPTS);
    ($work_dir, $anal_order, $save_name) = ($OPTS{d},$OPTS{l},$OPTS{s});
    for ($work_dir, $anal_order) {
	die "$usage" if (! defined($_));
    }

    print "Initializing...";
    $anal_order = Trim($anal_order);
    die "ERROR: Expected integer for level, got: $anal_order\n"
	if (! IsInteger($anal_order));
    die "ERROR: level has to be between 1 and 5\n"
	if ($anal_order < 1 or $anal_order > 5);

    while ($work_dir =~ /(\S+)/g) {
	$loc = $1;
	if (-e $loc) {
	    if (! -d $loc) {
		$bgffiles[0] = $loc;
	    } else {
		opendir (WORKDIR, $loc) or die "Cannot access specified directory: $!\n";
		@bgffiles = grep { /\w+\.bgf/ } map {"$loc/$_"} readdir WORKDIR;
		closedir WORKDIR;
	    }
	}
	
    }

    for $curr_file (@bgffiles) {
	$dat_file = $curr_file;
	$dat_file =~ s/\.bgf$/\.dat/;
	if (-e $dat_file) {
	    push @valid_files, (
				{
				    "BGFFILE" => $curr_file,
				    "DATFILE" => $dat_file,
				}
				);
	}
    }

    die "ERROR: No valid files found in $work_dir Aborting\n" if ($#valid_files < 0);
    
    print "Done\n";
    return \@valid_files;
}

sub ParsePDBFile(@) {
    my ($pdb_file) = $_[0];
    my ($Atom_Info) = GetPDBFileInfo($pdb_file, 0, " ", " ");
    my ($counter, $curr_res_name, %Res_Info, $curr_res_id, $strand_indicator, $strand_length);

    $curr_res_id = $strand_length = 0;
    $strand_indicator = "";
    for $counter (sort Numerically keys %{ $Atom_Info }) {
	if ($Atom_Info->{$counter}{"RES_ID"} > $curr_res_id) {
	    $curr_res_id = $Atom_Info->{$counter}{"RES_ID"};
	    $curr_res_name = $Atom_Info->{$counter}{"RES_NAME"};
	    if ($curr_res_name =~ /(\d+)/) {
		if ($strand_indicator eq "") {
		    $strand_indicator = $1;
		} elsif ($strand_indicator ne $1 and ! $strand_length) {
		    $strand_length = $curr_res_id;
		}
	    }
	    $Res_Info{$curr_res_id} = $curr_res_name;
	}
    }

    if (! $strand_length) {
	print "ERROR processing $pdb_file: CANNOT find strand_length\n";
    }
    return ($strand_length, \%Res_Info);

}

sub GetHelixLayout(@) {
    my ($length_strand) = @_;
    my (%Helix, $i, $j);

    for $j (1 ... 2) {
	for $i (1 ... $length_strand) {
	    $Helix{1}{$i} =  $i;
	    $Helix{2}{$i} = ($length_strand * 2) - $i + 1;
	}
    }
    return \%Helix;
}

sub CreateLinearSystem(@) {
    my ($Matrix, $Helix_Info, $Energy_Info, $level, $fileC, $VARS) = @_;
    my ($helix_counter, $unit_counter, $strand_length, @unit_holder, $weight_index);
    my ($index, $curr_unit, $matrix_ele, @Weights, $isEnd, $chain_counter, $crossed_middle);
    my ($sequence_energy, $sequence_name, $norm_factor, $unit_1, $unit_2, $nextunit, $eng_total);

    $level--;
    $eng_total = 0.0;
    ($norm_factor, @Weights) = GetWeights($level);

    @unit_holder = keys %{ $Helix_Info->{1} };
    $strand_length = ($#unit_holder + 1);
    @unit_holder = ();
    for $unit_counter (1 .. $strand_length) {
	$sequence_energy = $isEnd = 0;
	$sequence_name = "";
	if ($unit_counter == 1 or $unit_counter == $strand_length) {
	    $isEnd = 1;
	} else {
	    $isEnd = 0;
	}

	$curr_unit = $Helix_Info->{1}{$unit_counter};
	@unit_holder = ();
	for $index ((-$level + $unit_counter).. ($level + $unit_counter)) {
	    if (defined($Helix_Info->{1}{$index})) {
		push @unit_holder, $Helix_Info->{1}{$index};
	    }
	}
	
#	Internal Energy
	if ($isEnd) {
	    $matrix_ele = $PDB_Info->{$curr_unit};
	    if ($matrix_ele =~ /^D?(\w)(\d)/) {
		$matrix_ele = $2 ."'" . $1;
		$Matrix->{$fileC}{$unit_counter}{"UNITS"}{$matrix_ele} = 1;;
		$VARS->{$matrix_ele} += 1;
	    } else {
		die "UNKNOWN ELEMENT: $matrix_ele\n";
	    } 
	}

	for $chain_counter (1 .. 2) {
	    $curr_unit = $Helix_Info->{$chain_counter}{$unit_counter};
	    $matrix_ele = $PDB_Info->{$curr_unit};
	    if ($matrix_ele =~ /^D?(\w)/) {
		$matrix_ele = $1;
		if ($chain_counter == 1) {
		    $Matrix->{$fileC}{$unit_counter}{"UNITS"}{$matrix_ele} += 1;;
		    $VARS->{$matrix_ele} += 1;
		}
		$sequence_energy += $norm_factor * $Energy_Info->{$curr_unit}->{"Bond"};
	    } else {
		die "UNKNOWN ELEMENT: $matrix_ele\n";
	    }

	}

#	print "ME: $matrix_ele\n";
#	NonBond Energy Contribution
	$crossed_middle = 0;
	for $index (0 .. $#unit_holder) {
	    for $nextunit ($index  .. $#unit_holder) {
		$sequence_name = "";
		if (($index != $nextunit) or ($level == 0)) {
		    if ($unit_holder[$index] >= $unit_counter) {
			$crossed_middle = 1;
		    }
		    last
			if (! $crossed_middle and $unit_holder[$nextunit] > $unit_counter);
		    $weight_index = $#Weights/2 - ($nextunit - $index);
		    for $chain_counter (1 .. 2) {
			$unit_1 = $Helix_Info->{$chain_counter}{$unit_holder[$index]};
			$unit_2 = $Helix_Info->{$chain_counter}{$unit_holder[$nextunit]};
			
			$sequence_energy += GetSequenceEnergy($Weights[$weight_index], \%{ $Energy_Info->{$unit_1} }, 
							      \%{ $Energy_Info->{$unit_2} }, $unit_counter, $level, $isEnd);
			if ($level > 0) {
			    $sequence_name .= "$unit_1 $unit_2 ";
			} else {
			    $sequence_name .= "$unit_1 ";
			}
		    }
		    chop $sequence_name;
		    if ($level > 0) {
			$matrix_ele = DetermineSequence($sequence_name);
#			$Matrix->{$fileC}{$unit_counter}{"UNITS"}{$matrix_ele} += $Weights[$weight_index];
			$Matrix->{$fileC}{$unit_counter}{"UNITS"}{$matrix_ele} += 1;
			$VARS->{$matrix_ele} += 1;
#		    print "ME: $matrix_ele\n";
		    }
		}
		
	    }
	}
	$Matrix->{$fileC}{$unit_counter}{"ENERGY"} = $sequence_energy;
	$eng_total += $sequence_energy;
    }

    print "TOTAL ENERGY: $eng_total kcal/mol...";
}

sub GetSequenceEnergy(@) {
    # Employs a switching function too handle the ends: give them 3 times the weight for the nonbonds
    # Conserves the total energy of the system

    my ($curr_weight, $unit_1, $unit_2, $unit_counter, $level, $isEnd) = @_;
    my ($weight_1, $weight_2, $sequence_energy);

    $sequence_energy = 0.0;

    # Check too see if it's an end
    if ($level > 0) {
	if ($isEnd) {
	    if ($unit_counter == 1) {
		# Assign 3 times the weight too unit 1's nonbonds 5' End
		$weight_1 = 3 * $curr_weight;
		$weight_2 = $curr_weight;
	    } else {
		# 3' End
		$weight_1 = $curr_weight;
		$weight_2 = 3 * $curr_weight;

	    }
	} else {
	    $weight_1 = $weight_2 = $curr_weight;
	}
	    
    } else {
	$weight_1 = $weight_2 = $curr_weight;
    }

    $sequence_energy = $weight_1 * ($unit_1->{"VDW"} + $unit_1->{"Coulomb"}) +
	$weight_2 * ($unit_2->{"VDW"} + $unit_2->{"Coulomb"});

    return $sequence_energy;

}


sub GetWeights(@) {
    my ($num_pts) = $_[0];
    my (@weights_array, $denominator, $counter, $index, $running_sum, $a_val);
    
    $running_sum = 0.00;
    if ($num_pts == 0) {
	$weights_array[0] = 0.5;
	$denominator = 1;
    } else {
	for $counter (0 .. $num_pts) {
	    for $index (($counter + 1) .. $num_pts) {
		$running_sum += (1/($index - $counter));
	    }
	}

	$a_val = 0.25/$running_sum;
	$denominator = $num_pts/$a_val;
	
	for $index ((-$num_pts) .. $num_pts) {
	    if ($index != 0) {
		$counter = $a_val / (abs($index));
		push @weights_array, $counter;
	    } else {
		push @weights_array, 0;
	    }
	}
    }

    $denominator = 1;
    return ($denominator, @weights_array);

}

sub DetermineSequence (@) {
    my ($sequence) = @_;
    my ($return_str, @sequence_index, @sequence_name, $temp1, $temp2, $counter);

    @sequence_index = split /\s+/, $sequence;

    for $counter (@sequence_index) {
	push @sequence_name, GetUnitLetter($PDB_Info->{$counter});
    }
    if ($#sequence_index < 3) {
	$return_str = $sequence_name[0];
    } else {
	if (Priority($sequence_name[3]) > Priority($sequence_name[0])) {
	    $return_str = $sequence_name[3] .  $sequence_name[2] . "/" .
		$sequence_name[1] . $sequence_name[0];
	} else {
	    $return_str = $sequence_name[0] . $sequence_name[1] . "/" .
		$sequence_name[2] . $sequence_name[3];
	}
    }
    return $return_str;
}

sub Priority(@) {
#    This will assign weights according to: G > C > A > T
    my ($in_letter) = $_[0];

    if (uc($in_letter) eq "G") {
	return 4;
    } elsif (uc($in_letter) eq "C") {
	return 3;
    } elsif (uc($in_letter) eq "A") {
	return 2;
    } else {
	return 1;
    }
}

sub SolveLinearSystem(@) {
    my ($Matrix, $Variables) = @_;
    my ($link, $matrix_string, @var_holder, $counter, $vals_string, $cmd, $result);
    my (@holder, $err, $index, %STATS, %ERR, $EQNS, $file_counter, $cum, $total);

    $link = new Math::ematica '-linklaunch', '-linkname', 'math -mathlink';

    @var_holder = keys %{ $Variables };
    open ERRFILE, "> nn_dump" or die "Cannot write to file nn_dump: $!\n";
    $vals_string = "v = {";
    $matrix_string = "m = {";

    $total = 0;
    for $file_counter (sort Numerically keys %{ $Matrix }) {
	for $counter (sort Numerically keys %{ $Matrix->{$file_counter} }) {
	    $matrix_string .= "{";
	    for $index (@var_holder) {
		if (defined($Matrix->{$file_counter}{$counter}{"UNITS"}{$index})) {
		    $matrix_string .= sprintf("%12.10f", $Matrix->{$file_counter}{$counter}{"UNITS"}{$index}) . ",";
		} else {
		    $matrix_string .= "0.,";
		}
	    }
	    $total++;
	    chop $matrix_string;
	    $matrix_string .= "},";
	    $vals_string .= $Matrix->{$file_counter}{$counter}{"ENERGY"} . ",";
	}
	$EQNS->{$file_counter} = $total;
    }
    chop $matrix_string;
    chop $vals_string;
    $matrix_string .= "}";
    $vals_string .= "}";
    MathLinkCmd($matrix_string, 1, $link);
    MathLinkCmd($vals_string, 1, $link);
#    print ERRFILE "$matrix_string\n";
#    print ERRFILE "$vals_string\n";
#    Get Timing Info
    print "SVD Decomposition...";
    $cmd = "Timing[SingularValues[m];]";
    ($err, $result) = MathLinkCmd($cmd, 0, $link);
    ShowError($err, $link);
    if ($result =~ /(\d+\.*\d*)/) {
	print "Done in $1 seconds\n";
    } else {
	print "Done: $result\n";
    }
    
#    Do Actual Decomposition
    $cmd = "{u,sig,c} = SingularValues[m];";
    MathLinkCmd($cmd, 1, $link);
    $cmd = "eta = DiagonalMatrix[sig];";
    ($err, $result) = MathLinkCmd($cmd, 1, $link);
    ShowError($err, $link);
    $cmd = "sig[[1]]/Last[sig]";
    ($err, $result) = MathLinkCmd($cmd, 0, $link);
    ShowError($err, $link);
    if ($result =~ /(\d+\.\d+)/) {
	print "Spectral Condition: $1\n";
    } else {
	print "Spectral Condition: $result\n";
    }
    
    $cmd = "{Ur, Cr} = Map[Transpose, {u,c}];";
    MathLinkCmd($cmd, 1, $link);
    $cmd = "q = Cr.Inverse[eta].Transpose[Ur];";
    MathLinkCmd($cmd, 1, $link);
  #  $cmd = "Norm[m-m.q.m]";
  #  ($err, $result) = MathLinkCmd($cmd, 0, $link);
  #  ShowError($err, $link);
  #  if ($result =~ /(\d+\.\d+)/) {
#	print "Rounding Error: $1\n";
#    } else {
#	print "The Rounding Error: $result\n";
#    }
    $cmd = "x = q.v";
    ($err, $result) = MathLinkCmd($cmd, 0, $link);
    ShowError($err, $link);
    $counter = 0;
    while ($result =~ /(-?\d+\.?\d*)/g) {
	$STATS{$var_holder[$counter]}{"VALS"} = $1;
#	print "$var_holder[$counter]: $1\n";
	$counter++;
    }
#   Error
    $cmd = "soln = m.x;";
    MathLinkCmd($cmd, 1, $link);
    $cmd = "err = soln - v";
    ($err, @holder) = MathLinkCmd($cmd, 2, $link);
    ShowError($err, $link);
    $link->PutFunction("Exit", 0);
    $link->EndPacket;
    $link = ();

    $cum = $file_counter = $total = $result = $counter = 0;
    
    for $index (0 .. $#holder) {
	if (! defined($EQNS->{$file_counter})) {
	    $file_counter++;
	} elsif (($index + 1) > $EQNS->{$file_counter}) {
	    $result = sqrt($result);
	    $ERR{$file_counter}{"NORM"} = $result;
	    $ERR{$file_counter}{"TOTAL"} = $total;
	    $result = $counter = $total = 0;
	    $file_counter++;
	}
	$result += $holder[$index]**2;
	$cum += $holder[$index]**2;
	$total += $holder[$index];
	$ERR{$file_counter}{"UNITS"}{$counter} = $holder[$index];
	$counter++;
    }
    $result = sqrt($result);
    $ERR{$file_counter}{"NORM"} = $result;
    $ERR{$file_counter}{"TOTAL"} = $total;
#
    close ERRFILE;	
    $ERR{"TOTAL"} = sqrt($cum);
    return (\%STATS, \%ERR);
}

sub ShowError(@) {
    my ($err, $link) = @_;
    if ($err) {
	$link->PutFunction("Exit", 0);
	$link->EndPacket;
	$link = ();
	die "ERROR: $err\n";
    }
}
    
sub SolveLinearSystem_Defunct(@) {
    my ($Matrix, $Variables, $STATS) = @_;
    my (@var_holder, $counter, $matrix_string, $vals_string, $index, @holder);
    my ($return_str, $stop, @Eqns, @Energies, @Unknowns, %SubMatrix);

    print ".";
    @var_holder = keys %{ $Variables };
    @Eqns = @Energies = @Unknowns = ();

    for $counter (keys %{ $Matrix }) {
	$matrix_string = "{";
	for $index (@var_holder) {
	    if (defined($Matrix->{$counter}{"UNITS"}{$index})) {
		$matrix_string .= $Matrix->{$counter}{"UNITS"}{$index} . ",";
		if (! defined($Matrix->{$counter}{"USED"})) {
		    $Matrix->{$counter}{"USED"} = 0;
		}
#		$Matrix->{$counter}{$index} = 1;
#		$Unknowns[$counter]{$index} = 1;
	    } else {
		$matrix_string .= "0,";
	    }
	}
	chop $matrix_string;
	$matrix_string .= "}";
	$Eqns[$counter] = $matrix_string;
	$Energies[$counter] = $Matrix->{$counter}{"ENERGY"};
    }

  LOOP: while ($#Eqns >= $#var_holder) {
      for $counter (@var_holder) {
	  $stop = 0;
	  @Unknowns = keys %{ $Matrix };
	  %SubMatrix = ();
	  for $index (1 .. $#Unknowns) {
	      if (defined($Matrix->{$index}{"UNITS"}{$counter}) and ! $Matrix->{$index}{"USED"}) {
		  $SubMatrix{$index} = 1;
		  $stop = 1;
		  last;
	      }
	  }
	  if (! $stop) {
	      delete $Variables->{$counter};
	      @holder = keys %{ $Variables };
	      if ($#holder > 0) {
		  print "Decended: Removed $counter\n";
		  SolveLinearSystem(\%{ $Matrix }, \%{ $Variables }, \%{ $STATS });
	      }
	      last LOOP;
	  }
      }
      
      for $index (keys %SubMatrix) {
	  $matrix_string = $vals_string = "";

	  $matrix_string .= $Eqns[$index] . ","; 
	  $vals_string .= $Energies[$index] . ",";
	  $Matrix->{$index}{"USED"} = 1;
      }

      chop $matrix_string;
      chop $vals_string;
      $matrix_string = "m = {" . $matrix_string . "}";
      $vals_string = "v = {" . $vals_string . "}";
 
      $return_str = SolveMatrix($matrix_string, $vals_string);
      for $index (0 .. $#{ $return_str }) {
	  $STATS->{$var_holder[$index]} .= $return_str->[$index] . " ";
      }
      print "Got STATS\n";
  }
    return $STATS;

}

sub SolveMatrix(@) {
    my ($matrix_string, $vals_string) = @_;
    my ($MathLinkResult, @results_array, $link, $error);

    $link = new Math::ematica '-linklaunch', '-linkname', 'math -mathlink';
    MathLinkCmd($matrix_string, 1, $link);
    MathLinkCmd($vals_string, 1, $link);
    ($error, $MathLinkResult) = MathLinkCmd("LinearSolve[m,v]",0, $link);
    if ($MathLinkResult =~ /LinearSolve/) {
	$matrix_string =~ s/\},/\}\n/g;
	$matrix_string =~ s/\{\{/\{\n\{/g;
	$vals_string =~ s/\},/\}\n/g;
	print "ERROR in solving system!\n";
	print ERRFILE "$MathLinkResult\n";
	print ERRFILE "MathLinkError: " . $error . "\n";
	print ERRFILE "Coefficient Matrix==\n" . $matrix_string . "\n";
	print ERRFILE "Energy Vector==\n" . $vals_string . "\n";
	$link->PutFunction("Exit", 0);
	$link->EndPacket;
	$link = ();
	print "Aborting\n";
    } else {
	$link->PutFunction("Exit", 0);
	$link->EndPacket;
	$link = ();
	while ($MathLinkResult =~ /(-?\d+\.?\d*)/g) {
	    push @results_array, $1;
	}
    }

    return \@results_array;
}    

sub WriteOutput(@) {
    my ($fle_name, $STATS, $Err, $Variables) = @_;
    my ($counter, $out_string,$avg, $stdev, $index, @holder, $total);

    $out_string = "";
    for $counter (keys %{ $STATS }) {
	$avg = $STATS->{$counter}{"VALS"};
	$out_string .= sprintf("%-10s%12.3f%8d\n", $counter, $avg, $Variables->{$counter});
	    
    }

    open OUTFILE, "> $fle_name" or die "Cannot write too $save_name: $!\n";
    printf OUTFILE "%-10s%12s%8s\n", "Sequence", "del(H)", "Points";
    for $counter (0 ... 40) {
	print OUTFILE ".";
    }
    print OUTFILE "\n";
    print OUTFILE $out_string;
    
    printf OUTFILE "\n--==ERROR==--\n";
    printf OUTFILE "%-10s%12s\n", "ResID", "Error";
    $out_string = "";
    for $counter (0 ... 30) {
	print OUTFILE ".";
    }
    
    @holder = keys %{ $Err };
    $index = 0;
    while ($index <= $#holder) {
	if ($holder[$index] !~ /\d+/) {
	    splice @holder, $index, 1;
	} else {
	    $index++;
	}
    }

    $total = 1;
    for $index (sort Numerically @holder) {
	next
	    if ($index eq "TOTAL");
	$out_string .= "\nFILE $index\n";
	for $counter (sort Numerically keys %{ $Err->{$index}{"UNITS"} }) {
	    $out_string .= sprintf("%-10s%12.8f\n",($counter + 1), $Err->{$index}{"UNITS"}{$counter});
	    $total++;
	}
	
	$out_string .= sprintf("%-10s%12.8f\n", "TOTAL", $Err->{$index}{"TOTAL"});
	$out_string .= sprintf("%-10s%12.8f\n", "NORM", $Err->{$index}{"NORM"});
	
    }
    $out_string .= "\n" . sprintf("%-10s%12.8f\n", "E_NORM", $Err->{"TOTAL"});
	    
    printf OUTFILE $out_string;
    close OUTFILE;

}

sub MathLinkCmd(@) {

    my ($inExpression, $isTerminal, $link) = @_;
    my ($tmp, $err, @holder);

    print ERRFILE "$inExpression\n";
    $link->PutFunction("EvaluatePacket",1);
    if (! $isTerminal) {
	$link->PutFunction("ToString", 1);
    }
    $link->PutFunction("ToExpression", 1);
    $link->PutString("$inExpression");
    $link->EndPacket;

    while ($link->NextPacket != RETURNPKT) {
	$link->NewPacket;
    } 
    if ($isTerminal == 1) {
	$link->NewPacket;
	if (! $link) {
	    print "ERROR: " . $link->ErrorMessage . "\n";
	}
    } else {
	$err = $link->ErrorMessage;
	if ($isTerminal == 2) {
	    @holder = $link->GetRealList;
	} else  {
	    $tmp = $link->GetString;
	}
	if ($err =~ /ok so far/) {
	    $err = "";
	}
	if ($isTerminal == 2) {
	    return ($err, @holder);
	} else {
	    return ($err, $tmp);
	}	    
    }
}

sub GetUnitLetter(@) {
    my ($instring) = $_[0];

    if ($instring =~ /D(\w)/) {
	$instring = $1;
    }

    return $instring;

}

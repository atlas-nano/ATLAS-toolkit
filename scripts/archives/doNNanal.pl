#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}

use strict;
use File::Basename;
use Packages::General;
use Packages::FileFormats;
use Packages::HelixLayout;
use Packages::GetParms;

sub GetValidFiles(@);
sub ParsePDBFile(@);
sub CalculateNN(@);
sub WriteOutput(@);
sub GetHelixLayout(@);
sub DetermineSequence(@);
sub GetUnitLetter(@);
sub Numerically;
sub GetReferenceEnergy(@);
sub GetWeight(@);
sub PrintVio();

die "usage: $0 directory|file ref_energies level\n"
    if (! @ARGV or $#ARGV < 2);

my ($work_dir, $ref_energy, $level) = @ARGV;
my ($file_list, $dna, %NN_Data, $Helix_Info, $Energy_Info, $Energy_Ref);
my ($PDB_Info, $strand_length, $single_helix, %Tracker, %Vio_Counter);

{
    $file_list = GetValidFiles($work_dir);

    print "--=== START ===--\n\nSucessfully recognized " . ($#{ $file_list } + 1) . " files\n";

    $Energy_Ref = GetReferenceEnergy($ref_energy);
    
    for $dna (@{ $file_list }) {
	print "Processing ". $dna->{"PDBFILE"} . "....";
	%NN_Data = ();
	%Tracker = ();
	%Vio_Counter = ();
	($strand_length, $PDB_Info) = ParsePDBFile($dna->{"PDBFILE"}, 0);
	($Helix_Info, $single_helix) = GetHelixLayout($strand_length, $dna->{"PARMFILE"});
	if ($single_helix) {
	    ($Energy_Info) = GetMPSimEnergy($dna->{"DATFILE"}, ($strand_length * 2));
	}else {
	    ($Energy_Info) = GetMPSimEnergy($dna->{"DATFILE"}, ($strand_length * 4));
	}  
	
	CalculateNN($strand_length, $level);
	WriteOutput($dna->{"OUTFILE"});
    }

    PrintVio();
    
    print "\n---==== END ===--\n";
}

sub GetValidFiles(@) {
    my ($loc) = $_[0];
    my (@pdbfiles, $curr_file, $dat_file, $parm_file, @holder, $counter, $rec, $saveName);

    if (-e $loc) {
	if (! -d $loc) {
	    $holder[0] = $loc;
	} else {
	    opendir (WORKDIR, $loc) or die "Cannot access specified directory: $!\n";
	    @holder = grep { /\w+\.pdb$/ } map {"$loc/$_"} readdir WORKDIR;
	    closedir WORKDIR;
	}
    }

    $level = Trim($level);
    die "ERROR: Expected integer for level got: $level\n"
	if (! IsInteger($level));
    die "ERROR: level has to be between 1 and 5\n"
	if ($level < 1 or $level > 5);
    die "ERROR: Cannot find energy reference file: $ref_energy\n"
	if (! -e $ref_energy);
    die "ERROR: Cannot find any valid pdbfile here: $loc\n"
	if (! @holder);

    for $curr_file (@holder) {
	$curr_file =~ /^(.+)\.pdb$/;
	$dat_file = $1 . ".dat";
	$parm_file = $1 . ".parm";
	$saveName = $1 . ".anal";
	if (-e $dat_file and -e $parm_file) {
	    $rec = (
		    {
			"PDBFILE"  => $curr_file,
			"DATFILE"  => $dat_file,
			"PARMFILE" => $parm_file,
			"OUTFILE"  => $saveName,
		    }
		    );
	    push @pdbfiles, $rec;
	}
	
    }

    die "ERROR: No valid files found in $work_dir. Aborting\n"
	if (! @pdbfiles);
    return \@pdbfiles;

}

sub GetVals(@) {
    my ($in_string) = $_[0];
    my (@result);

    if ($in_string =~ /^(.+)\_(\d+)bp\_([1|0])/) {
	@result = ($1, $2, $3);
    } else {
	print "Error extracting information from $in_string\n";
    }

    return @result;
}

sub ParsePDBFile(@) {
    my ($pdb_file) = $_[0];
    my ($Atom_Info) = GetPDBFileInfo($pdb_file, 0, " ", " ",0);
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
    my ($length_strand, $parmfile) = @_;
    my (%Helix, $rec, $i, $j, $is_single_helix);
    my ($P_File, $h1, @h_info, $helix_counter, $counter); 
    my ($start_res, $end_res, $curr_strand);

    $P_File = Packages::GetParms->new();
    die "Error in Paramater file\n"
	if (! $P_File->IsValidParams($parmfile));
    if ($P_File->{"Molecule"}{"crossovers"} and ($#{ $P_File->{"Molecule"}{"crossovers"} } > 0 or
						 $P_File->{"Molecule"}{"crossovers"}[0] != -1)) {
	$is_single_helix = 0;
    } else {
	$is_single_helix = 1;
    }

    if ($is_single_helix) {
	for $j (1 ... 2) {
	    for $i (1 ... $length_strand) {
		$Helix{1}{1}{$i} =  $i;
		$Helix{1}{2}{$i} = ($length_strand * 2) - $i + 1;
	    }
	}
    } else {

	die "ERROR: Parmfile is required for double crossover molecules"
	    if (! $parmfile);
	die "ERROR: Cannot locate Paramfile: $parmfile: $!\n"
	    if (! -e $parmfile);
	$h1 = Packages::HelixLayout->spawn();
	$h1->DetermineHelixLayout(
				  $P_File->{"Molecule"}->{"major_groove"}, 
				  $P_File->{"Molecule"}->{"minor_groove"}, 
				  $P_File->{"Molecule"}->{"is3PrimeIn"}, 
				  $P_File->{"Molecule"}->{"bases_at_end"}, 
				  ($P_File->{"Molecule"}->{"total_bases"}/4), 
				  @{ $P_File->{"Molecule"}->{"crossovers"} }
				  );
	
	@h_info = $h1->GetHelixInfo();

	for $helix_counter (0 .. $#h_info) {
	    for $j (0 .. $#{ $h_info[$helix_counter] } ) {
		$counter = 1;
		for $i (0 .. $#{ $h_info[$helix_counter][$j] }) {
		    $start_res = $h_info[$helix_counter][$j][$i]->{"StartUnit"};
		    $end_res = $h_info[$helix_counter][$j][$i]->{"EndUnit"};
		    $curr_strand = $h_info[$helix_counter][$j][$i]->{"Strand"};
		    if (($curr_strand % 2) == 0) {
			for (reverse $start_res ... $end_res) {
			    $Helix{$helix_counter}{$j}{$counter} = ( ($curr_strand - 1) * $length_strand ) + $_;
			    $counter++;
			}
		    } else {
			for ($start_res ... $end_res) {
			    $Helix{$helix_counter}{$j}{$counter} = ( ($curr_strand - 1) * $length_strand ) + $_;
			    $counter++;
			}
		    }
		}
	    }
	}
			
    }

    return (\%Helix, $is_single_helix);
	    
}

sub CalculateNN(@) {
    my ($strand_l, $level) = @_;
    my ($helix_counter, $unit_counter, $curr_unit, $previous_unit, $next_unit, $holder, $Ref, $isEnd);
    my ($sequence_energy, $sequence_name, $sequence_value, $chain_counter, @keys_array, $unit_tracker, $indicator);
    my ($base_name, $prev_name, $next_name, $base_val, $factorC, $factorN, $factorP, $factor);

    %NN_Data = ();

    $factor = 0.25;
    for $helix_counter (sort Numerically keys %{ $Helix_Info }) {
	@keys_array = keys %{ $Helix_Info->{$helix_counter}{1} };
	@keys_array = ();
	for $unit_counter (1 .. $strand_l) {
	    $sequence_energy = $isEnd = 0;
	    $sequence_name = $sequence_value = $holder = "";
	    if ($unit_counter == 1 or $unit_counter == $strand_l) {
		$isEnd = 1;
	    } else {
		$isEnd = 0;
	    }

	    $prev_name = $next_name = "";
	    $Ref = $sequence_energy = 0;
	    $chain_counter = 1;
	    for $chain_counter (sort Numerically keys %{ $Helix_Info->{$helix_counter} }) {
		$curr_unit = $Helix_Info->{$helix_counter}{$chain_counter}{$unit_counter};
		$sequence_energy += $Energy_Info->{$curr_unit}->{"Bond"};		
#		Internal Energy
		$base_name = $PDB_Info->{$curr_unit};
		if ($isEnd and $chain_counter == 1) {
		    $base_name = $PDB_Info->{$Helix_Info->{$helix_counter}{"1"}{$unit_counter}};
		    # Add the end corrections
		    if ($base_name =~ /D?(\w)(\d)/) {
			$base_name = $2 . "'" . $1;
			$base_val = $2 . ord($1);

			if (! defined($Energy_Ref->{$base_val})) {
			    print "NONE: $base_name\n";
			    $Ref += ($Energy_Info->{$curr_unit}->{"VDW"} +
				     $Energy_Info->{$curr_unit}->{"Coulomb"});
			} else {
			    $Ref += $Energy_Ref->{$base_val};
			}
		    }
		}

		if ($chain_counter == 1) {
		# Add the bond energy too the reference
		    $base_name = $PDB_Info->{$curr_unit};
		    $base_name = GetUnitLetter($base_name);
		    $base_val = ord($base_name);
		    if (! defined($Energy_Ref->{$base_val})) {
			print "NONE: $base_name $base_val\n";
			$Ref += ($Energy_Info->{$curr_unit}->{"VDW"} +
				 $Energy_Info->{$curr_unit}->{"Coulomb"});
		    } else {
			$Ref += $Energy_Ref->{$base_val};
		    }
		}

#		NonBond Energy
		if ($level == 1) {
		    $sequence_energy += ($Energy_Info->{$curr_unit}->{"VDW"} + 
					 $Energy_Info->{$curr_unit}->{"Coulomb"});
		    $factorC = 1;
		} else {
		    if ($unit_counter == 1) {
			$factorC = 0.75;
			$factorN = 0.25;
		    } elsif ($unit_counter == $strand_length) {
			$factorC = 0.75;
			$factorP = 0.25;
		    } else {
			$factorC = $factorP = $factorN = 0.25;
		    }
		
		    if (defined($Helix_Info->{$helix_counter}{$chain_counter}{$unit_counter - 1})) {
			$next_unit = $Helix_Info->{$helix_counter}{$chain_counter}{$unit_counter - 1};
			$prev_name .= $next_unit . " " . $curr_unit . " ";

			$sequence_energy += $factorC * ($Energy_Info->{$curr_unit}->{"VDW"} + 
							$Energy_Info->{$curr_unit}->{"Coulomb"});
			$sequence_energy += $factorP * ($Energy_Info->{$next_unit}->{"VDW"} + 
							$Energy_Info->{$next_unit}->{"Coulomb"});
		    }
		    
		    if (defined($Helix_Info->{$helix_counter}{$chain_counter}{$unit_counter + 1})) {
			$next_unit = $Helix_Info->{$helix_counter}{$chain_counter}{$unit_counter + 1};
			$next_name .= $curr_unit . " " . $next_unit . " ";

			$sequence_energy += $factorC * ($Energy_Info->{$curr_unit}->{"VDW"} + 
							$Energy_Info->{$curr_unit}->{"Coulomb"});
			$sequence_energy += $factorN * ($Energy_Info->{$next_unit}->{"VDW"} + 
							$Energy_Info->{$next_unit}->{"Coulomb"});
		    }
		}
		
	    }

	    if ($prev_name ne "") {
		($sequence_name, $sequence_value) = DetermineSequence($prev_name);
		if (! defined($Energy_Ref->{$sequence_value})) {
		    print "NOT FOUND: $sequence_name, $sequence_value\n";
		} else {
		    $Ref += $Energy_Ref->{$sequence_value};
		    $base_name = $sequence_name . "-" . $base_name;
		    $base_val = $sequence_value . $base_val;
		}
	    }

	    if ($next_name ne "") {
		($sequence_name, $sequence_value) = DetermineSequence($next_name);
		if (! defined($Energy_Ref->{$sequence_value})) {
		    print "NOT FOUND: $sequence_name, $sequence_value\n";
		} else {
		    $Ref += $Energy_Ref->{$sequence_value};
		    $base_name = $base_name. "-" . $sequence_name;		
		    $base_val = $base_val . $sequence_value;
		}
	    }
		
	    
	    $unit_tracker = $helix_counter * 100 + $unit_counter;
	    $Tracker{$unit_tracker} = (
				       {
					   "SEQUENCE NAME" => $base_name,
					   "SEQUENCE ENERGY" => $Ref,
					   "REFERENCE ENERGY" => $sequence_energy,
					   "DIFFERENCE" => $sequence_energy - $Ref,
					   
				       }
				       );
	    
	    $NN_Data{$base_val}{"ENERGY"} .= "$sequence_energy ";
	    $NN_Data{$base_val}{"NAME"} = $sequence_name;

	}
    }
}

sub DetermineSequence (@) {
    my ($sequence) = @_;
    my (@return_array, @sequence_index, @sequence_name, $temp1, $temp2, $counter);

    @sequence_index = split /\s+/, $sequence;

    for $counter (@sequence_index) {
	push @sequence_name, GetUnitLetter($PDB_Info->{$counter});
    }

    if (GetWeight($sequence_name[3]) > GetWeight($sequence_name[0])) {
	$return_array[0] = $sequence_name[3] . "-" . $sequence_name[2] . ":" .
	    $sequence_name[1] . "-" . $sequence_name[0];
	$return_array[1] = ord($sequence_name[3]) . ord($sequence_name[2]) .
	    ord($sequence_name[1]) . ord($sequence_name[0]);
    } else {
	$return_array[0] = $sequence_name[0] . "-" . $sequence_name[1] . ":" .
	    $sequence_name[2] . "-" . $sequence_name[3];
	$return_array[1] = ord($sequence_name[0]) . ord($sequence_name[1]) .
	    ord($sequence_name[2]) . ord($sequence_name[3]);
    }
    
    $return_array[0] =~ s/:/\//g;
    $return_array[0] =~ s/\-//g;
    return @return_array;
}


sub GetWeight(@) {
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

sub WriteOutput(@) {
    my ($out_name) = @_;
    my ($tracker, $running_total, $predicted_total, $strain_total, $out_string);
    my ($unit_counter, $seq_name, $counter);

    $tracker = "tracker_" . basename($out_name);

    $running_total = $predicted_total = 0;
    $strain_total = 0;
    $out_string = sprintf("%-12s%12s%8s%8s\n", "SEQUENCE", "AVG", "+/-", "Total");

    $predicted_total = $strain_total = $running_total = 0;
    print "Creating $tracker....";
    open OUTFILE, "> $tracker" or die "Cannot create $tracker: $!\n";
    printf OUTFILE "%-8s%15s%12s%12s%12s\n", "#", "SEQ NAME", "SEQ ENG", "REF ENG", "DIFF";
    $counter = 1;
    for $unit_counter (sort Numerically keys %Tracker) {
	$predicted_total +=  $Tracker{$unit_counter}{"REFERENCE ENERGY"};
	$running_total += $Tracker{$unit_counter}{"SEQUENCE ENERGY"};
	$strain_total += $Tracker{$unit_counter}{"DIFFERENCE"};
	printf OUTFILE "%8d%15s%12.3f%12.3f%12.3f\n", $counter, $Tracker{$unit_counter}{"SEQUENCE NAME"}, 
	$Tracker{$unit_counter}{"SEQUENCE ENERGY"},
	$Tracker{$unit_counter}{"REFERENCE ENERGY"}, $Tracker{$unit_counter}{"DIFFERENCE"};
	$seq_name = $Tracker{$unit_counter}{"SEQUENCE NAME"};
	if ($Tracker{$unit_counter}{"DIFFERENCE"} > 0) {
	    $Vio_Counter{$seq_name}{"POSITIVE"} += 1;
	}
	$Vio_Counter{$seq_name}{"TOTAL"} += 1;
	$counter++;
	
    }
    printf OUTFILE "%-8s%15s%12.3f%12.3f%12.3f\n", "#", "TOTAL", $running_total, $predicted_total, $strain_total;
    close OUTFILE;
    print "Done\n";
}

sub PrintVio() {
    my ($out_name) = "nn_stats.txt";
    my ($seq_name, $pos);

    $out_name =~ s/\.anal$/nn_stats.anal/g;
    print "Creating Stats...";
    open OUTFILE, "> $out_name" or die "Cannot create $out_name: $!\n";
    printf OUTFILE "%-12s%5s%5s%5s\n", "SEQ NAME", "TOT", "POS", "NEG";
    for $seq_name (keys %Vio_Counter) {
	$pos = 0;
	if (defined($Vio_Counter{$seq_name}{"POSITIVE"})) {
	    $pos = ($Vio_Counter{$seq_name}{"POSITIVE"});
	}
	printf OUTFILE "%-12s%5d%5d%5d\n", $seq_name, $Vio_Counter{$seq_name}{"TOTAL"}, $pos,
	($Vio_Counter{$seq_name}{"TOTAL"} - $pos);
    }
    close OUTFILE;
    print "Done\n";

}
	
sub Numerically {
    ($a <=> $b);
}

sub GetUnitLetter(@) {
    my ($instring) = $_[0];

    if ($instring =~ /D(\w)/) {
	$instring = $1;
    }

    return $instring;

}

sub GetReferenceEnergy(@) {
    my ($infile) = $_[0];
    my ($sequence_key, %Ref_Energy, $has_ref);

    $has_ref = 0;
    if (open(INFILE, $infile) or die "Cannot open reference file: $infile\n") {
	while (<INFILE>) {
	    chomp;
	    if ($_ =~ /^(\w)(\w)\/(\w)(\w)\s+(\-?\d+\.\d+)/) {
		$sequence_key = ord($1) . ord($2) . ord($3) . ord($4);
		$Ref_Energy{$sequence_key} = $5;
		$has_ref = 1;
	    }elsif ($_ =~ /^(\d).(\w)\s+(\-?\d+\.\d+)/) {
		$sequence_key = $1 . ord($2);
		$Ref_Energy{$sequence_key} = $3;
#		print "$_ $sequence_key $3\n";
		$has_ref = 1;
	    }elsif ($_ =~ /^([a-zA-Z])\s+(\-?\d+\.\d+)\s+\d+$/) {
		$sequence_key = ord($1);
		$Ref_Energy{$sequence_key} = $2;
		$has_ref = 1;
	    }
	}
	close INFILE;
    }
    die "ERROR: No Valid data found in reference file: $infile\n"
	if (! $has_ref);
    return \%Ref_Energy;
}

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

sub GetValidFiles();
sub GetVals(@);
sub ParsePDBFile(@);
sub CalculateNN();
sub WriteOutput(@);
sub GetHelixLayout(@);
sub DetermineSequence(@);
sub GetUnitLetter(@);
sub Numerically;
sub GetReferenceEnergy(@);
sub IsValidAtom(@);

die "usage: $0 directory [parm_file] [ref_energies]\n"
    if (! @ARGV);

my ($work_dir, $parm_file, $ref_energy) = @ARGV;
my (%Composite, @file_list, $pdbfile, $energyfile);
my (%NN_Data, $PDB_Info, $Helix_Info, $Energy_Info, $Energy_Ref, $has_ref);
my ($out_name, $strand_length, $single_helix, %Tracker);

GetValidFiles();

$has_ref = 0;
print "--=== START ===--\n\nSucessfully recognized " . ($#file_list + 1) . " files\n";

if ($ref_energy && -e $ref_energy) {
    ($Energy_Ref) = GetReferenceEnergy($ref_energy);
}

for $pdbfile (@file_list) {
    print "Processing $pdbfile....";
    ($out_name, $strand_length, $single_helix) = GetVals($pdbfile);

    if ($out_name) {
	$out_name .= ".anal";
	$energyfile = $pdbfile;
	$energyfile =~ s/\.pdb/\.dat/g;
	$energyfile = $work_dir . "/" . $energyfile;
	(%NN_Data) = ();
#	(%Atoms_Info) = GetPDBFileInfo($pdbfile, 0, " ", " ");
	($PDB_Info) = ParsePDBFile($pdbfile, 0);
	($Helix_Info) = GetHelixLayout($strand_length, $single_helix, $parm_file);
	if ($single_helix) {
	    ($Energy_Info) = GetMPSimEnergy($energyfile, ($strand_length * 2));
	}else {
	    ($Energy_Info) = GetMPSimEnergy($energyfile, ($strand_length * 4));
	}  
#	}	    
	
	CalculateNN();
	WriteOutput($out_name, 0);
    }
}

if ($#file_list > 0) {
    print "Creating composite Nearest Neighbour Analysis File...";
    WriteOutput("nn_composite.anal", 1);
}
print "\n---==== END ===--\n";

sub GetValidFiles() {
    my (@pdbfiles, $curr_file, $dat_file);

    opendir (WORKDIR, $work_dir) or die "Cannot access specified directory: $!\n";
    @pdbfiles = grep { /^\w+\.pdb$/ } readdir WORKDIR;
    closedir WORKDIR;

    for $curr_file (@pdbfiles) {
	if ($curr_file =~ /^(.+)\.pdb$/) {
	    $dat_file = $work_dir . "/" . $1 . ".dat";
	    if (-e $dat_file) {
		push @file_list, $work_dir . "/" . $curr_file;
	    }
	}
    }

    die "ERROR: No valid files found in $work_dir. Aborting\n"
	if ($#file_list < 0);

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
    my ($pdb_file, $include_solvent) = @_;
    my (%Atom_Info) = GetPDBFileInfo($pdb_file, $include_solvent, " ", " ");
    my ($counter, $curr_res_name, %Res_Info, $curr_res_id);

    $curr_res_id = 0;
    for $counter (sort Numerically keys %Atom_Info) {
	if ($Atom_Info{$counter}{"RES_ID"} > $curr_res_id) {
	    $curr_res_id = $Atom_Info{$counter}{"RES_ID"};
	    $curr_res_name = $Atom_Info{$counter}{"RES_NAME"};
#	    if (IsValidAtom($Atom_Info{$counter}{"LABEL"})) {
		$Res_Info{$curr_res_id} = $curr_res_name;
#	    }
	}
    }

    return \%Res_Info;

}

sub GetHelixLayout(@) {
    my ($length_strand, $is_single_helix, $parmfile) = @_;
    my (%Helix, $rec, $i, $j);
    my ($P_File, $h1, @h_info, $helix_counter, $counter); 
    my ($start_res, $end_res, $curr_strand);

    if ($is_single_helix) {
	for $j (1 ... 2) {
	    for $i (1 ... $length_strand) {
		$Helix{1}{$j}{$i} = (($length_strand * ($j -1)) + $i);
	    }
	}
    } else {

	die "ERROR: Parmfile is required for double crossover molecules"
	    if (! $parmfile);
	die "ERROR: Cannot locate Paramfile: $parmfile: $!\n"
	    if (! -e $parmfile);
	$P_File = Packages::GetParms->new();
	die "Error in Paramater file\n"
	    if (! $P_File->IsValidParams($parmfile));
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
		    ($start_res, $end_res) = ($end_res, $start_res) 
			if ($start_res > $end_res);
		    for ($start_res ... $end_res) {
			$Helix{$helix_counter}{$j}{$counter} = ( ($curr_strand - 1) * $length_strand ) + $_;
			$counter++;
		    }
		}
	    }
	}
			
    }

    return \%Helix;
	    
}

sub CalculateNN() {
    my ($helix_counter, $unit_counter, $curr_unit, $previous_unit, $next_unit);
    my ($sequence_energy, $sequence_name, $sequence_value, $chain_counter, @keys_array, $unit_tracker);

    %NN_Data = ();

    for $helix_counter (sort Numerically keys %{ $Helix_Info }) {
	for $chain_counter (sort Numerically keys %{ $Helix_Info->{$helix_counter} }) {
	    $unit_counter = 2;
	    while ($unit_counter <= ($strand_length - 1)) {
		$curr_unit = $Helix_Info->{$helix_counter}{$chain_counter}{$unit_counter};
#		$previous_unit = $Helix_Info->{$helix_counter}{$chain_counter}{$unit_counter - 1};
		$next_unit = $Helix_Info->{$helix_counter}{$chain_counter}{$unit_counter + 1};
		$sequence_energy = $Energy_Info->{$curr_unit}->{"Bond"} + $Energy_Info->{$next_unit}->{"Bond"} +
		    0.5 * ($Energy_Info->{$curr_unit}->{"VDW"} + $Energy_Info->{$next_unit}->{"VDW"} +
			   $Energy_Info->{$curr_unit}->{"Coulomb"} + $Energy_Info->{$next_unit}->{"Coulomb"});
		($sequence_name, $sequence_value) = DetermineSequence($curr_unit, $next_unit);
		if ($has_ref) {
		    $unit_tracker = $helix_counter * 10000 + $chain_counter * 100 + $unit_counter;
		    $Tracker{$unit_tracker} = (
					       {
						   "SEQUENCE NAME" => $sequence_name,
						   "SEQUENCE ENERGY" => $sequence_energy,
						   "REFERENCE ENERGY" => $Energy_Ref->{$sequence_value},
						   "DIFFERENCE" => $sequence_energy - $Energy_Ref->{$sequence_value}
						   
					       }
					       );
		}
		if (! $NN_Data{$sequence_value}) {
		    $NN_Data{$sequence_value}{"ENERGY"} = $sequence_energy ;
		    $NN_Data{$sequence_value}{"NAME"} = $sequence_name;
		} else {
		    $NN_Data{$sequence_value}{"ENERGY"} .= " $sequence_energy";
		    $NN_Data{$sequence_value}{"NAME"} = $sequence_name;
		}
		
#	    Composite
		if (! $Composite{$sequence_value}) {
		    $Composite{$sequence_value}{"ENERGY"} = $sequence_energy ;
		    $Composite{$sequence_value}{"NAME"} = $sequence_name;
		} else {
		    $Composite{$sequence_value}{"ENERGY"} .= " $sequence_energy";
		    $Composite{$sequence_value}{"NAME"} = $sequence_name;
		}
		$unit_counter++;
	    }
	}
    }
}

sub DetermineSequence (@) {
    my ($curr_u, $next_u) = @_;
    my ($curr_name, $next_name, $curr_alt, $next_alt);
    my (@return_array);

    $curr_name = GetUnitLetter($PDB_Info->{$curr_u});
#    $previous_name = GetUnitLetter($PDB_Info->{$previous_u});
    $next_name = GetUnitLetter($PDB_Info->{$next_u});

    $return_array[0] = $curr_name . "-" . $next_name;
    $return_array[1] = ord($curr_name) . ord($next_name);

    return @return_array;
}

sub WriteOutput(@) {
    my ($avg, $total, $stdev, $out_string, $unit_counter, $seq_name, $Out_Data, $running_total);
    my ($out_name, $is_composite) = @_;
    my ($strain_total, $predicted_total, %Vio_Counter, $pos, $tracker, $counter);
    if ($is_composite) {
	$Out_Data = \%Composite;
    } else {
	$Out_Data = \%NN_Data;
    }

    if ($has_ref) {
	if ($is_composite) {
	    $tracker = "nn_tracker.txt";
	} else {
	    $tracker = "tracker_" . basename($out_name);
	}
    }

    $running_total = $predicted_total = 0;
    $strain_total = "";
    $out_string = sprintf("%-12s%12s%8s%8s\n", "SEQUENCE", "AVG", "+/-", "Total");
    for $unit_counter (sort Numerically keys %{ $Out_Data }) {
	($avg, $stdev, $total) = STDev($Out_Data->{$unit_counter}{"ENERGY"});
	$seq_name = $Out_Data->{$unit_counter}{"NAME"};
	$out_string .= sprintf("%-12s%12.3f%8.3f%8d\n", $seq_name, $avg, $stdev, ($total/$avg));
	$running_total += $total;
	if ($has_ref) {
	    $predicted_total += ($Energy_Ref->{$unit_counter} * (($total/$avg)));
	    $strain_total .= ($total - ($Energy_Ref->{$unit_counter} * (($total/$avg)))) . " ";
	}
    }	    
    if (! $is_composite) {
	$out_string .= sprintf("%-12s%12.3f\n", "#Total", $running_total);
	if ($has_ref) {
	    $out_string .= sprintf("%-12s%12.3f\n", "#Predicted", $predicted_total);
	    chop $strain_total;
	    ($avg, $stdev, $total) = STDev($strain_total);
	    $out_string .= sprintf("%-12s%12.3f%8.3f\n", "#Difference", $total, $stdev);
	}
    }
    if (open OUTFILE, "> $out_name") { 
	print OUTFILE $out_string;
	close OUTFILE;
	
	print "Done\nCreated file: $out_name\n";
    } else {
	print "ERROR: Cannot create file $out_name: $!\n";
    }
    if ($has_ref) {
	$predicted_total = $strain_total = $running_total = $counter = 0;
	print "Creating NN tracker file....";
	open OUTFILE, "> $tracker" or die "Cannot create $tracker: $!\n";
	printf OUTFILE "%-12s%12s%12s%8s\n", "#SEQ NAME", "SEQ ENG", "REF ENG", "DIFF";
	for $unit_counter (sort Numerically keys %Tracker) {
	    $counter++;
	    $predicted_total +=  $Tracker{$unit_counter}{"REFERENCE ENERGY"};
	    $running_total += $Tracker{$unit_counter}{"SEQUENCE ENERGY"};
	    $strain_total += $Tracker{$unit_counter}{"DIFFERENCE"};
	    printf OUTFILE "%8d%12s%12.3f%12.3f%8.3f\n", $counter, $Tracker{$unit_counter}{"SEQUENCE NAME"}, $Tracker{$unit_counter}{"SEQUENCE ENERGY"},
	    $Tracker{$unit_counter}{"REFERENCE ENERGY"}, $Tracker{$unit_counter}{"DIFFERENCE"};
	    $seq_name = $Tracker{$unit_counter}{"SEQUENCE NAME"};
	    if ($Tracker{$unit_counter}{"DIFFERENCE"} > 0) {
		$Vio_Counter{$seq_name}{"POSITIVE"} += 1;
	    }
	    $Vio_Counter{$seq_name}{"TOTAL"} += 1;

	}
	printf OUTFILE "%-12s%12.3f%12.3f%8.3f\n", "#TOTAL", $predicted_total, $running_total, $strain_total;
	close OUTFILE;
	print "Done\n";
	print "Creating Stats...";
	open OUTFILE, "> nn_stats.txt" or die "Cannot create nn_stats.txt: $!\n";
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
    my ($sequence_key, $is_valid);
    my (%Ref_Energy);

    if (open(INFILE, $infile)) {
	while (<>) {
	    chomp;
	    if ($_ =~ /^(\w)\-(\w)\s+(\-?\d+\.\d+)/) {
		$sequence_key = ord($1) . ord($2);
		$Ref_Energy{$sequence_key} = $3;
		$has_ref = 1;
	    }
	}
	close INFILE;
    }
    return \%Ref_Energy;
}

sub IsValidAtom(@) {
    my (@invalid_list); #list of atoms labels that are invalid
    my ($curr_atom_name) = $_[0];
    my ($return_val) = 1;

    @invalid_list = (
		     "C1'", "H1'", "C2'", "1H2'", "2H2'", "C3'", "H3'", "O3'",
		     "C4'", "H4'", "O4'", "C5'", "1H5'", "2H5'", "O5'", "P", 
		     "O1P", "O2P");
    
    for (@invalid_list) {
	if ($_ eq $curr_atom_name) {
	    $return_val = 0;
	    last;
	}
    }
    return $return_val;
}

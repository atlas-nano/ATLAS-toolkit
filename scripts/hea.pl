#!/usr/bin/perl -w
use strict;
use warnings;

my ($paramfile, $start_time, $end_time, $interval) = @ARGV;
my ($topfile, $trajfile, $mol_name, $total_bases);
my (%EngDatStructure, @unit_array);
my ($global_cuttoff, $run_ptraj);

$global_cuttoff = 0.0;
$run_ptraj = 1;

sub GetParameters();
sub ValidateInput();
sub CreateRstFiles();
sub DoEngAnal();
sub ProcessFile(@);
sub WriteOutput();
sub CalculateTotals();
sub CollectStats(@);
sub numerically { ($a<=>$b); }

# -- Start --
{
    if (! @ARGV or $#ARGV < 2) {
	die "usage: hea.pl paramfile start_time end_time interval\n";
    }

    print "\n-==Start Program ==-\n\n";
    GetParameters();
    ValidateInput();
    print "\n\nStep 1. Executing Ptraj...\n";
    print "-----------------------------------------------------------------\n\n";
    if ($run_ptraj) {
	CreateRstFiles();
    } else {
	print "Assuming restart files already exist\n";
    }
    print "-----------------------------------------------------------------\n\n";
    print "Step 2. Analyzing files...\n";
    print "-----------------------------------------------------------------\n\n";
    DoEngAnal();
    print "-----------------------------------------------------------------\n\n";
    
    @unit_array = sort numerically keys %EngDatStructure;
    
    WriteOutput();
    print "-==End Program==-\n\n";

}


sub GetParameters() {
    open PARAMFILE, $paramfile or die "Cannot open $paramfile: $!\n";
    my ($in_data, $c_list);

    if (! -d "outfiles") {
	mkdir "outfiles";
    }

   while (<PARAMFILE>) {
	chomp;
        $in_data = $_;
     
	if ($in_data =~ /Mol:\s(\d+).(\d+)/) {
	    $mol_name = $1 . $2;
	} elsif ($in_data =~ /Total bases:\s(\d+)/) {
	    $total_bases = $1;
	} elsif ($in_data =~ /Topology file:\s(.+)$/) {
	    $topfile = $1;
	} elsif ($in_data =~ /Trajectory file:\s(.+)$/) {
	    $trajfile = $1;
	} elsif ($in_data =~ /Cuttoff:\s(\d+\.\d+)$/) {
	    $global_cuttoff = $1;
	} elsif ($in_data =~ /Run Ptraj:\s([1|0])$/) {
	    $run_ptraj = $1;
	}
    }

    if ( ! $mol_name or ! $total_bases ) {
	die "Invalid Parameter file\n";
    }

}

sub ValidateInput() {
    -e $topfile or die "Cannot locate $topfile: $!\n";
    -e $trajfile or die "Cannot locate $trajfile: $!\n";
    
    if ($start_time =~ /(\d+)/) {
	$start_time = $1;
    } else {
	die "Invalid starting number. Expected integer\n";
    }
    
    if ($end_time =~ /(\d+)/) {
	$end_time = $1;
    } else {
	die "Invalid end number. Expected integer\n";
    }
    if ($interval =~ /(\d+)/) {
	$interval = $1;
    } else {
	die "Invalid interval. Expected integer\n";
    }
    
    if ($start_time > $end_time) {
	($start_time, $end_time) = Swap($start_time, $end_time);
    }
}

sub CreateRstFiles() {
    my ($out_cmd, $i, $j, $strand_ln);
    my ($starttm, $num_precision);

    if (! -d "rstfiles") {
	mkdir "rstfiles";
    }

    if (! -d "pdbfile") {
	mkdir "pdbfiles";
    }

    $num_precision = $end_time . "";
    $num_precision = length($num_precision);

    $starttm = sprintf("%0" . $num_precision . "d", $start_time);
#    $starttm = $start_time;


    $strand_ln = int($total_bases/2);

# Create restart files----------------------
    $out_cmd = "> tmp_gen_rst_file";

    open OUTCMD, $out_cmd or die "Cannot create file tmp_gen_rst_file: $!\n";
    print OUTCMD "trajin $trajfile $start_time $end_time $interval\n";

    $i = 1;
    while (($i * $strand_ln) < $total_bases) {

	print OUTCMD "center :1-" . ($i * $strand_ln) . " mass origin\n";
	print OUTCMD "image origin center\n";
	$i ++;
    }
    print OUTCMD "trajout rstfiles/" . $mol_name . "_rst restrat\n";

    close OUTCMD;

    $out_cmd = "ptraj $topfile < tmp_gen_rst_file > junk";
    
    system $out_cmd;
    system "rm -f junk tmp_gen_rst_file";

}

sub DoEngAnal() {

    my ($i, $curr_rst_file, $num_precision, $out_cmd, $j, $curr_eng_out_file);
    my ($my_data, $a_val, $b_val, $c_val, $d_val);
    my ($hash_key, $cutt_off);

    $cutt_off = 0;

    $num_precision = $end_time . "";
    $num_precision = length($num_precision);

    $i = $start_time;
    while ($i <= $end_time) {
	$curr_rst_file = "rstfiles/" . $mol_name . "_rst." . $i;

	print "Analyzing $curr_rst_file...";

	$out_cmd = "tail -1 $curr_rst_file |";
	open OUTCMD, $out_cmd or die "Cannot tail $curr_rst_file: $!\n";
	while (<OUTCMD>) {
	    $my_data = $_;
	    $my_data =~ /(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/;
	    $a_val = $1;
	    $b_val = $2;
	    $c_val = $3;
	    $d_val = $4;
	}
	close OUTCMD;

	$out_cmd = "> tmp_anal_res.in";

	if ($global_cuttoff > 0) {
	    $cutt_off = $global_cuttoff;
	} else {
	    $cutt_off = $a_val;
	    if ($b_val > $cutt_off) {
		$cutt_off = $b_val;
	    }
	    
	    if ($c_val > $cutt_off) {
		$cutt_off = $c_val;
	    }
	    
	    if ($d_val > $cutt_off) {
		$cutt_off = $d_val;
	    }
	}

	$cutt_off = sprintf("%1.1f", $cutt_off);

	print "Running anal using $cutt_off cutoff...";


	open OUTCMD, $out_cmd or die "Cannot create tmp_anal_res.in: $!\n";

	print OUTCMD "TITLE \'Residue Energy\'\n";
	print OUTCMD "\t1 0 0 0 " . ($total_bases + 1) . " 1\n";
	print OUTCMD "\t0 $a_val $b_val $c_val $d_val\n";
	print OUTCMD "\t1 0 1 0 50 1\n";
	print OUTCMD "\t$cutt_off 2.0 1.2 1.0\n";
	print OUTCMD "\t0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
	print OUTCMD "ENERGY\n";
	for $j (1 .. $total_bases) {
	    print OUTCMD "Title ResidueC\n";
	    print OUTCMD "RES $j\n";
	    print OUTCMD "STOP\n";
	    print OUTCMD "END\n";
	}
	print OUTCMD "END\nEND\nSTOP\n";
	close OUTCMD;

	$curr_eng_out_file = "outfiles/eng_" . $i . ".out";
	$out_cmd = "anal -O -i tmp_anal_res.in -o $curr_eng_out_file -p $topfile -c $curr_rst_file\n";
	system $out_cmd;
	system "rm -f tmp_anal_res.in";

	$hash_key = sprintf("%0" . $num_precision . "d", $i);
	ProcessFile($curr_eng_out_file, $hash_key, 0);
	$i += $interval;
    }
}

sub ProcessFile(@) {
    my ($curr_fle, $hash_key, $which_structure) = @_;
    my ($data_start, $in_data, $unit_nm);

    my (%structure_data);

    $data_start = 0;

    if (open INPFILE, $curr_fle) {
	while (<INPFILE>) {
	    chomp;

	    if ($_ =~ /ENERGY CONTRIBUTION BY GROUPS/) {
		$data_start = 1;
#		print "Storing data..";
	    }

	    if ($data_start) {
		$in_data = $_;
		if ($in_data = /^\s+(\d+)\s+(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)/) {
		    if ($1 == ($total_bases +1)) {
			$data_start = 0;
			last;
			last;
		    } else {
			$unit_nm = sprintf("%4d", $1);
			$EngDatStructure{$unit_nm}->{$hash_key} = (
							      {
								  "BOND" => $2,
								  "ANGLE" => $3,
								  "DIHEDRAL" => $4,
								  "VDW14" => $5,
								  "EEL14" => $6,
								  "VDWNB" => $7,
								  "EELNB" => $8,
								  "HBOND" => $9,
								  "CONSTRAINT" => $10
								  }
							      );
		    }
		}
	    }
	    
	} 
	print "Done\n";
	close INPFILE;
    }else {
	print "\nWARNING: Invalid file: $curr_fle\n";
    }
}
 
sub WriteOutput() {

    my ($unit_key, $eng_key, $in_energy);
    my ($avg_total_energy, $component_total_energy, $component_counter);
    my ($out_line, $avg_structure_out);

    CalculateTotals();

    if (! -d "avg_structures") {
	mkdir "avg_structures";
    }
    open AVGSTRU, "> avg_structures/Total_all_Energies.dat" or die "$!\n"; 
    $avg_total_energy = 0.0;
	print "     GROUP   ANGLE      BOND  CONSTRAINT DIHEDRAL    EEL14     EELNB     HBOND     VDW14     VDWNB      TOTAL\n";
    print AVGSTRU "\tENERGY CONTRIBUTION BY GROUPS\n\n";
	print AVGSTRU"     GROUP   ANGLE      BOND  CONSTRAINT DIHEDRAL    EEL14     EELNB     HBOND     VDW14     VDWNB      TOTAL\n";

    $avg_structure_out = "";
    for $unit_key ( sort numerically keys %EngDatStructure) {
	$component_total_energy = 0.0;
	$component_counter = 0;
	$out_line = sprintf("%8s", $unit_key);

	for $eng_key (sort numerically keys %{ $EngDatStructure{$unit_key}->{"Statistics"} } ) {
	    $in_energy = $EngDatStructure{$unit_key}->{"Statistics"}->{$eng_key}->{"Average"};
	    $component_total_energy += $in_energy;
	    $out_line .= sprintf("%10.2f", $in_energy);
	    $component_counter++;
	}
	$out_line .= sprintf("%12.3f", $component_total_energy);
	$avg_structure_out .= sprintf("%5d%12.3f\n", $unit_key, $component_total_energy);
	print "$out_line\n";
	print AVGSTRU "$out_line\n";
	$avg_total_energy += ($component_total_energy/$component_counter);
	$component_total_energy = 0.0;
	$component_counter = 0;
    }

    close AVGSTRU;

    if (open OUTFILE, "> avg_structures/Total_" . $start_time . "ps_" . $end_time . "ps.dat") {
	print OUTFILE "$avg_structure_out";
	close OUTFILE;
    }

    print "\nAverage Total Energy: $avg_total_energy KCal\n";

}

sub CalculateTotals() {

    my ($unit_key, $time_key, $eng_key, $unit_val, $time_val, $eng_val, $counter);
    my ($avg, $StDev, $total);
    $avg = $StDev = $total = 0.0;

    while ( ($unit_key, $unit_val) = each %EngDatStructure) {
	$counter = 0;
	while ( ($time_key, $time_val) = each %{ $unit_val } ) {
	    if ($time_key ne "Statistics") {
		while ( ($eng_key, $eng_val) = each %{ $time_val } ) {
		    $unit_val->{"Statistics"}->{$eng_key}->{"DataVals"} .= $eng_val;
		    $unit_val->{"Statistics"}->{$eng_key}->{"DataVals"} .= " ";
		    $counter++;
		}
	    }
	}
    }

    while ( ($unit_key, $unit_val) = each %EngDatStructure) {
	while ( ($eng_key, $eng_val) = each %{ $unit_val->{"Statistics"} }) {
	    ($avg, $StDev, $total) = CollectStats($eng_val->{"DataVals"});
	    $eng_val->{"Average"} = $avg;
	    $eng_val->{"StDev"} = $StDev;
	    $eng_val->{"Total"} = $total;
	}
    }

}

sub CollectStats(@) {

    my (@datavalues, $n_total, $avg, $result, $i);

    @datavalues = split / /, $_[0];
    $n_total = $#datavalues;
    $avg = 0.0;
    $result = 0.0;

    foreach $i (@datavalues) {
        $avg += $i;
    }

    $avg = $avg/($n_total + 1);

    foreach (@datavalues) {
        $result += ($_ - $avg) **2;
    }

    if ($n_total ==0) {
        $n_total = 1;
    }


    $result = sqrt($result/$n_total);
    return ($avg, $result, ($avg * $n_total));

}

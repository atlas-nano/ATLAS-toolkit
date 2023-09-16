#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename;
use General;
use GetParms;
use HelixLayout;
use FileFormats qw(GetBGFFileInfo);

# This program will perform energy analysis on a md trajectory file using anal
# It attempts to extract the components of the energy matrix (vdw, electrostatics,
# etc.) and create time dependant graphs additionally, if the energy_files for 
# the individual helices are specified, it will provide the appropritate strain 
# energy plots.
#
# usage: energy_anal_nonamot.pl param_file start_time end_time interval 
#        [helix1.out helix2.out]

sub ValidateInput();
sub GetParameters();
sub CreateRstFiles();
sub CreatePDBFiles();
sub DoEngAnal();
sub DoMPSimEngAnal();
sub WriteOutput();
sub ProcessFile(@);
sub Process_File_New(@);
sub CalculateTotals();
sub CreateGraphs();
sub DetermineHelixLayout();
sub ReCalculateEnergy();
sub GetUnitMap();
sub FilterAtom(@);
sub GetReferenceEnergy(@);
sub GenConnects();
sub RunCmd(@);

my ($paramfile, $start_time, $end_time, $interval, $reference, $filter_atm) = @ARGV;
my (%helix_info, %EngDatStructure, %Unit_Holder, $Ref_Energy, $printStr);
my (@eng_header, @unit_array, @unit_array_hash, @crossovers, @time_array, @helix);
my ($calc_strain, $spec_pic, $has_valid_files, $P_File);
my ($Ptraj_cmd) = "/ul/maiti/ptraj-6.3/linux/ptraj ";
my ($ambPDB_cmd) = "/exec/amber8/exe/ambpdb ";
my ($anal_cmd) = "/exec/amber8/exe/anal ";
my ($mpsim_cmd) = "/ul/maiti/src/non_pme/build/linux/mpsim.20030918 ";
my ($ATOMS, $BONDS);

$|++;

if (!@ARGV or $#ARGV < 3) {
    die "usage: $0 param_file start_time end_time ",
    "interval [reference file] [filter atom]\n";
}

$has_valid_files = 0;
# -- Start --
{
    print "\n-==Start Program ==-\n\n";
    open DUMPFILE, "> dumpfile.txt" or die "Cannot create dumpfile.txt: $!\n";

    $P_File = GetParms->new();
    if (! $P_File->IsValidParams($paramfile)) {
	die "Error in Paramater file\n";
    }

    DetermineHelixLayout();
    ValidateInput();
    print "\n\nStep 1. Executing Ptraj...";
    if ($P_File->{"Amber_Options"}->{"run_ptraj"}) {
	#CreatePDBFiles();
	CreateRstFiles();
	print "Done\n";
    } else {
	print "Assuming pdb files already exist..Done\n";
    }
   
    GetUnitMap();

    $printStr = "Step 3. Analyzing files...";
    DoMPSimEngAnal();

    @unit_array = sort { $a <=>$b } keys %EngDatStructure;
    if (! $has_valid_files) {
	close DUMPFILE;
	die "Cannot find any valid files. Check the dumpfile for more info\n";
    }
  
    if ($calc_strain) {
	($Ref_Energy, $calc_strain) = GetReferenceEnergy($reference);
    }

    if ($calc_strain) {
	$spec_pic = "/ul/tpascal/scripts/spectrum.png";
	if (! -d "total_energy") {
	    mkdir "total_energy", 0777;
	}
	chdir "total_energy";
    }	

    print "Step 4: Calculating Total Energy...";
    WriteOutput();

    if ($calc_strain) {
	$spec_pic = "/ul/tpascal/scripts/strain_spectrum.png";
	chdir "../";
	if (! -d "strain_energy") {
	    mkdir "strain_energy", 0777;
	}
	chdir "strain_energy";
	ReCalculateEnergy();
	print "Step 5: Calculating Strain Energy...";
	WriteOutput();
    }
    
    print "\nAll Tasks Completed\n\n-==End Program==-\n\n";
    close DUMPFILE;
}

sub ValidateInput() {
    
    if ($reference) {
	$calc_strain = 1;
	if (! -e $reference) { 
	    print "Warning: Cannot locate $reference: $!\n";
	    $calc_strain = 0;
	}
    }
    
    if (! $calc_strain) {
	print "NOTE: Will not perform strain energy calculations\n";
	$spec_pic = "/ul/tpascal/scripts/spectrum.png";
    } else {
	print "NOTE: Performing strain energy calculations using reference data: $reference\n";
	$spec_pic = "/ul/tpascal/scripts/strain_spectrum.png";
    }
    
    die "Invalid starting number. Expected integer\n"
	if (! IsInteger($start_time) );
    
    die "Invalid end number. Expected integer\n"
	if (! IsInteger($end_time) );

    die "Invalid interval. Expected integer\n"
	if (! IsInteger($interval) );
    
    if ($start_time > $end_time) {
	($start_time, $end_time) = ($end_time, $start_time);
    }
}


sub CreateRstFiles() {
    my ($out_cmd, $i, $j, $strand_ln);

    mkdir "rstfiles", 0777 if (! -d "rstfiles");

    $strand_ln = int($P_File->{"Molecule"}->{"total_bases"}/4);

# Create restart & pdb files----------------------
    $out_cmd = "> tmp_gen_rst_file";

    $i = $j = 1;

    open OUTCMD, $out_cmd or die "Cannot create file tmp_gen_rst_file: $!\n";
    print OUTCMD "trajin " .
	 $P_File->{"Files"}->{"trajectory"} . " $start_time $end_time $interval\n";

    while ($i <  $P_File->{"Molecule"}->{"total_bases"}) {
	print OUTCMD "center :1-" . ($j * $strand_ln) . " mass origin\n";
	print OUTCMD "image origin center\n";
	$i +=$strand_ln;
	$j++;
    }

    print OUTCMD "trajout rstfiles/" . $P_File->{"Molecule"}->{"name"} . "_rst restart\n";
    print OUTCMD "\n";
    close OUTCMD;

    $out_cmd = $Ptraj_cmd . " " . $P_File->{"Files"}->{"topology"}; 
    $out_cmd .= " < tmp_gen_rst_file >& junk";
    RunCmd($out_cmd, "");

    system "rm -f junk tmp_gen_rst_file";
}

sub DoEngAnal() {

    my ($i, $curr_rst_file, $num_precision, $out_cmd, $j, $curr_eng_out_file);
    my ($my_data, @vals, $cutt_off, $hash_key);

    $num_precision = $end_time . "";
    $num_precision = length($num_precision);

    system "mkdir -p outfiles";
    $i = $start_time;
    $cutt_off = 0;

    if (! $P_File->{"Amber_Options"}->{"run_anal"}) {
	print "Assuming the anal files already exist...";
    }
    while ($i <= $end_time) {
	$curr_eng_out_file = "outfiles/eng_" . $i . ".out";
	if ($P_File->{"Amber_Options"}->{"run_anal"}) {
	    $curr_rst_file = "rstfiles/" . $P_File->{"Molecule"}->{"name"} . "_rst." . $i;
	    print DUMPFILE "Analyzing $curr_rst_file...";
	    $out_cmd = "tail -1 $curr_rst_file |";
	    
	    open OUTCMD, $out_cmd or die "Cannot tail $curr_rst_file: $!\n";
	    while (<OUTCMD>) {
		@vals = /(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/;
	    }
	    close OUTCMD;
	    
	    if ($P_File->{"Amber_Options"}->{"cutoff"} > 0) {
		$cutt_off = $P_File->{"Amber_Options"}->{"cutoff"};
	    } else {
		for (@vals) {
		    $cutt_off = $_ if ($_ > $cutt_off);
		}
	    }
	    
	    $cutt_off = sprintf("%4.1f", $cutt_off);
	    
	    print DUMPFILE "Running anal using $cutt_off cutoff...";
	    
	    $out_cmd = "> tmp_anal_res.in";
	    open OUTCMD, $out_cmd or die "Cannot create tmp_anal_res.in: $!\n";
	    
	    print OUTCMD "TITLE \'Residue Energy\'\n";
	    print OUTCMD "\t1 0 0 0 " . ($P_File->{"Molecule"}->{"total_bases"} + 1) . " 1\n";
	    print OUTCMD "\t0 ";
	    for (@vals) {
		print OUTCMD "$_ ";
	    }
	    print OUTCMD "\n";
	    print OUTCMD "\t1 0 1 0 50 1\n";
	    print OUTCMD "\t$cutt_off 2.0 1.2 1.0\n";
	    print OUTCMD "\t0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
	    print OUTCMD "ENERGY\n";
	    for (1 .. $P_File->{"Molecule"}->{"total_bases"}) {
		print OUTCMD "Title ResidueC\n";
		print OUTCMD "RES $_\n";
		print OUTCMD "STOP\n";
		print OUTCMD "END\n";
	    }
	    print OUTCMD "END\nEND\nSTOP\n";
	    close OUTCMD;
	    
	    $out_cmd = "$anal_cmd -O -i tmp_anal_res.in -o $curr_eng_out_file " . 
		"-p " . $P_File->{"Files"}->{"topology"} . " -c $curr_rst_file\n";
	    RunCmd($out_cmd);
	    system "rm -f tmp_anal_res.in";
	}
 
	push @time_array, $i;
	$hash_key = $i;
	ProcessFile($curr_eng_out_file, $hash_key, 0);
	$i += $interval;
    }
}

sub ProcessFile(@) {
    my ($curr_fle, $hash_key, $which_structure) = @_;
    my ($data_start, $in_data, $unit_nm, $counter, $eng);
    my ($curr_val, $curr_total);

    $data_start = 0;
    if (open INPFILE, $curr_fle) {
	while (<INPFILE>) {
	    chomp;
	    if ($_ =~ /ENERGY CONTRIBUTION BY GROUPS/) {
		$data_start = 1;
	    }elsif ($data_start) {
		$in_data = $_;
		if ($in_data =~ /^\s+GROUP/ && $#eng_header <= 0) {
		    while ($in_data =~ /\s+(\w+)/g) {
			if ($1 ne "GROUP" && $1 ne "TOTAL") {
			    push @eng_header, $1;
#			    print DUMPFILE "GOT TOO GROUP HEADINGS: $1\n";
			}
		    }
		} elsif ($in_data =~ /^\s+(\d+)\s+/ && $#eng_header > 0) {
		    $unit_nm = $1;
#		    print DUMPFILE "unit: $unit_nm, time: $hash_key\n";
		    if ($unit_nm == ($P_File->{"Molecule"}->{"total_bases"} +1)) {
			$data_start = 0;
			last;
			last;
		    } elsif ($which_structure && 
			     $unit_nm == ($P_File->{"Molecule"}->{"total_bases"}/2 +1) ) {
			$data_start = 0;
			last;
			last;
		    } else {
			if (! $which_structure) {
			    $counter = 0;
			    $curr_val = $curr_total = 0;
			    while ($in_data =~ /(\-?\d+\.\d+)\s*/g && $counter <= $#eng_header) {
				$eng = $eng_header[$counter];
				$curr_val = $1;
				$curr_total += $curr_val;
#				print "counter: $counter eng: $eng\n";
				$EngDatStructure{$unit_nm}->{$hash_key}->{$eng} = $curr_val; 
				$counter++;
			    }
			    $EngDatStructure{$unit_nm}->{$hash_key}->{"TOTAL"} = $curr_total;
			} else {
			    $unit_nm = GetCorrectUnit($unit_nm, $which_structure);
			    $counter = 0;
			    while ($in_data =~ /(\-?\d+\.\d+)\s*/g) {
				$eng = $eng_header[$counter];
				$helix_info{$unit_nm}->{"0"}->{$eng} += $1;
				$counter++;
			    }
			}
		    }
		}
	    }
	    
	} 
	if ($#eng_header > 0) {
	    $has_valid_files = 1;
	    print DUMPFILE "Done\n";
	}else {
	    print DUMPFILE "Failure. $curr_fle is invalid\n";
	}
	close INPFILE;
    }else {
	print DUMPFILE "\nWARNING: Invalid file: $curr_fle\n";
    }
}
 
sub DetermineHelixLayout() {
    my $hl = HelixLayout->spawn();
    $hl->DetermineHelixLayout(
			      $P_File->{"Molecule"}->{"major_groove"}, 
			      $P_File->{"Molecule"}->{"minor_groove"}, 
			      $P_File->{"Molecule"}->{"is3PrimeIn"}, 
			      $P_File->{"Molecule"}->{"bases_at_end"}, 
			      ($P_File->{"Molecule"}->{"total_bases"}/4), 
			      @{ $P_File->{"Molecule"}->{"crossovers"} }
			      );

    @helix = $hl->GetHelixInfo();
}

sub WriteOutput() {

    my ($unit_key, $eng_key, $in_energy, $time_avg, $time_stdev);
    my ($avg_total_energy, $component_total_energy, $component_counter);
    my ($out_line, $avg_structure_out, $outfile, $time_total);

    CalculateTotals();

    $avg_total_energy = 0.0;

    mkdir "avg_structures", 0777 if (! -d "avg_structures");

    $out_line = sprintf("%8s", "GROUP");

    for (@eng_header) {
	$out_line .= sprintf("%10s", $_);
    }

    $out_line .= sprintf("%10s\n", "Total");
    $component_counter = 0;
    $avg_structure_out = "";
    for $unit_key ( @unit_array ) {
	$component_total_energy = 0.0;
	$out_line .= sprintf("%8s", $unit_key);
	for $eng_key (@eng_header) {
	    $in_energy = $EngDatStructure{$unit_key}->{"Statistics"}->{$eng_key}->{"Average"};
	    $component_total_energy += $in_energy;
	    $out_line .= sprintf("%10.2f", $in_energy);
	    $component_counter++;
	}
	$out_line .= sprintf("%12.3f", $component_total_energy) . "\n";
	$time_avg = $EngDatStructure{$unit_key}->{"Statistics"}->{"TIMEINFO"}->{"Average"};
	$time_stdev = $EngDatStructure{$unit_key}->{"Statistics"}->{"TIMEINFO"}->{"StDev"};
	$time_total = $EngDatStructure{$unit_key}->{"Statistics"}->{"TIMEINFO"}->{"Total"};
	$avg_structure_out .= sprintf("%5d%12.3f%10.3f%5d\n", 
				      $unit_key, $time_avg, $time_stdev, $time_total);
	$avg_total_energy += $time_avg;
	$component_total_energy = 0.0;
    }

    $outfile = "avg_structures/Total_all_Energies.dat";
    open AVGSTRU, "> $outfile" or die "Cannot create $outfile: $!\n"; 
    print AVGSTRU $out_line;
    close AVGSTRU;

    $outfile = "avg_structures/Total_" . $start_time . "ps_" . $end_time . "ps.dat";
    open OUTFILE, "> $outfile" or die "Cannot create $outfile: $!\n";
    print OUTFILE "$avg_structure_out";
    close OUTFILE;

    $avg_total_energy = sprintf("%12.3f", $avg_total_energy);
    print "Done\n\tAverage Total Energy: $avg_total_energy KCal\n";

    print "\tCreating graphs...";
    CreateGraphs();
    print "Done\n";

}

sub CreateGraphs() {

    my (@eng_comp_array, $this_unit, $total_energies);
    my ($i, $j, $outfile, $eng_total, $eng_data, $running_total);
    my ($graph_2d_out, $graph_3d_out, $my_total, $my_stdev, $my_num);

#   populate the arrays with the keys for the time and energy contributions


    mkdir "2d_graphs", 0777 if (! -d "2d_graphs");
    mkdir "3d_graphs", 0777 if (! -d "3d_graphs");

#    Create Per_Energy Graphs
    for $i (@eng_header) {
	$graph_3d_out = "";
	$graph_2d_out = "";
	for $j (@time_array) {
	    $eng_total = $eng_data = "";
	    for $this_unit ( @unit_array ) {
		$eng_data .= sprintf("%5d%5d%10.2f\n", $this_unit, $j, 
				     $EngDatStructure{$this_unit}->{$j}->{$i});
		$eng_total .= "$EngDatStructure{$this_unit}->{$j}->{$i} ";
	    }
	    ($my_total, $my_stdev, $my_num) = STDev($eng_total);
	    $graph_2d_out .= sprintf("%5d%10.2f %6.2f\n", $j,$my_total, $my_stdev);
	    $graph_3d_out .= $eng_data . "\n";
	}

	$outfile = "2d_graphs/" . $i . "_" . $start_time . "ps_" . $end_time . "ps.dat";
	open OUTFILE, "> $outfile" || die "Cannot create $outfile: $!\n";
	print OUTFILE "$graph_2d_out";
	close OUTFILE;
	
	
	$outfile = "3d_graphs/" . $i . "_" . $start_time . "ps_" . $end_time . "ps.dat";
	open OUTFILE, "> $outfile" || die "Cannot create $outfile: $!\n";
	print OUTFILE "$graph_3d_out";
	close OUTFILE;
    }

    $graph_2d_out = $graph_3d_out = "";

#    Create Total Energy per unit graphs
    for $j (@time_array) {
	$running_total = $eng_data = "";
	for $this_unit ( @unit_array ) {
	    $eng_total = 0.0;
	    for $i (@eng_header) {
		$eng_total += $EngDatStructure{$this_unit}->{$j}->{$i};
	    }
	    $eng_data .= sprintf("%5d%5d%10.2f\n", $this_unit, $j, $eng_total);
	    $running_total .= "$eng_total ";
	}
	chop $running_total;
	($my_total, $my_stdev, $my_num) = STDev($running_total);
	$graph_2d_out .= sprintf("%5d%10.2f\n", $j,$my_num);
	$total_energies .= "$my_num ";
	$eng_data .= sprintf("%-5s%5d%10.2f%6.2f\n", "#Avg", $j, $my_total, $my_stdev);
	$eng_data .= sprintf("%-5s%5d%10.2f\n", "#Tot", $j, $my_num);
	$graph_3d_out .= ($eng_data . "\n");
    }

    chop $total_energies;
    ($my_total, $my_stdev, $my_num) = STDev($total_energies);
    $graph_2d_out .= sprintf("%5s%10.2f %6.2f\n", "Avg", $my_total, $my_stdev);

    $outfile = "2d_graphs/TOTAL_" . $start_time . "ps_" . $end_time . "ps.dat";
    open OUTFILE, "> $outfile" || die "Cannot create $outfile: $!\n";
    print OUTFILE "$graph_2d_out";
    close OUTFILE;
        
    $outfile = "3d_graphs/TOTAL_" . $start_time . "ps_" . $end_time . "ps.dat";
    open OUTFILE, "> $outfile" || die "Cannot create $outfile: $!\n";
    print OUTFILE "$graph_3d_out";
    close OUTFILE;
    
#    create average structures
    for $j (@eng_header) {
	$eng_data = "";
	for $this_unit ( @unit_array ) {
	    $eng_data .= sprintf("%5d", $this_unit);
	    $eng_data .= sprintf("%10.2f", 
				 $EngDatStructure{$this_unit}->{"Statistics"}->{$j}->{"Average"});
	    $eng_data .= sprintf("%10.2f", 
				 $EngDatStructure{$this_unit}->{"Statistics"}->{$j}->{"StDev"});
	    $eng_data .="\n";
	}
	$outfile = "avg_structures/" . $j . "_" . $start_time . "ps_" . $end_time . "ps.dat";
	open OUTFILE, "> $outfile" || die "Cannot create $outfile: $!\n";
	print OUTFILE "$eng_data";
	close OUTFILE;
    }
}

sub CalculateTotals() {

    my ($unit_key, $time_key, $eng_key, $unit_val, $time_val, $eng_val, $counter);
    my ($avg, $StDev, $total);
    $avg = $StDev = $total = "";

# Reset the Statistics information
    for $unit_key (@unit_array) {
	if ($EngDatStructure{$unit_key}->{"Statistics"}) {
	    delete $EngDatStructure{$unit_key}->{"Statistics"};
	}
    }

# Calc the Energy per unit per Contributor
    for $time_key (@time_array) {
	for $unit_key (@unit_array) {
	    $total = 0;
	    for $eng_key (@eng_header) {
		$eng_val = $EngDatStructure{$unit_key}->{$time_key}->{$eng_key};
		$EngDatStructure{$unit_key}->{"Statistics"}->{$eng_key}->{"DataVals"} .= "$eng_val ";
		$total += $eng_val;
	    }
	    $EngDatStructure{$unit_key}->{"Statistics"}->{"TIMEINFO"}->{"DataVals"} .= "$total ";
	    $total = 0;
	}
    }
		
    $total = 0;
# Store the Energy per unit per Contributor
    for $unit_key (@unit_array) {
	while ( ($eng_key, $eng_val) = each %{ $EngDatStructure{$unit_key}->{"Statistics"} }) {
	    chop $eng_val->{"DataVals"};
	    ($avg, $StDev, $total) = STDev($eng_val->{"DataVals"});
	    $eng_val->{"Average"} = $avg;
	    $eng_val->{"StDev"} = $StDev;
	    if ($avg < 0 || $avg > 0) {
		$eng_val->{"Total"} = ($total/$avg);
	    } else {
		$eng_val->{"Total"} = 0;
	    }
	    
	}
    }

}

sub ReCalculateEnergy() {
    my ($this_unit, $i, $j);
    
#   populate the arrays with the keys for the time and energy contributions
    
    for $this_unit (@unit_array) {
	for $i (@time_array) {
	    for $j (@eng_header) {
		$EngDatStructure{$this_unit}->{$i}->{$j} -= $Ref_Energy->{$this_unit}->{$j};
	    }
	}
    }
}

sub GetCorrectUnit(@) {
    my ($curr_unit, $helix_no) = @_;
    my ($strand_len) = $P_File->{"Molecule"}->{"total_bases"}/4;
    my ($result, $curr_chain);

    $result = 0;
    $helix_no -= 1;
    if ($curr_unit <= $strand_len) {
	$curr_chain = 0;
    } else {
	$curr_chain = 1;
	$curr_unit -= $strand_len;
    }

    for my $regions (@{ $helix[$helix_no][$curr_chain] }) {
	if ($curr_unit >= $regions->{"StartUnit"} && $curr_unit <= $regions->{"EndUnit"}) {
	    $result = ( ($regions->{"Strand"} -1) * $strand_len ) + $curr_unit;
	    last;
	}
    }

    return $result;

}

sub DoMPSimEngAnal() {
# This sub will execute MPSim too do one-point energy calculations
    my ($i, $curr_pdb_file, $num_precision, $out_cmd, $j, $curr_eng_out_file, $strLen);
    my ($my_data, @vals, $cutt_off, $hash_key, $curr_rst_file, $isvalid, $tot);
    my ($curr_bgf_file, $outstring, $curr_project, $counter, $connect_file, $start);

    $num_precision = $end_time . "";
    $num_precision = length($num_precision);

    system "mkdir -p outfiles";
    system "mkdir -p bgffiles";
    $i = $start_time;
    $cutt_off = $counter = $start = 0;
    $tot = int(($end_time - $start_time)/$interval);

    if (! $P_File->{"Amber_Options"}->{"run_anal"}) {
	print "Assuming the mpsim output files already exist...";
    } else {
	#$connect_file = GenConnects();
    }
    
    print "$printStr" . "calculating time remaining\r";
    while ($i <= $end_time) {
	$isvalid = 1;
	$curr_eng_out_file = "outfiles/" . $P_File->{"Molecule"}->{"name"} . "_" . $i . ".out";
	if ($P_File->{"Amber_Options"}->{"run_anal"} && (! -e $curr_eng_out_file)) {
	    $curr_pdb_file = "pdbfiles/" . $P_File->{"Molecule"}->{"name"} . "_pdb." . $i;
	    $curr_rst_file = "rstfiles/" . $P_File->{"Molecule"}->{"name"} . "_rst." . $i;
	    $curr_bgf_file = "bgffiles/" . $P_File->{"Molecule"}->{"name"} . "_${i}.bgf";

	    print DUMPFILE "Analyzing $curr_bgf_file...";
	    
	    $out_cmd = "/ul/tpascal/scripts/amber2bgf.pl " . $P_File->{"Files"}{"topology"};
	    $out_cmd .= " $curr_rst_file $curr_bgf_file /ul/tpascal/ff/AMBER95.cnv >> junk";
	    RunCmd($out_cmd);
	    $isvalid = 1;

	    if ($isvalid) {
		$outstring = "";
		$curr_project = $P_File->{"Molecule"}->{"name"} . "_" . $i;
		system "cp /ul/tpascal/mpsim/ewald.par .";
		open CTL, "/ul/tpascal/mpsim/1-pme.ctl" || die "Cannot open 1-pme.ctl: $!\n";
		while (<CTL>) {
		    chomp;
		    if ($_ =~ /PROJECT/) {
			$outstring .= "PROJECT               $curr_project";
		    } elsif ($_ =~ /STRUCTURE/) {
			$outstring .= "STRUCTURE             $curr_bgf_file";
		    } else {
			$outstring .= $_;
		    }
		    $outstring .= "\n";
		}
		close CTL;
		open CTL, "> 1-pme.ctl" || die "Cannot modify control file: $!\n";
		print CTL $outstring;
		close CTL;
		$out_cmd = "$mpsim_cmd 1-pme.ctl >& mpsim.out";
		if (system($out_cmd)) {
		    print "Mpsim Error!\n";
		} else {
		    system "rm -f 1-pme.ctl ewald.par mpsim.out";
		    system "cp $curr_project" . ".fin.ener outfiles/$curr_project" . ".out";
		    system "rm -f $curr_project" . ".*";
		    $curr_project = "outfiles/" . $curr_project . ".out";
		    Process_File_New($curr_project, $i, 0);
		    $counter ++;
		    push @time_array, $i;
		}
	    }
	}else {
	    $curr_project = $P_File->{"Molecule"}->{"name"} . "_" . $i;
	    $curr_project = "outfiles/" . $curr_project . ".out";
	    if (-e $curr_project) {
		Process_File_New($curr_project, $i);
		$counter ++;
		push @time_array, $i;
	    }
	}
	$i += $interval;
	$start = time() if (! $start);
	$strLen = PrintProgress($counter, $tot, $start, $printStr);
    }

    if ($counter ==0) {
	die "No valid energy files found!\n";
    }
    system "rm -fr _restart.1 _connects";

    printf "$printStr%-${strLen}s\n", "Done";
}

sub CreatePDBFiles() {
    my ($out_cmd, $i, $j, $strand_ln);

    mkdir "pdbfiles", 0777 if (! -d "pdbfiles");

    $strand_ln = int($P_File->{"Molecule"}->{"total_bases"}/4);

# Create restart & pdb files----------------------
    $out_cmd = "> tmp_gen_rst_file";

    $i = $j = 1;

    open OUTCMD, $out_cmd or die "Cannot create file tmp_gen_rst_file: $!\n";
    print OUTCMD "trajin " .
	 $P_File->{"Files"}->{"trajectory"} . " $start_time $end_time $interval\n";

    while ($i <  $P_File->{"Molecule"}->{"total_bases"}) {
	print OUTCMD "center :1-" . ($j * $strand_ln) . " mass origin\n";
	print OUTCMD "image origin center\n";
	$i +=$strand_ln;
	$j++;
    }

    print OUTCMD "trajout pdbfiles/" . $P_File->{"Molecule"}->{"name"} . "_pdb pdb dumpq\n";
    print OUTCMD "\n";
    close OUTCMD;

    $out_cmd = $Ptraj_cmd . " " . $P_File->{"Files"}->{"topology"}; 
    $out_cmd .= " < tmp_gen_rst_file >& junk";
    RunCmd($out_cmd);

    system "rm -f junk tmp_gen_rst_file";
}

sub Process_File_New(@) {
    my ($curr_fle, $hash_key) = @_;
    my ($data_start, $in_data, $unit_nm, $counter, $eng);
    my ($curr_val, $curr_total, $unit_label);

    $data_start = 0;
    if (open INPFILE, $curr_fle) {
	while (<INPFILE>) {
	    chomp;
	    $in_data = $_;
	    if ($in_data =~ /^\s*[Atom|GROUP]/ && $#eng_header <= 0) {
		while ($in_data =~ /\s+(\w+)/g) {
		    if ( ($1 !~ /GROUP/i) && ($1 !~ /TOTAL/i) ) {
			push @eng_header, $1;
#			    print "GOT TOO GROUP HEADINGS: $1\n";
		    }
		}
	    } elsif ($in_data =~ /^\s*(\d+)\s+/ && $#eng_header > 0) {
		$unit_nm = $ATOMS->{$1}{RESNUM};
		$unit_label = $ATOMS->{$1}{ATMNAME};
		    
#		print DUMPFILE "unit: $unit_nm, time: $hash_key\n";
		if ($unit_nm == ($P_File->{"Molecule"}->{"total_bases"} +1)) {
		    $data_start = 0;
		    last;
		} elsif (! FilterAtom($unit_label)) {
		    $counter = 0;
		    $curr_val = $curr_total = 0;
		    while ($in_data =~ /(\-?\d+\.\d+)\s*/g && $counter <= $#eng_header) {
			$eng = $eng_header[$counter];
			$curr_val = $1;
			$curr_total += $curr_val;
#				print "counter: $counter eng: $eng\n";
			if ($EngDatStructure{$unit_nm}->{$hash_key}->{$eng}) {
			    $EngDatStructure{$unit_nm}->{$hash_key}->{$eng} += $curr_val; 
			} else {
			    $EngDatStructure{$unit_nm}->{$hash_key}->{$eng} = $curr_val; 
			}				
			$counter++;
		    }
		    if ($EngDatStructure{$unit_nm}->{$hash_key}->{"TOTAL"}) {
			$EngDatStructure{$unit_nm}->{$hash_key}->{"TOTAL"} += $curr_total;
		    }else {
			$EngDatStructure{$unit_nm}->{$hash_key}->{"TOTAL"} = $curr_total;
		    }
		    $curr_val = $curr_total = 0;
		}
	    }
	}
	
	if ($#eng_header > 0) {
	    $has_valid_files = 1;
	    print DUMPFILE "Done\n";
	}else {
	    print DUMPFILE "Failure. $curr_fle is invalid\n";
	}
	close INPFILE;
    }else {
	print DUMPFILE "\nWARNING: Invalid file: $curr_fle\n";
    }
}

sub GetUnitMap() {
    my ($currRstFile) = "rstfiles/" . $P_File->{"Molecule"}->{"name"} . "_rst." . $start_time;
    my ($bgfFile) = $P_File->{Molecule}{name} . ".bgf";
    print "Step 2. Getting unit map...";
    my ($bgfCmd) = "/ul/tpascal/scripts/amber2bgf.pl $P_File->{Files}{topology} $currRstFile $bgfFile >& junk";

    die "ERROR: Cannot locate $currRstFile!\n" if (! -e $currRstFile);
    die "ERROR: Command terminated abruptly: $bgfCmd\n" if (system($bgfCmd));

    ($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
    print "Done\n";
}

sub FilterAtom(@) {
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
    if (! $filter_atm) {
	return 0;
    }else {
        return $return_val;
    }
}

sub GetReferenceEnergy(@) {
    my ($infile) = $_[0];
    my ($sequence_key, $is_valid);
    my (%Ref_Energy, $has_ref);

    if (open(INFILE, $infile) or die "Cannot open reference file: $infile\n") {
	while (<INFILE>) {
	    chomp;
	    if ($_ =~ /^(\w)(\w)\/(\w)(\w)\s+(\-?\d+\.\d+)/) {
		$sequence_key = ord($1) . ord($2) . ord($3) . ord($4);
		$Ref_Energy{$sequence_key}{"TOTAL"} = $5;
		$has_ref = 1;
	    }elsif ($_ =~ /^(\d).\/(\w)(\w)\s+(\-?\d+\.\d+)/) {
		$sequence_key = $1 . ord($2) . ord($3);
		$Ref_Energy{$sequence_key}{"TOTAL"} = $4;
#		print "$_ $sequence_key $4\n";
		$has_ref = 1;
	    }		
	}
	close INFILE;
    }
    return (\%Ref_Energy, $has_ref);
}

sub GenConnects() {
# This sub will generate a restart file from ptraj, then use it to generate a connectivity file from ambpdb

# create ptraj_in file
    my ($out_cmd);

    open OUTFILE, "> ptraj_in" or die "Cannot create temporary file ptraj_in: $!\n";
    print OUTFILE "trajin " . $P_File->{"Files"}{"trajectory"} . " 1 1 1\n";
    print OUTFILE "trajout _restart restart\n";
    close OUTFILE;
    
    $out_cmd = $Ptraj_cmd . " " . $P_File->{"Files"}->{"topology"} . " < ptraj_in >& junk";
    RunCmd($out_cmd, "_restart.1");
    
    $out_cmd = $ambPDB_cmd . " -p " . $P_File->{"Files"}->{"topology"} . " -bnd < _restart.1 >& _connects";
    RunCmd($out_cmd, "_connects");
    
    return "_connects";

}

sub RunCmd(@) {
    my ($cmd, $mfile) = @_;

    if (system("$cmd")) {
	die "Error while running cmd $cmd\n";
    }

    if (defined($mfile) and $mfile ne "" and ! -e $mfile) {
	die "Expected file $mfile was not created from command $cmd\n";
    }
}

#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}

use strict;
use File::Basename;
use Packages::General;
use Packages::GetParms;
use Packages::HelixLayout;
use Packages::Namot;
use Packages::MovieGen;
use Packages::ManipImages;

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
sub CreatePtrajFiles(@);
sub DoEngAnal();
sub WriteOutput();
sub ProcessFile(@);
sub CalculateTotals();
sub CreateGraphs();
sub DetermineHelixLayout();
sub ReCalculateEnergy();
sub CreateNamot2Files();
sub MovieGen(@);

my ($paramfile, $start_time, $end_time, $interval, $helix1, $helix2) = @ARGV;
my (%helix_info, %EngDatStructure);
my (@eng_header, @unit_array, @unit_array_hash, @crossovers, @time_array, @helix);
my ($calc_strain, $spec_pic, $has_valid_files, $P_File);

if (!@ARGV or $#ARGV < 3) {
    die "usage: $0 param_file start_time end_time ",
    "interval [helix1.out helix2.out]\n";
}

$has_valid_files = 0;
# -- Start --
{
    print "\n-==Start Program ==-\n\n";
    open DUMPFILE, "> dumpfile.txt" or die "Cannot create dumpfile.txt: $!\n";

    $P_File = Packages::GetParms->new();
    if (! $P_File->IsValidParams($paramfile)) {
	die "Error in Paramater file\n";
    }

    DetermineHelixLayout();
    ValidateInput();
    print "\n\nStep 1. Creating Restart Files...";
    if ($P_File->{"Amber_Options"}->{"run_ptraj"}) {
	CreatePtrajFiles("restart", "rstfiles", "rst");
	print "Done\n";
    } else {
	print "Assuming restart files already exist..Done\n";
    }
   
    print "Step 2. Analyzing files...";
    DoEngAnal();

    @unit_array = sort { $a <=>$b } keys %EngDatStructure;
    if ($has_valid_files) {
	print "Done\n";
    } else {
	close DUMPFILE;
	die "Cannot find any valid files. Check the dumpfile for more info\n";
    }
  
    if ($calc_strain) {
	$spec_pic = "/ul/tpascal/scripts/spectrum.png";
	if (! -d "total_energy") {
	    mkdir "total_energy", 0777;
	}
	chdir "total_energy";
    }	

    print "Step 3: Calculating Total Energy...";
    WriteOutput();

    if ($calc_strain) {
	$spec_pic = "/ul/tpascal/scripts/strain_spectrum.png";
	chdir "../";
	if (! -d "strain_energy") {
	    mkdir "strain_energy", 0777;
	}
	chdir "strain_energy";
	DetermineStrandLayout();
	ReCalculateEnergy();
	print "Step 4: Calculating Strain Energy...";
	WriteOutput();
    }
    
    print "\nAll Tasks Completed\n\n-==End Program==-\n\n";
    close DUMPFILE;
}

sub ValidateInput() {
    
    if ($helix1 && $helix2) {
	$calc_strain = 1;
	if (! -e $helix1) { 
	    print "Warning: Cannot locate $helix1: $!\n";
	    $calc_strain = 0;
	}
	
	if (! -e $helix2) {
	    print "Warning: Cannot locate $helix2: $!\n";
	    $calc_strain = 0;
	}
    }
    
    if (! $calc_strain) {
	print "NOTE: Will not perform strain energy calculations\n";
	$spec_pic = "/ul/tpascal/scripts/spectrum.png";
    } else {
	print "NOTE: Performing strain energy calculations using helix1 data: ",
	"$helix1 and helix2 data: $helix2\n";
	$spec_pic = "/ul/tpascal/scripts/strain_spectrum.png";
    }
    
    die "Invalid starting number. Expected integer\n"
	if (! IsInteger($start_time) );
    
    die "Invalid end number. Expected integer\n"
	if (! IsInteger($end_time) );

    die "Invalid interval. Expected integer\n"
	if (! IsInteger($interval) );
    
    if ($start_time > $end_time) {
	($start_time, $end_time) = Swap($start_time, $end_time);
    }
}


sub CreatePtrajFiles(@) {
    my ($file_type, $file_dir, $file_ext) = @_;
    my ($out_cmd, $i, $j, $strand_ln);

    system "mkdir -p $file_dir";

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

    print OUTCMD "trajout $file_dir/" . $P_File->{"Molecule"}->{"name"} . "_$file_ext $file_type\n";
    print OUTCMD "\n";
    close OUTCMD;

    $out_cmd = "ptraj " . $P_File->{"Files"}->{"topology"} . " < tmp_gen_rst_file >& junk";
    
    system $out_cmd;
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
	    
	    $out_cmd = "/ul/maiti/amber7/src/anal/anal -O -i tmp_anal_res.in -o $curr_eng_out_file " . 
		"-p " . $P_File->{"Files"}->{"topology"} . " -c $curr_rst_file\n";
	    system $out_cmd;
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
				$helix_info{$unit_nm}->{"0"}->{$eng} = $1;
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
    my $hl = Packages::HelixLayout->spawn();
    $hl->DetermineHelixLayout(
			      $P_File->{"Molecule"}->{"major_groove"}, 
			      $P_File->{"Molecule"}->{"minor_groove"}, 
			      $P_File->{"Molecule"}->{"is3PrimeIn"}, 
			      $P_File->{"Molecule"}->{"bases_at_end"}, 
			      $P_File->{"Molecule"}->{"total_bases"}, 
			      $P_File->{"Molecule"}->{"crossovers"}
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
    my (@helix_keys, $this_unit);
    my ($i, $j, $outfile);

    ProcessFile($helix1, "0", 1);
    ProcessFile($helix2, "0", 2);

#   populate the arrays with the keys for the time and energy contributions
    
    @helix_keys = keys %helix_info;

    if ($#helix_keys == $#unit_array) {
	
	for $i (@time_array) {
	    for $this_unit (@unit_array) {
		for $j (@eng_header) {
		    $EngDatStructure{$this_unit}->{$i}->{$j} -= $helix_info{$this_unit}->{"0"}->{$j};
		}
	    }
	}
    } else {
	die "The files $helix1 or $helix2 contains invalid data\n";
    }
}

sub GetCorrectUnit(@) {
    my ($curr_unit, $helix_no) = @_;
    my ($strand_len) = $P_File->{"Molecules"}->{"total_bases"}/4;
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
	    $result = ( ($regions->{"Chain"} - 1) * $strand_len ) + $curr_unit;
	    last;
	}
    }

    return $result;

}

sub CreateNamot2Files() {
# Assign an rgb color too each base depending on the strain energy of that base
    my ($i, $j, $eng_max, $eng_min, $temp, $pic_file, $counter);
    my ($curr_pdb_file, @namot_string, $text_length, $jpeg_file, @valid_pics);

    $eng_max = 0;
    $eng_min = -99999;

    system "mkdir -p pics";

    $counter = $time_array[0];
    for $i (@time_array) {
	for $j (@unit_array) {
	    $temp = $EngDatStructure{$j}->{$i}->{"TOTAL"};
	    $eng_min = $temp if $temp > $eng_min;
	    $eng_max = $temp if $temp < $eng_max;
	}
    }
    
    for $i (@time_array) {
	@namot_string = [];
	$curr_pdb_file = "pdbfiles/" . $P_File->{"Molecule"}->{"name"} . 
	    "_pdb." . $i;
	system "/ul/tpascal/scripts/stripNaH20.pl $curr_pdb_file";
	system "/ul/maiti/src/util/scripts/fixcurvepdb.pl $curr_pdb_file > junk";
	$pic_file = "pics/" . $P_File->{"Molecule"}->{"name"} . "_" . 
	    $i . ".png";
	if (-e $curr_pdb_file) {
	    push @namot_string, "load pdb na $curr_pdb_file";
	    for $j (@unit_array) {
		$temp = $EngDatStructure{$j}->{$i}->{"TOTAL"};
		$temp = DetermineRGBVal($temp, $eng_max, $eng_min);
		push @namot_string, "set color m1:1:$j $temp";
	    }
	    push @namot_string, "set background black";
	    push @namot_string, "write png $pic_file";
	    DoNamotCmd(\@namot_string);
	    if (-e $spec_pic) {
		system "cp $spec_pic . ";
	    }
	    if (-e $pic_file) {
		AnnotatePic($pic_file, "Time $i ps", 20, 20, 0, 3, 2);
		ColasceImages($pic_file, $spec_pic, $pic_file, 0, 50, 2, 3);
		$text_length = length("$eng_max KCal") + 65;
		AnnotatePic(
			    $pic_file, 
			    "$eng_max KCal", 
			    13, 30, 
			    $text_length, 
			    2, 3
			    );
		$text_length = length("$eng_min KCal") + 65;
		AnnotatePic(
			    $pic_file, 
			    "$eng_min KCal", 
			    13, 325, 
			    $text_length, 2, 3
			    );
		$jpeg_file = "pics/" . 
		    $P_File->{"Molecule"}->{"name"} . "_" . $counter . ".jpg";
		ConvertPicture($pic_file, $jpeg_file);
#		system "rm -f $pic_file";
		push @valid_pics, $jpeg_file;
		$counter++;
	    }
		
	}
    }

    return @valid_pics;
}

sub MovieGen(@) {
    my ($files) = $_[0];
    my ($start) = $files->[0];

    ($start =~ /(\d+)\.jpg$/) ?
	$start = $1 :
	$start = 1;

    my $rec = (
	       {
		   "start" => $start,
		   "end"   => ($start + $#{$files}) ,
		   "filebase" => "pics/" . $P_File->{"Molecule"}->{"name"} . "_",
		   "extension" => "jpg",
		   "host"      => $P_File->{"Movie"}->{"host"},
		   "user"      => $P_File->{"Movie"}->{"user_name"},
		   "savename"  => $P_File->{"Molecule"}->{"name"} . "_movie.avi",
		   "framerate" => $P_File->{"Movie"}->{"length"}/($#{$files}),
		   "moviesize" => $P_File->{"Movie"}->{"size"},
	       }
	       );

    CreateMovie($rec) ?
	return 1 :
	return 0;
    

}

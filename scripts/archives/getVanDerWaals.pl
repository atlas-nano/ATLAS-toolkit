#!/usr/bin/perl -w

# getVanDerWaals.pl: This script will take an input PDB file, load two instances,
# and apply random rotations to both. It will then apply a random translation to the
# second unit, create an amber topology and coordinate file and perform MD minimization.
# It will then analyze the simulation, record the total nonbond interactions and store it
# as a function of the center of mass distance between the two units.
# Finally, it will use Mathematica to fit the resulting energy curve to a morse potential
BEGIN {
    push @INC, "/ul/tpascal/.libs/Math-ematica-1.108";
    push (@INC, "/ul/tpascal/scripts/");
}

use strict;
use File::Basename;
use Packages::General;
use Packages::PDBReader;
use Math::ematica qw(:PACKET :TYPE :FUNC);

sub VerifyFiles(@);
sub GetParms(@);
sub GetPDBInfo(@);
sub GetRandomNumber(@);
sub ApplyRotation(@);
sub ApplyTranslation(@);
sub WriteFile(@);
sub CreateAmberParms(@);
sub PerformMinimization(@);
sub AnalyzeFile(@);
sub GetComDist(@);
sub GetLeastSquaresFit(@);
sub WriteOutput(@);
sub Numerically;
sub RecordStat(@);
sub GetMasses(@);
sub MathLinkCmd(@);
sub ShowError(@);
sub DoCmd(@);
sub CreateMinimizedStructure(@);
sub SmoothCurve(@);
sub CreateCtlFile(@);
sub GexMax(@);
sub RunMpSim(@);
sub DoMpSimAnal(@);

die "usage: $0 pdbfile parm_file save_name\n"
    if (! @ARGV or $#ARGV < 2);

my ($pdb_file, $parm_file, $save_name) = @ARGV;
my ($PARMS, $PDBFILE, $Angles, $CURR_FILE, $Trans, $out_file, $counter, $mpsim_in);
my ($topology, $coordinate, $isConverged, $Eng, $com_dist, $data_file, $DATA);
my ($d0, $alpha, $r0, %STATS, $shouldProceed, $MASSES, $num_atms, $min_file, $p_name);

{
    system "mkdir -p SANDER PDBFILES";
    print "Initializing...";
    $shouldProceed = VerifyFiles($ARGV[0], $ARGV[1]);
    die "\n"
	if (! $shouldProceed);

    print "Done\nGetting Parameters...";
    $PARMS = GetParms($parm_file);
    print "Done\nGetting PDB File info...";
    ($PDBFILE, $num_atms) = GetPDBInfo($pdb_file, $PARMS);
    $MASSES = GetMasses($PARMS->{"MassesFile"});
    print "Done\n";
    
    $counter = $PARMS->{"Start"};
    $isConverged = 0;
    while ($counter < ($PARMS->{"Num_Structures"} + $PARMS->{"Start"})) {
	if ($isConverged == 0) {
	    print "Generating Structure $counter...";
	}
	$Angles = GetRandomNumber(360);
	$CURR_FILE->{1} = ApplyRotation($PDBFILE, $Angles);
	RecordStat(\%{ $STATS{$counter} }, $Angles, "ANGLES", 1);
	$Angles = GetRandomNumber(360);
	$CURR_FILE->{2} = ApplyRotation($PDBFILE, $Angles);
	RecordStat(\%{ $STATS{$counter} }, $Angles, "ANGLES", 2);
	
	$Trans = GetRandomNumber(($PARMS->{"CuttOff"}/sqrt(3)));
	ApplyTranslation(\%{ $CURR_FILE->{2} }, $Trans);
	RecordStat(\%{ $STATS{$counter} }, $Trans, "TRANS", 2);
	
	$out_file = WriteFile($CURR_FILE, $counter);
	
#	($mpsim_in, $p_name) = CreateCtlFile($out_file);
	($topology, $coordinate) = CreateAmberParms($out_file, $PARMS->{"AmberLeapFile"}, $PARMS->{"AmberOffFile"});
	
	if ($isConverged == 0) {
	    print "Performing Minimization...";
	}
	($data_file, $coordinate) = PerformMinimization($topology, $coordinate, \%{ $PARMS->{"Sims"} });
#	($data_file) = RunMpSim($mpsim_in, $p_name);
	if ($data_file eq "") {
	    $isConverged = 2;
	} else {
	    ($isConverged, $Eng) = AnalyzeFile($data_file, \%{ $PARMS->{"Anal"} });
#	    ($isConverged, $Eng) = DoMpSimAnal($data_file);
	}
	if ($isConverged == 1) {
#	    $min_file = CreateMinimizedStructure($topology, $coordinate, $out_file);
	    $min_file = $out_file;
	    $com_dist = GetComDist($min_file, $MASSES, $num_atms);
	    $STATS{$counter}{"ENERGY"} = $Eng;
	    $STATS{$counter}{"COM"} = $com_dist;
	    print "Converged!\n";
	    $counter++;
#	    DoCmd("rm -fr $topology $coordinate *.crd");
	    $isConverged = 0;
	} else {
	    print ".";
	    $isConverged = 2;
	}
 
    }

    open ERRFILE, "> dumpfile.txt" or die "Cannot create dumpfile.txt: $!\n";
    print "Obtaining Least SquareFit...";
    $DATA = SmoothCurve(\%STATS);
    %STATS = ();
#    GetLeastSquaresFit($DATA);
    close ERRFILE;
    print "Done\n--===RESULTS==--\n";
    WriteOutput($save_name, $DATA);
}
    
sub VerifyFiles(@) {
    my ($counter, $isValid);

    $isValid = 1;
    for $counter (@_) {
	if (! -e $counter or ! -r $counter or ! -T $counter) {
	    print "Cannot access $counter: $!\n";
	    $isValid = 0;
	    last;
	}
    }

    return $isValid;
}

sub GetParms(@) {
    my ($in_file) = $_[0];
    my ($counter, %Parms);

    $counter = 0;
    open INFILE, $in_file or die "Cannot open $in_file: $!\n";
    $Parms{"Start"} = 1;

    while (<INFILE>) {
	chomp;
	if ($_ =~ /^Total Structures:\s+(\d+)/) {
	    $Parms{"Num_Structures"} = $1;
	    $counter++;
	} elsif ($_ =~ /^NonBond Cuttoff:\s+(\d+)/) {
	    $Parms{"CuttOff"} = $1;
	    $counter++;
	} elsif ($_ =~ /^xLeap template:\s+(\S+)/) {
	    if (VerifyFiles($1)) {
		$Parms{"AmberLeapFile"} = $1;
		$counter++;
	    }
	} elsif ($_ =~ /^xLeap Off File:\s+(\S+)/) {
	    if (VerifyFiles($1)) {
		$Parms{"AmberOffFile"} = $1;
		$counter++;
	    }
	} elsif ($_ =~ /^Minimization File:\s+(\S+)/) {
	    if (VerifyFiles($1)) {
		$Parms{"Sims"}{"MinFile"} = $1;
		$counter++;
	    }
	} elsif ($_ =~ /^Residue Name:\s+(\w+)/) {
	    $Parms{"Res_Name"} = $1;
	    $counter++;
	} elsif ($_ =~ /^Sim Executable:\s+(\S+)/) {
	    if (-e $1) {
		$Parms{"Sims"}{"EXE"} = $1;
		$counter++;
	    }
	} elsif ($_ =~ /^Sim Processors:\s+(\d+)/) {
	    $Parms{"Sims"}{"Processors"} = $1;
	    $counter++;
	} elsif ($_ =~ /^Sim MpiRun:\s+(\S+)/) {
	    if (VerifyFiles($1)) {
		$Parms{"Sims"}{"MPIRun"} = $1;
		$counter++;
	    }
	} elsif ($_ =~ /^Sim MachineFile:\s+(\S+)/) {
	    $Parms{"Sims"}{"MachineFile"} = $1;
	    $counter++;
	} elsif ($_ =~ /^Masses File:\s+(\S+)/) {
	    if (VerifyFiles($1)) {
		$Parms{"MassesFile"} = $1;
		$counter++;
	    }
	} elsif ($_ =~ /^Do Analysis Only:\s+(\d+)/) {
	    $Parms{"Anal_Only"} = $1;
	} elsif ($_ =~ /^Structures Start:\s+(\d+)/) {
	    $Parms{"Start"} = $1;
	}
    }
    close INFILE;

    if (defined($ENV{PBS_NODEFILE})) {
	$counter++;
	$Parms{"Sims"}{"MachineFile"} = $ENV{PBS_NODEFILE};
	print "Using value in PBS_NODEFILE $ENV{PBS_NODEFILE}...";
    }

    die "ERROR: No or more parameters are missing. Check $in_file\n"
	if ($counter < 11);

    return \%Parms;
 
	    
}

sub GetPDBInfo(@) {
    my ($pdb, $Parms) = @_;
    my ($Atom_Info) = GetPDBFileInfo($pdb, 0, " ", " ");
    my ($counter, $unit_counter);

    $unit_counter = 0;
    for $counter (keys %{ $Atom_Info }) {
	$Atom_Info->{$counter}{"RES_NAME"} = $Parms->{"Res_Name"};
	$unit_counter++;
    }
    return ($Atom_Info, $unit_counter);
}

sub GetRandomNumber(@) {
    my ($max) = $_[0];
    my (@output_array, $rand_angle, $counter);

    for $counter (0 .. 2) {
	$rand_angle = rand $max;
	$output_array[$counter] = sprintf("%.2f", $rand_angle);
    }

    return \@output_array;
    
}

sub ApplyRotation(@) {
   my ($in_file, $rot_array) = @_;
   my (@rotation_matrix, $counter, $index, $matrix_keys, %F_Data);

   for $counter (keys %{ $in_file }) {
       for $index (keys %{ $in_file->{$counter} }) {
	   $F_Data{$counter}{$index} = $in_file->{$counter}{$index};
       }
   }

   @rotation_matrix = (
		       {
			   "XCOORD" => [1, 0, 0],
			   "YCOORD" => [0, cos($rot_array->[0]), sin($rot_array->[0])],
			   "ZCOORD" => [0, -sin($rot_array->[0]), cos($rot_array->[0])],
			},
		       
		       {
			   "XCOORD" => [cos($rot_array->[1]), 0, -sin($rot_array->[1])],
			   "YCOORD" => [0, 1, 0],
			   "ZCOORD" => [sin($rot_array->[1]), 0, cos($rot_array->[1])],
		       },
		       {
			   "XCOORD" => [cos($rot_array->[2]), sin($rot_array->[2]), 0],
			   "YCOORD" => [-sin($rot_array->[2]), cos($rot_array->[2]), 0],
			   "ZCOORD" => [0, 0, 1],
		       },
		       );

   for $matrix_keys (keys %F_Data) {
       $counter = \%{ $rotation_matrix[0] };
       for $index ("XCOORD", "YCOORD", "ZCOORD") {
	   $F_Data{$matrix_keys}{$index} = ($in_file->{$matrix_keys}{"XCOORD"} * $counter->{$index}[0] +
					    $in_file->{$matrix_keys}{"YCOORD"} * $counter->{$index}[1] +
					    $in_file->{$matrix_keys}{"ZCOORD"} * $counter->{$index}[2]);
       }
   }
   @rotation_matrix = ();
   return \%F_Data;
	       
}

sub ApplyTranslation(@) {
    my ($in_file, $trans_array) = @_;
    my ($matrix_keys, $counter, $index);

   for $matrix_keys (keys %{ $in_file }) {
       $counter = 0;
       for $index ("XCOORD", "YCOORD", "ZCOORD") {
	   $in_file->{$matrix_keys}{$index} += $trans_array->[$counter];
	   $counter++;
       }
   }
}

sub WriteFile(@) {
    my ($F_Data, $F_counter) = @_;
    my ($out_name, $counter, $fmt, $out_text, $f_counter, $atom, $unit_counter);

    $fmt = "%-6s%5d%4s%5s%6d%12.3f%8.3f%8.3f\n";
    $out_name = "orient" . $F_counter . ".pdb";
    $unit_counter = 0;

    if (defined($PARMS->{"Anal_Only"}) and $PARMS->{"Anal_Only"} == 1) {
	return $out_name;
    }

    for $f_counter (1, 2) {
	for $counter (sort Numerically keys %{ $F_Data->{$f_counter} }) {
	    $atom = \%{ $F_Data->{$f_counter}{$counter} };
	    $unit_counter++;
	    $out_text .= sprintf($fmt, $atom->{"HEADER"}, $unit_counter, $atom->{"LABEL"}, $atom->{"RES_NAME"},
				 $f_counter, $atom->{"XCOORD"}, $atom->{"YCOORD"}, $atom->{"ZCOORD"});
	}
	$out_text .= "TER\n";
    } 
	    
    if (defined($PARMS->{"Anal_Only"}) and $PARMS->{"Anal_Only"} == 1) {
	return $out_name;
    }

    open OUTFILE, "> PDBFILES/$out_name" or die "Cannot write to $out_name: $!\n";
    print OUTFILE $out_text;
    close OUTFILE;

    return $out_name;
}

sub Numerically {
    ($a<=>$b);
}

sub RecordStat(@) {
    my ($curr_stat, $in_data, $title, $which_structure) = @_;
    my (@dim_data, $counter);

    @dim_data = ("X", "Y", "Z");

    for $counter (0 .. 2) {
	$curr_stat->{$which_structure}{$title}{$dim_data[$counter]} = sprintf("%.2f", $in_data->[$counter]);
    }
}

sub CreateAmberParms(@) {
    my ($file_nm, $leap_src, $off_src) = @_;
    my ($mycmd, $top_file, $crd_file);

    die "\n"
	if (! VerifyFiles($leap_src) or ! VerifyFiles($off_src));
    
    $mycmd = "cp $leap_src ./leaprc";
    if (! DoCmd($mycmd)) {
	die "ERROR while attempting to copy $leap_src to curr directory: $!\n";
    }
    
    $top_file = $crd_file = "SANDER/" . $file_nm;
    $top_file =~ s/\.pdb$/\.top/;
    $crd_file =~ s/\.pdb/\.crd/;

    open MYLEAPRC, ">> leaprc" or die "Cannot write to leaprc: $!\n";
    print MYLEAPRC "loadoff $off_src\n";
    print MYLEAPRC "dna = loadpdb PDBFILES/$file_nm\n";
    print MYLEAPRC "saveamberparm dna $top_file $crd_file\n";
    
    print MYLEAPRC "quit\n";
    
    close MYLEAPRC;
    
    if (! DoCmd("tleap >& xleap_out")) {
	die "ERROR executing tleap\n";
    }

#    DoCmd("rm -fr leaprc leap.log");
    return ($top_file, $crd_file);

}

sub PerformMinimization(@) {
    my ($top_file, $crd_file, $Sim_opts) = @_;
    my ($parallel_cmd, $out_text, $mdout_file, $restart_file, $curr_dir);


    $curr_dir = $ENV{PWD};

    $out_text = "SANDER";
    
    $mdout_file = $restart_file = basename($top_file);
    $mdout_file =~ s/\.top$/\.out/;
    $restart_file =~ s/\.top$/\.restart/;

    $mdout_file = $out_text . "/" . $mdout_file;
    $restart_file = $out_text . "/" . $restart_file;
    if ($Sim_opts->{"Processors"} > 1) {
	$parallel_cmd = $Sim_opts->{"MPIRun"} . " -np " . $Sim_opts->{"Processors"} . " -machinefile " .
	    $Sim_opts->{"MachineFile"} . " ";
    }
    
    $parallel_cmd .= $Sim_opts->{"EXE"} .  " -O -i " . $Sim_opts->{"MinFile"} . " -p $top_file" .
	" -c $crd_file -o $mdout_file -r $restart_file -ref $crd_file >& sander_output"; 
    
    DoCmd($parallel_cmd);
    if (! -e $mdout_file or ! -e $restart_file ) {
	return ("", "");
    }

    return  ($mdout_file, $restart_file);


}

sub AnalyzeFile(@) {
    my ($eng_file, $Anal_opts) = @_;
    my ($isvalid, $converged, %Eng, $tstep, $use_next_line);
    my ($in_string, $curr_eng, $isEnd, @holder, $counter);
    $tstep = $use_next_line = $isEnd = 0;
    $isConverged = 1;

    open ENGFILE, $eng_file or die "Cannot open simulation output file $eng_file: $!\n";
    while (<ENGFILE>) {
	if ($_ =~ /REPEATED LINMIN FAILURE/) {
#	    $isConverged = 0;
	    last;
	} elsif ($_ =~ /TIMINGS/) {
	    last;
	} elsif ($_ =~ /^\s+FINAL RESULTS/) {
	    last;
	}elsif ($_ =~ /^\s+NSTEP\s+ENERGY\s+RMS\s+GMAX\s+NAME\s+NUMBER/) {
	    $use_next_line = 1;
	} elsif ($use_next_line and $_ =~ /^\s+(\d+)\s+(\-?\d+\.\d+\S*)\s+(\-?\d+\.\d+\S*)\s+(\-?\d+\.\d+\S*)\s+(\w+)\s+(\d+)/) {
	    $tstep = $1;
	    $Eng{"Etot"} = $2;
	    $Eng{"RMS"}   = $3;
	    $Eng{"GMAX"}   = $4;
	    $Eng{"NAME"}   = $5;
	    $Eng{"NUMBER"}   = $6;
	    
	} elsif ($tstep > 0 and $_ =~ /^\s[BOND|1-4 VDW|VDWAALS]/) {
	    $in_string = $_;
	    $in_string =~ s/=//g;
	    $in_string = Trim($in_string);
	    @holder = split /\s\s+/, $in_string;
	    $counter = 0;
	    while ($counter <= ($#holder - 1)) {
		$Eng{"ENG"}{$holder[$counter]} = $holder[$counter + 1];
		$counter += 2;
	    }
	}
    }

    close ENGFILE;
    
    if (! %Eng) {
	return (0, 0);
    }

    for $counter ("VDWAALS", "EEL") { #, "EEL", "1-4 VDW", "1-4 EEL") {
	if ($Eng{"ENG"}{$counter} =~ /^\-?\d+/) {
	    $curr_eng += $Eng{"ENG"}{$counter};
	} else {
	    $curr_eng = "ERROR";
	    $isConverged = 0;
	    last;
	}
    }
    return ($isConverged, $curr_eng);
}

sub GetMasses(@) {
    my ($in_file) = $_[0];
    my (%Masses);

    open MASSFILE, $in_file or die "Cannot open masses file $in_file: $!\n";
    while (<MASSFILE>) {
	if ($_ =~ /^MASS_(\w+)\s+(\d+\.*\d*)/) {
	    $Masses{$1} = $2;
	}
    }
    close MASSFILE;

    die "ERROR: $in_file does not contain any mass information\n"
	if (! %Masses);
    
    return (\%Masses);
}

sub GetComDist(@) {
    my ($file_nm, $masses, $total_atms) = @_;
    my ($comDist, $counter, $f_counter, $atom, $total_mass, $sum, %POS);
    my ($currMass, $atmLabel, $index, $unit_counter, $structures, %MASS);


    $unit_counter = $f_counter = 1;
    $total_mass = 0;
    $file_nm = "PDBFILES/" . $file_nm;
    $structures = GetPDBFileInfo($file_nm, 0, " ", " ");
    for $counter (sort Numerically keys %{ $structures } ) {
	if ($unit_counter == ($total_atms + 1)) {
	    $MASS{1} = $total_mass;
	    $total_mass = 0;
	    $f_counter = 2;
	}
	$atom = \%{ $structures->{$counter} };
	$atmLabel = $atom->{"LABEL"};
	if ($atmLabel =~ /([A-Z])/i) {
	    $atmLabel = $1;
	}
	$currMass = $masses->{$atmLabel};
	$total_mass += $currMass;
	for $index ("XCOORD", "YCOORD", "ZCOORD") {
	    $POS{$f_counter}{$index} += $atom->{$index} * $currMass;
	}
	$unit_counter++;
    }
    $MASS{2} = $total_mass;
    
    for $f_counter (1 .. 2) {
	for $index ("XCOORD", "YCOORD", "ZCOORD") {
	    $POS{$f_counter}{$index} = $POS{$f_counter}{$index} / $MASS{$f_counter};
	}
    }

    $comDist = 0;
    for $index ("XCOORD", "YCOORD", "ZCOORD") {
	$comDist += ($POS{1}{$index} - $POS{2}{$index})**2;
    }

    $comDist = sqrt($comDist);
    return $comDist;
}

sub GetLeastSquaresFit(@) {
    my ($stats) = $_[0];
    my ($counter, $points);
    my ($MathLinkResult, $holder, $link, $error, $math_cmd);
    $link = new Math::ematica '-linklaunch', '-linkname', 'math -mathlink';

    $math_cmd = "points = {";
    for $counter (keys %{ $stats }) {
	$math_cmd .= "{" . $counter . "," . $stats->{$counter}{"ENG"} . "}, ";
    }

    $math_cmd = substr($math_cmd, 0, length($math_cmd) - 2) . "};";
    MathLinkCmd($math_cmd, 1, $link);
    $math_cmd = "model = a * (Exp[-2*b*(x -c)] - 2*Exp[-b*(x - c)]);";
    MathLinkCmd($math_cmd, 1, $link);
    $math_cmd = "f1 = FindFit[points, model, {a, b, c}, x]";
    ($error, $holder) = MathLinkCmd($math_cmd, 0, $link);
    ShowError($error, $link);
    $math_cmd = "Plot[model/.f1,{x, 0, 20}, AxesOrigin->{0,0},PlotRange->{Automatic, {0, -20}}," . 
	"Epilog->Prepend[Point/@ points,PointSize[0.02]]];";
    MathLinkCmd($math_cmd, 1, $link);
    $link->PutFunction("Exit", 0);
    $link->EndPacket;
    $link = ();
    $holder =~ s/\->//g;
    $holder =~ s/,//g;
    if ($holder =~ /^\{a\s+(\-?\d+\.\d+)\s+b\s+(\-?\d+\.\d+)\s+c\s+(\-?\d+\.\d+)\}/) {
	$stats->{"Depth"} = $1;
	$stats->{"alpha"} = $2;
	$stats->{"Req"} = $3;
    } else {
	die "ERROR: Got junk from GetLeastSquaresFit: $holder\n";
    }

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

sub SmoothCurve(@) {
    my ($stats) = @_;
    my ($counter, $com, $energy, %P_Out, $avg, $std, $tot);
    
    for $counter (keys %{ $stats }) {
	if ($counter =~ /\d+/) {
	    $com = $stats->{$counter}{"COM"};
	    $com = sprintf("%.1f", $com);
	    $energy = $stats->{$counter}{"ENERGY"};
	    if ($energy =~ /^\-?\d+/) {
		$P_Out{$com}{"ENG"} .= "$energy ";
		$P_Out{$com}{"FILE"} .= "$counter ";
	    }
	}
    }

    for $counter (sort Numerically keys %P_Out) {
	chop $P_Out{$counter}{"ENG"};
	($avg, $std, $tot) = STDev($P_Out{$counter}{"ENG"});
	$P_Out{$counter}{"ENG"} = $avg;
	$P_Out{$counter}{"STDEV"} = $std;
	$P_Out{$counter}{"TOTAL"} = $tot;
    }
    return \%P_Out;
}

sub WriteOutput(@) {
    my ($out_name, $stats) = @_;
    my ($counter, $out_string, $avg, $std, $tot, @holder);

    $out_string = "";

    for $counter (keys %{ $stats }) {
	if ($counter =~ /^\d/) {
	    push @holder, $counter;
	}
    }

    for $counter (sort Numerically @holder) {
	($avg, $std, $tot) = ($stats->{$counter}{"ENG"}, $stats->{$counter}{"STDEV"}, $stats->{$counter}{"TOTAL"});
#	if (abs($std/$avg) < 0.01) {
	    $out_string .= sprintf("%12.2f %12.3f %12.5G", $counter, $avg, $std);
	    $out_string .= "  " . $stats->{$counter}{"FILE"} . "\n";
#	}
#	print "HERE\n";
    }

    for $counter ("Depth", "alpha", "Req") {
#	$out_string .= "$counter: " . $stats->{$counter} . "\n";
    }

    open OUTFILE, "> $out_name" or die "Cannot write to $out_name: $!\n";
    print OUTFILE $out_string;
    close OUTFILE;

}

sub DoCmd(@) {
    my ($cmd_str) = $_[0];
    my ($returnval) ;

    $returnval = 1;
    if (defined($PARMS->{"Anal_Only"}) and $PARMS->{"Anal_Only"} == 1) {
	return 1;
    }
    open STDERR, ">> errors";
    if (system($cmd_str)) {
	print STDERR "ERROR executing $cmd_str\n";
	$returnval = 0;
    }
    close STDERR;

    return $returnval;
}

sub CreateMinimizedStructure(@) {
    my ($top, $crd, $old_name) = @_;
    my ($out_cmd, $old_file, $save_name);

    $save_name = $old_name;
    $save_name =~ s/\.pdb$/_min\.pdb/g;
    
    $out_cmd = "ambpdb -p $top < $crd > PDBFILES/$save_name";
    DoCmd($out_cmd);
    return $save_name;
}

sub CreateCtlFile(@) {
    my ($pdb_file) = $_[0];
    my ($out_cmd, @vals, $bgf_file, $project_name, $outstring);

    $pdb_file = "PDBFILES/" . $pdb_file;
    $project_name = basename($pdb_file);
    $project_name =~ s/\.pdb$//;
    $bgf_file = "BGFFILES/" . $project_name . ".bgf";

    @vals = GetMax($pdb_file);
    
    $out_cmd = "/ul/tpascal/scripts/amber2bgf.pl ";
    $out_cmd .= "$pdb_file /ul/tpascal/ff/AMBER95.cnv ";
    $out_cmd .= "$bgf_file $vals[0] $vals[1] $vals[2]";

    if (! DoCmd($out_cmd)) {
	print "ERROR Creating bgf file $bgf_file\n";
	return 0;
    }

    $out_cmd = "cp /ul/tpascal/mpsim/ewald.par .";
    if (! DoCmd($out_cmd)) {
	print "ERROR copying /ul/tpascal/mpsim/ewald.par\n";
    }

    open CTL, "/ul/tpascal/mpsim/1-pme.ctl" or die "Cannot open 1-pme.ctl: $!\n";
    while (<CTL>) {
	chomp;
	if ($_ =~ /\[p_name\]/) {
	    $outstring .= "PROJECT               $project_name";
	} elsif ($_ =~ /\[bgf_file\]/) {
	    $outstring .= "STRUCTURE             $bgf_file";
	} else {
	    $outstring .= $_;
	}
	$outstring .= "\n";
    }
    close CTL;
    open CTL, "> 1-pme.ctl" || die "Cannot modify control file: $!\n";
    print CTL $outstring;
    close CTL;
    
    return ("1-pme.ctl", $project_name);

}

sub GetMax(@) {
    my ($in_file) = $_[0];
    my ($PDB_INFO) = GetPDBFileInfo($in_file, 0, " ", " ");
    my ($counter, $BBOX, $j, $k);

     for $counter (keys %{ $PDB_INFO }) {
	for $j ("XCOORD", "YCOORD", "ZCOORD") {
	    for $k  ("MAX", "MIN") {
		if (!defined($BBOX->{$j}{"MAX"}) or $PDB_INFO->{$counter}{$j} > $BBOX->{$j}{"MAX"}) {
		    $BBOX->{$j}{"MAX"} = $PDB_INFO->{$counter}{$j};
		}
		if (!defined($BBOX->{$j}{"MIN"}) or $PDB_INFO->{$counter}{$j} < $BBOX->{$j}{"MIN"}) {
		    $BBOX->{$j}{"MIN"} = $PDB_INFO->{$counter}{$j};
		}
	    }
	}
    }
    
    return ($BBOX->{"XCOORD"}{"MAX"} - $BBOX->{"XCOORD"}{"MIN"}, $BBOX->{"YCOORD"}{"MAX"} - $BBOX->{"YCOORD"}{"MIN"},
		     $BBOX->{"ZCOORD"}{"MAX"} - $BBOX->{"ZCOORD"}{"MIN"});
} 
		
sub RunMpSim(@) {
    my ($ctl_file, $project_name) = @_;
    my ($out_cmd);

    $out_cmd = "/ul/maiti/src/non_pme/build/linux/mpsim.20030918 ";
    $out_cmd .= "$ctl_file >& mpsim.out";
    if (! DoCmd($out_cmd)) {
	print "Mpsim Error!\n";
	return "";
    }
    system "rm -f 1-pme.ctl ewald.par mpsim.out";
    system "cp $project_name" . ".fin.ener outfiles/$project_name" . ".out";
    system "rm -f $project_name" . ".*";
    return "outfiles/" . $project_name . ".out";
}

sub DoMpSimAnal(@) {
    my ($curr_fle) = $_[0];
    my ($in_data, $counter, @eng_header, $eng_label);
    my ($curr_val, $curr_total, $eng_total);

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
		$eng_total = 0;
		while ($in_data =~ /(\-?\d+\.\d+)\s*/g && $counter <= $#eng_header) {
		    $eng_label = $eng_header[$counter];
		    if ($eng_label eq "VDW" or $eng_label eq "Coulomb") {
			$eng_total += $1;
		    }
		}

	    }
	}
	
	CLOSE INPFILE;
	
	if ($#eng_header > 0) {
	    return (1, $eng_total);
	}else {
	    print "Failure. $curr_fle is invalid";
	    return (0, "");
	}
	close INPFILE;
    } else {
	print "WARNING: Invalid file: $curr_fle";
	return (0, "");
    }

}

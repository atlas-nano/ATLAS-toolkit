#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}
use strict;
use File::Basename;
use Packages::GetParms;
use Cwd;

# This script will take a dna topology file, and coordinate file
# And perform a simulation on them. It will then perform analysis

sub GetParms();
sub DoSim();
sub DoAnalysis();
sub CreateNewDirectory(@);
sub DoCmd(@);
sub CreateSanderCmd(@);
sub FixInFile;
sub CreateNewParmFile();
sub EngStructureAnal(@);

die "Usage: $0 parameter_file sim_cmd [parallel_cmd]\n"
    if (! @ARGV or $#ARGV < 1);

my($parmfile, $pmemd_cmd, $parallel_cmd) = @ARGV;

print "\n---START---\n\n";
my $Parms;
GetParms();

if (! $parallel_cmd) {
    $Parms->{"Parallel_Cmd"} = $ENV{DO_PARALLEL}
	if ($ENV{DO_PARALLEL});
} else {
    $Parms->{"Parallel_Cmd"} = $parallel_cmd;
}

DoSim()
    if ($Parms->{"Perform_Sim"});
DoAnalysis()
    if ($Parms->{"Perform_Anal"});
print "\n---END---\n";

sub GetParms() {
    my ($in_data);
    die "Cannot locate parameter file $parmfile: $!\n"
	if (! -e $parmfile);

    $Parms = Packages::GetParms->new();
    if (! $Parms->IsValidParams($parmfile)) {
        die "Error in Paramater file\n";
    }
 
    $Parms->{"StartDir"} = cwd();
    $Parms->{"Perform_Anal"} = 1;
    $Parms->{"Perform_Sim"} = 1;
    $Parms->{"Multi_Helix"} = 1;
    $Parms->{"Do_Vels"} = 0;

    die "Cannot open $parmfile: $!\n"
	if (! open(PARMFILE, $parmfile));
    while (<PARMFILE>) {
	chomp;
	$in_data = $_;
	if ($in_data =~ /^Directory Name: (\w+)/) {
	    $Parms->{"SaveName"} = $1;
	} elsif ($in_data =~ /^Production Steps: (\d+)/) {
	    $Parms->{"Prod_Steps"} = $1;
	} elsif ($in_data =~ /^Perform Analysis: 0/) {
	    $Parms->{"Perform_Anal"} = 0;
	    print "NOTE: Will not Perform Analysis\n";
	} elsif ($in_data =~ /^Perform Simulation: 0/) {
	    $Parms->{"Perform_Sim"} = 0;
	    print "NOTE: Will not perform Simulation\n";
	} elsif ($in_data =~ /^Parallel Command: (.+)$/) {
	    $Parms->{"Parallel_Cmd"} = $1;
	} elsif ($in_data =~ /^Do Equil: 0/) {
	    $Parms->{"Exclude_Equil"} = 1;
	    print "NOTE: Skipping Equilibration\n";
	} elsif ($in_data =~ /^Do Dynamics: 0/) {
	    print "NOTE: Skipping NVT Dynamics\n";
	    $Parms->{"Exclude_Dynamics"} = 1;
	} elsif ($in_data =~ /^Do Production: 0/) {
	    print "NOTE: Skipping Production Run\n";
	    $Parms->{"Exclude_Production"} = 1;
	} elsif ($in_data =~ /^Helix1 Info: (\w+)/) {
	    $Parms->{"Helix1"} = $1;
	} elsif ($in_data =~ /^Helix2 Info: (\w+)/) {
	    $Parms->{"Helix2"} = 2;
	} elsif ($in_data =~ /^Multi Helix: 0/) {
	    $Parms->{"Multi_Helix"} = 0;
	} elsif ($in_data =~ /^Is Restart: 1/) {
	    $Parms->{"Restart_Run"} = 1;
	} elsif($in_data =~ /^Analysis Start: (\d+)/) {
	    $Parms->{"Anal"}{"Start"} = $1;
	} elsif($in_data =~ /^Analysis End: (\d+)/) {
            $Parms->{"Anal"}{"End"} = $1;
        } elsif($in_data =~ /^Analysis Interval: (\d+)/) {
            $Parms->{"Anal"}{"Interval"} = $1;
        } elsif ($in_data =~ /^Production Traj Dump: (\d+)/) {
	    $Parms->{"Traj_Dump"} = $1;
	} elsif ($in_data =~ /^Analysis Reference: (\w+)/) {
	    if (-e $1) {
		$Parms->{"Anal"}{"Reference"} = $1;
	    }
	} elsif ($in_data =~ /^Skip Structureanal/) {
	    $Parms->{"Skip_Structure"} = 1;
	} elsif ($in_data =~ /^Skip Energyanal/) {
	    $Parms->{"Skip_Energy"} = 1;
	} elsif ($in_data =~ /^Velocity Step: (\d+)/) {
	    print "NOTE: Will dump velocities for $1 steps for entropy calc\n";
	    $Parms->{"Do_Vels"} = 1;
	    $Parms->{"Vel_Steps"} = $1;
	}
    }
    close PARMFILE;

    $Parms->{"Perform_Sim"} = 0
	if ($Parms->{"Exclude_Equil"} and $Parms->{"Exclude_Dynamics"} and 
	    $Parms->{"Exclude_Production"});
    if (! $Parms->{"Anal"}{"Start"} or ! $Parms->{"Anal"}{"End"} or ! $Parms->{"Anal"}{"Interval"}) {
	$Parms->{"Perform_Anal"} = 0;
	print "Will not Perform Analysis\n";
    }
 
    die "Nothing to do! Terminating Execution\n"
	if (! $Parms->{"Perform_Anal"} and ! $Parms->{"Perform_Sim"});
    die "ERROR: Name of Temporary Directory not specified!\n"
	if (! $Parms->{"SaveName"});
    die "ERROR: Number of Steps for Productions Run not specified!\n"
	if (! $Parms->{"Prod_Steps"});
    die "ERROR: Numer of steps of Production has to exceed 1000!\n"
	if ($Parms->{"Prod_Steps"} < 1000);
  
}

sub DoSim() {

    my ($out_text);

    print "\n====START SIMULATION===\n\n";
#    $Parms->{"SaveDir"} = CreateNewDirectory($Parms->{"StartDir"}, $Parms->{"SaveName"});
#    ($Parms->{"SaveDir"}) ?
#	print "NOTE: Results will be stored in " . $Parms->{"SaveDir"} . "\n" :
#	die "Unable to create directory to store results. Aborting\n";

    $out_text = "/temp1/" . $ENV{USER};
    $Parms->{"TmpDirName"} = $out_text . "/" . $Parms->{"SaveName"};
    system ("mkdir -p " . $Parms->{"TmpDirName"})
	if (! -d $Parms->{"TmpDirName"});
    print "NOTE: Simulation will be run in: " . $Parms->{"TmpDirName"} . "\n";
	    
    chdir $Parms->{"TmpDirName"} or die "Cannot access " . $Parms->{"TmpDirName"} . ": $!\n";

    $out_text = "cp -r /ul/tpascal/md_simulations/amber8_parallel/* .";
    DoCmd($out_text);

    $out_text = "cp " . $Parms->{"Files"}->{"topology"} . " .";
    DoCmd($out_text);

    $out_text = "cp " . $Parms->{"Files"}->{"trajectory"} . " .";
    DoCmd($out_text);

    $Parms->{"Files"}->{"topology"} = $Parms->{"TmpDirName"} . "/" . 
	basename($Parms->{"Files"}->{"topology"});
    $Parms->{"Files"}->{"trajectory"} = $Parms->{"TmpDirName"} . "/" . 
	basename($Parms->{"Files"}->{"trajectory"});

#   ----APPEND number of bases to mininize/step1.in and dynamics/dyn.in---
    
    $out_text = "./addBaseNum.pl " . $Parms->{"Molecule"}->{"total_bases"};
    DoCmd($out_text);

#   ----EQUILIBRATION----
    if (! $Parms->{"Exclude_Equil"}) {
	print "\nSTARTING EQUILIBRATION PROCEDURE\n";
	chdir "./minimize" or die "ERROR: Cannot access the minimize sub-directory\n";
	
	system "cp " . $Parms->{"Files"}->{"trajectory"} . " ./start.restart";
	print "STEP 1: WATER/IONS EQUILIBRATION: 500kcal SOLUTE RESTRAINT, 4000 STEPS MIN\n";
	$out_text = CreateSanderCmd("step1","start", 1);
	DoCmd($out_text);

	print "STEP 2: NO SOLUTE RESTRAINT, 6000 STEPS MIN\n";
	$out_text = CreateSanderCmd("step2","step1", 1);
	DoCmd($out_text);

	print "\n--END EQUILIBRATION\n\n";
	$Parms->{"Files"}->{"trajectory"} = $Parms->{"TmpDirName"} . "/minimize/step2.trj";
	chdir "../";
    }
 
#-----NVT DYNAMICS
    if (! $Parms->{"Exclude_Dynamics"}) {
	print "\nSTARTING NVT DYNAMICS PROCEDURE\n";
	chdir "./dynamics" or die "ERROR: Cannot access the dynamics sub-directory\n";

	print "STEP 3:  10Kcal SOLUTE RESTRAINT, 25000 STEPS NVT dynamics (50 ps)\n";
	$out_text = CreateSanderCmd("dyn", "../minimize/step2", 2);
	DoCmd($out_text);

	print "\n--END NVT DYNAMICS\n\n";
	$Parms->{"Files"}->{"trajectory"} = $Parms->{"TmpDirName"} . "/dynamics/dyn.trj";
	chdir "../";
    }

#----PRODUCTION (NPT) DYNAMICS
    if (! $Parms->{"Exclude_Production"}) {
	print "\nSTARTING PRODUCTION (NPT) PROCEDURE\n";
	chdir "./production" or die "ERROR: Cannot access the production sub-directory\n";
	
	print "STEP 4:  " . $Parms->{"Prod_Steps"} . " STEPS Production Dynamics\n";
	if ($Parms->{"Restart_Run"}) {
            system "cp " . $Parms->{"Files"}->{"trajectory"} . " ./prod_rst.restart";
	    FixInFile("Prod", "prod.in");
	    $out_text = CreateSanderCmd("prod", "prod_rst", 3);
	} elsif ($Parms->{"Exclude_Dynamics"}) {
	    $out_text = "cp ../dynamics/dyn.restart ./prod_rst.restart";
	    system($out_text);
            FixInFile("Prod", "prod.in");
            $out_text = CreateSanderCmd("prod", "prod_rst", 3);
	} else {
 	    FixInFile("Prod", "prod.in");
 	    $out_text = CreateSanderCmd("prod", "../dynamics/dyn", 3);
	}
	DoCmd($out_text);

	print "\n--END PRODUCTION\n\n";
	$Parms->{"Files"}->{"trajectory"} = $Parms->{"TmpDirName"} . "/production/prod.trj";
	chdir "../";
    }

# ----VELOCITY (NPT) DYNAMICS
    if ($Parms->{"Do_Vels"}) {
	print "\nSTARTING NPT PROCEDURE DUMPING VELOCITIES\n";
	chdir "./velocity" or die "ERROR: Cannot access the velocity sub-directory\n";
	
	print "STEP 5: " . $Parms->{"Vel_Steps"} . " STEP NPT DYNAMICS (DUMPING VELOCITY EVERY 2fs)\n";
	FixInFile("Vel", "vel.in");
	$out_text = CreateSanderCmd("vel", "../production/prod", 4);
	DoCmd($out_text);

	print "\n--END VELOCITY PROCEDURE\n\n";
	chdir "../";
    }

    print "=====END SIMULATION====\n";	
}

sub CreateNewDirectory(@) {
    my ($curr_dir, $save_name) = @_;
    my ($return_name, $i);

    $i = 1;

    if (-d $curr_dir) {
	$return_name = $curr_dir . "/" . $save_name;

	if (! -d $return_name) {
	    system "mkdir $return_name";
	} else {
	    while ($i > 0) {
		$return_name = $curr_dir . "/" . $save_name . "_" . $i;
		if (! -d $return_name) {
		    system "mkdir $return_name";
		    last;
		}
		$i++;
	    }
	}
    }

    return $return_name;
}

sub DoCmd(@) {
    my ($sys_cmd) = $_[0];

    if (system($sys_cmd)) {
	die "Error executing command: $sys_cmd\n";
    }
}

sub CreateSanderCmd(@) {
    my ($curr_prefix, $previous_prefix, $which_sim) = @_;
    my ($return_string);

    if ($Parms->{"Parallel_Cmd"})  {
	$return_string = $Parms->{"Parallel_Cmd"} . " ";
#       print "PARALLE: $return_string\n";
    }

    $return_string .= "$pmemd_cmd -O -i " .
	$curr_prefix . ".in -o " . $curr_prefix . ".out -p " .
	$Parms->{"Files"}->{"topology"} . " -c " . $previous_prefix .
	".restart -r " . $curr_prefix . ".restart";
    if ($which_sim > 1) {
	$return_string .= " -x " . $curr_prefix . ".trj";
    }
    if ($which_sim  < 3) {
	$return_string .= " -ref " . $previous_prefix . ".restart";
    }
    if ($which_sim == 4) {
	$return_string .= " -ref " . $previous_prefix . ".restart -v " . $curr_prefix . ".vel";
    }

    return $return_string;
    
}

sub FixInFile {
    my ($field, $myFile) = @_;
    my ($out_text);

    $field .= "_Steps";
    $out_text = "";
    open PRODFILE, $myFile or die "Cannot open sander input file $myFile: $!\n";
    while (<PRODFILE>) {
	chomp;
	if ($_ =~ /^(\s+)nstlim = (\d+)/) {
	    $out_text .= $1 . "nstlim = " . $Parms->{$field} . ",\n";
	} elsif ($_ =~ /^(.+)ntwx\s+\=\s+\d+(.+)/ and $Parms->{"Traj_Dump"}) {
	    $out_text .= $1 . "ntwx = " . $Parms->{"Traj_Dump"} . $2 . "\n";
	}else {
	    $out_text .= $_ . "\n";
	}
    }
    close PRODFILE;

    open PRODFILE, "> $myFile" or die "Cannot modify sander input file $myFile: $!\n";
    print PRODFILE $out_text;
    close PRODFILE;
}

sub DoAnalysis() {
    my ($start, $end, $interval, $out_cmd, $curr_dir, $new_parm, $reference);

    $start = $Parms->{"Anal"}{"Start"};
    $end = $Parms->{"Anal"}{"End"};
    $interval = $Parms->{"Anal"}{"Interval"};
    $reference = "";

    if ($Parms->{"Anal"}{"Reference"}) {
	$reference = $Parms->{"Anal"}{"Reference"};
    }

    print "\n\n========ANALYSIS========\n\n";
    system "mkdir -p analysis";
    $curr_dir = "../";
    chdir $curr_dir;
    $new_parm = CreateNewParmFile();

    if (! $Parms->{"Exclude_Production"}) {
# --- Energy Analysis
	if (! $Parms->{"Skip_Energy"}) {
        print "\n---===ENERGY ANALYSIS===--\n";
 	$out_cmd = CreateNewDirectory($Parms->{"TmpDirName"}, "energy");
	if ($out_cmd) {
	    print "NOTE: Results will be stored in $out_cmd\n";
	    chdir $out_cmd;
	} else {
	    die "ERROR: Cannot create an analysis directory to store results. Aborting\n";
	}
	$curr_dir = $out_cmd;
	print "\n---Running Energy Analysis for frames $start - $end every $interval frames---\n";
	$out_cmd = "/ul/tpascal/scripts/energy_anal_nonamot.pl $new_parm " .
	    "$start $end $interval $reference";
	DoCmd($out_cmd);
	print "\n----DONE ENERGY ANALYSIS\n\n";
	}

# ---- Structure Analysis
	if (! $Parms->{"Skip_Structure"}) {
	print "\n---===STRUCTURE ANALYSIS===--\n";
	$out_cmd = CreateNewDirectory($curr_dir, "structure");
	if ($out_cmd) {
	    print "\nStoring Results in $out_cmd\n";
	    chdir $out_cmd;
	} else {
	    die "ERROR: Cannot create an directory to store results. Aborting\n";
	}

	$curr_dir = $out_cmd;
	$out_cmd = "/ul/tpascal/scripts/structure_anal.pl $new_parm " . $Parms->{"SaveName"} .
	           " $start $end $interval " . $Parms->{"Multi_Helix"};
	DoCmd($out_cmd);
	print "\n--===DONE STRUCTURE ANALYSIS===--\n";
	}

	$out_cmd = "mv structure energy " . $Parms->{"SaveDir"};
#	DoCmd($out_cmd);
	$curr_dir = "../";
	chdir $curr_dir;
    }
 
    print "\n=====END ANALYSIS\n\n";
}

sub CreateNewParmFile() {
    my ($return_name);

    $return_name = $Parms->{"TmpDirName"} . "/" . $Parms->{"SaveName"} . "_parm";
    open OUTFILE, "> $return_name" or die "Cannot write to $return_name: $!\n";
    print OUTFILE "Mol: " . $Parms->{"Molecule"}->{"major_groove"} . ":" .
	$Parms->{"Molecule"}->{"minor_groove"} . "\n";
    print OUTFILE "Total bases: " . $Parms->{"Molecule"}->{"total_bases"} . "\n";
    print OUTFILE "Bases at end: " . $Parms->{"Molecule"}->{"bases_at_end"} . "\n";
    print OUTFILE "Crossovers: ";
    for (@{$Parms->{"Molecule"}->{"crossovers"}}) {
	print OUTFILE "$_ ";
    }
    print OUTFILE "\n";
    print OUTFILE "Topology file: " . $Parms->{"Files"}->{"topology"} . "\n";
    print OUTFILE "Trajectory file: " . $Parms->{"Files"}->{"trajectory"} . "\n";
    print OUTFILE "3 prime in: " . $Parms->{"Molecule"}->{"is3PrimeIn"} . "\n";
    close OUTFILE;

    return $return_name;
}


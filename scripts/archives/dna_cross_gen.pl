#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}
                                                                                                                   
# dna_cross_gen.pl - automates the process of generating dna crossover molecules
#  allows creation of two independant helixes of b-dna and allows them to be
#  crossed over a specific points.
# note: it is very important when specifying the crossover points to be
#       very careful. Bad specifications will lead to bad structures
# process: utilizes a number of small perl scripts
# usage: dna_cross_gen.pl crossed-dna-name.pdb [scripts-dir]

# Load the Perl5 Namot2 module


use Packages::Namot;

p5namot::Cmd("set hush ERROR off");
p5namot::Cmd("set hush INFO off");
p5namot::Cmd("set hush REQUESTED off");
p5namot::Cmd("set hush WARNING off");
# Constants Declaration Section
$home_dir = "/ul/tpascal/scripts"; # set's the home directory of the scripts

# Variable declaration Section
$can_fixpdb = 1; # specifies whether the fixpdb.pl file is present
$output_file = "";
$mol_name = "";
$crossover_spec_file = ""; 
$is5prime = 0;

@helix_info = (
{
    "Crossovers" => 0,
    "Major_Groove" => 0,
    "Minor_Groove" => 0,
    "Base_Pairs" => 0,
    "Transalation" => 0,
    "5PrimeInside" =>0
}
);

sub GetYesNo($) {
# This sub will ask a yes/no question and get the correct answer

    my $valid_ans = 0;
    my $my_reply = 0;

    while (! $valid_ans) {
	my $ans = input($_[0]);
	$ans = lc($ans);
	if ($ans =~ /^y/) {
	    $my_reply =  1;
	    $valid_ans = 1;
	} else {
	    $my_reply = 0;
	    $valid_ans = 1;
	}
    }

    return $my_reply;
}

sub CreateAHelix($) {

    $whichstrand = $_[0];

    $output_string = $home_dir . "/createdna.pl ";

    $has_sequence = input ("Do you have a sequence file for this helix? [y/n]: ");
    if ($has_sequence =~ /^y|yes$/i) {
	$sequence1 = input ("Enter File Name: ");
	$sequence1 =~ s/\s+//;
	$file_not_found = 1;
	while ($file_not_found) {
	    if (! -e $sequence1) {
		print "Cannot locate $sequence1: $!\n";
		$sequence1 = input ("Enter File Name: ");
		$sequence1 =~ s/\s+//;
	    } else {
		last;
	    }
	}
	$output_string .= $sequence1;
    }
    $return_Str =  $mol_name . "_" . $whichstrand . ".pdb";
    $output_string .= " " . $return_Str;#. " 1";

    $build_status = system($output_string);
    
    if ($build_status != 0) {
	die "\nAn error has occured or you may have terminated the program.\n";
    }else {
	$execute_Str = "execute " . $return_Str . ".script";
	p5namot::Cmd($execute_Str);
	print "Done\n";
    }

# check for the parameter file - it holds the helix info
    $param_fle = $return_Str . ".parm";
    if (!open(PARAMFILE, $param_fle)) {
	die "Cannot open molecule paramater file. $!\nExecution Terminated.\n";
    } else {
	while (<PARAMFILE>) {
	    $incoming_line = $_;
	    if ($incoming_line =~ /PX (\d+):(\d+)/) { # first line
		$helix_info[0]{"Major_Groove"} = $1;
		$helix_info[0]{"Minor_Groove"} = $2;
	    }elsif ($incoming_line =~ /(\d+),(\w+)/) { #second line
		$helix_info[0]{"Base_Pairs"} = $1;
		$helix_info[0]{"Transalation"} = $2;
	    }
	}
	close PARAMFILE;
    }

# remove the paramater file
	$kill_str = "rm " . $param_fle;
	system($kill_str);


    return $return_Str;
}


sub Check_for_Scripts() {
# checks for the small perl scripts that this program uses to power itself
# only one is critical - createdna.pl

    $file_loc = $home_dir . "/createdna.pl";
    -e $file_loc or die "Critical error! Unable to locate $file_loc: $!.\nExecution Terminated.\n";

    $file_loc = $home_dir . "/makecrossovers.pl";
    -e $file_loc or die "Unable to locate $file_loc: $!.\nExecution Terminated\n";
    
#non critical files

    $file_loc = $home_dir . "/fixpdb4amber.pl";
-e $file_loc or warn "Unable to locate $file_loc: $!.\nWill not be able to fix pdb files\n";
	$can_fixpdb = 0;

}

sub GetCrossoverPoints() {
# GetCrossoverPoints - allows the user to input the crossover points

    #print "\nPlease ensure that the base numbers you are about to enter are correct.\n";
    #print "Incorrect base number will lead to bad crossover structures\n\n";

    print "\nEnter crossover points from bottom up (-99 to stop)\n";
    print "---------------------------------------------------\n";

    $line_counter = 1;

    $always_true = 1;
    while ($always_true) {
	print "$line_counter : ";
	$line_input = <STDIN>;

	chomp($line_input);

	if ($line_input =~ /-99/ and $line_counter > 1) {
	    $always_true = 0;
	} else {

	    if (! ($line_input =~ /\D+/ or length($line_input)==0)) {
		$line_input =~ /(\d+)/;
		$line_input = $1;

		if ($line_input > 0 and $line_input < $helix_info[0]{"Base_Pairs"}) {
		    push @crossover_pts, $line_input;
		    $line_counter += 1;
		}
	    }
	}
    }
    
    $helix_info[0]{"Crossovers"} = $#crossover_pts+1;
    CreateSpecsFile();

}

sub CreateSpecsFile() {
# CreateSpecsFile - creates the specs file

    $base_pair_ratio = $helix_info[0]{"Major_Groove"} . ":" . $helix_info[0]{"Minor_Groove"};
    $total_crossovers =  $helix_info[0]{"Crossovers"};
    $total_bases =  $helix_info[0]{"Base_Pairs"};
    $helix_transalation =  20;

if (!open(SPECSFLE, "> " . $crossover_spec_file)) {
    die "Cannot create file $crossover_spec_file. $!\n";
} else {
    print SPECSFLE "PX " . $base_pair_ratio . "\n";
    print SPECSFLE "Total Crossovers:" . $total_crossovers;
    print SPECSFLE " in " . $total_bases . " bases";
    print SPECSFLE ", Transalation:" . $helix_transalation . "\n";
    print SPECSFLE "CROSSOVERS\n";

    if ($helix_info[0]{"Crossovers"} > 1) {
	for ($i=0;$i<$helix_info[0]{"Crossovers"};$i++) {
	    $curr_base = $crossover_pts[$i];
	    print SPECSFLE "1:" . $curr_base . "->2:" . $curr_base . "\n";
	}
    } else {
	$curr_base = $crossover_pts[0];
	print SPECSFLE "1:" . $curr_base . "->2:" . $curr_base . "\n";
    }
}
    close SPECSFLE;
}


sub input($) {
    $printstring = $_[0];
    print "$printstring";
    $returnval = <STDIN>;
    return $returnval;
}



# -==Start==-

# validate cmd line argument
if (!@ARGV) {
    die "usage: dna_cross_gen.pl crossed-dna-name.pdb [scripts-dir]\n";
} else {
    $file_exists = 1;
    $output_file = $ARGV[0];

    if ($ARGV[1]) {
	$home_dir = $ARGV[1];
    }

    while ($file_exists) {
	if (open(PDBFILE, $output_file)) {
	    $overwrite = input ("The file $output_file already exists. Overwrite? [y/n]: ");
	    if ($overwrite =~ /^y|yes$/i) {
		last;
	    } else {
		$invalidname = 1;
		while ($invalidname) {
		    $output_file = input ("New Name: ");
		    if ($output_file =~ /w+/) {
			$output_file =~ s/\s+$//;
			last;
		    }
		}
	    }
	    close PDBFILE;
	} else {
	    $file_exists = 0;
	    last;
	}
    }
}

$mol_name = $output_file;

if ($output_file =~ /(\w+)\.pdb$/) {
    $mol_name = $1;
}

# Check for scripts
Check_for_Scripts();

print "\n-===Starting Program===-\n";
print "------------------------\n\n";


# Create helixes

$invalid_ans = 1;

while ($invalid_ans) {
    my $answer = input("\nAre the helixes arranged from 5' to 3'? [Y/N]: ");
    $answer = lc($answer);
    if ($answer =~ /^y$/) {
	$is5prime = "y";
	$invalid_ans = 0;
    }else {
	if ($answer =~/^n$/) {
	    $is5prime = "n";
	    $invalid_ans = 0;
	}
    }
}

print "\nCreating helix 1.\n";
print "-----------------\n";

$crossover_Str = CreateAHelix ("Strand1");

print "\nCreating helix 2.\n";
print "-----------------\n";

$crossover_Str = $crossover_Str ." " . CreateAHelix ("Strand2");

#Crossover helixes
$pdb_file_nm = $mol_name . ".pdb";

# --> 1. create crossover specification file

$crossover_spec_file = $mol_name . "_cross-spec.txt"; 

GetCrossoverPoints();

#rotate the helixes and get the appropriate angles
$helix1 = $mol_name . "_Strand1.pdb";
$helix2 = $mol_name . "_Strand2.pdb";

$script_fle = $home_dir . "/find_o_angle.pl $helix1 $helix2 $crossover_pts[0] $is5prime $home_dir";
system $script_fle;


# --> 2. execute program
$output_str = $home_dir . "/makecrossovers.pl $crossover_Str $pdb_file_nm $crossover_spec_file $is5prime";

$script_fle = "execute " . $pdb_file_nm . ".script";

system $output_str;

p5namot::Cmd($script_fle);
print "Done\n";

print "\nFixing Base Names...";

# Fix Pdb File to be converted
$output_str = $home_dir . "/fixpdb4amber.pl amber-" . $pdb_file_nm;
system $output_str;

#Prepare the structure for sander simulation

$output_str = $home_dir . "/prep_sander.pl amber-" . $pdb_file_nm . " $mol_name $home_dir";
system $output_str;

print "\n     Results\n";
print "--------------------\n\n";

$topo_file = $mol_name . ".top";
$coord_file = $mol_name . ".crd";

print "The pdb file: $pdb_file_nm was created\n";
print "The following files were created for simulation with sander:\n";
print "Topology File: $topo_file\n";
print "Coordinate File: $coord_file\n";
print "The following pictures were created: pic_" . $mol_name .".png,  pic_cpk_" . $mol_name . ".png\n";

$delete_Str = "rm -fr " . $mol_name. "*script " . $mol_name. "*_Strand*.pdb "  . $mol_name. "_cross-spec.txt";
system $delete_Str;

print "\n\n-===End Program===-\n\n";

if (-e $home_dir . "/submit_job.pl") {
    $run_now = GetYesNo("Would you like to start sander simulation now? [y/n]: ");
    if ($run_now) {
	$spec_machine = 0;
	$spec_machine  = GetYesNo("Do you want to run sander on a specific machine? [y/n]: ");
	$machine_name = "";
	if ($spec_machine) {
	    $valid_name =0;
	    while (! $valid_name) {
		$machine_name = input("Enter valid hostname [wolf or node]: ");
		$machine_name = lc($machine_name);
		if ($machine_name =~ /[wolf|node]\d+/) {
		    $valid_name = 1;
		}
	    }
	}
	$my_cmd = $home_dir . "/submit_job.pl $topo_file $coord_file ";
	$my_cmd .= ($helix_info[0]{"Base_Pairs"} * 2) . " $home_dir $machine_name &";
	system $my_cmd;
    }
}

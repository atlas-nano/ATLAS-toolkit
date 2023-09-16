#!/usr/bin/perl -w

# makecrossovers.pl - this program will create crossovers between two
#     strands of DNA molecules
# usage: makecrossovers.pl helix-1.pdb helix-2.pdb outfile.pdb crossover_list_file is5prime
# process: writes a script for namot2. namot2 creates the crossovers

# Variable Declaration Section

$helix1_loc = "";
$helix2_loc = "";
$output_fle = "";
$crossovers_file ="";
$periodicity = 0;
$is5prime = 0;

@crossover_data = (
		   {
		       "Molecule_1" => 0,
		       "Chain_1" => 0,
		       "Group_1" => 0,
		       "Molecule_2" => 0,
		       "Chain_2" => 0,
		       "Group_2" => 0
		   }
		   ); #hash containing all of the info for the crossover

@molecule_info = (
		  {
		      "Major_Groove" => 0,
		      "Minor_Groove" => 0,
		      "crossover_no" => 0, #number of crossovers
		      "basepair_no" => 0, #number of base pairs
		      "trans_vec" =>0.0 #distance apart that the molecules should be 
		  }
		  ); #hash that holds all the info about building the crossovers


sub file_test() {
# file_test - test to see if  all of the command line files are valid

   -e $helix1_loc or die "Error trying to access $helix1_loc\n";
   -e $helix2_loc or die "Error trying to access $helix2_loc\n";
   -e $crossovers_file or die "Error trying to access the crossover specification file $crossovers_file\n";

   $invalid_output_file = 0;
   -e $output_fle or $invalid_output_file = 1;
    if ($invalid_output_file) {
	print "The file $output_fle already exists.\n";
    }

    while ($invalid_output_file) {
	$overwrite = input ("Overwrite? [y, n]: ");
	$overwrite =~ s/\s+$//; #remove the final newline
	if ($overwrite =~ /^y|Y$/) {
	    last;
	}elsif ($overwrite =~ /^n|N$/) {
	    $output_fle = "";
	    while (!($output_fle =~ /\w+/)) {
		$output_fle = input ("New file name: ");
		chomp($output_fle);
		$output_fle =~ s/\s+$//;
	    }
	    $invalid_output_file = 0;
	    -e $output_fle or $invalid_output_file = 1;
	    if ($invalid_output_file) {
		print "The file $output_fle already exists.\n";
	    }
	}
     }
}

sub GetCrossovers() {
# GetCrossovers - opens the crossover confiruration file and extracts the required info
    
#Variable Declaration Section
    $line_in = "";
    $got_groove = 0; # specifies whether the major:minor groove data was obtained
    $got_cross_no = 0; # specifies whether the crossover number was obtained
    $valid_format = 0; # specifies whether the input file is a valid file
    $has_crossovers = 0; #specifies whether the input file has crossover data

# ---Start==-

    open CROSSFILE, $crossovers_file;
    while (<CROSSFILE>) {
	$line_in = $_;
	if (($line_in =~ /PX (\d+):(\d+)/) and !$got_groove) { # major groove: minor groove
	    $molecule_info[0]{"Major_Groove"} = $1;
	    $molecule_info[0]{"Minor_Groove"} = $2;
	    $got_groove = 1;
	}elsif (($line_in =~ /Total Crossovers:(\d+) in (\d+) bases, Transalation:(\d+)/) and !$got_cross_no) { 
# get the number of crossovers
	    $molecule_info[0]{"crossover_no"} = $1;
	    $molecule_info[0]{"basepair_no"} = $2;
	    $molecule_info[0]{"trans_vec"} = $3;
	    $periodicity = $molecule_info[0]{"Major_Groove"} + $molecule_info[0]{"Minor_Groove"};
	    $got_cross_no =1; 
	}elsif ($line_in =~ /CROSSOVERS/) {
	    $valid_format = 1; # file must have word "CROSSOVERS"
	}elsif ($line_in =~ /^(\d+):(\d+)->(\d+):(\d+)$/) { # valid crossover spec
	    # since this is a crossover vector, send the results to be processed
	    RegisterCrossover($1,$2,$3,$4);
	    $has_crossovers = 1; 
	}
    }
    close CROSSFILE;
    if (!$got_groove or !$got_cross_no or !$valid_format or !$has_crossovers) {
	die "Invalid file format in $crossovers_file. Terminating execution..\n";
    }

# -==End==-
}

sub GetChain($) {
# GetChain - gets the chain of a molecule

    my $periodicity = $molecule_info[0]{"Minor_Groove"} + $molecule_info[0]{"Major_Groove"};

    $curr_group = $_[0] - ($molecule_info[0]{"Minor_Groove"}/2);

    while ($curr_group > $periodicity) {
	$curr_group -= $periodicity;
    }

    #print "5Prime " . $is5prime . "\n";
    if ($curr_group <= $molecule_info[0]{"Minor_Groove"}) {
	if ($is5prime) {
	    $returnval = 1;
	} else {
	    $returnval = 2;
	}
    }else {
	if ($is5prime) {
	    $returnval = 2;
	} else {
	    $returnval = 1;
	}
    }
    
    return $returnval;
}

sub RegisterCrossover(@) {
# RegisterCrossover - pushed the link info into the crossover data hash

($mol1, $group1, $mol2, $group2) = @_;



# now to determine which the chains for the corresponding group

$chain1_no = 0;
$chain2_no = 0;

$complete_array = {
	"Molecule_1" => $mol1,
	"Chain_1" => $chain1_no,
	"Group_1" => $group1,
	"Molecule_2" => $mol2,
	"Chain_2" => $chain2_no,
	"Group_2" => $group2
	};

push @crossover_data, $complete_array;
}

sub WriteNamot2Script() {
# WriteNamot2Script - writes the contents of the molecule_data and crossover_info hashes into
# a compatible namot2 script
    $transvector = 20;
    $script_fle = $output_fle . ".script";
    if (!open(SCRIPTFLE, "> $script_fle")) {
	die "Cannot create script file $script_fle, $!\n";
    } else {
	print SCRIPTFLE "set WCad on\n";
	print SCRIPTFLE "set WCadFirst on\n";
	print SCRIPTFLE "load pdb na $helix1_loc\n";
	print SCRIPTFLE "load pdb na $helix2_loc\n";
	print SCRIPTFLE "trans 2  0  " . $transvector . "  0.0\n";
	print SCRIPTFLE "render\n";
	
	for ($i=1;$i<=$#crossover_data; $i++) {
	    $chain1_group = $crossover_data[$i]{"Group_1"};
	    $chain2_group = $crossover_data[$i]{"Group_2"};
	    
	    $ischain1 = GetChain($chain1_group);
	    #print "Called GetChain with $chain1_group, got $ischain1\n";
	    if ($ischain1 == 1) {
		print SCRIPTFLE "\n# Crossover #$i\n";
		print SCRIPTFLE "nick h" . $crossover_data[$i]{"Molecule_1"} . ":" .  $chain1_group .  ":1\n";
		print SCRIPTFLE "render\n";
		print SCRIPTFLE "nick h" . $crossover_data[$i]{"Molecule_2"} . ":" .  $chain2_group . ":1\n";
		
		print SCRIPTFLE "render\n";
		print SCRIPTFLE "link h1:" .  $chain1_group . ":1 h2:" .  ($chain2_group + 1) . ":1\n";
		print SCRIPTFLE "render\n";
		print SCRIPTFLE "link h2:" .  $chain2_group . ":1 h1:" .  ($chain1_group + 1) . ":1\n";
		print SCRIPTFLE "render\n";
	    } else {
		print SCRIPTFLE "\n# Crossover #$i\n";
		    $chain1_group += 1;
		$chain2_group += 1;
		
		print SCRIPTFLE "nick h" . $crossover_data[$i]{"Molecule_1"} . ":" .  $chain1_group .  ":2\n";
		print SCRIPTFLE "render\n";
		print SCRIPTFLE "nick h" . $crossover_data[$i]{"Molecule_2"} . ":" .  $chain2_group . ":2\n";
		
		print SCRIPTFLE "render\n";
		print SCRIPTFLE "link h1:" .  $chain1_group . ":2 h2:" .  ($chain2_group -1) . ":2\n";
		print SCRIPTFLE "render\n";
		print SCRIPTFLE "link h2:" .  $chain2_group . ":2 h1:" .  ($chain1_group -1). ":2\n";
		print SCRIPTFLE "render\n";
	    }
	}
	$pic_fle = $output_fle . ".png";
	print SCRIPTFLE "\nwrite pdb $output_fle\n";
	print SCRIPTFLE "write amber amber-" . $output_fle . "\n";
	print SCRIPTFLE "# This will color the different strand of the molecule\n";
	print SCRIPTFLE "set depth_cueing on\n";
	print SCRIPTFLE "render\n";
	print SCRIPTFLE "set color m1:1:*:* blue\n";
	print SCRIPTFLE "render\n";
	print SCRIPTFLE "set color m1:2:*:* yellow\n";
	print SCRIPTFLE "render\n";
	print SCRIPTFLE "set color m2:1:*:* red\n";
	print SCRIPTFLE "render\n";
	print SCRIPTFLE "set color m2:2:*:* green\n";
	print SCRIPTFLE "render\n";
	print SCRIPTFLE "set color m2:2:*:* green\n";
	print SCRIPTFLE "render\n";
	print SCRIPTFLE "set background white\n";
	print SCRIPTFLE "render\n";
	print SCRIPTFLE "# Create Pictues...\n";
	print SCRIPTFLE "write png pic_" . $pic_fle . "\n";
	print SCRIPTFLE "set space *\n";
	print SCRIPTFLE "render\n";
	print SCRIPTFLE "write png pic_cpk_" . $pic_fle . "\n";
	print SCRIPTFLE "close\n";
    }
    close SCRIPTFLE;
    print "\nGenerating crossovers...";
}

sub input($) {
    $printstring = $_[0];
    print "$printstring";
    $returnval = <STDIN>;
    return $returnval;
}


# -== Start Program ==-


# first ensure that the right arguments were passed
if (!@ARGV or !$ARGV[1] or !$ARGV[2] or !$ARGV[3] or $ARGV[4] eq "") {
    die "usage: makecrossovers.pl helix-1.pdb helix-2.pdb outfile.pdb crossover-spec-file\n";
}

# extract command line parameters
($helix1_loc, $helix2_loc, $output_fle, $crossovers_file, $is5prime) = @ARGV;


if (! $is5prime =~ /^[y|n]/) {
    die "Invalid value for is5prime, expected y or n\n";
} else {
    if ($is5prime =~ /y/) {
	$is5prime = 1;
    } else {
	$is5prime = 0;
    }
}

# test the files
file_test();

# open the crossover specification file and read the get the values
GetCrossovers();

# Write out the results
WriteNamot2Script();

# -== End Program ==-

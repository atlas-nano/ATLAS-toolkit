#!/usr/bin/perl -w
use strict;

# createdna.pl
# Perl Script used to create helixes of the PX molecules of Ned Seaman
# application: create dna with varying numbers in the major and minor groves
# process: utilizes namot2 to create these structures
# use: createdna.pl [dnasequencefile] outputfile.pdb [writevals]

sub GetBaseNames();
sub CreateNamotScript();
sub input;
sub Writefile();
sub Writefile();
sub GetUserInput();

# Variable Declaration Section
my ($major_no, $minor_no, @base_names, @bases_output, $writevals, @twist_angle_data);
my ($invalidfile, $isoverwrite, $temp1);

$writevals = $major_no = $minor_no = 0;

#-==Start-==
# First check to see if the outputfile was specified
if (!@ARGV) {
    die "usage: $0 outputfile.pdb [scripts_dir] [dnasequencefile]\n";
}

# Then check to see if a file was specified, if so then open it
my ($outputfile, $scripts_dir, $sequencefile) = @ARGV;

if (defined($sequencefile)) {
    GetBaseNames();
} else { # if no file was specified, start to read in the bases from cmd line
    GetUserInput();
}

if ($minor_no > $major_no) { # if the user swithed the major and minor grooves
    $temp1 = $minor_no;
    $major_no = $temp1;
    $minor_no = $major_no; # swap then
}

# Determine if the writevals parameter was passed

#if ($ARGV[2] or ($ARGV[1] eq "1")) {
    $writevals =1;
#}

# Then create the entries for the output file
    CreateNamotScript();

# Validate the outputfile and Write the Script

$invalidfile =1;
while($invalidfile){
    if (open(OUTPUTFILE, $outputfile)) {
	$isoverwrite = input("The file: $outputfile already exists. Overwrite? [y/n] ");
	if ($isoverwrite =~ /^y|Y$/) {
	    $invalidfile = 0;
	    Writefile();
	} elsif($isoverwrite =~ /^n|N$/)  {
	    $outputfile = input("Enter new filename: ");
	    $invalidfile =1;
	} 
    }elsif ($outputfile =~ /^\w+/) {
	$invalidfile = 0;
	Writefile();
    } else {
	while (!($outputfile =~ /^\w+/)) {
	    $outputfile = input ("Enter new filename: ");
	}
	$invalidfile = 1;
    }
    close OUTPUTFILE;
}


sub GetBaseNames() {
# opens the base sequence file, reads in the bases
# then places the sequence in the $bases_output array
# format of base sequence file:
# major groove: minor groove ratio e.g. 6:5
# name of neucleotide e.g. at
# name of next neucleotide
# ....

    my ($majorgroove, $minorgroove);
    $minorgroove = $majorgroove = 0;

    open SEQUENCE, $sequencefile or die "Cannot open $sequencefile. $!\n";
    print "Reading Sequences from $sequencefile...";
    while (<SEQUENCE>) {
	if ($_ =~ /^(\d+):(\d+)/) { # first line so get vals for grooves
	    $majorgroove = $1;
	    $minorgroove = $2;
	    if (4>$majorgroove or 10<$majorgroove) {
		die "Error in file: invalid major groove. Got $majorgroove, expected 5 - 9\n";
	    }
	    if (4>$minorgroove or 10<$minorgroove) {
		die "Error in file: invalid minor groove. Got $minorgroove, expected 5 - 9\n";
	    }
	} else { # assume that it is a sequence
	    $_ = lc($_);
	    if ($_ =~ /^(at|ta|gc|cg)/) {
		push @base_names, $1;
	    } elsif ($_ =~ /\w+/) {
		print "Skipping Invalid record: $_ ...\n";
	    }
	}
    }
    
    close SEQUENCE;
    die "Invalid file format. No record found for major: minor grooves\n"    
	if ($majorgroove eq "" or $minorgroove eq "");
    print "Done.\nFound " . ($#base_names + 1) . " base-pairs for $majorgroove:$minorgroove PX molecule.\n";

    $major_no = $majorgroove;
    $minor_no = $minorgroove;
}

sub GetYesNo($) {
# This sub will ask a yes/no question and get the correct answer
    my ($valid_ans, $ans, $my_reply);

    while (! $valid_ans) {
	$ans = input($_[0]);
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
sub GetUserInput() {
# GetUserInput - allows the user to input the base sequence from the command line
    my ($majorgroove, $minorgroove, $numbervalid, $legit_bases, $bp_sequence);
    my (@temp_names, $count, $basepair, $i, $tempholder, $my_ans);

    $majorgroove = $minorgroove = $numbervalid = 0;

     while (!$numbervalid) {
	$majorgroove = input("No. of bases in major groove [4-9]: ");
	chomp($majorgroove);
	if ($majorgroove =~ /(\d+)/) {
	    if ($1 < 4 or $1> 9) {
		print "Invalid value\n";
		$numbervalid = 0;
	    } else {
		$numbervalid = 1;
	    }
	} else {
	    print "Invalid value. Expected Integer\n";
	    $numbervalid = 0;
	}
    }

    $numbervalid = 0;
    while (!$numbervalid) {
	$minorgroove = input("No. of bases in minor groove [4-9]: ");
	chomp($minorgroove);
	if ($minorgroove =~ /(\d+)/) {
	    if ($1 < 4 or $1> 9) {
		print "Invalid value\n";
		$numbervalid = 0;
	    } else {
		$numbervalid = 1;
	    }
	} else {
	    print "Invalid value\n";
	    $numbervalid = 0;
	}
    }
   
    $major_no = $majorgroove;
    $minor_no = $minorgroove;

    $legit_bases = 0;

    while ($legit_bases == 0) {
	print "\nEnter two-letter base pairs sequence below delimited by a space\n";
	$bp_sequence = <STDIN>;
	chomp($bp_sequence);

	$bp_sequence = lc($bp_sequence);

	@temp_names = split " ", $bp_sequence;
	$count = $#temp_names +1;

	for ($i=0; $i<$count; $i++) {
	    $basepair = $temp_names[$i];
	    if ($basepair =~ /(at|ta|gc|cg)/) {
		$tempholder = $1;
		push @base_names, $tempholder;
		$legit_bases += 1;
	    }
	}
	if ($legit_bases > 0 and $legit_bases > $minorgroove) {
	    $my_ans = GetYesNo("Recognised $legit_bases base-pairs. Continue [y/n]? ");
	    if (! $my_ans) {
		$legit_bases = 0;
	    }
	} else {
	    if ($legit_bases > 0) {
		print "Please enter at least $minorgroove legit base-pairs\n";
		$legit_bases = 0;
	    }
	}
    }
}

sub CreateNamotScript() {
# CreateNamotScript - creates the lines to be written in the script file to be executed by namot2
# populates the @basesoutput array with the info in @bases_names
# the number of atoms before the first turn is dependant on the major groove
# i.e. 4 in the major groove will have 1 atom before the first turn
# 5 & 6 will have 2, 7 & 8 will have 3 and 9 will have 4

# Variable Declaration Section
    my ($j, $scriptloc, $atoms_b4_turn, $totalsequences, $curr_sequence_no, @sequence_twist);
    my ($invalid_input, $next_groove_end, $i, $transalation);

    $atoms_b4_turn = $totalsequences = $curr_sequence_no = 0;

    $totalsequences = $#base_names +1;
    $sequence_twist[0] = "45";
    $sequence_twist[1] = "36";
    $sequence_twist[2] = "30";
    $sequence_twist[3] = "25.714286";
    $sequence_twist[4] = "22.5";
    $sequence_twist[5] = "20";
#    $major_angle = $sequence_twist[$major_no-4]; #twist angle for the major groove
#    $minor_angle = $sequence_twist[$minor_no-4]; #twist angle for the minor groove

# now determine the number of atoms before the first turn
    if (! defined($scripts_dir)) {
	$scripts_dir = $ENV{PWD};
    }
    
    $scriptloc = $scripts_dir . "/function.pl";
    -e $scriptloc or die "Cannot locate file $scriptloc, so cannot calculate twist angles\n";

    open TWISTANGLE, "$scriptloc $minor_no $major_no |" or die "Cannot execute function.pl";
    while (<TWISTANGLE>) {
	chomp;
	push @twist_angle_data, $_;
    }
    close TWISTANGLE;

    $invalid_input = 1;

    while ($invalid_input) {
	$atoms_b4_turn = input("\nNo. of base pairs in first minor groove: ");
	chomp($atoms_b4_turn);
	if ($atoms_b4_turn =~ /^[4-9]$/) {
	    $invalid_input = 0;
	    print "\n";
	}
    }

    if ($atoms_b4_turn > $totalsequences) {
	$atoms_b4_turn = $totalsequences;
    }

    $j = 0;
    for $i (0 .. $atoms_b4_turn) { # place the atoms before the first turn
	CreateLine($base_names[$i], $twist_angle_data[$j]);
	$j++;
    }
    $curr_sequence_no = $atoms_b4_turn +1;
    while ($curr_sequence_no <= $totalsequences) {

	$next_groove_end = $curr_sequence_no + $major_no;
	if ($next_groove_end > $totalsequences) {
	    $next_groove_end = $totalsequences;
	}
	$j = $minor_no;
	for ($i=$curr_sequence_no; $i<($next_groove_end); $i++) {
	    CreateLine($base_names[$i], $twist_angle_data[$j]); # create major groove
	    $j++;
	}

        # print "Finished minor groove: $minor_no. CurrPos: $next_groove_end of $totalsequences\n";

	if ($next_groove_end < $totalsequences) {
	    $curr_sequence_no += $major_no;
	    $next_groove_end = $curr_sequence_no + $minor_no;
	    if ($next_groove_end > $totalsequences) {
		$next_groove_end = $totalsequences;
	    }
	    
	    $j = 0;
	    for ($i=$curr_sequence_no; $i<($next_groove_end); $i++) {
		CreateLine($base_names[$i], $twist_angle_data[$j]); # create minor groove
		$j++;
	    }
	    
	    $curr_sequence_no += $minor_no;
	    
	   # print "Finished major groove: $major_no. CurrPos: $next_groove_end of $totalsequences\n";
	} else {
	    last;
	}
    }
# if writevals was called then write the parameters of this helix to a file


if ($writevals) {
    $transalation = DetermineTransalation();
    $totalsequences -= 1;
    PrintParameters($totalsequences, $atoms_b4_turn, $transalation);
}

}

sub CreateLine(@) {
# CreateLine - create a line compatible to that of a namot2 script for adding a unit
# format:
# add unit [base_name]-wc B 0 0 0 0 0.000000 0.000000 [rise] 0.000000 [twist] 0.000000
# eg:
# add unit at-wc B 0 0 0 0 0.000000 0.000000 3.400000 0.000000 30.000000 0.000000

    my ($basenm, $twistangle) = @_;
    my ($helix_rise, $returnstr);

#print "$basenm, $twistangle\n";
    $helix_rise = "3.380000";
    $returnstr = "add unit " . $basenm . "-wc B 0 0 0 0 0.000000 0.000000 ". $helix_rise . " 0.000000 " . $twistangle . ".000000 0.000000";
#print "$returnstr\n";
    push @bases_output, $returnstr;
}

sub Writefile() {
# writes the lines in $bases_output into the file $outputfile
    my ($scriptfile, $j);

    $scriptfile = $outputfile . ".script";
    if (!open(OUTFILE, "> " . $scriptfile)) {
	die "Error when trying to write to $outputfile. $!\n";
    } else {
	print OUTFILE "# script to create PX $major_no:$minor_no DNA Double Crossover Molecule\n";
	# print "$#bases_output\n";
	for $j (0 .. $#bases_output) {
	    print OUTFILE $bases_output[$j] . "\n";
	    print OUTFILE "render\n";
	   # print "$j = $bases_output[$j]\n";
	}
	print OUTFILE "write pdb " . $outputfile . "\n";
	#print OUTFILE "write amber amber-" . $outputfile ."\n";
	print OUTFILE "close\n";
	#print OUTFILE "quit\n";
	print "Generating helix....";
    }
    close OUTFILE;

}


sub PrintParameters(@) {
# Prints out the parameters of the helix to a file
    my ($para_fle);

    $para_fle = $outputfile . ".parm";
    if (open(PARAFILE, "> " . $para_fle)) {
	print PARAFILE "PX " . $major_no . ":" . $minor_no . "\n";
	print PARAFILE $_[0] . "," . $_[1] . "," . $_[2] . "\n";
	close PARAFILE;
    }else {
	print "\nError! Could not create parameter file. $!\n";
    }
}
		
sub DetermineTransalation() {
# determines how far apart the helixes should be
# this is determined by the number of base pairs in the major groove

    my (@transalation_vector) = (
				 {
				     "4" => 23.5,
				     "5" => 23.5,
				     "6" => 23.5,
				     "7" => 23.5,
				     "8" => 23.5,
				     "9" => 23.5
				     }
				 );
    return $transalation_vector[0]{$major_no};
    
}
sub input {
    my ($printstring) = $_[0];
    my ($returnval);
                                                                                                                  
    print "$printstring";
    $returnval = <STDIN>;
    return $returnval;
}
                                                                                                                  


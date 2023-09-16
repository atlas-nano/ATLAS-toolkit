#!/usr/bin/perl -w

# generateCrossovers.pl - automates the process of generating dna crossover molecules
#  allows creation of two independant helixes of b-dna and allows them to be
#  crossed over a specific points.
# note: it is very important when specifying the crossover points to be
#	   very careful. Bad specifications will lead to bad structures
# process: utilizes a number of small perl scripts
# usage: dna_cross_gen.pl crossed-dna-name.pdb [scripts-dir]

# Load the Perl5 Namot2 module
BEGIN {
	unshift @INC, "/home/tpascal/scripts/";
}
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use Packages::Namot2::Namot;
use File::Basename;
use Cwd qw(realpath);
use strict;

p5namot::Cmd("set hush ERROR off");
p5namot::Cmd("set hush INFO off");
p5namot::Cmd("set hush REQUESTED off");
p5namot::Cmd("set hush WARNING off");

#prototype subroutines
sub CreateAHelix;
sub Check_for_Scripts;
sub GetCrossoverPoints;
sub CreateSpecsFile;
sub input;
sub Validate;
sub RunScript;
sub SubmitJob;
sub GetMotif;
sub GetAns;
sub Trim;

# Variable declaration Section
my ($can_fixpdb, $numHelix, $myCrossovers, $Helix, $crossType, $sLen, $counter, $bp);
my ($mol_name, $mol1, $mol2, $tmp, $tmpFiles, $topo_file, $coord_file, $crossovers);

die "usage: $0 saveName [scripts-dir]\n"
	if (!@ARGV);
my ($output_file, $home_dir) = @ARGV;

$can_fixpdb = 1; # specifies whether the fixpdb.pl file is present
$bp = 0;

# validate cmd line argument
Validate;
print "\nWelcome to the DNA Crossover builider.\nThis script will allow you to create super structures " .
	"of DNA double helices.\nHelices are arranged from left to right and can be of either 5'->3' or " .
	"3' -> 5' topology\n(from the bottom) and of any crossover type (PX or DX)\n\n";

print "\n-===Starting Program===-\n";
print "------------------------\n\n";

# Get options
($numHelix, $crossType, $Helix) = GetMotif;

# Create helixes
print "\nCreating helices.\n-----------------";
for $counter (1 .. $numHelix) {
	print "\nHelix $counter\n";
	$Helix->[($counter - 1)]{"INFO"} = CreateAHelix("Strand" . $counter);
}

#Crossover helixes
print "\nEnter crossover points from bottom up (-99 to stop) for each helical pair\n";
print "----------------------------------------------------------------------------\n";
if ($crossType == 3) {
	print "Enter the crossover type after the crossover [1 - 2]\n";
}

for $counter (1 .. ($numHelix  - 1)) {
	$sLen = $Helix->[($counter - 1)]{"INFO"}{"TotalBases"};
	if ($Helix->[$counter]{"INFO"}{"TotalBases"} < $sLen) {
		$sLen = $Helix->[$counter]{"INFO"}{"TotalBases"};
	}
	$crossovers->{$counter} = GetCrossoverPoints($counter, $crossType, $sLen);
}

#rotate the helixes and get the appropriate angles
print "\nOptimizing Helical distances\n";
print "--------------------------------\n";
for $counter (1 .. ($numHelix  - 1)) {
	print "  HELIX $counter <-> " . ($counter + 1) . "..";
	$mol1 = \%{ $Helix->[($counter - 1)] };
	$mol2 = \%{ $Helix->[($counter)] };
	
	$tmp = CreateSpecsFile($crossovers->{$counter}, $mol1, $mol2);
	$tmpFiles .= $mol1->{"INFO"}{"Name"} . " " . $mol2->{"INFO"}{"Name"} . " $tmp";
	
	RunScript("find_o_angle.pl " . $mol1->{"INFO"}{"Name"} . " " . $mol2->{"INFO"}{"Name"} . 
		  " " . $tmp . " " . ($counter - 1) . " $home_dir");
}

#Create Crossovers
print "\nCreating Crossovers\n";
print "---------------------\n";
for $counter (1 .. ($numHelix  - 1)) {
	print "  HELIX $counter <-> " . ($counter + 1) . "..";
	$mol1 = \%{ $Helix->[($counter - 1)] };
	$mol2 = \%{ $Helix->[($counter)] };
	
	$tmpFiles .= " $tmp";
	$tmp = CreateSpecsFile($crossovers->{$counter}, $mol1, $mol2);
	RunScript("makecrossovers.pl " . $mol1->{"INFO"}{"Name"} . " " . $mol2->{"INFO"}{"Name"} . 
		  " " . $mol_name . ".pdb " . $tmp);

	$tmp = "execute " . $mol_name . ".pdb.script";
	p5namot::Cmd($tmp);
	print "Done\n";
}

print "\n	 Results\n--------------------\n\n";
print "The pdb file: $mol_name" . ".pdb was created\n";

print "The following pictures were created: " . $mol_name . ".png,  " . $mol_name . "_cpk.png\n";

#RunScript("rm -fr " . $mol_name. "*script " . $mol_name. "*Strand*.pdb "  . $mol_name. "_cross-spec.txt");
print "\n\n-===End Program===-\n\n";

if ($can_fixpdb) {
	print "\n---Creating Amber Files----\n";
# Fix Pdb File to be converted
	print "\nFixing Base Names...";
	RunScript("fixpdb4amber.pl $mol_name"  . "_amber.pdb");

#Prepare the structure for sander simulation
	RunScript("prep_sander.pl $mol_name" . "_amber.pdb $mol_name $home_dir");
	$topo_file = $mol_name . ".top";
	$coord_file = $mol_name . ".crd";
	print "The following files were created for simulation with sander:\n";
	print "Topology File: $topo_file\n";
	print "Coordinate File: $coord_file\n";

	SubmitJob;
}


sub SubmitJob {
	my ($run_now, $spec_machine, $machine_name, $valid_name, $my_cmd);

	if (-e $home_dir . "/submit_job.pl") {
	$run_now = GetAns("Would you like to start sander simulation now? [y/n]: ", "Y", "N");
	if ($run_now) {
		$spec_machine = 0;
		$spec_machine  = GetAns("Do you want to run sander on a specific machine? [y/n]: ", "Y", "N");
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
		$my_cmd = "submit_job.pl $topo_file $coord_file ";
		$my_cmd .= $bp . " $home_dir $machine_name &";
		RunScript($my_cmd);
	}
	}
}

sub CreateAHelix {
	my ($whichstrand) = $_[0];
	my ($output_string, $has_sequence, $file_not_found, $sequence1);
	my ($incoming_line, $kill_str, $execute_Str, $return_Str, $param_fle, $rec);

	$return_Str = $mol_name . "_" . $whichstrand . ".pdb";
	$output_string = "createdna.pl " . $return_Str . " "; 

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
	$output_string .= " $home_dir " . $sequence1;
	}

	RunScript($output_string);
	
	$execute_Str = "execute " . $return_Str . ".script";
	p5namot::Cmd($execute_Str);
	print "Done\n";


# check for the parameter file - it holds the helix info
	$param_fle = $return_Str . ".parm";
	if (!open(PARAMFILE, $param_fle)) {
	die "Cannot open molecule paramater file. $!\nExecution Terminated.\n";
	} else {
	while (<PARAMFILE>) {
		$incoming_line = $_;
		if ($incoming_line =~ /PX (\d+):(\d+)/) { # first line
		$rec->{"Major_Groove"} = $1;
		$rec->{"Minor_Groove"} = $2;
		}elsif ($incoming_line =~ /(\d+),(\w+)/) { #second line
		$rec->{"TotalBases"} = $1;
		$bp += $1;
		$rec->{"Minor1"} = $2;
		}
	}
	close PARAMFILE;
	}

# remove the paramater file
	$kill_str = "rm " . $param_fle;
	RunScript($kill_str);

	$rec->{"Name"} = $return_Str;
	return $rec;
}


sub Check_for_Scripts {
	my ($file_loc);

# critical files	
	for $file_loc ("createdna.pl", "makecrossovers.pl", "find_o_angle.pl", "rotatehelix.pl", "function.pl") {
	die "Critical error! Unable to locate $file_loc: $!.\nExecution Terminated.\n"
		if (! -e ($home_dir . "/" . $file_loc));
	}
	
#non critical files
	
	for $file_loc ("fixpdb4amber.pl", "prep_sander.pl") {
	if (! -e ($home_dir . "/" . $file_loc)) {
		warn "WARNING: Unable to locate $file_loc: $!.\nWill not run convert structures to AMBER.\n";
		$can_fixpdb = 0;
		last;
	}
	}

}

sub GetCrossoverPoints {
	my ($helixPair, $ct, $maxVal) = @_;
	my ($isValid, $lineInput, @cross, $rec, @tmp, $lineCounter);

	print "  HELIX $helixPair <-> " . ($helixPair+1) . "\n";

	$isValid = $lineCounter = 1;
	while ($isValid) {
	$lineInput = Trim(input("\t$lineCounter : "));
	next
		if ($lineInput !~ /\d+/);
	last
		if ($lineInput eq "-99" and $lineCounter > 1);
	@tmp = split /\s+/, $lineInput;
	if ($ct < 3 or $#tmp == 1) {
		if ($tmp[0] =~ /(\d+)/ ) {
		if ($1 > 1 and $1 < $maxVal) {
			$rec = (
				{
				"TYPE" => $ct,
				"VAL"  => $1,
				}
				);
			push @cross, $rec;
			$lineCounter++;

		}
		}
	}
	}
	return \@cross;

}

sub CreateSpecsFile {
	my ($crossInfo, $helix1, $helix2) = @_;
	my ($outStr, $trans, $isParallel, $is5prime, $specFile);

	$specFile = $mol_name . "_crossSpec.txt";

	$outStr = "Helix 1 " .
	$helix1->{"INFO"}{"Major_Groove"} . ":" . $helix1->{"INFO"}{"Minor_Groove"} .
	" " . $helix1->{"INFO"}{"TotalBases"} . " bases 1st Minor: " . $helix1->{"INFO"}{"Minor1"} . "\n";
	$outStr .= "Helix 2 " .
	$helix2->{"INFO"}{"Major_Groove"} . ":" . $helix2->{"INFO"}{"Minor_Groove"} .
	" " . $helix2->{"INFO"}{"TotalBases"} . " bases 1st Minor: " . $helix2->{"INFO"}{"Minor1"} . "\n";
	
	$trans = "0 18.5 0";
	$outStr .= "Transalation: " . $trans . "\n";
	$isParallel = $is5prime = 0;
	$isParallel = 1
	if ($helix1->{"TOPO"} != $helix2->{"TOPO"});
	$is5prime = 1
	if ($helix1->{"TOPO"} == 1);
	$outStr .= "Parallel: $isParallel\n5PrimeToLeft: $is5prime\nCROSSOVERS\n";
	
	for $counter (0 .. $#{ $crossInfo }) {
	$outStr .= "1:" . $crossInfo->[$counter]{"VAL"} . "->2:" . $crossInfo->[$counter]{"VAL"} .
		" TYPE:" . $crossInfo->[$counter]{"TYPE"} . "\n";
	}
	open SPECSFLE, "> $specFile" or die "Cannot create crossover specification file $specFile: $!\n";
	print SPECSFLE $outStr;
	close SPECSFLE;

	return $specFile;
}


sub input {
	my ($printstring) = $_[0];
	my ($returnval);

	print "$printstring";
	$returnval = <STDIN>;
	return $returnval;
}

sub Validate {
	my ($file_exists, $overwrite, $invalidname);

	$file_exists = 0;
	$file_exists = 1
	if (-e $output_file);
	if (! defined($home_dir)) {
	$home_dir = realpath($0);
		$home_dir =~ s/\/\w+\.pl$//;
		
	}

	$home_dir =~ s/\/$//;

	while ($file_exists) {
	if (-e $output_file) {
		$overwrite = GetAns("The file $output_file already exists. Overwrite? [y/n]: ", "Y","N");
		if ($overwrite) {
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
	} else {
		$file_exists = 0;
		last;
	}
	}

	$mol_name = basename($output_file);
	
	if ($output_file =~ /\.pdb$/) {
	$mol_name =~ s/\.pdb$//;
	}

# Check for scripts
	Check_for_Scripts();
	
}

sub RunScript {
	my ($myCmd) = $_[0];
	my ($location);

	if ($myCmd !~ /\.pl/) {
	$location = "";
	} elsif ($home_dir eq $ENV{PWD}) {
	$location = "./";
	} else {
	$location = $home_dir . "/";
	}

	$myCmd = $location . $myCmd;
#	print $myCmd;
	die "ERROR while running $myCmd... Terminating\n" 
	if (system($myCmd));

}

sub GetMotif {
	my ($totHelix, @topo, $typeCross, $counter);

	$totHelix = GetAns("Number of Double helical structures [2 - 10]: ", 2, 10);

	print "\nHelical topology (at bottom, from left to right):\n" .
		"\t1. 5' -> 3'\n" .
		"\t2. 3' -> 5'\n";

	for $counter (1 .. $totHelix) {
	$topo[($counter - 1)]{"TOPO"} = GetAns("Helix #$counter: ", 1, 2); 
	}

	print "\nType of Crossovers:\n" .
		"\t1. All PX like (\"X\" pattern)\n" .
		"\t2. All DX like (loop pattern)\n" .
		"\t3. Mixed\n";

	$typeCross = GetAns(": ", 1, 3);
	
	return ($totHelix, $typeCross, \@topo);
}

sub GetAns {
	my ($question, $start, $end) = @_;
	my ($isChar, $answer, $isValid);

	$isChar = $isValid = 0;
	
	if (lc($start) eq "y") {
	$start = 1;
	$end = 0;
	$isChar = 1;
	} else {
	$start--;
	$end++;
	}


	while (! $isValid) {
	$answer = Trim(input($question));
	if ($isChar) {
		if (lc($answer) eq "y") {
		$answer = 1;
		$isValid = 1;
		} elsif (lc($answer) eq "n") {
		$answer = 0;
		$isValid = 1;
		}
	} else {
		if ($answer =~ /(\d+)/) {
		if ($1 > $start and $1 < $end) {
			$isValid = 1;
			$answer = $1;
		}
		}
	}
	}
	
	return $answer;
}

sub Trim {
	my ($inStr) = $_[0];

	$inStr =~ s/^\s+//;
	$inStr =~ s/\s+$//;

	return $inStr;
}

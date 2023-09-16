#!/usr/bin/perl -w

# This program will rotate both angles and find the optimum angle

# usage: find_o_angle.pl helix1 helix2 crossoverpoint is5prime scriptdir

# check for arguments
$h1_angle = 0;
$h2_angle =0;

($helix1, $helix2, $c_point, $is5prime, $script_dir) = @ARGV;

if (! $helix1 or ! $helix2 or ! $c_point or ! $script_dir or ! $is5prime) {
    die "usage: find_o_angle.pl helix1 helix2 crossoverpoint 5primein scriptdir\n";
}

-e $helix1 or die "Error: Cannot find $helix1\n";
-e $helix2 or die "Error: Cannot find $helix2\n";
-d $script_dir or die "Error: Cannot find $script_dir, $!\n";

if (! $c_point =~ /(\d+)/) {
    die("Invalid crossover point, exepcted integer\n");
}

$c_point =~ /(\d+)/;
$c_point = $1;

if (! $is5prime =~ /^[y|n]/) {
    die "Invalid value for is5prime, expected y or n\n";
} else {
    if ($is5prime =~ /y/) {
	$is5prime = 1;
    } else {
	$is5prime = 0;
    }
}

if ($is5prime) {
    $c_point +=1;
}

$shell_cmd = "";

sub ExtractInfo($) {

    my $which_helix = $_[0];
    my $is_second = 0;
    my $min_dist = 9999;

    if (! open INFILE, "results" . $which_helix . ".txt") {
	die "Cannot open output file\n";
    } else {
	while (<INFILE>) {
	    $instring = $_;
	    chomp($instring);
	    if ($instring =~ /^(\d+)\sREQUESTED:Distance\sis\s(\d+\.\d+)/) {
		$curr_angle = $1;
		$h_dist = $2;
		if ($h_dist < $min_dist) {
		    if ($which_helix ==1) {
			$h1_angle = $curr_angle;
		    } else {
			$h2_angle = $curr_angle;
		    }
		    $min_dist = $h_dist;
		}
	    }
	}
	close INFILE;
    }
}
print "\nFinding optimum rotation angles...\n";
print "----------------------------------\n";

$shell_cmd = $script_dir . "/getdist.pl $helix1 $helix2 $c_point 1 $is5prime > results1.txt";

print "Rotating Helix 1....";
system $shell_cmd;
print "Done\n";


ExtractInfo(1);


#second helix

$shell_cmd = $script_dir . "/getdist.pl $helix1 $helix2 $c_point 2 $is5prime > results2.txt";

print "Rotating Helix 2....";
system $shell_cmd;
print "Done\n";

ExtractInfo(2);

#rotate the first helix

$shell_cmd = $script_dir . "/rotatehelix.pl $helix1 $h1_angle > junk";
system $shell_cmd;

#now rotate the second helix

$shell_cmd = $script_dir . "/rotatehelix.pl $helix2 $h2_angle > junk";
system $shell_cmd;

system "rm -f results1.txt results2.txt junk";

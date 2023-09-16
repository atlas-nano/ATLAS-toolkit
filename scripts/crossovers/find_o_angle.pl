#!/usr/bin/perl -w
# This program will rotate both angles and find the optimum angle
# usage: find_o_angle.pl helix1 helix2 crossoverpoint is5prime scriptdir

use strict;

sub Initialize;
sub ExtractInfo;
sub RunScript;
sub getResults;

# check for arguments
die "usage: find_o_angle.pl helix1 helix2 specsFile rotateOnlyOne [scriptsDirectory]\n"
    if (! @ARGV or $#ARGV < 3);

my ($helix1, $helix2, $specsFile, $rotateOne, $script_dir) = @ARGV;
my ($dist, $h1_angle, $h2_angle, $converged, $new_dist, $tolerance);

Initialize;
$h1_angle = $h2_angle = $converged = $new_dist = 0;
$tolerance = 1;
$dist = 99999;

while (! $converged) {
# first helix
    if (! $rotateOne) {
	RunScript("getdist.pl $helix1 $helix2 1 $specsFile > results1.txt");
	($new_dist, $h1_angle) = ExtractInfo(1);
	$new_dist = $new_dist/4;
	RunScript("rotatehelix.pl $helix1 $h1_angle > junk");
#	print "distance: $new_dist angle: $h1_angle...";
    }
    
    $converged = 1
	if (($dist - $new_dist) < $tolerance);
    $dist = $new_dist;
    last if ($converged);
#second helix
    RunScript("getdist.pl $helix1 $helix2 2 $specsFile > results2.txt");
    ($new_dist, $h2_angle) = ExtractInfo(2);
    $new_dist = $new_dist/4;
    RunScript("rotatehelix.pl $helix2 $h2_angle > junk");
#    print "distance: $new_dist angle: $h2_angle...";
    $converged = 1
	if (($dist - $new_dist) < $tolerance);
    $dist = $new_dist;
    
    last
	if ($rotateOne);
}
print "Shortest distance: $new_dist...Done\n";

#system "rm -f results1.txt results2.txt junk";

sub Initialize {
    my ($tmp);

    if (! defined($script_dir)) {
	$script_dir = $ENV{PWD};
    }

    for $tmp ($helix1, $helix2, $specsFile) {
	die "Error accessing $tmp\n"
	    if (! -e $tmp and ! -d $tmp);
    }
}

sub ExtractInfo {

    my $which_helix = $_[0];
    my ($counter, $angle, $min_dist, $instring, $curr_dist);
    $counter = $angle = 0;
    $min_dist = 99999;

    open INFILE, "results" . $which_helix . ".txt" or die "Cannot open output file\n";
    while (<INFILE>) {
	$instring = $_;
	chomp($instring);
	if ($instring =~ /^REQUESTED:Distance\sis\s(\d+\.\d+)/) {
	    $curr_dist += $1;
	} elsif ($instring =~ /^REQUESTED:TER/) {
	    $counter++;
	    if ($curr_dist < $min_dist) {
		$angle = $counter;
		$min_dist = $curr_dist;
	    }
	    $curr_dist = 0;
	}
    }
    close INFILE;
    
    return ($min_dist, $angle);
}

sub RunScript {
    my ($myCmd) = $_[0];
                                                                                                               
    $myCmd = $script_dir . "/" . $myCmd;
    die "ERROR while running $myCmd... Terminating\n"
        if (system($myCmd));
                                                                                                               
}

sub getResults {
    my ($myFile) = $_[0];
    my ($dist, $angle, $isValid);

    $isValid = 0;
    $myFile = "results" . $myFile . ".txt";

    open RFILE, $myFile or die "Cannot open $myFile: $!\n";
    while (<RFILE>) {
	chomp;
	if ($_ =~ /Smallest distance: (\d+\.\d+)\s+angle: (\d+)/) {
	    $isValid = 1;
	    ($dist, $angle) = ($1, $2);
	}
    }
    close RFILE;

    die "Error parsing $myFile. No valid data found\n"
	if (! $isValid);

    return ($dist, $angle);

}


#!/usr/bin/perl -w
use strict;
use warnings;

sub TestForInt(@);
sub ValidateFiles();
sub CreateAvgStructure();
sub CalcRms;

if (! @ARGV or $#ARGV < 6) {
    die "usage: $0 filebase traj_file top_file starttime endtime interval strandlen [start_crd]\n";
}

my ($filebase, $traj_file, $top_file, $start_tm, $end_tm, $interval, $strandlength, $crd_file) = @ARGV;
my (@validfiles, $avg_structure);

-e $traj_file or die "Cannot find $traj_file: $!\n";
-e $top_file or die "Cannot find $top_file: $!\n";

TestForInt($start_tm, "start time");
TestForInt($end_tm, "end time");
TestForInt($interval, "interval");
TestForInt($strandlength, "strand length");

print "Step 1. Create average structure file\n";
print "-------------------------------------\n";
CreateAvgStructure();

print "\n";
print "Step 2. Calculating RMS using average structure as reference\n";
print "-----------------------------------------------------------\n";
CalcRms(1);
print "\n";

if (defined($crd_file) and -e $crd_file) {
print "\n"; 
print "Step 3. Calculating RMS using starting structure as reference\n";
print "-----------------------------------------------------------\n";
CalcRms(0); 
print "\n"; 
}

print "-==End==-\n";

sub CreateAvgStructure() {

    my ($out_cmd, $counter, $start_time);

    # use last 20 frames as average
    $start_time = ($end_tm - ($interval * 20));
    if ($start_time < 0) {
	$start_time = $start_tm;
    }
    print "Creating average file from frames $start_time - $end_tm (every $interval step(s))\n";
    $avg_structure = $filebase . "_avg.pdb";
    $out_cmd = "> tmp_gen_rst_file";
    open OUTCMD, $out_cmd or die "Cannot create file tmp_gen_rst_file: $!\n";
    print OUTCMD "trajin $traj_file $start_time $end_tm $interval\n";
    print OUTCMD "strip :Na+\nstrip :WAT\n";

    for $counter (1 .. 4) {
        print OUTCMD "center :1-" . ($counter * $strandlength) . " mass origin\n";
        print OUTCMD "image origin center\n";
    }

    print OUTCMD "average $avg_structure pdb\n";

    close OUTCMD;

    $out_cmd = "ptraj $top_file < tmp_gen_rst_file >& junk";

    print "Executing Ptraj....";

    if ( system "$out_cmd") {
	die "Failure!! Error executing ptraj\n";
    }

    print "Sucess.\nCreated $avg_structure as the average structure\n";
}

sub CalcRms {
    my ($isAVG) = $_[0];
    my ($out_cmd, $counter);

    $out_cmd = "> tmp_gen_rst_file";
    open OUTCMD, $out_cmd or die "Cannot create file tmp_gen_rst_file: $!\n";
    print OUTCMD "trajin $traj_file\n";
    print OUTCMD "strip :Na+\nstrip :WAT\n";

    if (! $isAVG) {
	print OUTCMD "reference $crd_file\n";
    } else {
        print OUTCMD "reference $avg_structure\n";
    }

    for $counter (1 .. 4) {
        print OUTCMD "center :1-" . ($counter * $strandlength) . " mass origin\n";
        print OUTCMD "image origin center\n";
    }

    print OUTCMD "rms reference out $filebase" . "_rms";
    if (! $isAVG) {
	print OUTCMD "_start";
    }
    print OUTCMD ".dat time 1 nofit\n";

    close OUTCMD;

    $out_cmd = "ptraj $top_file < tmp_gen_rst_file >& junk";

    print "Executing Ptraj...";

    if ( system "$out_cmd") {
	die "Failure!! Error executing ptraj\n";
    }

    print "Sucess.\nCreated $filebase" . "_rms";

    if (! $isAVG) { 
        print "_start"; 
    } 
    print ".dat as the rms datafile\n";

    system "rm -f tmp_gen_rst_file junk";
}

sub TestForInt(@) {
    my ($testval, $errortext) = @_;
    if (! $start_tm =~ /^\d+$/) {
	die "Invalid $errortext. Expected integer. Got $testval\n";
    }
}

 

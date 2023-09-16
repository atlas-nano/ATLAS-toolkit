#!/usr/bin/perl -w
use strict;
use constant pi => 3.14159265358979;

sub CRadDegrees(@);
my (@twist_angles);
my ($global_avg);

if (!@ARGV or $#ARGV < 1) {
    die "usage: function.pl minor_groove major_groove\n";
}

my ($minor_groove, $major_groove) = @ARGV;

my ($minor_avg) = 180/($minor_groove + 1); 
my ($major_avg) = 180/($major_groove + 1);
my ($in_angle, $curr_angle, $i);
my ($groove_1_avg, $groove_2_avg);

$global_avg = 360/($major_groove + $minor_groove);

#print "global_avg: $global_avg\n";

my ($major_peak) = abs((180/$major_groove) - $global_avg)/0.636;
my ($minor_peak) = abs((180/$minor_groove) - $global_avg)/0.636;

#print "major: $major_peak, minor: $minor_peak\n";
$major_peak = CRadDegrees($major_peak, 0);
$minor_peak = CRadDegrees($minor_peak, 0);


for $i (1 .. $minor_groove) {
    $in_angle = CRadDegrees(($minor_avg * $i), 0);
#    print "in: " . CRadDegrees(($in_angle), 1) . "\n";
    $curr_angle = $global_avg + CRadDegrees(($minor_peak * sin($in_angle)),1);
#    print " out: $curr_angle\n";
    $twist_angles[$i -1] = $curr_angle;
    $groove_1_avg += $curr_angle;
}
#print "Total: $groove_1_avg\n";

for $i (1 .. $major_groove) {
    $in_angle = CRadDegrees((($major_avg * $i) + 180), 0);
#    print "in: " . CRadDegrees(($in_angle), 1) . "\n";
    $curr_angle = $global_avg + CRadDegrees(($major_peak * sin($in_angle)),1);
#    print " out: $curr_angle\n";
    $twist_angles[$i + $minor_groove -1] = abs($curr_angle);
    $groove_2_avg += abs($curr_angle);
}
#print "Total: $groove_2_avg\n";

for $i (0 .. $#twist_angles) {
    if ($i == 0) {
#	print "Twist angle for groove with \"$minor_groove\"\n";
    } else {
	if ($i == $minor_groove) {
#	    print "Total: $groove_1_avg\n\n";
#	    print "Twist angle for groove with \"$major_groove\"\n";
	}
    }
    print "$twist_angles[$i]\n";
}
#print "Total: $groove_2_avg\n\n";
#print "Turn total: " .  ($groove_2_avg +  $groove_1_avg) . "\n";
sub CRadDegrees(@) {

    my ($convertToDegrees) = $_[1];
    my ($inital_angle) = $_[0];
    my ($resulting_angle) = 0.0;

    if ($convertToDegrees) { $resulting_angle = $inital_angle * 180 / pi; }
    else { $resulting_angle = $inital_angle * pi/180; }

    return $resulting_angle;
}

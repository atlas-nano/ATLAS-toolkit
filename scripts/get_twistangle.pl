#!/usr/bin/perl -w

use strict;

sub asin { atan2($_[0], sqrt(1 - $_[0] * $_[0])) }
sub DetermineTwistAngle(@);
sub CRadDegrees(@);

if (! @ARGV) {
    die "get_twistangle.pl datafile\n";
}

-e $ARGV[0] or die "Cannot find $ARGV[0]: $!\n";

open INFILE, $ARGV[0] or die "Cannot open $ARGV[0]: $!\n";

while (<INFILE>) {
    chomp;
    my $instr = $_;
    if ($instr =~ /^\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
	printf("%4d%9.2f\n", $1, DetermineTwistAngle($2, $3));
    }
}

close INFILE;

sub CRadDegrees(@) {

    my ($convertToDegrees) = $_[1];
    my ($inital_angle) = $_[0];
    my ($resulting_angle) = 0.0;

    my ($pi) = atan2(1,1) *4;

    if ($convertToDegrees) { $resulting_angle = $inital_angle * 180 / $pi; }
    else { $resulting_angle = $inital_angle * $pi/180; }

    return $resulting_angle;
}

sub DetermineTwistAngle(@) {

    my ($result) = 0.0;
    my ($twist) = CRadDegrees($_[0], 0);
    my ($rize) = $_[1];

    $result = CRadDegrees(asin(3.34/($rize/sin($twist))), 1);

    return $result;
}

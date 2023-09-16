#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Math::SphericalHarmonics qw(plgndr ylm);

my ($val, $val2, $l, $m, $x, $theta, $phi);

$l = 0;
$m = 0;
$x = 0;
$theta=0.5;
$phi=0.5;
$val = plgndr($l,$m,$x);
$val2 = ylm($l,$m,$theta,$phi);
print "l: $l m: $m x: $x Plm(x): $val theta: $theta phi: $phi Ylm(theta,phi): $val2->{re} -i $val2->{im}\n";

$l = 1; $m = 0; $x = 0.5;
$val = plgndr($l,$m,$x);
$val2 = ylm($l,$m,$theta,$phi);
print "l: $l m: $m x: $x Plm(x): $val theta: $theta phi: $phi Ylm(theta,phi): $val2->{re} -i $val2->{im}\n";

$l = 2; $m = 0; $x = 0.5;
$val = plgndr($l,$m,$x);
$val2 = ylm($l,$m,$theta,$phi);
print "l: $l m: $m x: $x Plm(x): $val theta: $theta phi: $phi Ylm(theta,phi): $val2->{re} -i $val2->{im}\n";

$l = 2; 
$theta=2.90981; $phi=2.68447;
for $m(-$l..$l) {
	$val2 = ylm($l,$m,$theta,$phi);
	print "l: $l m: $m x: $x Plm(x): $val theta: $theta phi: $phi Ylm(theta,phi): $val2->{re} -i $val2->{im}\n";
}

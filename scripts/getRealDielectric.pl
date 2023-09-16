#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use Math::Erf qw(erf);
use General qw(IsInteger IsDecimal);
use constant PI => atan2(1,1) * 4;

sub init;

$|++;
my ($gewaldFactor, $cuttoff, $dielectric) = init;
my ($qFactor);

$qFactor = erf($gewaldFactor * $cuttoff) - 2 * $gewaldFactor * $cuttoff * 
	   exp(-1 * $gewaldFactor**2 * $cuttoff**2)/sqrt(PI);

$dielectric = (($qFactor + 2)*($dielectric - 1) + 3)/(($qFactor - 1)*($dielectric - 1) + 3);

print "Q factor: $qFactor. Real dielectric: $dielectric\n";

sub init {
    my (%OPTS);

    getopt('ecd',\%OPTS);
    
    for ("e", "c", "d") {
	die "usage: $0 -e ewald factor -c real space cuttoff -d dielectric\n" if (! defined($OPTS{$_}));
    }
    
    print "Initializing...";
    for ("e", "c", "d") {
	die "ERROR: Expected number. Got \"$OPTS{$_}'\"\n" if (! IsInteger($OPTS{$_}) and ! IsDecimal($OPTS{$_}));
    }
    
    print "Done\n";

    return ($OPTS{e}, $OPTS{c}, $OPTS{d});
}


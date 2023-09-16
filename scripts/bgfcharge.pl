#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(ParseStructFile GetSystemCharge);
use General qw(FileTester);

die "usage: $0 bgf_file\n"
    if (! @ARGV);

my ($bgfFile) = $ARGV[0];
FileTester($bgfFile);

my ($ATOMS, $CONS) =ParseStructFile($bgfFile, 0,0);
delete $ATOMS->{HEADER} if (exists($ATOMS->{HEADER})); 
my ($totCharge) = GetSystemCharge($ATOMS);

printf "Total Charge: %-11.7f", $totCharge;
print "+" if ($totCharge <= 0);
printf "%-11.7f\n", (-1*$totCharge);

#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats;
use Packages::General;

die "usage: $0 bgf_file\n"
    if (! @ARGV);

my ($bgfFile) = $ARGV[0];
FileTester($bgfFile);

my ($ATOMS, $CONS) = GetBGFFileInfo($bgfFile, 0);

my ($totCharge) = GetSystemCharge($ATOMS);

printf "Total Charge: %-11.5f\n", $totCharge;

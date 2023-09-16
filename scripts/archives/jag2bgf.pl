#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General;

sub parseJagOut;
sub writeBGF;
sub initialize;

die "usage: $0 jaguar_output_file [save_name]\n"
    if (! @ARGV);

my ($jagFile, $saveName) = @ARGV;
initialize;

print "Parsing Jaguar Output file $jagOut...";
my ($AtmData) = parseJagOut($jagOut);
print "Done\nWriting BGF file $save_name...";
writeBGF($AtmData);
print "Done\n";

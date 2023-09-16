#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;

sub initialize;

die "usage: $0 topologyFile coordinateFile [save_name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($topFile, $crdFile, $save_name) = @ARGV;

initialize;

my ($PARMS, $CONN) = parseTopFile($topFile);
my ($ATOMS) = parseCoordFile($crdFile);

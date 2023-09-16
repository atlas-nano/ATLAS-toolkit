#!/usr/bin/perl -w

use strict;

die "usage: $0 bgf_file cerius_ff [save_name] [force_field type] [use_element]\n"
    if (! @ARGV or $#ARGV < 1);

my ($bgf_file, $ff_file, $save_name, $ff_type, $use_ele) = @ARGV;
my ($newScript);

$newScript = "/ul/tpascal/scripts/createLammpsInput.pl";
print "NOTE: Script is now obsolete. Use $newScript\n";
$newScript .= " -f $ff_file -b $bgf_file -s $save_name";
print "EXECUTING: $newScript\n\n";

system($newScript);

#!/usr/bin/perl -w
use strict;

die "usage: $0 amberLibs\n"
    if (! @ARGV);

my ($libFile) = $ARGV[0];
my (@LIBS, $outStr);

die "Error accessing regular file $libFile: $!\n"
   if (! -e $libFile or ! -r $libFile or ! -T $libFile);

open INDATA, $libFile or die "ERROR: Cannot open file $libFile: $!\n";
while (<INDATA>) {
    chomp;
    while ($_ =~ /(\S+)/g) {
	push @LIBS, $1;
    }
}

close INDATA;

$outStr = "parm = {";
for (@LIBS) {
    $outStr .= "$_ ";
}

$outStr .= "}\n";
$outStr .= "saveamberparm parm parm06.top parm06.crd\n";

open OUTDATA, "> leapadd" or die "ERROR: Cannot create file leapadd: $!\n";
print OUTDATA $outStr;
close OUTDATA;


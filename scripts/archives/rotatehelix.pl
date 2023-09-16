#!/usr/bin/perl -w

BEGIN {
    unshift (@INC, "/ul/tpascal/.libs/blib/arch/auto/p5namot");
}

use p5namot;

my ($pdbfile, $rotangle) = @ARGV;

if (! $pdbfile) {
    die "usage: rotatehelix.pl pdbfile rotationangle\n";
}

if (! $rotangle =~ /^\d+$/) {
    die "Invalid rotation angle, expected integer\n";
}

-e $pdbfile or die "Cannot open $pdbfile, $!\n";

p5namot::Cmd("set hush");
p5namot::Cmd("render");
p5namot::Cmd("load pdb na $pdbfile");
p5namot::Cmd("render");
p5namot::Cmd("rotate 1 3 $rotangle");
p5namot::Cmd("render");
p5namot::Cmd("write pdb $pdbfile");
p5namot::Cmd("render");
p5namot::Cmd("close");
p5namot::Cmd("quit");

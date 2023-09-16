#!/usr/bin/perl -w
BEGIN {
	unshift @INC, "/home/tpascal/scripts/";
}
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use Packages::Namot2::Namot;
use strict;

my ($pdbfile, $rotangle) = @ARGV;

die "usage: $0 pdbfile rotationangle\n"
    if (! @ARGV or $#ARGV == 0);

-e $pdbfile or die "Cannot open $pdbfile, $!\n";
$rotangle =~ /^\-?\d+/ or die "Integer expected for rotationangle. got $rotangle\n";
print "Rotating by $rotangle degrees\n";
p5namot::Cmd("set hush");
p5namot::Cmd("render");
p5namot::Cmd("load pdb na $pdbfile");
p5namot::Cmd("render");
p5namot::Cmd("rotate 1 3 $rotangle");
p5namot::Cmd("render");
p5namot::Cmd("write pdb $pdbfile");
p5namot::Cmd("render");
p5namot::Cmd("close");
#p5namot::Cmd("quit");

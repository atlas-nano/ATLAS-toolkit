#!/usr/bin/perl -w

# namot_cmd.pl - executes a namot2 script
# usage: namot_cmd.pl script

# Load the Perl5 Namot2 module
#
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Namot2::Namot;

p5namot::Cmd("set hush ERROR off");
p5namot::Cmd("set hush INFO off");
p5namot::Cmd("set hush REQUESTED off");
p5namot::Cmd("set hush WARNING off");

if (!@ARGV) {
    die "usage: namot_cmd.pl script\n";
}

-e $ARGV[0] or die "Cannot find $ARGV[0]: $!\n";

my ($instr, @outarray, $i);

open INFILE, $ARGV[0] or die "Cannot open $ARGV[0]: $!\n";
while (<INFILE>) {
    $instr = $_;
    chomp $instr;
    push @outarray, $instr;
}

close INFILE;

for $i (@outarray) {
p5namot::Cmd("$i");
p5namot::Cmd("render");
}

p5namot::Cmd("quit");

#!/usr/bin/perl -w

BEGIN {
    unshift (@INC, "/ul/tpascal/.libs/blib/arch/auto/p5namot");
}

($helix1, $helix2, $c_point, $whichHelix, $is5prime) = @ARGV;

if (! $helix1 or ! $helix2 or ! $c_point or ! $whichHelix) {
    die "usage getdist.pl helix1 helix2 crossover_point which_helix 5_prime_inside?\n";
}

-e $helix1 or die "Cannot find $helix1\n";
-e $helix2 or die "Cannot find $helix2\n";

if (! $c_point =~ /(\d+)/) {
    die "Invalid crossover point, expected integer\n";
}

$c_point =~ /(\d+)/;
$c_point = $1;

if (! $whichHelix =~/\d/) {
    die "Invalid helix, expected 1 or 2\n";
}

$whichHelix =~ /(\d)/;

if ($1 > 2 || $1 <1) {
    die "Invalid helix, expected 1 or 2\n";
}

$whichHelix = $1;

print "WhichHelix: $whichHelix\n";
if ($is5prime <0 or $is5prime >1) {
    die "find_o_angle: Invalid input for 5primein, expected 0 or 1\n";
}

use p5namot;
#p5namot::Cmd("set hush ERROR off");
p5namot::Cmd("render");
p5namot::Cmd("set hush INFO off");
p5namot::Cmd("render");
p5namot::Cmd("set hush WARNING off");
p5namot::Cmd("render");

p5namot::Cmd("load pdb na $helix1");
p5namot::Cmd("render"); 
p5namot::Cmd("load pdb na $helix2");
p5namot::Cmd("render");
p5namot::Cmd("trans 1 0 20 0.0");
p5namot::Cmd("render");

if ($whichHelix ==2 and $is5prime) {
    $c_point -= 1;
} else {
    if ($whichHelix==2) {
	$c_point += 1;
    }
}

for ($i = 0; $i<=361; $i++) {
  p5namot::Cmd("rotate $whichHelix 3 1");
  print "$i ";
  p5namot::Cmd("render");
  if ($is5prime) {
      p5namot::Cmd("query dist 1:$c_point:2:P 2:$c_point:2:P");
  } else {
      p5namot::Cmd("query dist 1:$c_point:1:P 2:$c_point:1:P");
  }
  p5namot::Cmd("render");
}


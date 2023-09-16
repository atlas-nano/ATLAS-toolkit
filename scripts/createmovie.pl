#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use p5namot;

p5namot::Cmd("set hush ERROR off");
p5namot::Cmd("set hush INFO off");
p5namot::Cmd("set hush REQUESTED off");
p5namot::Cmd("set hush WARNING off");

open MYSCRIPT, "run.script";
while (<MYSCRIPT>) {
    chomp;
    push @outarray, $_;
}

close MYSCRIPT;

for ($i=0;$i<=$#outarray;$i++) {
    p5namot::Cmd("$outarray[$i]");
      $addext = $i * 2;
      $addext1 = ($i *2) +1;
      if ($i<50) {
	  $addext = "0" . $addext;
	  $addext1 = "0" . $addext1;
      }
      if ($i <5) {
	  $addext = "0" . $addext;
	  $addext1 = "0" . $addext1;
      }
      p5namot::Cmd("set background black");
      p5namot::Cmd("render");
      p5namot::Cmd("write png dump" . $addext . ".png");
      $out_cmd = "/temp1/tpascal/projects/ImageMagick/bin/convert dump" . $addext . ".png dump" . $addext . ".jpg";
      system $out_cmd;
      $my_cmd = "cp dump" . $addext . ".jpg  dump" . $addext1 . ".jpg"; 
      system $my_cmd;
}

$out_cmd = "rm -f dump*.png 65*";
system $out_cmd;
p5namot::Cmd("quit");

#!/usr/bin/perl -w

BEGIN {
    unshift (@INC, "./blib/lib/");
	unshift (@INC, "./blib/arch/auto/p5namot/p5namot");
}

use p5namot;

$i=0;

p5namot::Cmd("generate s d b aaa");
p5namot::Cmd("render");
p5namot::Cmd("set background black");
p5namot::Cmd("render");
p5namot::Cmd("write png pho.$i.png");
p5namot::Cmd("close");

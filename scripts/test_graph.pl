#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$Bin/Packages";
use Graph;
use Graph::Directed;
use Graph::ChuLiuEdmonds;

my $graph = Graph::Directed->new(vertices=>[qw(a b c d)]);
$graph->add_weighted_edges(qw(a b 3 c d 7 d a 2 d b 1 c a 2));
my $msts = $graph->MST_ChuLiuEdmonds($graph);
print "cycles: " . $graph->has_a_cycle . "\n";

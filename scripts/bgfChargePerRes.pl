#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use FileFormats qw(GetBGFFileInfo sortByRes);
use General qw(FileTester);

my ($bgfFile);
my ($ATOMS, $RES);

&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, undef, undef) = GetBGFFileInfo($bgfFile);
$RES = sortByRes($ATOMS);
print "Done\n";
&printChargePerRes($ATOMS, $RES);

sub printChargePerRes {
    my ($atoms, $res) = @_;
    my ($i, $j, @reslist, $charge);

    @reslist = sort {$a<=>$b} keys %{ $res };
    for $i (@reslist) {
	$charge = 0;
	for $j (keys %{$res->{$i}{ATOMS}}) {
	    $charge += $atoms->{$j}{CHARGE};
	}
	printf "%10d%12.3f\n",$i,$charge;
    }
}

sub init {
    my (%OPTS);

    getopt('b',\%OPTS);
    die "usage: $0 -b bgffile\n" if (! exists($OPTS{b}));
    print "Initializing...";
    $bgfFile = $OPTS{b};
    FileTester($bgfFile);
    print "Done\n";
}

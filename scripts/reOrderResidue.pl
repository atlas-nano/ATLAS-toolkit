#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF sortByRes);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub reOrderResNumbers;
sub numerically {($a<=>$b); }

my ($bgfFile, $saveFile, $startRes);
my ($ATOMS, $BONDS, $HEADERS, $RES);
$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nReordering residue numbers...";
$RES = sortByRes($ATOMS);
&reOrderResNumbers($ATOMS, $RES, $startRes);
print "Done\nCreating BGF file $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub reOrderResNumbers {
    my ($atoms, $res, $resNum) = @_;
    my ($i, $j, $oldRes, $currRes);

    $oldRes = $currRes = -1;
    $resNum--;
    for $i (sort numerically keys %{ $atoms }) {
	$currRes = $atoms->{$i}{RESNUM};
	if ($currRes != $oldRes) {
	    $oldRes = $currRes;
	    $resNum++;
	}
	$atoms->{$i}{RESNUM} = $resNum;
    }
}

sub init {
    my (%OPTS);
    getopt('bsr',\%OPTS);
    die "usage: $0 -b bgf file -r (res # to start) -s (save file name)\n" 
	if (! exists($OPTS{b}));
    ($bgfFile, $saveFile, $startRes) = ($OPTS{b}, $OPTS{s}, $OPTS{r});
    print "Initializing...";
    FileTester($bgfFile);
    $startRes = 1 if (! defined($startRes) or $startRes !~ /^\d+$/);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_resSorted.bgf";
    }
    print "Done\n";
}

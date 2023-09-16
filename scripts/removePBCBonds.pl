#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use General qw(FileTester GetBondLength);
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use BOX qw(GetBox);

sub init;
sub removePBCBonds;

$|++;
&init;

my ($ATOMS, $BONDS, $HEADERS, $saveFile, $bgfFile, $BOX);

print "Parsing bgf file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$BOX = GetBox($ATOMS, undef, $HEADERS);
print "Done\nRemoving any bonds spanning periodic boundaries...";
&removePBCBonds($ATOMS, $BONDS, $BOX);
print "Done\nCreate bgf file $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub removePBCBonds {
    my ($atoms, $bonds, $box) = @_;
    my ($i, $j, $dist, $atom1, $atom2, $minBoxLen);

    $minBoxLen = 0;
    for $i ("X", "Y","Z") {
		$minBoxLen = $box->{$i}{len}/2 if (($box->{$i}{len}/2) > $minBoxLen);
    }

    for $i (keys %{ $bonds }) {
		$atom1 = \%{ $atoms->{$i} };
		$j = 0;
		while ($j <= $#{ $bonds->{$i} }) {
		    $atom2 = \%{ $atoms->{$bonds->{$i}[$j]} };
		    $dist = GetBondLength($atom1, $atom2);
			if ($dist > $minBoxLen) { 
				splice @{ $bonds->{$i} }, $j, 1;
			} else {
				$j++;
			}
		}
    }
}

sub init {
    my (%OPTS);

    getopt('bs',\%OPTS);
    die "usage: $0 -b bgf file -s [save name]\n" if (! defined($OPTS{b}));
    print "Initializing...";
    ($bgfFile, $saveFile) = ($OPTS{b}, $OPTS{s});

    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_nopbc.bgf";
    }

    print "Done\n";
}

#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";

use General qw(FileTester);
use FileFormats qw(GetBGFFileInfo addHeader createBGF createHeaders);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub updateBGFBox;
		     
my ($bgfFile, $saveName, $cell);
my ($ATOMS, $BONDS, $HEADERS, $BOX);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS,$HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nUpdating box info...";
&updateBGFBox($HEADERS,$cell);
&addHeader($ATOMS, $HEADERS);
print "Done\nSaving BGF file $saveName...";
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub updateBGFBox {
    my ($headers, $cell) = @_;
    my ($i, $bInfo);

    for $i ("PERIOD","CRYSTX","AXES","SGNAME","CELLS") {
        for (0 .. $#{ $HEADERS }) {
	    delete $HEADERS->[$_] if ($HEADERS->[$_] and $HEADERS->[$_] =~ /^$i/);
        }
    }

    $bInfo = "PERIOD 111\nAXES   ZYX\nSGNAME P 1                  1    1\n" .
    sprintf("CRYSTX %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n", @{ $cell }) .
    "CELLS    -1    1   -1    1   -1    1";
    push @{ $headers }, $bInfo;
}

sub init {
    my (%OPTS, @tmp);
  
    getopt('bsc',\%OPTS);

    die "usage: $0 -b bgf_file -c 'a b c alpha beta gamma' -s (savename)\n" if (! exists($OPTS{b}) or ! exists($OPTS{c}));
    print "Initializing...";
    ($bgfFile, $saveName) = ($OPTS{b}, $OPTS{s});
    FileTester($bgfFile);
    die "ERROR: Expected xx.xx yy.yy zz.zz from cell. Got \"$OPTS{c}\"!\n"
        if($OPTS{c} !~ /(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s*(.*)/);
    @{ $cell } = ($1,$2,$3);
    if ($4) {
	@tmp = split /\s+/, $4;
	for (@tmp) { push @{ $cell }, $_; }
    }
    for (3 .. 5) {
	push @{ $cell }, 90 if ($_ > $#{ $cell });
    }
    if (! $saveName) {
	$saveName = basename($bgfFile);
	$saveName =~ s/\.\w+$//;
	$saveName .= "_mod.bgf";
    }
    print "Done\n";
}


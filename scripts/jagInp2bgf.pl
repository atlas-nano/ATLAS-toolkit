#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo createBGF addHeader);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub ParseJagInput;
sub updateCoords;

my ($jagInput, $bgfFile, $saveFile);
my ($ATOMS, $BONDS, $HEADERS, $JAG);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
print "Done\nParsing Jaguar Input File $jagInput...";
$JAG = ParseJagInput($jagInput);
print "Done\nUpdating BGF coordinates...";
&updateCoords($ATOMS, $JAG);
print "Done\nCreating BGF file $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub updateCoords {
    my ($bgf, $jag) = @_;
    my ($i, $dim, $atmName);

    for $i (keys %{ $bgf }) {
	$atmName = $bgf->{$i}{ATMNAME};
	die "ERROR: Atom $i ($atmName) not found in Jaguar Input file!\n"
	    if (! exists($jag->{$atmName}));
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    $bgf->{$i}{$dim} = $jag->{$atmName}{$dim};
	}
    }
}

sub ParseJagInput {
    my ($inputFile) = $_[0];
    my (%DATA);

    open JAGINPUT, $inputFile or die "ERROR: Cannot open Jaguar Input file $inputFile: $!\n";
    while (<JAGINPUT>) {
	chomp;
	if ($_ =~ /^\s*(\w+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
	    $DATA{$1}{XCOORD} = $2;
	    $DATA{$1}{YCOORD} = $3;
	    $DATA{$1}{ZCOORD} = $4;
	}
    }
    close JAGINPUT;
    
    die "ERROR: No valid data found while reading $inputFile!\n"
	if (! %DATA);
    return \%DATA;
}

sub init {
    my (%OPTS);
    
    getopt('bjs',\%OPTS);
    die "usage: $0 -b bgf file -j jaguar input file -s (save bgf name)\n"
	if (! exists($OPTS{b}) or ! exists($OPTS{j}));
    print "Initializing...";
    ($bgfFile, $jagInput, $saveFile) = ($OPTS{b}, $OPTS{j}, $OPTS{s});
    FileTester($bgfFile);
    FileTester($jagInput);
    
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_mod.bgf";
    }
    print "Done\n";
}

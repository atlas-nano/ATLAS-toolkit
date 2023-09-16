#!/usr/bin/perl -w
#TODO - add impropers!!
#
#
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";

use Getopt::Std qw(getopt);
use File::Basename qw(basename dirname);

use FileFormats qw(GetBGFFileInfo addHeader createHeaders createBGF);
use General qw(FileTester LoadFFs LoadElements ReadFFs);
use CERIUS2 qw(saveCeriusFF);
use BOX qw(GetBox);
use GROMACS qw(ParseGroFF parseGroFile UpdateGroFFData);
use constant PI => 4*atan2(1.0,1.0);

sub init;

my ($topFile, $savePrefix, $groFile);
my ($PARMS, $ELEMENTS);
my ($ATOMS, $BONDS, $BOX, $HEADERS, $FF);

$|++;
&init;
print "Parsing GROMACS topology file $topFile...";
$PARMS = ParseGroFF($topFile,0);
print "Done\nGenerating atom and bond information...";
($ATOMS, $BONDS, $BOX) = parseGroFile($PARMS,$groFile);
print "Done\nGenerating forcefield information...";
$FF = UpdateGroFFData($PARMS,$ELEMENTS,$ATOMS);
print "Done\nCreating BGF file ${savePrefix}.bgf...";
&addHeader($ATOMS, createHeaders($BOX, $savePrefix));
&createBGF($ATOMS,$BONDS, "${savePrefix}.bgf");
print "Done\nCreating CERIUS2 formatted FF ${savePrefix}.ff...";
&saveCeriusFF($FF, "${savePrefix}.ff", $ELEMENTS);
print "Done\n";

sub init {
	my (%OPTS);

	getopt('tgs',\%OPTS);
	($topFile, $savePrefix, $groFile) = ($OPTS{t}, $OPTS{s}, $OPTS{g});
	die "usage: $0 -t gromacs topology file -g gromacs gro file -s (save prefix)\n"
		if (! exists($OPTS{t}) or ! exists($OPTS{g}));

	print "Initializing...";
	FileTester($topFile);
	FileTester($groFile);

	if (! defined($savePrefix)) {
		$savePrefix = basename($topFile);
		$savePrefix =~ s/\.\w+$//;
	}
	$ELEMENTS = LoadElements();
	print "Done\n";
}

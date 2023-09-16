#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF);
use Packages::General qw(FileTester);

sub init;
sub updateBGF;
sub showHelp;
sub numerically { ($a<=>$b); }

my ($ATOMS, $BONDS, $HEADERS, $bgfFile, $saveFile);

$|++;
&init;
print "Parsing MESO BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nSplitting residues...";
updateBGF($ATOMS);
print "Done\nCreating BGF $saveFile...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";


sub updateBGF {
    my ($atoms) = $_[0];
    my ($i, $resCount);

    $resCount = 1;
    for $i (sort numerically keys %{ $atoms }) {
	if ($atoms->{$i}{ATMNAME} =~ /PHO|PHE|ADE|CYT|GUA|THY/) {
	    $resCount++;
	}
	$atoms->{$i}{RESNUM} = $resCount;
    }
}
 
sub init {
    my (%OPTS);
    
    getopt('bs',\%OPTS);
    ($bgfFile, $saveFile) = ($OPTS{b},$OPTS{s});
    &showUsage if (! defined($bgfFile));
    print "Initializing...";
    FileTester($bgfFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$/_mod\.bgf/;
    }
    print "Done\n";
}

sub showUsage {
    die "usage: $0 -b bgf file -s [output file]\n";
}


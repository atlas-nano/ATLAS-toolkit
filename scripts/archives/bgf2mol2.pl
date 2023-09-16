#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use Packages::General qw(FileTester);
use Packages::FileFormats qw(GetBGFFileInfo createMOL2);
use Getopt::Std;
use File::Basename;

my ($ATOMS, $BONDS, $HEADERS, $bgfFile, $mol2File);
$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nCreating MOL2 file $mol2File...";
createMOL2($ATOMS, $BONDS, $mol2File);
print "Done\n";

sub init {
    my (%OPTS);
    getopt('bm',\%OPTS);
    ($bgfFile, $mol2File) = ($OPTS{b}, $OPTS{m});
    die "usage: $0 -b bgfFile -m [mol2file]\n" if (! defined($bgfFile));
    FileTester($bgfFile);
    if (! defined($mol2File)) {
	$mol2File = basename($bgfFile);
	$mol2File =~ s/\.\w+$/\.mol2/;
    }
}

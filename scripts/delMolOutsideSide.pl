#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols);

sub init;
sub removeMolsOutsideCell;
sub getCellBox;

my ($bgfFile, $saveFile);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $select);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
&GetMols($ATOMS, $BONDS);
$BOX = getCellBox($HEADERS);
print "Done\nRemoving any atoms outside unit cell...";
$select = removeMolsOutsideCell($ATOMS, $BOX);
($ATOMS, $BONDS, undef) = GetBGFAtoms($select, $ATOMS, $BONDS);
print "Done\nCreating BGF file $saveFile...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub removeMolsOutsideCell {
    my ($atoms, $boxInfo) = @_;
    my ($i, $dim, $j, $selected);

    %{ $selected } = %{ $atoms };
    for $i (keys %{ $atoms }) {
	next if (! exists($selected->{$i}));
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    if ($atoms->{$i}{$dim} < 0 or $atoms->{$i}{$dim} > $boxInfo->{$dim}) {
		for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
		    delete $selected->{$j};
		}
		last;
	    }
	}
    }
    return $selected;
}

sub getCellBox {
    my ($headerInfo) = $_[0];
    my ($i, %CellBOX);

    for $i (@{ $headerInfo }) {
	if ($i =~ /^CRYSTX\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
	    $CellBOX{XCOORD} = $1;
	    $CellBOX{YCOORD} = $2;
	    $CellBOX{ZCOORD} = $3;
	}
    }

    die "ERROR: No Cell information found!\n" if (! %CellBOX);
    return \%CellBOX;
}

sub init {
    my (%OPTS);
    
    getopt('bs', \%OPTS);
    die "usage: $0 -b bgf file -s (save name)\n" if (! defined($OPTS{b}));
    
    print "Initializing...";
    ($bgfFile, $saveFile) = ($OPTS{b}, $OPTS{s});
    FileTester($bgfFile);
    
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_mod.bgf";
    }
    print "Done\n";
}

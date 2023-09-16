#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General qw(FileTester);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub updateCoordinates;
sub GetSeqGeomFileInfo;
sub showUsage;

my ($bgfFile, $seqGeomFile, $saveFile);
my ($ATOMS, $BONDS, $HEADERS, $COORDS);

$|++;

&init;
print "Parsing BGF file $bgfFile...";
($ATOMS,$BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nParsing Seqquest Geom File $seqGeomFile...";
($COORDS) = GetSeqGeomFileInfo($seqGeomFile);
print "Done\nUpdating Coordinates...";
updateCoordinates($ATOMS, $COORDS);
print "Done\nCreating BGF file $saveFile...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub updateCoordinates {
    my ($atoms, $seqData)  = @_;
    my ($atomC);
    
    for $atomC (keys %{ $atoms }) {
	die "ERROR: Atom $atomC not found in Seqquest file!\n" if (! exists($seqData->{$atomC}));
	for ("XCOORD", "YCOORD", "ZCOORD") {
	    $atoms->{$atomC}{$_} = $seqData->{$atomC}{$_} * 0.529177;
	}
    }
}

sub GetSeqGeomFileInfo {
    my ($inFile) = $_[0];
    my (%SEQDATA, $rec);

    open SEQFILE, $inFile or die "ERROR: Cannot open Seqquest file $inFile: $!\n";
    while (<SEQFILE>) {
	chomp;
	if ($_ =~ /^\s*(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
	    $rec = (
		    {
			"INDEX" => $1,
			"TYPE"  => $2,
			"XCOORD" => $3,
			"YCOORD" => $4,
			"ZCOORD" => $5,
		    }
		    );
	    $SEQDATA{$1} = $rec;
	}
    }
    close SEQFILE;

    die "ERROR: File $inFile is not a valid Seqqest Geometry file!\n" if (! keys %SEQDATA);

    return \%SEQDATA;
}

sub init {
    my (%OPTS, $usage);

    getopt('bso',\%OPTS);
    ($bgfFile, $seqGeomFile, $saveFile) = ($OPTS{b}, $OPTS{s}, $OPTS{o});
    $usage = &showUsage;
    for ($bgfFile, $seqGeomFile) {
	die "$usage\n" if (! defined($_));
    }

    print "Initializing...";
    FileTester($bgfFile);
    FileTester($seqGeomFile);
    
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$/_mod\.bgf/;
    }
    print "Done\n";
}

sub showUsage {
    my ($usage) = "usage: $0 -b bgf file -s seqquest geom file -o [output name]\n" . 
	"Options:\n\t-b bgf file: location of the bgf file\n" .
	"\t-t seqquest geom file: location of the Seqquest geometry output file\n" .
	"\t-o [output name]: (Optional) Name of the output file.\n\t\t" . 
	"\t  Will be {bgfname}_mod.bgf if not specified\n";
    return $usage;
}

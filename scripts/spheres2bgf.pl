#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(createBGF createHeaders addHeader);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub parseSpheresFile;

my ($spheresFile, $bgfFile);
my ($ATOMS, $HEADERS);

$|++;
&init;
print "Parsing Spheres file $spheresFile...";
$ATOMS = parseSpheresFile($spheresFile);
print "Done\nCreate BGF file $bgfFile...";
$HEADERS = createHeaders(undef, $bgfFile);
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, undef, $bgfFile);
print "Done\n";

sub parseSpheresFile {
    my ($inFile) = $_[0];
    my (%ATOMS, $atomCount);

    open SPHERESFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
    $atomCount = 0;
    while (<SPHERESFILE>) {
	chomp;
	if ($_ =~ /^\s+(\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)/) {
	    $atomCount++;
	    $ATOMS{$atomCount} = (
			  {
			      "LABEL"     => "ATOM",
			      "ATMNAME"   => "H${1}",
			      "RESNAME"   => "P1",
			      "RESNUM"    => "$1",
			      "CHARGE"    => "0",
			      "FFTYPE"    => "He",
			      "NUMBONDS"  => 0,
			      "LONEPAIRS" => 0,
			  }
			  );
	    $ATOMS{$atomCount}{XCOORD} = $2;
	    $ATOMS{$atomCount}{YCOORD} = $3;
	    $ATOMS{$atomCount}{ZCOORD} = $4;
	}
    }
    close SPHERESFILE;

    die "ERROR: Spheres file $inFile is not valid!\n" if (! keys %ATOMS );
    
    return \%ATOMS;
}

sub init {
    my (%OPTS);
    
    getopt('sb',\%OPTS);
    ($spheresFile, $bgfFile) = ($OPTS{s},$OPTS{b});
    
    die "usage: $0 -s spheres file -b [bgf file]\n" if (! defined($spheresFile));

    print "Initializing...";
    FileTester($spheresFile);
    
    if (! defined($bgfFile)) {
	$bgfFile = basename($spheresFile);
	$bgfFile =~ s/\.\w+$/\.bgf/;
    }
}

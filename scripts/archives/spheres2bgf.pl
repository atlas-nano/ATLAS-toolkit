#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(createBGF createHeaders addHeader);
use Packages::General qw(FileTester);
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
    my (%ATOMS);

    open SPHERESFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
    while (<SPHERESFILE>) {
	chomp;
	if ($_ =~ /^\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
	    $ATOMS{$1} = (
			  {
			      "LABEL"     => "ATOM",
			      "ATMNAME"   => "H${1}",
			      "RESNAME"   => "P1",
			      "RESNUM"    => "900",
			      "CHARGE"    => "0",
			      "FFTYPE"    => "He",
			      "NUMBONDS"  => 0,
			      "LONEPAIRS" => 0,
			  }
			  );
	    $ATOMS{$1}{XCOORD} = $2;
	    $ATOMS{$1}{YCOORD} = $3;
	    $ATOMS{$1}{ZCOORD} = $4;
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

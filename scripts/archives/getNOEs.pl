#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use Packages::AMBER qw(ParseAmberTrj getTopInfo getOpts GetAmberByteOffset);
use Packages::FileFormats qw(GetBGFFileInfo sortByRes);
use Packages::General qw(FileTester TrjSelections STDev LoadElements);
use Packages::NOEs qw(GetNOEs SaveNOEs);
use File::Basename;
use Packages::CERIUS2 qw (parseCerius2FF);
use strict;

sub init;
sub writeData;
sub getElement;
sub calcNOEs;

die "usage: $0 amberTop amberTrajectory cerius2FF \"trajectorySelection\" [saveName]\n"
    if (! @ARGV || $#ARGV < 3);
my ($topFile, $trjFile, $ff, $selection, $saveName) = @ARGV;

my ($DATA, $OPTS, $SELECT, $ATOMS, $NOEs, $RES); 
my ($printStr, $PARMS, $HYDROGENS, $totAtms);

$|++;
print "Initializing...";
&init;
print "Done\nParsing AMBER topology file $topFile...";
($DATA, $totAtms) = getTopInfo($topFile, $OPTS);
print "Done\nParsing Cerius2 FF $ff...";
$PARMS = parseCerius2FF($ff, 0);
print "Done\nGetting list of Hydrogens...";
$HYDROGENS = getElement($DATA->{ATOMS}, $PARMS);
print "Done\n";
$printStr = "Parsing AMBER trajectory $trjFile...";
&GetAmberByteOffset($SELECT, $trjFile, $totAtms);
ParseAmberTrj($DATA->{ATOMS}, $trjFile, $SELECT, $totAtms, \&calcNOEs, $printStr, $HYDROGENS);
print "Creating data file $saveName...";
writeData($NOEs, $saveName);
print "Done\n";

sub init {
    FileTester($topFile);
    FileTester($trjFile);
    FileTester($ff);

    $SELECT = TrjSelections($selection);

    $saveName = basename($topFile) if (! defined($saveName));
    $saveName =~ s/\.\w+$/\.dat/;
    $OPTS = &getOpts;
    
    $RES->{28} = 1;
    $RES->{29} = 1;
    $RES->{27} = 1;
    $RES->{18} = 1;
    $RES->{17} = 1;
    $RES->{8} = 1;
    $RES->{7} = 1;
    $RES->{30} = 1;
}

sub getElement {
    my ($atoms, $ffdata) = @_;
    my (%HLIST, $i, $element);

    for $i (keys %{ $atoms }) {
	$element = $ffdata->{ATOMTYPES}{ $atoms->{$i}{FFTYPE} }{ATOM};
	next if ! (defined($element));
	if (lc($element) eq "h") {
	    $HLIST{$i} = 1;
	}
    }

    return \%HLIST;
}

sub calcNOEs {
    my ($ATOMS, $BOX, $frameNum, $HLIST) = @_;
    my ($i, $j, $dim, $dist, $noeName);

    for $i (keys %{ $HLIST }) {
	next if (! exists($RES->{ $ATOMS->{$i}{RESNUM} }));
	$ATOMS->{$i}{RESNAME} =~ s/D//;
	for $j (keys %{ $HLIST }) {
	    next if ($i <= $j);
	    next if ($ATOMS->{$i}{RESNUM} == $ATOMS->{$j}{RESNUM});
	    next if (! exists($RES->{ $ATOMS->{$j}{RESNUM} }));
	    $ATOMS->{$j}{RESNAME} =~ s/D//;
	    $dist = 0;
	    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
		$dist += ($ATOMS->{$i}{$dim} - $ATOMS->{$j}{$dim})**2;
	    }
	    $dist = sqrt($dist);
	    if ($dist <= 5) {
		if ($ATOMS->{$i}{RESNUM} > $ATOMS->{$i}{RESNUM}) {
		    $noeName = $ATOMS->{$i}{RESNAME} . $ATOMS->{$i}{RESNUM} . " " . $ATOMS->{$i}{ATMNAME} . " - ";
		    $noeName .= $ATOMS->{$j}{RESNAME} . $ATOMS->{$j}{RESNUM} . " " . $ATOMS->{$j}{ATMNAME};
		} else {
		    $noeName = $ATOMS->{$j}{RESNAME} . $ATOMS->{$j}{RESNUM} . " " . $ATOMS->{$j}{ATMNAME} . " - ";
		    $noeName .= $ATOMS->{$i}{RESNAME} . $ATOMS->{$i}{RESNUM} . " " . $ATOMS->{$i}{ATMNAME};
		}
		$NOEs->{$noeName} = 0 if (! exists($NOEs->{$noeName}));
		$NOEs->{$noeName}++;
	    }
	}
    }
		
}

sub writeData {
    my ($noes, $save) = @_;
    my ($noeName, $counter, @tmp, $tot, $percent);

    $counter = 0;
    $tot = 26;

    open NOEDATA, "> $save" || die "ERROR: Cannot write to $save: $!\n";
    for $noeName (keys %{ $noes }) {
	$percent = $noes->{$noeName} * 100/$tot;
	next if ($percent < 60);
	$counter++;
	printf NOEDATA "%8d%8.3f # $noeName\n", $counter, $percent;
	}
    close NOEDATA;
}

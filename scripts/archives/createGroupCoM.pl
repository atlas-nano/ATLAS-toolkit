#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";	
}

use strict;
use Packages::AMBER qw(getTopInfo ParseAmberTrj getOpts GetAmberByteOffset);
use Packages::General qw(FileTester TrjSelections CoM);
use Packages::FileFormats qw (sortByRes createBGF);
use File::Basename;

sub doAnal;
sub CreateHelixStructureFile;
sub init;
sub numerically;

die "usage: $0 topFile trjFile parmFile \"trajectory selection\"\n"
    if (! @ARGV || $#ARGV < 3);

my ($topFile, $trjFile, $parmFile, $selection) = @ARGV;
my ($OPTS, $SELECT, $DATA, $totAtms, $RESDATA); 
my ($printStr, $VEC, $link, $atmCounter, $trjName);

print "Initializing...";
&init;
print "Done\nParsing AMBER topology file $topFile...";
($DATA, $totAtms) = getTopInfo($topFile, $OPTS);
$RESDATA = sortByRes($DATA->{ATOMS});
print "Done\n"; 
$printStr = "Parsing AMBER trajectory $trjFile...";
open AMBERTRJ, "> $trjName" or die "ERROR: Cannot create file $trjName: $!\n";
print AMBERTRJ "AMBER Trajectory of 3WAY helix COM\n";
$atmCounter = 0;
&GetAmberByteOffset($SELECT, $trjFile, $totAtms);
ParseAmberTrj($DATA->{ATOMS}, $trjFile, $SELECT, $totAtms, \&doAnal, $printStr, \*AMBERTRJ);
close AMBERTRJ;
print "Saved helix information to $trjName\n";

sub doAnal {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($helixNum, $group, $res, $groupAtms, $atom, $i, $com, $counter);
    my ($HELIX, $RES, @tmp);
    
    $HELIX = $OPTS->{HELIX};
    @tmp = sort numerically keys %{ $HELIX };
    $RES = $RESDATA;

    for $helixNum (@tmp) {
	$counter = 0;
	for $group (0 .. $#{ $HELIX->{$helixNum}{COM} }) {
	    $groupAtms = ();
	    $counter++;
	    for $res (values %{ $HELIX->{$helixNum}{COM}[$group] }) {
		for $atom (keys %{ $RES->{$res}{ATOMS} }) {
		    for $i (keys %{ $ATOMS->{$atom} }) {
			$groupAtms->{$atom}{$i} = $ATOMS->{$atom}{$i};
		    }
		}
	    }
	    $com = CoM($groupAtms);    # calc the CoM of the group
	    for $i ("XCOORD", "YCOORD", "ZCOORD") { #now print the atom info
		if ($atmCounter >= 10) {
		    print $fileHandle "\n";
		    $atmCounter = 0;
		}
		printf $fileHandle "%8.3f", $com->{$i};
		$atmCounter++;
	    }
	   
	}
    }

    #print the box info
    printf $fileHandle "\n%8.3f%8.3f%8.3f\n", $BOX->{2}{DATA}, $BOX->{3}{DATA}, $BOX->{4}{DATA};
    $atmCounter = 0;
}

sub CreateHelixStructureFile {
    my ($HELIX, $helixNum, $offset) = @_;
    my ($i, $j, $rec, $k);
    
    $HELIX->{OFFSET} = $offset;

    $rec = (
	    {
		"LABEL"     => "ATOM",
		"ATMNAME"   => "BP",
		"RESNAME"   => "CoM",
		"RESNUM"    => $helixNum,
		"FFTYPE"    => "Ne",
		"NUMBONDS"  => 0,
		"LONEPAIRS" => 0,
		"CHARGE"    => 0,
		"RADII"     => 5,
	    }
	    );
    
    for $i (0 .. $#{ $HELIX->{COM} }) {
	$j = $i + 1 + $offset;
	for $k (keys %{ $rec }) {
	    $HELIX->{ATOMS}{$j}{$k} = $rec->{$k};
	}
	$HELIX->{ATOMS}{$j}{XCOORD} = $j;
	$HELIX->{ATOMS}{$j}{YCOORD} = $j;
	$HELIX->{ATOMS}{$j}{ZCOORD} = $j;
	$HELIX->{ATOMS}{$j}{ATMNAME} .= $j;
	if ($i == 0) {
	    $HELIX->{BONDS}{$j}[0] = $j + 1;
	} elsif ($i == $#{ $HELIX->{COM} }) {
	    $HELIX->{BONDS}{$j}[0] = $j - 1;
	} else {
	    $HELIX->{BONDS}{$j}[0] = $j - 1;
	    $HELIX->{BONDS}{$j}[1] = $j + 1;
	}
    }
    
    return $j;
}

sub init {
    my ($inStr, @resInfo, $helixNum, $rec, $i, $isValid);
    my ($math_cmd, $offset, $ATOMS, $BONDS, $j, $k, $bgfName);

    $|++;
    FileTester($topFile);
    FileTester($trjFile);
    FileTester($parmFile);

    $OPTS = &getOpts;

    $isValid = $offset = 0;
    open PARMFILE, $parmFile || die "ERROR: Cannot open parameter file $parmFile: $!\n";
    while (<PARMFILE>) {
	chomp;
	if ($_ =~ /^Helix (\d): (.+)/) {
	    $helixNum = $1;
	    $inStr = $2;
	    while ($inStr =~ /(\S+)/g) {
		$isValid = 1;
		@resInfo = split /\,/, $1;
		$rec = ();
		for $i (0 .. $#resInfo) {
		    $rec->{$i} = $resInfo[$i];
		}
		push @{ $OPTS->{HELIX}{$helixNum}{COM} }, $rec;
	    }
	    $offset = CreateHelixStructureFile($OPTS->{HELIX}{$helixNum}, $helixNum, $offset);
	}
    }
    close PARMFILE;

    for $j (keys %{ $OPTS->{HELIX} }) {
	for $i (keys %{ $OPTS->{HELIX}{$j}{ATOMS} }) {
	    %{ $ATOMS->{$i} } = %{ $OPTS->{HELIX}{$j}{ATOMS}{$i} };
	    @{ $BONDS->{$i} } = @{ $OPTS->{HELIX}{$j}{BONDS}{$i} };
	}
    }
    
    $bgfName = basename($topFile);
    $bgfName =~ s/\.\w+$//g;
    $bgfName .= "_com.bgf";
    createBGF($ATOMS, $BONDS, $bgfName);

    $trjName = basename($trjFile);
    if ($trjName =~ /^(.+)\.(\w+)$/) {
	$trjName = $1 . "_com." . $2;
    } else {
	$trjName .= "_com.crdbox";
    }

    die "ERROR: Parameter file $parmFile is invalid!\n" if (! $isValid);

    $SELECT = TrjSelections($selection);
}

sub numerically {
    ($a<=>$b);
}

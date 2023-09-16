#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use Packages::AMBER qw(ParseAmberTrj getTopInfo getOpts GetAmberByteOffset);
use Packages::FileFormats qw(GetBGFFileInfo sortByRes);
use Packages::General qw(FileTester TrjSelections STDev);
use Packages::NOEs qw(GetNOEs SaveNOEs);
use File::Basename;
use strict;

sub init;
sub writeData;
sub getFNames;
sub getMode;
sub numerically;

die "usage: $0 amberTop amberTrajectory noeList \"trajectorySelection\" [saveName]\n"
    if (! @ARGV || $#ARGV < 3);
my ($topFile, $trjFile, $noeFile, $selection, $saveName) = @ARGV;

my ($DATA, $OPTS, $SELECT, $RES, $ATOMS, $NOEs, $totAtms, $printStr);

$|++;
print "Initializing...";
&init;
print "Done\nParsing AMBER topology file $topFile...";
($DATA, $totAtms) = getTopInfo($topFile, $OPTS);
$RES = sortByRes($DATA->{ATOMS});
print "Done\nReading NOEs from $noeFile...";
$NOEs = GetNOEs($noeFile, $DATA->{ATOMS}, $RES, $ATOMS);
print "Done\n";
$printStr = "Parsing AMBER trajectory $trjFile...";
&GetAmberByteOffset($SELECT, $trjFile, $totAtms);
ParseAmberTrj(undef, $trjFile, $SELECT, $totAtms, \&SaveNOEs, $printStr, $NOEs);
print "Creating data file $saveName...";
writeData($NOEs, $saveName);
print "Done\n";

sub init {
    FileTester($topFile);
    FileTester($trjFile);
    FileTester($noeFile);

    $SELECT = TrjSelections($selection);
    $saveName = basename($topFile) if (! defined($saveName));
    $saveName =~ s/\.\w+$/\.dat/;
    $OPTS = &getOpts;
}

sub writeData {
    my ($noes, $save) = @_;
    my ($i, $vals, $counter, $j, $fName, $listName, $tot);
    my (%TRACKER, $k, $tName, $aName, $noeName, $mode, %NOEList);

    ($listName, $aName, $tName) = getFNames($save);
    ($mode, $tot) = getMode($noes);
    system "mkdir -p data" if ($mode == 1);

    open NOEFILE, "> $save" || die "ERROR: Cannot create file $save: $!\n";
    open NOELIST, "> $listName" || die "ERROR: Cannot create file $listName: $!\n"; 
    open NOEAVG, "> $aName" || die "ERROR: Cannot create file $aName: $!\n";
    $counter = 0;
    for $i (@{ $noes }) {
	printf NOELIST "%8d %8d # $i->{NAME}\n", $i->{ATOM1}, $i->{ATOM2};
	$counter++;
	$noeName = $i->{NAME};
	$noeName =~ s/\s+/ \- /g;
	$noeName =~ s/_/ /g;
	$NOEList{ $i->{NAME} } = "";
	$vals = "";
	for $j (0 .. $#{ $i->{VALS} }) {
	    $vals .= "$i->{VALS}[$j] ";
	    $NOEList{ $i->{NAME} } .= sprintf("%5d%8.3f\n",($j + 1), $i->{VALS}[$j]);
	    if ($mode == 2) {
		for $k (1 .. $tot) {
		    if (! exists($TRACKER{$j}{$k})) {
			$TRACKER{$j}{$k} = 0;
		    }
		    if (exists($i->{AVG}{$k}) && (abs($i->{AVG}{$k} - $i->{VALS}[$j]) <= 0.5)) { # with range of measured NOE
			$TRACKER{$j}{$k}++;
		    }
		}
	    }
	}
	($i->{STATS}{AVG}, $i->{STATS}{STDEV}, $i->{STATS}{TOTAL}) = STDev($vals);
	$i->{STATS}{TOTAL} /= $i->{STATS}{AVG};
	$i->{STATS}{TOTAL}++;
	printf NOEFILE "%5d%8.3f%8.3f%8d # $i->{NAME}\n", $counter,
	    $i->{STATS}{AVG}, $i->{STATS}{STDEV}, $i->{STATS}{TOTAL};
	printf NOEAVG "$noeName %8.3f \n", $i->{STATS}{AVG};
    }
    close NOEFILE;
    close NOELIST;
    close NOEAVG;

    if ($mode == 1) {
	for $i (keys %NOEList) {
	    $counter++;
	    $fName = $i;
	    $fName =~ s/\s+/_/g;
	    $fName =~ s/\W+//g;
	    $fName = "data/${counter}_${fName}.dat";
	    open NOEOUT, "> $fName" || die "ERROR: Cannot create $fName: $!\n";
	    print NOEOUT "$NOEList{$i}";
	    close NOEOUT;	
	}
    } else {
	open TRAKERFILE, "> $tName" || die "ERROR: Cannot create $tName: $!\n";
	for $i (sort numerically keys %TRACKER) {
	    printf TRAKERFILE "%5d", ($i + 1);
	    for $j (1 .. $tot) {
		printf TRAKERFILE "%8.3f", ($TRACKER{$i}{$j}/($#{ $noes } + 1));
	    }
	    print TRAKERFILE "\n";
	}
	close TRAKERFILE;
    }
}

sub getFNames {
    my ($save) = $_[0];
    my ($listName, $aName, $tName);
    
    if ($save =~ /(.+)\.(\w+)$/) {
	$listName = "${1}_noelist.${2}";
	$aName = "${1}_avg.${2}";
	$tName = "${1}_noe_traker.${2}";
    } else {
	$listName = "${save}_noelist.dat";
	$aName = "${save}_avg.dat";
	$tName = "${save}_noe_traker.dat";
    }

    return ($listName, $aName, $tName);
}

sub getMode {
    my ($noes) = $_[0];
    my ($mode, $tot) = (1,0);
    my ($i, $j);

    for $i (@{ $noes }) {
	if (exists($i->{AVG})) {
	    $mode = 2;
	    for $j (keys %{ $i->{AVG} }) {
		$tot = $j if ($j > $tot);
	    }
	}
    }
    
    return ($mode, $tot);
}

sub numerically {
    ($a<=>$b);
}

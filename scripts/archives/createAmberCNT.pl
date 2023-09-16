#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo createBGF addHeader);
use Packages::General qw(FileTester Trim IsInteger);
use File::Basename;
use Getopt::Std;

sub init;
sub updateBGF;
sub createLEAPFile;
sub getBondList;

my ($bgfFile, $numAtms, $suffix, $cntRES);
my ($ATOMS, $BONDS, $HEADERS, $BOX);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nUpdate BGF file info...";
$cntRES = updateBGF($ATOMS);
print "Done\nCreating AMBER compatible BGF file ${suffix}.bgf...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, "${suffix}.bgf");
print "Done\nCreating LEAP input ${suffix}_leap...";
createLEAPFile($cntRES, "${suffix}");
print "Done\n";

sub createLEAPFile {
    my ($resInfo, $fName) = @_;
    my ($i, $bondList, $atmList, $atom);

    open LEAPFILE, "> ${fName}.leap" || die "ERROR: Cannot create LEAP file ${fName}.leap: $!\n";
    for $i (sort {$a<=>$b} keys %{ $resInfo }) {
	$atom = $resInfo->{$i};
	printf LEAPFILE "%-4s = createAtom%5s \"".  $atom->{FFTYPE} . "\"%8.3f\n",
	lc($atom->{ATMNAME}), $atom->{"ATMNAME"}, $atom->{"CHARGE"};
	$atmList .= "add r " . lc($atom->{ATMNAME}) . "\n";
	$bondList .= getBondList($i);
    }
    print LEAPFILE "\nr = createResidue CNT\n$atmList\n$bondList\n";
    print LEAPFILE "CNT = createUnit CNT\nadd CNT r\nset CNT.1 restype undefined\n";
    print LEAPFILE "saveOff CNT ${fName}.lib\nquit\n";
    close LEAPFILE;
}

sub getBondList {
    my ($atom) = $_[0];
    my ($atomName, $returnStr, $bond, @bondAtms, $line);

    $returnStr = "";
    $atomName = $ATOMS->{$atom}{"ATMNAME"};

    if (! exists($BONDS->{$atom}) or $#{ $BONDS->{$atom} } == -1) {
        return "";
    }
    @bondAtms = @{ $BONDS->{$atom} };
    $line = sprintf("bond %4s", lc($atomName));
    for $bond (@bondAtms) {
        next
            if ($atom < $bond );
        $returnStr  .= $line . sprintf("%4s\n", lc($ATOMS->{$bond}{"ATMNAME"}));
    }

    return $returnStr;
}

sub updateBGF {
    my ($atoms) = $_[0];
    my ($i, @tmp, $counter, $resID, %NEWRES);

    @tmp = sort {$a<=>$b} keys %{ $atoms };
    $counter = 1;
    $resID = 1;
    for $i (@tmp) {
	$atoms->{$i}{ATMNAME} = "C${counter}";
	$atoms->{$i}{RESNAME} = "CNT";
	$atoms->{$i}{RESNUM} = $resID;
	$atoms->{$i}{FFTYPE} = "CA";
	$counter++;
	if ($resID == 1) {
	    $NEWRES{$i} = $atoms->{$i};
	}
	if ($counter > $numAtms) {
	    $resID++;
	    $counter = 1;
	}
    }
    
    return \%NEWRES;
}

sub init {
    my (%OPTS);
    
    getopt('bns',\%OPTS);
    ($bgfFile, $numAtms, $suffix) = ($OPTS{b}, $OPTS{n}, $OPTS{s});
    die "usage: $0 -b bgfFile -n num_atms -s [saveName]\n"
	if (! defined($bgfFile) || ! defined($numAtms));
    print "Initializing...";

    FileTester($bgfFile);
    $numAtms = Trim($numAtms);
    die "ERROR: Expected integer for $numAtms. Got \"$numAtms\"\n" if (! IsInteger($numAtms));
    if (! defined($suffix)) {
	$suffix = basename($bgfFile);
	$suffix =~ s/\.\w+$/_mod/;
    } else {
	$suffix =~ s/\.\w+$//;
    }
    print "Done\n";
}

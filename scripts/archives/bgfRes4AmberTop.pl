#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo sortByRes);
use Packages::General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub numerically { ($a<=>$b); }
sub getResData;
sub writeResData;

my ($ATOMS, $BONDS, $bgfFile, $RES, $saveFile, $resInfo);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
print "Done\nCreating Residue list...";
$RES = sortByRes($ATOMS);
$resInfo = getResData($ATOMS, $RES);
print "Done\nCreating data file $saveFile...";
&writeResData($resInfo, $saveFile);
print "Done\n";

sub writeResData {
    my ($RESDATA, $saveName) = @_;
    my ($i, $resPointerData, $j, $k);

    open OUTDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    print OUTDATA "%FLAG RESIDUE_LABEL\n%FORMAT(20a4)\n";

    $j = $k = 0;
    for $i (@{ $RESDATA }) {
	printf OUTDATA "%-4s", $i->{RESNAME};
	$resPointerData .= sprintf("%8d", $i->{START_ATOM});
	$j++;
	$k++;
	if ($j == 20) {
	    $j = 0;
	    print OUTDATA "\n";
	}
	if ($k == 10) {
	    $k = 0;
	    $resPointerData .= "\n";
	}
    }
    print OUTDATA "\n" if ($j > 0);
    print OUTDATA "%FLAG RESIDUE_POINTER\n%FORMAT(10I8)\n$resPointerData";
    close OUTDATA;
}

sub getResData {
    my ($atoms, $res) = @_;
    my ($i, $atomID, @tmp, @RESDATA, $rec);

    for $i (sort numerically keys %{ $res }) {
	@tmp = sort numerically keys %{ $res->{$i}{ATOMS} };
	$atomID = $tmp[0];
	$rec = (
		{
		    "START_ATOM" => $atomID,
		    "RESNAME"    => $atoms->{$atomID}{RESNAME},
		}
		);
	push @RESDATA, $rec;
    }

    return \@RESDATA;
}

sub init {
    my (%OPTS);
    
    getopt('bs',\%OPTS);
    die "usage: $0 -b bgffile -s [data file]\n" if (! defined($OPTS{b}));
    ($bgfFile, $saveFile) = ($OPTS{b}, $OPTS{s});
    print "Initializing...";
    FileTester($bgfFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$/\.dat/;
    }
    print "Done\n";
}

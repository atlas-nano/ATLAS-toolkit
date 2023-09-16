#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use FileFormats qw(GetBGFFileInfo);
use AMBER qw(CreateAmberTrj);
use BOX qw(GetBox);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub numerically { ($a<=>$b); }
sub CreateAmberRst;
sub getAtomVel;

my ($bgfFile, $velFile, $rstFile);
my ($ATOMS, $BONDS, $BOX, $HEADERS);

$|++;
&init;
print "Parsing BGF File $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\n";
if (defined($velFile)) {
    print "Reading velocity file $velFile...";
    &getAtomVel($velFile, $ATOMS);
    print "Done\n";
}
print "Creating AMBER Restart file ";
print "(with velocities) $velFile..." if (defined($velFile));
open OUTFILE, "> $rstFile" or die "ERROR: Cannot create Amber Trj file $rstFile: $!\n";
print OUTFILE "TITLE: AMBER Restart file created by $0 on " .
    scalar(localtime(time())) . "\n";
print OUTFILE scalar(keys %{ $ATOMS }) . "  0.2505000E+05\n";
$BOX = GetBox($ATOMS, undef, $HEADERS);
&CreateAmberRst($ATOMS, $BOX, \*OUTFILE);
print "Done\n";
sub CreateAmberRst {
    my ($ATOMS, $BBOX, $OUTFILE) = @_;
    my ($atomC, $counter, @tmp, $dim, $BOX);

    $counter = 0;
    @tmp = sort numerically keys %{ $ATOMS };
    if (defined($BBOX)) {
        for $dim ("X", "Y", "Z") {
            $BOX->{$dim} = $BBOX->{$dim}{"hi"} - $BBOX->{$dim}{"lo"};
        }
    }

    for $atomC (@tmp) {
        for $dim ("XCOORD", "YCOORD", "ZCOORD") {
            $counter++;
            printf $OUTFILE "%12.7f", $ATOMS->{$atomC}{$dim};
            if ($counter == 6) {
                print $OUTFILE "\n";
                $counter = 0;
            }
        }
    }
    print $OUTFILE "\n" if ($counter > 0);
    $counter = 0;
    for $atomC (@tmp) {
        for $dim ("XVEL", "YVEL", "ZVEL") {
	    next if (! exists( $ATOMS->{$atomC}{$dim} ));
            $counter++;
            printf $OUTFILE "%12.7f", $ATOMS->{$atomC}{$dim};
            if ($counter == 6) {
                print $OUTFILE "\n";
                $counter = 0;
            }
        }
    }
    
    print $OUTFILE "\n" if ($counter > 0);
    if (defined($BBOX)) {
        for $dim ("X", "Y", "Z") {
            printf $OUTFILE "%12.7f", $BOX->{$dim};
        }
	for $dim (1 .. 3) {
	    printf $OUTFILE "%12.7f", 90.0;
	}
        print $OUTFILE "\n";
    }
}

sub getAtomVel {
    my ($infile, $atoms) = @_;
    my ($VELS, $start, $atomC, $counter, @dims, $dimC, $inStr);
    my ($atomCount);

    open VELFILE, $infile or die "ERROR: Cannot read from velocity file $infile: $!\n";
    $start = $atomC = 0;
    $counter = 1;
    @dims = ("XVEL", "YVEL", "ZVEL");
    $atomCount = scalar(keys %{ $atoms });

    while (<VELFILE>) {
	chomp;
	$inStr = $_;
	if ($inStr =~ /^\s*\-?\d+\.\d+\s+/) {
	    $start = 1;
	}
	if ($start) {
	    while ($inStr =~ /(\-?\d+\.\d+)/g) {
		if (($counter % 3) == 1) {
		    $dimC = 0;
		    $atomC++;
		} elsif (($counter % 3) == 2) {
		    $dimC = 1;
		} else {
		    $dimC = 2;
		}
		$VELS->{$atomC}{$dims[$dimC]} = $1;
		$counter++;
	    }
	}
    }
    close VELFILE;
    if (scalar(keys %{ $VELS }) == scalar(keys %{ $atoms })) {
	for $atomC (keys %{ $VELS }) {
	    for $dimC ("XVEL", "YVEL", "ZVEL") {
		$atoms->{$atomC}{$dimC} = $VELS->{$atomC}{$dimC};
	    }
	}
    } else {
	print "inconsistent atom count: " .  scalar(keys %{ $atoms }) . 
	    " in bgf but " . scalar(keys %{ $VELS }) . " in velocity...";
    }
}

sub init {
    my (%OPTS);
    getopt('bsv',\%OPTS);
    die "usage: $0 -b bgf file -v (velocity file) -s (save name)\n" if (! exists($OPTS{b}));
    print "Initializing...";
    ($bgfFile, $velFile, $rstFile) = ($OPTS{b}, $OPTS{v}, $OPTS{s});
    FileTester($bgfFile);
    undef($velFile) if (! defined($velFile) or ! -e $velFile or ! -r $velFile);
    if (! defined($rstFile)) {
	$rstFile = basename($bgfFile);
	$rstFile =~ s/\.\w+$//;
	$rstFile .= "_rst.inpcrd";
    }
    print "Done\n";
}

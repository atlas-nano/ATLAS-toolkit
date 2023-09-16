#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester STDev GetStats GetEquilPoint);

sub init;
sub GetEquilPoint;
sub parseDataFile;

my ($DATA, $dataFile, $STATS, $tolerance);
my ($equilDATA, $isFound);

$|++;
&init;
print "Getting data from $dataFile...";
$DATA = parseDataFile($dataFile);
print "Done\nDetermining Equilibration point...";
($isFound, $equilDATA) = GetEquilPoint($DATA, $tolerance);
print "no found...will use all data..." if (! $isFound);
print "Done\nCollecting Stats...";
$STATS = GetStats($equilDATA);
print "Done\n";
print "AVG: $STATS->{AVG}\nSTDEV: $STATS->{STDEV}\n# data vals: $STATS->{NUM}\n";

		
sub parseDataFile {
    my ($inFile) = $_[0];
    my (%data, $counter);

    $counter = 0;
    open DATFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
    while (<DATFILE>) {
	chomp;
	if ($_ =~ /^\s*(\-?\d+\.?\d*)\s+(\-?\d+\.?\d*)/) {
	    $data{$1} = $2;
	} elsif ($_ =~ /^\s*(\-?\d+\.?\d*)\s*$/) {
	    $data{$counter} = $1;
	    $counter++;
	}
    }
    close DATFILE;
    
    die "ERROR: $inFile does not contain any valid information!\n" if (! %data);
    return \%data;
}

sub init {
    my (%OPTS);
    getopt('dt',\%OPTS);

    die "usage: $0 -d data file -t (tolerance %)\n" if (! exists($OPTS{d}));
    print "Initializing...";
    ($dataFile, $tolerance) = ($OPTS{d}, $OPTS{t});
    FileTester($dataFile);
    $tolerance = 1 if (! defined($tolerance));
    die "ERROR: Expected integer/decimal for tolerance. Got \"$tolerance\"\n"
	if ($tolerance !~ /^\d+\.?\d*/);
    print "warning.. tolerance > 100% of average!..." if ($tolerance > 100);
    $tolerance /= 100;
    print "Done\n";
}

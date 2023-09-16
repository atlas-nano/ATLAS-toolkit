#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use General qw(FileTester);


sub init;
sub getBestSnap;
sub createBestSnap;

my ($lammpsTrj, $field, $timeStep, $saveName, $dumpFreq, $lammpsLog, $target);
my ($SELECT, $scriptLoc, $bestTstep);

$|++;
$scriptLoc = &init;
print "Getting Best snapshot based on $field from $lammpsTrj...";
$bestTstep = &getBestSnap;
$bestTstep = $dumpFreq * int($bestTstep/$dumpFreq) + $dumpFreq;
print "Done\nCreating LAMMPS trajectory of timestep $bestTstep to $saveName...";
&createBestSnap($lammpsTrj, $bestTstep, $dumpFreq, $saveName);
print "Done\n";

sub createBestSnap {
    my ($lammpsFile, $searchTstep, $dumpFreq, $saveFile) = @_;
    my ($grepCmd, @lines, $i, $lineNum, $headCmd, $count);
    
    $count = `wc -l < $lammpsFile`;
    #$searchTstep += $dumpFreq;

    for $i ($searchTstep, ($searchTstep + $dumpFreq)) {
	$grepCmd = "grep -n '^${i}\$' ${lammpsFile}";
	open GREPCMD, "${grepCmd} |" or die "ERROR: Cannot execute $grepCmd: $!\n";
	$lineNum = -1;
	while (<GREPCMD>) {
	    if ($_ =~ /^(\d+)\:${i}/) {
		$lineNum = $1 - 2;
		push @lines, $lineNum;
	    }
	}
	close GREPCMD;
    }
    die "ERROR: Timestep $searchTstep not found in $lammpsFile\n" if (! @lines);
    $headCmd = "tail -" . ($count - $lines[0]) . " $lammpsFile";
    if ($#lines == 1) {
	$headCmd .= " | head -" . ($lines[1] - $lines[0]);
    }
    $headCmd .= "> ${saveFile}";
    die "ERROR while executing $headCmd\n" if (system($headCmd));
}

sub getBestSnap {
    my ($cmd, $bestVal);

    $cmd = "$scriptLoc -l $lammpsLog -f $field -s ${field}_profile.dat $timeStep > _get_best_snap.dat";
    die "ERROR occured. See _get_best_snap.dat\n" if (system($cmd));
    system("rm -fr _get_best_snap.dat");

    open INFILE, "${field}_profile.dat" or die "ERROR: Cannot open ${field}_profile.dat: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /BEST\s+(\d+)/ and ! defined($target)) {
	    $bestVal = $1;
	    last;
	} elsif ($_ =~ /^\s*(\d+\.?\d*)\s+(\-?\d+\.?\d*)/ and defined($target)) {
	    $bestVal = $1 if (! defined($bestVal) or (abs($2 - $target) < $bestVal));
        }
    }
    close INFILE;
    die "ERROR: No valid data found in ${field}_profile.dat!\n" if (! defined($bestVal));
    system("rm -fr ${field}_profile.dat");
    return $bestVal;
}

sub init {
    my (%OPTS, $SCRIPTLOC, $suffix, $junk);

    getopt('ltsfdav',\%OPTS);
    die "usage: $0 -l lammps trj -d dump frequency -a lammps log file\n" . 
	" -f (field: volume default) -t (timestep range: all default)-s (save name) -v (target val = avg)\n"
	if (! exists($OPTS{l}) or ! exists($OPTS{d}) or ! exists($OPTS{a}));
    print "Initializing...";
    ($lammpsTrj, $field, $timeStep, $saveName, $dumpFreq, $lammpsLog, $target) = 
	($OPTS{l}, $OPTS{f}, $OPTS{t}, $OPTS{s}, $OPTS{d}, $OPTS{a}, $OPTS{v});
    FileTester($lammpsTrj);
    FileTester($lammpsLog);
    die "ERROR: Expected integer for decimal fequency. Got \"$OPTS{d}\"\n"
	if ($dumpFreq !~ /^\d+$/);
    if (defined($timeStep)) {
	$timeStep = "-t ${timeStep}";
    } else {
	$timeStep = "";
    }
    if (! defined($field) ) {
	$field = "volume";
    }
    if (! defined($saveName)) {
	if ($lammpsTrj =~ /(\w+)\.(\w+)$/) {
	    $saveName = "${1}_best_${field}.${2}";
	} else {
	    $saveName = basename($saveName);
	    $saveName =~ s/\.\w+$//;
	    $saveName .= "_best_${field}.lammpstrj";
	}
    }
    $SCRIPTLOC = "/home/tpascal/scripts/getLammpsLogInfo.pl";
    FileTester($SCRIPTLOC);
    if (defined($target)) {
	undef($target) if ($target !~ /^\-?\d+\.?\d*$/);
    }
    print "Done\n";
    return $SCRIPTLOC;
}

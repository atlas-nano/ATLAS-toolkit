#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use FileFormats qw(GetBGFFileInfo sortByRes);

sub init;
sub matchFiles;
sub getData;
sub writeStats;

my ($FILES, $DATA, $bgfFile, $saveName,, $ATOMS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, undef) = GetBGFFileInfo($bgfFile, 0);
$RES = sortByRes($ATOMS);
print "Done\nGetting data...";
$DATA = getData($FILES, $ATOMS, $RES);
print "Done\nWriting stats to $saveName...";
&writeStats($DATA, $saveName);
print "Done\n";

sub getData {
    my ($datafiles, $atoms, $res) = @_;
    my ($i, $groups, $groupID, $headers, $thermo, $RES);

    for $i (keys %{ $datafile }) {
	$groups = ();
	open GRPFILE, $datafile->{$i}{GROUP} or die "ERROR: Cannot open group file $datafile->{$i}{GROUP}: $!\n";
	while (<GRPFILE>) {
	    chomp;
	    if ($_ =~ /^(\d+)\s+RES\s+(\d+)/) {
		$groupID = $1 - 1;
		@tmp = keys %{ $res->{$2}{ATOMS} };
		$name = $atoms->{$tmp[0]}{RESNAME} . "$2";
		$groups->{$groupID} = $name;
	    } elsif ($_ =~ /^(\d+)\s+(.+)/) {
		$groups->{$1 - 1} = $2;
	    }
	}
	close GRPFILE;
	die "ERROR: Not a valid file $datafile->{$i}{GROUP}\n" if (! $groups);
	$thermo = ();
	open THERMOFILE, $datafile->{$i}{THERMO} or die "ERROR: Cannot open thermo file $datafile->{$i}{THERMO}: $!\n";
	while (<THERMOFILE>) {
	    chomp;
	    if ($_ =~ /^Group/) {
		$index = 0;
		while ($_ =~ /(\w+)/g) {
		    $thermo->{GROUPS}{$index}{$1} = 1;
		    $index++;
		}
	    } elsif ($_ =~ /^(\d+)/) {
		$index = 0;
		$j = $1;
		while ($_ =~ /(\-?\d+\.?\d*)\s+\d+\.?\d*/g) {
		    $thermo->{DATA}{ $groups->{$j} }{ $thermo->{GROUPS}{$j}{$index} } = $1;
		}
	    }
	}
	close THERMOFILE;
	die "ERROR: Not a valid thermo file $datafile->{$i}{THERMO}: $!\n" if (! $thermo);
    }

    return $thermo;
}

sub init {
    my (%OPTS, $findCmd);
    my ($grpFiles, $thermoFiles, $saveName, $dir, $fileList, $i);

    getopt('gtsb',\%OPTS);
    for $i ("g", "t", "b") {
	die "usage: $0 -g \"group files\" -t \"thermo files\" -s (save name)\n"
	    if (! exists($OPTS{$i}));
    }
    print "Initializing...";
    $bgfFile = $OPTS{b};
    die "ERROR: Cannot access bgf file $bgfFile: $!\n" if (! -e $bgfFile || ! -r $bgfFile || ! -T $bgfFile);
    for $i ("g", "t") {
	$findCmd = "find " . dirname($OPTS{$i}) . " -name '" . basename($OPTS{$i}) . "' -print" if (! -e $OPTS{$i});
	$findCmd = "ls $OPTS{$i}" if (-e $OPTS{$i});
	die "ERROR: No valid file found while searching $OPTS{$i}!\n"
	    if (! open(DATFILES, "$findCmd |"));
	while (<DATFILES>) {
	    chomp;
	    if (-e $_ and -r $_ and -T $_) {
		$dir = dirname($_);
		#next if ($dir ne ".");
		$FILES->{$i}{$_} = 1;
	    }
	}
	close DATFILES;
	die "ERROR: No valid files found in path \"$fileList\"\n"
	    if (! defined($FILES->{$i}));
    }
    $FILES->{LIST} = matchFiles($FILES);
    $saveName = $OPTS{s};
    $saveName = basename($FILES->{LIST}{1}) if (! defined($saveName));
    $saveName =~ s/\.\w+$//;
    $saveName .= "_vac_thermo.dat";
    print "Done\n";
}

sub matchFiles {
    my ($files) = $_[0];
    my ($i, $lists, $j, $index, $jndex);
    
    for $i (keys %{ $files->{g} }) {
	next if ($i !~ /(\d+)/);
	$index = $1;
	for $j (keys %{ $files->{t} }) {
	    next if ($j !~ /(\d+)/);
	    $jndex = $1;
	    next if ($index != $jndex);
	    $lists->{$index} = {"GROUP" => $i, "THERMO" => $j };
	    delete $files->{t}{$j};
	}
    }
    die "ERROR: No matching files found!\n" if (! $lists);
    
    return $lists;
}

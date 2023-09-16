#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester GetSelections GetBondLength GetStats);
use ManipAtoms qw(GetAtmList);
use FileFormats qw(GetBGFFileInfo);

sub init;
sub getInteractions;
sub getMolData;
sub numerically { ($a cmp $b ); }

my ($FILES, $molInfo, $accum);
my ($ATOMS);

$|++;

&init;
print "Calculating interactions...";
&getInteractions($FILES, $molInfo);

sub getInteractions {
    my ($files, $molecule) = @_;
    my ($i, $ATOMS, $mol1, $mol2, $j, $k, $HEADERS); 
    my ($DATA, $t1, $t2, $STATS, $dist, $count);

    $count = 0;
    for $k (@{ $files }) {
	($ATOMS, undef) = GetBGFFileInfo($k, 0);
	($mol1, $mol2) = getMolData($ATOMS, $molecule);
	for $i (keys %{ $mol1 }) {
	    for $j (keys %{ $mol2 }) {
		next if ($i == $j);
		$dist = GetBondLength($ATOMS->{$i}, $ATOMS->{$j});
		$t1 = $ATOMS->{$i}{FFTYPE};
		$t2 = $ATOMS->{$j}{FFTYPE};
		($t1, $t2) = ($t2, $t1) if ($t1 lt $t2);
		push @{ $DATA->{$count}{VALS}{$t1}{$t2} }, $dist;
		$HEADERS->{$t1}{$t2} = 1;
	    }
	    $DATA->{$count}{FILE} = $k;
	}
	$count++;
    }
    print "Done\n";

    for $count (keys %{ $DATA }) {
	print $DATA->{$count}{FILE} . "\n" if (! $accum);
	for $i (keys %{ $DATA->{$count}{VALS} }) {
	    for $j (keys %{ $DATA->{$count}{VALS}{$i} }) {
		$STATS->{$count}{$i}{$j} = GetStats($DATA->{$count}{VALS}{$i}{$j});
		if (! $accum) {
		    printf "%-8s%-8s%12.3f%8.3f%8d\n", $i, $j, $STATS->{$count}{$i}{$j}{AVG}, 
		    $STATS->{$count}{$i}{$j}{STDEV}, $STATS->{$count}{$i}{$j}{NUM};
		}
	    }
	}
    }
    return "" if (! $accum);
    printf "%-30s", "FILE";
    for $i (sort numerically keys %{ $HEADERS }) {
	for $j (sort numerically keys %{ $HEADERS->{$i} }) {
	    printf "      %4s:%-4s      ", $i,$j;
	}
    }
    print "\n";
    for $k (sort {$a<=>$b} keys %{ $DATA }) {
	printf "%-30s", $DATA->{$k}{FILE};
	for $i (sort numerically keys %{ $HEADERS }) {
	    for $j (sort numerically keys %{ $HEADERS->{$i} }) {
		if (! exists($DATA->{$k}{VALS}{$i}{$j})) {
		    printf "%8s%5s%7s", "-", "-", "-";
		} else {
		    printf "%8.3f%5.1f%7d", $STATS->{$k}{$i}{$j}{AVG}, 
		    $STATS->{$k}{$i}{$j}{STDEV}, $STATS->{$k}{$i}{$j}{NUM};
		}
	    }
	}
	print "\n";
    }
}

sub getMolData {
    my ($ATOMS, $molecule) = @_;
    my ($mol1, $mol2, $i);

    if (defined($molecule)) {
	$mol1 = GetSelections($molecule, 0);
	$mol1 = GetAtmList($mol1, $ATOMS);
	if (! keys %{ $mol1 }) {
	    $mol1 = $mol2 = $ATOMS;
	} else {
	    for $i (keys %{ $ATOMS }) {
		$mol2->{$i} = 1 if (! exists($mol1->{$i}));
	    }
	}
    } else {
	$mol1 = $mol2 = $ATOMS;
    }
    
    return ($mol1, $mol2);
}

sub init {
    my (%OPTS, $mol1, $findCmd, $fileLoc);

    getopt('bac',\%OPTS);
    die "usage: $0 -b bgffile|location -a (mol1 atoms) -c (accumulate = no)\n" 
	if (! exists($OPTS{b}));
    print "Initializing...";
    ($fileLoc, $mol1, $accum) = ($OPTS{b}, $OPTS{a}, $OPTS{c});
    $accum = 0 if (! defined($accum) or $accum !~ /1|yes/i);
    $accum = 1 if ($accum =~ /1|yes/i);

    if (-e $fileLoc and ! -d $fileLoc) {
	FileTester($fileLoc);
	push @{ $FILES }, $fileLoc;
    } elsif (-d $fileLoc) {
	$findCmd = "find $fileLoc -name '*.bgf' -print";
	if (open(FINDCMD, "$findCmd |")) {
	    while (<FINDCMD>) {
		chomp;
		push @{ $FILES }, $_;
	    }
	    close FINDCMD;
	} else {
	    die "ERROR: No valid bgf files found while searching $fileLoc\n";
	}
    } else {
	$findCmd = "ls $fileLoc | grep '.bgf'";
	if (open(FINDCMD, "$findCmd |")) {
	    while (<FINDCMD>) {
		chomp;
		push @{ $FILES }, $_;
	    }
	    close FINDCMD;
	} else {
	    die "ERROR: No valid bgf files found while searching $fileLoc\n";
	}
    }

    if (defined($mol1)) {
	if ($mol1 =~ /\s+/) {
	    @{ $molInfo->[0] } = split /\s+/, $mol1;
	} else {
	    $molInfo->[0] = $mol1;
	}
    }
    print "Done\n";
}

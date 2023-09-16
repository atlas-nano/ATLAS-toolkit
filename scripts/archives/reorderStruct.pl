#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General qw(FileTester);
use Packages::FileFormats qw(createPDB createMOL2 createBGF sortByRes
			     GetBGFFileInfo GetMOL2FileInfo GetPDBFileInfo);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub numerically { ($a<=>$b); }
sub getAtomOrder;
sub reorderAtoms;

my ($structFile, $saveName, $ATOMS, $BONDS, $HEADERS, $ATMLIST); 
my ($readFunc, $writeFunc, $ORDER, $fileType, $orderType);

$|++;
&init;
print "Parsing $fileType file $structFile...";
($ATOMS, $BONDS, $HEADERS) = $readFunc->($structFile);
print "Done\n";
$ATMLIST = &getAtomOrder($ATOMS);
print "Reordering atoms...";
($ATOMS, $BONDS) = &reorderAtoms($ATOMS, $BONDS, $ORDER, $ATMLIST);
print "Done\nCreating $fileType file $saveName...";
$writeFunc->($ATOMS, $BONDS, $saveName);
print "Done\n";

sub reorderAtoms {
    my ($atoms, $bonds, $orderRequested, $orderCurrent) = @_;
    my ($atoms_mod, $bonds_mod, $i, $j, $atomCount, $resCount, $resID, $k);


    $atomCount = $resCount = $resID = 0;
    for $i (@{ $orderRequested }) {
	next if (! exists($orderCurrent->{$i}));
	for $j (sort numerically keys %{ $orderCurrent->{$i}{ATOMS} }) {
	    $atomCount++;
	    if ($resID != $atoms->{$j}{RESNUM}) {
		$resCount++;
		$resID = $atoms->{$j}{RESNUM};
	    }
	    %{ $atoms_mod->{$atomCount} } = %{ $atoms->{$j} };
	    $atoms_mod->{$atomCount}{OINDEX} = $j;
	    $atoms_mod->{$atomCount}{RESNUM} = $resCount;
	    $atoms->{$j}{NINDEX} = $atomCount;
	}
    }

    for $i (keys %{ $atoms_mod }) {
	$k = $atoms_mod->{$i}{OINDEX};
	$bonds_mod->{$i} = ();
	next if (! exists($bonds->{$k}));
	for $j (@{ $bonds->{$k} }) {
	    next if (! exists($atoms->{$j}{NINDEX}));
	    push @{ $bonds_mod->{$i} }, $atoms->{$j}{NINDEX};
	}
    }

    return ($atoms_mod, $bonds_mod);
}

sub getAtomOrder {
    my ($atoms) = $_[0];
    my ($INDEX, @tmp);

    if ($orderType == 1) {
	print "Sorting atoms by residue...";
	$INDEX = sortByRes($atoms);
    } else {
	@tmp = sort numerically keys %{ $atoms };
	for (@tmp) {
	    $INDEX->{$_}{ATOMS}{$_} = 1;
	}
    }

    return $INDEX;
}

sub init {
    my (%OPTS, $order, $tmp);

    getopt('fsotr',\%OPTS);
    
    for ("f", "o") {
	die "usage: $0 -f structure file -o order \n" . 
	    "\t-t [bgf(default)|pdb|mol2] -r [residue(default)|atom] -s [save name]\n" 
	    if (! defined($OPTS{$_}));
    }

    ($structFile, $saveName, $fileType, $order, $orderType) = ($OPTS{f}, $OPTS{s}, $OPTS{t}, $OPTS{o}, $OPTS{r});

    print "Initializing...";
    FileTester($structFile);

    while ($order =~ /(\S+)/g) {
	$tmp = $1;
	if ($tmp =~ /(\d+)\-(\d+)/) { # range 
	    if ($1 > $2) {
		for (reverse $2 .. $1) {
		    push @{ $ORDER }, $_;
		}
	    } else {
		for ($1 .. $2) {
		    push @{ $ORDER }, $_;
		}
	    }
	} elsif ($tmp =~ /(\d+)/) {
	    push @{ $ORDER }, $1;
	}
    }

    die "ERROR: No valid order found in $order\n" if (! @{ $ORDER });

    $fileType = "bgf" if (! defined($fileType));

    if ($fileType =~ /bgf|mol2|pdb/i) {
	$fileType = uc($fileType);
	$readFunc = eval ('\&Get' . $fileType . 'FileInfo');
	$writeFunc = eval('\&create' . $fileType);
    } else {
	die "ERROR: Expected bgf|mol2|pdb while parsing filetype. Got $fileType\n";
    }

    if (! defined($orderType) || $orderType !~ /atom/i) {
	$orderType = 1;
    } else {
	$orderType = 2;
    }

    if (! defined($saveName)) {
	$saveName = basename($saveName);
	if ($saveName = /^(.*)\.(\w+)/) {
	    $saveName = "${1}_reorder.${2}";
	}
    }
    print "Done\n";
}


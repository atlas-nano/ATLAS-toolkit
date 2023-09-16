#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use Storable qw(dclone);
use General qw(CoM FileTester);
use FileFormats qw(GetBGFFileInfo);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString GetAtmData);

sub init;
sub ParseJagInput;
sub updateCOORDs;
sub scaleCoM;
sub writeJagInput;
sub numerically { ($a<=>$b); }

my ($bgfFile, $jagInput, $savePrefix, $COM, $solu, $solv);
my ($ATOMS, $BONDS, $SOLVENT, $SOLUTE, $COORDS, $DATA);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
&GetMols($ATOMS, $BONDS);
$SOLUTE = SelectAtoms($solu, $ATOMS);
$SOLVENT = GetMols(GetAtmData($ATOMS, SelectAtoms($solv, $ATOMS)), $BONDS);
print "Done\nParsing Jaguar Input file $jagInput...";
$COORDS = ParseJagInput($jagInput);
&updateCOORDs($ATOMS, $COORDS);
print "Done\nGenerating " . scalar(@{ $COM }) . " scaled input files...";
$DATA = scaleCoM($ATOMS, $SOLUTE, $SOLVENT, $COM);
print "Done\nWriting Jaguar Input files...";
&writeJagInput($DATA, $jagInput, $savePrefix);
print "Done\n";

sub scaleCoM {
    my ($atoms, $solute, $solvent, $scaleF) = @_;
    my ($data, $i, $j, $com, $mol, $k, $l, $offset, $snap);

    for $i (@{ $scaleF }) {
	$data->{$i} = dclone($atoms);
	for $j (keys %{ $solvent }) {
	    $mol = GetAtmData($ATOMS, $solvent->{$j}{MEMBERS});
	    $com = CoM($mol);
	    for $k (keys %{ $com }) {
		$offset = ($i - 1) * $com->{$k};
		for $l (keys %{ $mol }) {
		    $data->{$i}{$l}{$k} += $offset;
		}
	    }
	}
    }
    return $data;
}

sub writeJagInput {
    my ($data, $jag, $prefix) = @_;
    my ($i, $j, $header);

    open JAGIN, "$jag" or die "ERROR: Cannot open $jag: $!\n";
    while (<JAGIN>) {
	last if ($_ =~ /&zmat/);
	$header .= "$_";
    }
    close JAGIN;

    for $i (keys %{ $data }) {
	open JAGFILE, "> ${prefix}.${i}.in" or die "ERROR: Cannot create ${prefix}.${i}.in: $!\n";
	print JAGFILE "${header}" . '&zmat' . "\n";
	for $j (sort numerically keys %{ $data->{$i} }) {
	    print JAGFILE "$data->{$i}{$j}{ATMNAME} $data->{$i}{$j}{XCOORD} $data->{$i}{$j}{YCOORD} $data->{$i}{$j}{ZCOORD}\n";
	}
 	print JAGFILE '&' . "\n";
	close JAGFILE;
    }
}

sub updateCOORDs {
    my ($bgf, $jag) = @_;
    my ($i, $dim, $atmName);

    for $i (keys %{ $bgf }) {
        $atmName = $bgf->{$i}{ATMNAME};
	$atmName =~ s/^\s+//;
	$atmName =~ s/\s+$//;
	delete $bgf->{$i}{MOLECULE};
	delete $bgf->{$i}{MOLECULEID};
	delete $bgf->{$i}{MOLSIZE};
        die "ERROR: Atom $i ($atmName) not found in Jaguar Input file!\n"
            if (! exists($jag->{$atmName}));
        for $dim ("XCOORD", "YCOORD", "ZCOORD") {
            $bgf->{$i}{$dim} = $jag->{$atmName}{$dim};
        }
    }
}

sub ParseJagInput {
    my ($inputFile) = $_[0];
    my (%DATA);

    open JAGINPUT, $inputFile or die "ERROR: Cannot open Jaguar Input file $inputFile: $!\n"; 
    while (<JAGINPUT>) {
        chomp; 
        if ($_ =~ /^\s*(\w+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
            $DATA{$1}{XCOORD} = $2;
            $DATA{$1}{YCOORD} = $3;
            $DATA{$1}{ZCOORD} = $4;
        }
    }
    close JAGINPUT;
    
    die "ERROR: No valid data found while reading $inputFile!\n"
        if (! %DATA);
    return \%DATA;
}

sub init {
    my (%OPTS, $comStr);
    
    getopt('bjuvsc',\%OPTS);
    ($bgfFile, $jagInput, $solu, $solv, $savePrefix, $comStr) = 
	($OPTS{b},$OPTS{j}, $OPTS{u}, $OPTS{v}, $OPTS{s}, $OPTS{c});
    
    die "usage: $0 -b bgf file -j jag input file -c com scale factor(s) -u (solute atoms) -v (solvent atoms) -s (save prefix)\n" 
	if (! defined($bgfFile) or ! defined($jagInput) or ! defined($comStr));

    print "Initializing...";
    FileTester($bgfFile);
    FileTester($jagInput);
    while ($comStr =~ /(\d+\.?\d*)/g) {
	push @{ $COM }, $1;
    }
    die "ERROR: Expected integer or decimal for com scale factor, got \"$comStr\"!\n"
	if (! $COM);
    if (! defined($savePrefix)) {
	$savePrefix = basename($bgfFile);
	$savePrefix =~ s/\.\w+$//;
    }
    $solu = "moleculeid==1" if (! defined($solu));
    $solv = "moleculeid>1" if (! defined($solv));
    $solu = BuildAtomSelectionString($solu);
    $solv = BuildAtomSelectionString($solv);
    print "Done\n";
}

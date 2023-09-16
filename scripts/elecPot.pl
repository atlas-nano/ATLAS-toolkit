#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use FileFormats qw(GetBGFFileInfo createPQR);
use General qw(ReadFFs FileTester LoadFFs);
use ManipAtoms qw(GetMols);

sub init;
sub computeElecEnergyPerAtom;
sub writeElecPot;
sub numerically { ($a<=>$b); }
sub addRadii;

my ($bgfFile, $ffFiles, $saveFile);
my ($ATOMS, $BONDS, $HEADERS, $elecPot, $PARMS);

$|++;
&init;
$PARMS = LoadFFs($ffFiles);
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
&GetMols($ATOMS, $BONDS);
for (keys %{ $ATOMS }) {
    &addRadii($ATOMS->{$_}, $PARMS);
}
print "Done\nCreating PQR file...";
&createPQR($ATOMS, $BONDS, "_tmp.pqr");
print "Done\nCalculating Electrostatic Potential...";
$elecPot = &computeElecEnergyPerAtom($ATOMS, "_tmp.pqr");
print "Done\nWriting results to $saveFile...";
&writeElecPot($ATOMS, $elecPot, $saveFile);
print "Done\n";

sub writeElecPot {
    my ($atoms, $molEng, $saveName) = @_;
    my ($i);

    open OUTFILE, "> $saveName" || die "ERROR: Cannot write to $saveName: $!\n";
    #write the per atom elec pot
    for $i (sort numerically keys %{ $atoms }) {
	printf OUTFILE "ATOMS %-8d %12.5f\n",$i,$atoms->{$i}{ENERGY};
    }
    #write the per molecule elec pot
    print OUTFILE "\n";
    for $i (sort numerically keys %{ $molEng }) {
        printf OUTFILE "MOL %-8d %12.5f\n",$i,$molEng->{$i};
    }
    close OUTFILE;
}

sub computeElecEnergyPerAtom {
    my ($atoms, $pqrFile) = @_;
    my ($apbsCmd, $isValid, %ENERGIES, $molID);

    $apbsCmd = "/ul/tpascal/programs/bin/coulomb -e $pqrFile";
    open COULCMD, "$apbsCmd |" or die "ERROR: Cannot run cmd $apbsCmd: $!\n";
    while (<COULCMD>) {
        chomp;
        if ($_ =~ /Atom\s+(\d+):\s+Energy\s+\=\s+(\-?\d+\.\d+E?.?\d*)/) {
            $molID = $atoms->{$1}{MOLECULEID};
            $ENERGIES{$molID} += $2;
            $atoms->{$1}{ENERGY} = $2;
        }
    }
    die "ERROR: No valid energies obtained from $apbsCmd!\n" if (! %ENERGIES);

    system("rm -fr _tmp.pqr io.mc");
    return \%ENERGIES;
}

sub addRadii {
    my ($atom, $parms) = @_;
    my ($i, $ffType);

    $ffType = $atom->{FFTYPE};
    die "ERROR: Atomtype $ffType not found in forcefield(s)!\n"
	if (! exists($parms->{VDW}{$ffType}{$ffType}));
    $atom->{RADII} = $parms->{VDW}{$ffType}{$ffType}{1}{VALS}[1];
}

sub init {
    my (%OPTS, $forceFields, $ffType);

    getopt('bfs',\%OPTS);
    die "usage: $0 -b bgf file -f force field -s (save file)\n" 
	if (! exists($OPTS{b}) || ! exists($OPTS{f}));

    print "Initializing...";
    ($bgfFile, $forceFields, $saveFile) = ($OPTS{b}, $OPTS{f}, $OPTS{s});

    ($ffFiles, $ffType) = ReadFFs($forceFields);

    FileTester($bgfFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_elecPot.dat";
    }
    print "Done\n";
} 

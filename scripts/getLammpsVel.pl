#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo  GetBGFAtoms);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsTrjType GetLammpsByteOffset);
use General qw(FileTester TrjSelections GetSelections);
use AMBER qw(CreateAmberTrj);
use ManipAtoms qw(GetAtmList GetMols SelectAtoms BuildAtomSelectionString);
use File::Basename;
use Getopt::Std qw(getopt);

my ($velFile, $saveName, $saveType, $SELECT, $pStr, $bgfFile, $selection);
my ($ATOMS, $BONDS, $SOLVENT, $velFactor, $isShake, $ATOMSELECT);

sub init;
sub saveVel;
sub getShakeDOF;

$|++;
&init;
$velFactor = 1;
if (defined($bgfFile) and -e $bgfFile and -r $bgfFile and -T $bgfFile) {
    print "Parsing bgf file $bgfFile...";
    ($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
    &GetMols($ATOMS, $BONDS);
    $velFactor = getShakeDOF($ATOMS) if ($isShake);
    print "Done\nParsing atom/residue selection...";
	$ATOMSELECT = SelectAtoms($selection, $ATOMS);
	print "Done\n";
}

$velFactor *= 1000/20.455 if ($saveType ne "lammps");

$pStr = "Parsing LAMMPS Velocity file $velFile...";
if ($saveType eq "lammps") {
    open OUTDATA, "> $saveName" || die "ERROR: Cannot create $saveName: $!\n";
    print OUTDATA "Velocities\n\n";
} elsif ($saveType eq "amber") {
    open OUTDATA, "> $saveName" || die "ERROR: Cannot create $saveName: $!\n";
    print OUTDATA "AMBER velocity file created by $ENV{USER} at " . 
	scalar(localtime) . "\n";
}
&GetLammpsByteOffset($SELECT, $velFile, scalar keys %{ $ATOMS });
#&GetLammpsTrjType($SELECT, $velFile, undef);
&ParseLAMMPSTrj($saveType, $velFile, $SELECT, "vel", \&saveVel, $pStr, \*OUTDATA);
if ($saveType ne "data") {
  close OUTDATA or die "ERROR: Cannot close $saveName: $!\n";
}

sub getShakeDOF {
    my ($atoms) = $_[0];
    my ($i, $dof, $factor, $atmCount);
    print "adjusting velocities for dof removed by shake...";
    $atmCount = scalar keys %{ $atoms };
    $dof = $atmCount * 3;

    for $i (keys %{ $atoms }) {
	if ($atoms->{$i}{FFTYPE} =~ /^H/) { #if this is a hydrogen
	    $dof--;
	    if (${ $atoms->{$i}{MOLSIZE} } == 3) {
		$dof -= 0.5;
	    }
	}
    }

    $factor = sqrt($atmCount * 3 / $dof);

    return $factor;
}

sub saveVel {
    my ($DATA, $type, $fileHandle) = @_;
    my ($i, $j, $BOX, $count, $datFile);

    #%{ $ATOMSELECT } = (keys %{ $DATA->{ATOMS} }) if (! defined($ATOMSELECT));
    $count = 1;
    if ($type eq "data") {
      $datFile = "${saveName}." . $DATA->{TIMESTEP}[0] . ".vel";
      open DATAFILE, "> $datFile" or die "ERROR: Cannot write to $datFile: $!\n";
      print DATAFILE "Velocities\n\n"; 
      $fileHandle = \*DATAFILE;
    }
    if ($type =~/lammps|data/) {
	for $i (sort {$a<=>$b} keys %{ $ATOMSELECT }) {
	    next if (! exists($DATA->{ATOMS}{$i}));
	    printf $fileHandle "%8d", $i;
	    for $j ("XVEL", "YVEL", "ZVEL") {
		printf $fileHandle "%12.8f", $DATA->{ATOMS}{$i}{$j} * $velFactor;
	    }
	    print $fileHandle "\n";
	    $count++;
	}
    } else {
	for $i (1 .. $DATA->{"NUMBER OF ATOMS"}[0]) {
	    for $j ("X", "Y", "Z") {
		$DATA->{ATOMS}{$i}{$j . "COORD"} = $DATA->{ATOMS}{$i}{$j . "VEL"} * $velFactor;
	    }
	}
	CreateAmberTrj($DATA->{ATOMS},$BOX, $fileHandle);
    }
    close DATAFILE if ($type eq "data");
}

sub init {
    my (%OPTS, $selections, $atomSelect);
    
    getopt('vtsfbia',\%OPTS);
    ($velFile, $selections, $saveType, $saveName, $bgfFile, $isShake, $atomSelect) = 
	($OPTS{v}, $OPTS{f}, $OPTS{t}, $OPTS{s}, $OPTS{b}, $OPTS{i}, $OPTS{a});
    die "usage:$0 -v lammps velocity file -f trajectory frames (* for all)\n -b [bgf file] -a (atom selection)" . 
	"\t-t ouput type [=lammps|amber] (optional) " . "-s [saveName] (optional) -i (is shake = no)\n"
	if (! defined($velFile) || ! defined($selections));

    print "Initializing...";
    FileTester($velFile);
    if (! defined($saveName)) {
	$saveName = basename($velFile);
	$saveName =~ s/\.\w+$//;
	$saveName .= "_lammps.vel";
    }
    $saveType = "lammps" if (! defined($saveType) || $saveType !~ /^(lammps|amber|data)/i);
    $saveType = lc($saveType);
    if($saveType eq "data") {
      my $tmp = $saveName;
      $tmp =~ s/\.\w+$//;
      $saveName = dirname($tmp) . "/" . "data." . basename($tmp);
    }
    $SELECT = TrjSelections($selections);
	$atomSelect = "index>0" if (! defined($atomSelect));
	$selection = BuildAtomSelectionString($atomSelect);

    $isShake = 0 if (! defined($isShake) or $isShake !~ /^(1|yes)$/i);
    print "Done\n";
}

#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use General qw(FileTester GetSelections);
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use ManipAtoms qw(GetAtmList);

sub init;
sub createScreamCmd;
sub runScream;
sub numerically { ($a<=>$b); }

my ($bgfFile, $selection, $saveName, %short);
my ($ATOMS, $BONDS, $HEADERS, $screamCmd, $SELECT, $RES);

$|++;
print "Initializing...";
&init;
print "Done\nParsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$SELECT = GetSelections($selection, 0);
$SELECT = GetAtmList($SELECT, $ATOMS);
print "Done\nRunning Scream...";
$screamCmd = createScreamCmd($ATOMS, $SELECT, \%short);
&runScream($bgfFile, $screamCmd);
print "Done\nSaving BGF file $saveName...";
($ATOMS, $BONDS, undef) = GetBGFFileInfo("best_1.bgf");
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
system("rm -fr scream.par Anneal-Energies.txt Field1.bgf timing.txt best_1.bgf Residue-E.txt scream.out");
print "Done\n";

sub runScream {
    my ($bgfName, $mutate) = @_;
    my ($screamExec) = "/project/Biogroup/Software/SCREAM/builds/current_build/SCREAM_wrap.py";
    my ($screamCmd);
    my ($prep) = "/project/Biogroup/Software/SCREAM/builds/current_build/utilities/Make_SCREAM_compatible_Bgf.py";

    $screamCmd = "$prep $bgfName _scream.bgf";
    die "ERROR: Cannot run $screamCmd!\n" if (system($screamCmd) or ! -e "_scream.bgf");

    $screamCmd = "$screamExec _scream.bgf 10 FLAT 0.4 $mutate >& _screamout";

    die "ERROR while executing SCREAM! See _screamout\n" if (system($screamCmd) or ! -e "best_1.bgf");
    system "rm -fr _screamout _scream.bgf";
}

sub createScreamCmd {
    my ($atomData, $indices, $resLetters) = @_;
    my (%RES, $i, $mutateRes, $resNum, $resName, $rL);

    for $i (keys %{ $indices }) {
	$resNum = $atomData->{$i}{RESNUM};
	$resName = $atomData->{$i}{RESNAME};
	die "ERROR: Residue \"$resName\" is not a valid amino acid residue!\n"
	    if (! exists($resLetters->{uc($resName)}));
	#$rL = $resLetters->{uc($resName)};
	#$rL = "A";
	next if exists($RES{$resNum});
	$RES{$resNum} = "A${resNum}_" . $atomData->{$i}{CHAIN};
    }

    for $i (sort numerically keys %RES) {
	$mutateRes .= "$RES{$i} ";
    }

    return $mutateRes;
}

sub init {
    my (%OPTS, $select);
    getopt('bsr',\%OPTS);

    for ("b", "r") {
	die "usage: $0 -b bgf file -r reslist -s [save name]\n" if (! exists($OPTS{$_}));
    }

    ($bgfFile, $select, $saveName) = ($OPTS{b}, $OPTS{r}, $OPTS{s});
    FileTester($bgfFile);
    if ($select =~ /\s+/) {
        @{ $selection } = split /\s+/, $select;
    } else {
        $selection->[0] = $select;
    }
    if (! defined($saveName)) {
	$saveName = basename($bgfFile);
	$saveName =~ s/\.\w+$//;
	$saveName .= "_ala.bgf";
    }

    %short = (
              "ALA" => "A",
              "CYS" => "C",
              "CYX" => "C",
              "ASP" => "D",
              "GLU" => "E",
              "PHE" => "F",
              "GLY" => "G",
              "HIS" => "H",
              "HIE" => "H",
              "HSD" => "H",
              "HSE" => "H",
              "HSP" => "H",
              "ILE" => "I",
              "LYS" => "K",
              "LEU" => "L",
              "MET" => "M",
              "ASN" => "N",
              "PRO" => "P",
              "GLN" => "Q",
              "ARG" => "R",
              "SER" => "S",
              "THR" => "T",
              "VAL" => "V",
              "TRP" => "W",
              "TYR" => "Y" );

}

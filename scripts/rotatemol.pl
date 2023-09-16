#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(ParseStructFile createBGF addHeader);
use General qw(FileTester Rotate CoM GetFileTypeStr);
use ManipAtoms qw(GetMols BuildAtomSelectionString SelectAtoms);
use constant PI => atan2(1,1) * 4;

my ($struct_file, $atomSel, $comSel, $rAngles, $saveName, $writeJag, $OFFSET);
my ($ATOMS, $BONDS, $HEADERS, $aSELECT, $cSELECT, $MOLS, $savePrefix);

$|++;
&init;
print "Parsing BGF file $struct_file...";
($ATOMS, $BONDS, $HEADERS) = ParseStructFile($struct_file, 1);
$MOLS = GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$aSELECT = SelectAtoms($atomSel, $ATOMS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $aSELECT });
$cSELECT = SelectAtoms($comSel, $ATOMS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $cSELECT });
print "Done\nRotating Atoms...";
&rotateAtoms($ATOMS, $aSELECT, $cSELECT, $rAngles, $OFFSET);
print "Done\nCreating $saveName...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";
&writeJagInputFile($ATOMS,$MOLS,$savePrefix) if ($writeJag);

sub writeJagInputFile {
    my ($atoms, $mols, $prefix) = @_;
    my ($i, $curr, $aPrefix, $header, $atmfield);

    $header->{a} = <<DATA;
&gen
basis=cc-pvdz++
dftname=m06
DATA

    $header->{b} = <<DATA;
mulken=1
icfit=1
iacc=2
ip172=2
DATA

    $header->{c} = <<DATA;
&
&zmat
DATA

    print "Writing Jaguar input file ${prefix}.jagin...";
    open JAGINPFILE, "> ${prefix}.jag.in" or die "ERROR: Cannot write to ${prefix}.jag.in: $!\n";
    open JAGACPFILE, "> ${prefix}_a_cp.jag.in" or die "ERROR: Cannot write to ${prefix}_a_cp.jag.in: $!\n";
    open JAGBCPFILE, "> ${prefix}_b_cp.jag.in" or die "ERROR: Cannot write to ${prefix}_b_cp.jag.in: $!\n";

    print JAGINPFILE $header->{a} . $header->{b} . $header->{c};
    print JAGACPFILE $header->{a} . $header->{c};
    print JAGBCPFILE $header->{a} . $header->{c};

    for $i (sort {$a<=>$b} keys %{ $atoms }) {
	$curr = $atoms->{$i};
	$aPrefix = $curr->{ATMNAME};
	$aPrefix =~ /([a-zA-Z]+)/;
	$aPrefix = $1;
	$aPrefix =~ s/^\s+//;
	$aPrefix =~ s/\s+$//;
	$aPrefix .= $i;

	printf JAGINPFILE " %-10s %20.10f %20.10f %20.10f\n",$aPrefix,$curr->{XCOORD},$curr->{YCOORD},$curr->{ZCOORD};
	if(!exists($MOLS->{1}{MEMBERS}{$i})) {
	    printf JAGACPFILE " %-10s %20.10f %20.10f %20.10f\n",$aPrefix . "@",$curr->{XCOORD},$curr->{YCOORD},$curr->{ZCOORD};
	    printf JAGBCPFILE " %-10s %20.10f %20.10f %20.10f\n",$aPrefix,$curr->{XCOORD},$curr->{YCOORD},$curr->{ZCOORD};
	} else {
	    printf JAGBCPFILE " %-10s %20.10f %20.10f %20.10f\n",$aPrefix . "@",$curr->{XCOORD},$curr->{YCOORD},$curr->{ZCOORD};
	    printf JAGACPFILE " %-10s %20.10f %20.10f %20.10f\n",$aPrefix,$curr->{XCOORD},$curr->{YCOORD},$curr->{ZCOORD};
	}
    }
    print JAGINPFILE "&\n";
    print JAGACPFILE "&\n";
    print JAGBCPFILE "&\n";

    close JAGINPFILE;
    close JAGACPFILE;
    close JAGBCPFILE;
    print "Done\n";
}

sub rotateAtoms {
    my ($atoms, $aSel, $cSel, $rAngles, $offset) = @_;
    my ($i, $j, $CoM, $sAtoms, $cAtoms, $angles, $del);

    @{ $angles } = (0,0,0);
    $angles->[0] = $rAngles->{X} if (exists($rAngles->{X}));
    $angles->[1] = $rAngles->{Y} if (exists($rAngles->{Y}));
    $angles->[2] = $rAngles->{Z} if (exists($rAngles->{Z}));
    for $i (0 .. 2) {
		$angles->[$i] *= PI/180; #radians
    }

	for $i (keys %{ $cSel }) {
		$cAtoms->{$i} = $atoms->{$i}
	}
    $CoM = CoM($cAtoms);

	for $i (keys %{ $aSel }) {
		$sAtoms->{$i} = \%{ $atoms->{$i} };
		for $j ("XCOORD", "YCOORD", "ZCOORD") {
		    $sAtoms->{$i}{$j} -= $CoM->{$j};
		}
    }

    for $i (0 .. $#{ $angles }) {
		next if ($angles->[$i] == 0);
		&Rotate($sAtoms, $angles, $i);
    }

    for $j ("XCOORD", "YCOORD", "ZCOORD") {
		$del = 0;
		$del = $offset->{$j} if (defined($offset) and exists($offset->{$j}));
		for $i (keys %{ $cAtoms }) {
		    $sAtoms->{$i}{$j} += $CoM->{$j} + $del;
		}
    }
}

sub init {
    my (%OPTS, $atomSelstr, $offsetStr, $rAngleStr, $rotStr, $comSelstr);

    getopt('barsjoc',\%OPTS);

    ($struct_file, $atomSelstr, $rAngleStr, $saveName, $writeJag, $offsetStr, $comSelstr) =
	($OPTS{b}, $OPTS{a},    $OPTS{r},   $OPTS{s},  $OPTS{j},  $OPTS{o},   $OPTS{c});

	&usage if (! defined($struct_file) or ! defined ($rAngleStr));

    print "Initializing...";
    FileTester($struct_file);
    
	while ($rAngleStr =~ /(x|y|z):\s*(\-?\d+\.?\d*e?\-?\d*)/ig) {
		$rAngles->{uc $1} = $2;
		$rotStr .= "$1_$2";
    }
    die "ERROR: Expected X|Y|Z: xx.xxx for rotation angles. Got \"$rAngleStr\"\n"
		if (! defined($rAngles));

    $atomSelstr = "index>0" if (! defined($atomSelstr));
    $atomSel = BuildAtomSelectionString($atomSelstr);
    $comSelstr = "index>0" if (! defined($comSelstr));
    $comSel = BuildAtomSelectionString($comSelstr);

	$writeJag = 0 if (! defined($writeJag) or $writeJag !~ /1|yes/i);
    $writeJag = 1 if ($writeJag =~ /1|yes/i);
    if (!defined($saveName)) {
		$struct_file =~ /^\s*(\S+)/;
		$saveName = basename($1);
		$saveName =~ s/\.\w+$//;
		$saveName .= "${rotStr}.bgf";
    }
    if ($writeJag) {
	$savePrefix = $saveName;
	$savePrefix =~ s/\.\w+$//;
    }
    $offsetStr = "" if (! defined($offsetStr));
    while ($offsetStr =~ /(x|y|z):\s*(\+|\-)\s*(\d+\.?\d*e?\-?\d*)/gi) {
	$OFFSET->{uc $1 . "COORD"} = $2 . $3;
    }
    print "Done\n";
}

sub usage {
	my ($fTypeStr) = GetFileTypeStr;	
    die <<DATA
usage: $0 -b struct_file -r 'rot_angle(s)' -a (atom selection) -c (com selection) -o '(offset)' -j (jaguar_flag) -s (save name)
arguments:
	struct_file: (Required) location of structure file
$fTypeStr	
	rot_angle(s): (Required) x|y|z: xx.xx
	atom selection: (Optional) Field selection for atoms to rotate. Default: all
	com selection: (Optional) Field selection for center of mass atom to place at origin. Default: all
	offset: (Optional) x|y|z:+/- xx.xx. Default: none
	jaguar_flag: (Optional) Flag to write Jaguar input file. Default: no
	save name: (Optional). Name of BGF file to create. Dafault: (input_file_name).mod.bgf
DATA

}

#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(ParseStructFile addHeader createBGF AddMass insertHeaderRemark GetBGFAtoms);
use General qw(FileTester CoM LoadFFs ReadFFs);
use BOX qw(GetBox MoveAtomsToOrigin);
use Getopt::Std qw(getopt);
use ManipAtoms qw(CenterSystem SelectAtoms BuildAtomSelectionString);

sub init;
sub showUsage;
sub getCenter;

my ($bgfFile, $cerius2FF, $saveName, $centerType, $dimOpts, $refAtomStr, $selection);
my (%Offset, $dim, $atomC, $radii, $CENTER, $PARMS, $SELECT);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
my ($ATOMS, $BONDS, $HEADERS) = ParseStructFile($bgfFile, 1);
print "Done\n";
$PARMS = ();
if (defined($cerius2FF)) {
	print "Parsing CERIUS2 forcefield $cerius2FF...";
	$PARMS = LoadFFs($cerius2FF);
	AddMass($ATOMS, $PARMS);
	print "Done\n";
}
print "Selecting reference atoms based on ${refAtomStr}...";
$SELECT = SelectAtoms($selection, $ATOMS);
print "\nGetting box information...";
my ($BOX) = GetBox($ATOMS, $PARMS, $HEADERS);
$CENTER = getCenter($ATOMS, $BONDS, $BOX, $centerType, $dimOpts, $SELECT);
print "Done\nCentering atoms...";
if ($centerType eq "box_origin") {
	&MoveAtomsToOrigin($ATOMS);
} else {
	CenterSystem($ATOMS, $CENTER, ());
}
&insertHeaderRemark($HEADERS, "REMARK $bgfFile centered $centerType");
&addHeader($ATOMS, $HEADERS);
print "Done\nCreating BGF file $saveName...";
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub getCenter {
	my ($atoms, $bonds, $box, $centerOpt, $dim, $refSel) = @_;
	my (@dim, %CENTER, $i, $atomCenter, $refAtoms);

	($refAtoms, undef, undef) = GetBGFAtoms($refSel, $atoms, $bonds);
	$atomCenter = CoM($refAtoms);
	@dim = keys %{ $dim };
	if ($centerOpt eq "box_origin") {
		$box = GetBox($ATOMS, $PARMS, undef);
		for $i (@dim) {
			$CENTER{"${i}COORD"} = $box->{$i}{lo};
		}
	} elsif ($centerOpt eq "com_origin") {
		for $i (@dim) {
			$CENTER{"${i}COORD"} = $atomCenter->{"${i}COORD"};
		}
	} elsif ($centerOpt eq "box_center") {
		for $i (@dim) {
			 $CENTER{"${i}COORD"} = -1*$box->{$i}{len}/2 - $box->{$i}{lo};
		}
	} elsif ($centerOpt eq "com_center") {
		for $i (@dim) {
			 $CENTER{"${i}COORD"} = -($box->{$i}{len}/2 - $box->{$i}{lo} - $atomCenter->{"${i}COORD"});
		}
	} else {
		for $i (@dim) {
			$CENTER{"${i}COORD"} = $box->{$i}{len}/2 - $atomCenter->{"${i}COORD"} + $box->{$i}{lo};
			#$CENTER{"${i}COORD"} = $atomCenter->{"${i}COORD"};
		}
	}

	return \%CENTER;
}

sub init {
	my (%OPTS, $usage, $centerOpt, $ffData, $dimStr, $tmp);
	getopt('bfscdr',\%OPTS);
	($bgfFile, $ffData, $saveName, $centerOpt, $dimStr, $refAtomStr) = ($OPTS{b},$OPTS{f},$OPTS{s},$OPTS{c},$OPTS{d},$OPTS{r});
	$usage = &showUsage;
	for ($bgfFile) {
		die "$usage\n" if (! defined($_));
	}

	print "Initializing...";
	FileTester($bgfFile);
	($cerius2FF, undef) = ReadFFs($ffData) if (defined($ffData));
	$centerOpt = "box_origin" if (! defined($centerOpt));
	
	if ($centerOpt =~ /(box_origin|com_center|com_origin|box_center)/i) {
		$centerType = lc($1);
	}
	if (! $saveName) {
		$saveName = $bgfFile;
		$saveName =~ s/\.\w+$/_${centerType}\.bgf/;
	}

	$dimStr = "xyz" if (! defined($dimStr));
	while ($dimStr =~ /(x|y|z)/ig) {
		$tmp = uc $1;
		$dimOpts->{$tmp} = 1;
	}
	$refAtomStr = "index>0" if (! defined($refAtomStr));
	$selection = BuildAtomSelectionString($refAtomStr);
	die "ERROR: Invalid centering option \"${centerOpt}\"\n$usage\n" if (! defined($centerType));
}

sub showUsage {
	my ($usage) = "usage: $0 -b bgf file -f cerius2 forcefield -s [save name] -c [centering type] -r [reference atom selection]\n" .
		"options:\n\t-b bgf file: location of the bgf file\n" .
		"\t-f cerius2 forcefield: location of CERIUS2 formatted forcefield.\n\t\tBGF file must be typed with this forcefield\n" .
		"\t-s [save name]: (optional) the name to save the centered bgf. \n\t\tIf not specified will be {bgf file}_centered.bgf\n" .
		"\t-c [centering type]: (optional) the type of centering to perform.\n" .
		"\t\tChoices:\n\t\tcom_origin - center the center of mass of the molecule to the origin\n" .
		"\t\tcom_center - place the center of mass of the molecule in the center of the box\n" .
		"\t\tbox_origin - (default) place the lowest extrema of the molecule at the origin\n";
	return $usage;
}

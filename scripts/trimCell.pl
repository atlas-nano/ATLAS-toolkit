#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use FileFormats qw(GetBGFFileInfo addHeader createHeaders createBGF addBoxToHeader 
                   GetBGFAtoms insertHeaderRemark);;
use BOX qw(GetBox InitBox Map2UnitCell Cart2Frac);
use REPLICATE qw (GetBoxDisplacementTensor);
use ManipAtoms qw(SplitAtomsByMol GetAtmList GetMols SelectAtoms BuildAtomSelectionString 
                  GetAtmData ReimageAtoms AddMolsToSelection);
use General qw(FileTester CoM GetSelections);

sub init;
sub getNewBox;
sub trimCell;
sub makeAtomsMols;

my ($bgfFile, $saveFile, $newCell, $startOrigin, $isMol, $cStr, $atomSel);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $MOLS, $selection, $SELECT);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$BOX = GetBox($ATOMS, undef, $HEADERS);
if ($isMol) {
	$MOLS = &GetMols($ATOMS, $BONDS);
} else {
	$MOLS = makeAtomsMols($ATOMS);
}
print "Done\nParsing atom/residue selection...";
$SELECT = GetSelections($selection, 0);
$SELECT = GetAtmList($SELECT, $ATOMS);
print "Done\nRemoving atoms outside new box...";
&centerSys($ATOMS, $BOX, $startOrigin, $isMol);
&InitBox($BOX, $ATOMS);
&ReimageAtoms($ATOMS, $BONDS, $MOLS, $BOX, $SELECT);
&InitBox($newCell, $ATOMS);
&trimCell($ATOMS, $BONDS, $MOLS, $newCell, $SELECT);
($ATOMS, $BONDS) = GetBGFAtoms($ATOMS, $ATOMS, $BONDS);
print "Done\nCreating BGF file $saveFile...";
&insertHeaderRemark($HEADERS, "REMARK $bgfFile trimed $cStr");
&addBoxToHeader($HEADERS, $newCell);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub makeAtomsMols {
	my ($atoms) = $_[0];
	my ($i, %atomsMols);

	for $i (keys %{ $atoms }) {
		$atomsMols{$i}{$i} = 1;
	}

	return \%atomsMols;
}

sub trimCell {
	my ($atoms, $bonds, $mols, $box, $select) = @_;
	my ($i, $molCenter, $j, $k, $currMol, $isOutside);

	for $i (keys %{ $mols }) {
		$currMol = ();
		for $j (keys %{ $mols->{$i}{MEMBERS} }) {
			%{ $currMol->{$j} } = %{ $atoms->{$j} };
		}
		$molCenter = CoM($currMol);
		#&Map2UnitCell($molCenter, $box);
		$currMol = ();
		$currMol->{1} = $molCenter;
		&Cart2Frac($currMol, $box);
		$isOutside = 0;
		for $j ("FA", "FB", "FC") {
			$isOutside++ if($currMol->{1}{$j}<0 or $currMol->{1}{$j}>1);
		}
		for $j (keys %{ $mols->{$i}{MEMBERS} }) {
			if ($isOutside>0) {
				next if (! exists($select->{$j}));
				delete $atoms->{$j};
				delete $bonds->{$j};
			} else {
				for $k ("XCOORD", "YCOORD", "ZCOORD") {
					#$atoms->{$j}{$k} += $molCenter->{SHIFT}{$k};
				}
			}

		}
	}
}

sub centerSys {
	my ($atoms, $oldBox, $startPoint, $molOpt) = @_;
	my ($i, $j, $k, $selection, $SELECTIONS, $centerAtoms, $atomsShift); 

	for $i (keys %{ $oldBox }) {
		next if ($i =~ /^(X|Y|Z)$/);
		delete $oldBox->{$i};
	}
	if ($startPoint == -1) { #atom selection
		$selection = BuildAtomSelectionString($atomSel);
		$SELECTIONS = SelectAtoms($selection, $ATOMS);
		die "ERROR: No valid atoms selected for centering!\n"
			if(! keys %{ $SELECTIONS });
		&AddMolsToSelection($SELECTIONS, $ATOMS) if ($molOpt);
		$centerAtoms = GetAtmData($ATOMS, $SELECTIONS);
		$atomsShift = CoM($centerAtoms);
	} elsif ($startPoint == 1) { #origin (0,0,0)
		$atomsShift = {XCOORD=>0,YCOORD=>0,ZCOORD=>0};
	} elsif ($startPoint == 2) { # box min
		$atomsShift = {XCOORD=>$oldBox->{X}{lo}, YCOORD=>$oldBox->{Y}{lo}, ZCOORD=>$oldBox->{Z}{lo}};
	} elsif ($startPoint == 3) { # box max
		$atomsShift = {XCOORD=>$oldBox->{X}{hi}, YCOORD=>$oldBox->{Y}{hi}, ZCOORD=>$oldBox->{Z}{hi}};
	} else { # symmetric cut from box center
		$atomsShift = {XCOORD=>$oldBox->{X}{len}/2, YCOORD=>$oldBox->{Y}{len}/2, ZCOORD=>$oldBox->{Z}{len}/2};
	}

	&GetBoxDisplacementTensor($oldBox);
	for $i (keys %{ $atoms }) {
		for $j ("X", "Y", "Z") {
			for $k ("X", "Y", "Z") {
				$atoms->{$i}{"${j}COORD"} -= $oldBox->{$j}{DISP_V}{$k}/2;
			}
		}
	}
}

sub init {
	my (%OPTS, $atomStr, $tmp, $angle_str);
	getopt('bcsoma',\%OPTS);
	for ("b", "c") {
		if (! exists($OPTS{$_})) {
			&usage;
			die "\n";
		}
	}
	print "Initializing...";
	($bgfFile, $cStr, $atomStr, $saveFile, $startOrigin, $isMol) = ($OPTS{b}, $OPTS{c}, $OPTS{a}, $OPTS{s}, $OPTS{o}, $OPTS{m});
	FileTester($bgfFile);
	$atomStr = "*" if (! defined($atomStr));
	if ($atomStr =~ /\s+/) {
		@{ $selection } = split /\s+/, $atomStr;
	} else {
		$selection->[0] = $atomStr;
	}

	$atomSel = $startOrigin;
	$startOrigin = 0 if (! defined($startOrigin));
	$startOrigin = 1 if ($startOrigin =~ /^(1|yes)/i);
	$startOrigin = -1 if($startOrigin !~ /^(1|0)$/);

	$isMol = 1 if (! defined($isMol) or $isMol !~ /0|no/i);
	if ($cStr !~ /(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)(.*)/) {
		die "ERROR: Expected integers for x,y and z cell length. Got \"$cStr\"\n";
	} else {
		$newCell = (
					{
						"X"	     => {hi=>$1,lo=>0,angle=>90,len=>$1},
						"Y"	     => {hi=>$2,lo=>0,angle=>90,len=>$2},
						"Z"	     => {hi=>$3,lo=>0,angle=>90,len=>$3},
						"STRING" => "${1}x${2}x${3}",
					}
					);
		if (defined($4)) {
			$angle_str = $4;
			$angle_str =~ s/^\s+//;
			@{ $tmp } = split /\s+/,$angle_str;
			$newCell->{a} = $tmp->[0] if ($#{ $tmp } > -1);
			$newCell->{b} = $tmp->[1] if ($#{ $tmp } > 0);
			$newCell->{c} = $tmp->[2] if ($#{ $tmp } > 1);
		}
	}
	if (! defined($saveFile)) {
		$saveFile = basename ($bgfFile);
		$saveFile =~ s/\.\w+$//;
		$saveFile .= "_trim_" . $newCell->{STRING} . ".bgf";
	}
	print "Done\n";
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf file -c "x y z (alpha beta gamma)"
		  -a (atom selection) -s (save name) -o (origin) -m (molopt)
options:
	-b bgf file: location of input bgf file
	-c "x y z": new cell parameters. Must be smaller than current cell
	-a atom selection: (optional)
		any valid bgf field expression. E.g. resname eq 'WAT' will select
		all the "WAT" residues while index > 10 will select all indices > 10.
		combine multiple expressions to make complicated selections: e.g.
		(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
		Default: 'all'
	-s save name: (Optional) name of new bgf file
	-m molopt: (Optional) Flag to keep molecules together. Expected 1|yes|0|no. Default 1
	-o origin: 
		1: (0, 0, 0) 
		2: box min 
		3: box max 
		0 (default): box center
		or specify an atom selection in quotes
DATA
}

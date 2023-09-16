#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF);
use Packages::General qw(GetSelections FileTester CombineMols);
use Packages::ManipAtoms qw(GetAtmList);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub showUsage;
sub init;
sub addGroupToStructure;
sub numerically { ($a<=>$b); }
sub originCenter;
sub getOffset;
sub getAltPos;

my ($ATOMS, $BONDS, $HEADERS, $FILES, $SELECT, $randomize);
my ($SELECTATMS, $selection, $DIM, $CENTER, $isAlternate);

$|++;
&init;
print "Parsing Structure file $FILES->{STRUCTURE}...";
($ATOMS->{STRUCTURE}, $BONDS->{STRUCTURE}, $HEADERS) = GetBGFFileInfo($FILES->{STRUCTURE},1);
print "Done\nParsing group file $FILES->{GROUP}...";
($ATOMS->{GROUP}, $BONDS->{GROUP}) = GetBGFFileInfo($FILES->{GROUP}, 0);
originCenter($ATOMS->{GROUP});
print "Done\nParsing atom/residue selection...";
$SELECT = GetSelections($selection, 0);
$SELECTATMS = GetAtmList($SELECT, $ATOMS->{STRUCTURE});
print "Done\nAdding Groups to structure...";
addGroupToStructure($ATOMS, $BONDS, $SELECTATMS, $DIM);
addHeader($ATOMS->{STRUCTURE}, $HEADERS);
print "Done\nCreating $FILES->{SAVE}...";
createBGF($ATOMS->{STRUCTURE}, $BONDS->{STRUCTURE}, $FILES->{SAVE});
print "Done\n";

sub originCenter {
    my ($atoms) = $_[0];
    my (@tmp) = sort numerically keys %{ $atoms };
    my ($FIRST, $i);

    %{ $FIRST } = %{ $atoms->{$tmp[0]} };
    
    for $i (keys %{ $atoms }) {
	for ("XCOORD", "YCOORD", "ZCOORD") {
	    $atoms->{$i}{$_} -= $FIRST->{$_};
	}
    }
}

sub addGroupToStructure {
    my ($atoms, $bonds, $atomSelection, $dimension) = @_;
    my ($i, @tmp, $atmCounter, $OFFSET, $currGroup, $j, $alternate, $k);

    @tmp = sort numerically keys %{ $atoms->{STRUCTURE} };
    $atmCounter = pop(@tmp) + 1;
    @tmp = sort numerically keys %{ $atoms->{GROUP} };
    for $i (sort numerically keys %{ $atomSelection }) {
	$alternate = getAltPos($atoms, $bonds, $i);
	$atoms->{STRUCTURE}{$i}{ALTPOS} = $alternate;
	$OFFSET = (); 
	$currGroup = ();
	for ("XCOORD", "YCOORD", "ZCOORD") {
	    $k = getOffset($dimension->{$_});
	    $k *= $alternate if ($isAlternate);
	    $OFFSET->{$_} = $atoms->{STRUCTURE}{$i}{$_} + $k;
	}
	for $j (@tmp) {
	    %{ $currGroup->{ATOMS}{$j} } = %{ $atoms->{GROUP}{$j} };
	    for ("RESNAME") {
		$currGroup->{ATOMS}{$j}{$_} = $atoms->{STRUCTURE}{$i}{$_};
	    }
	    if (defined($bonds->{GROUP}{$j})) {
		@{ $currGroup->{BONDS}{$j} } = @{ $bonds->{GROUP}{$j} };
	    } else {
		$currGroup->{BONDS}{$j} = ();
	    }
	    for ("XCOORD", "YCOORD", "ZCOORD") {
		$currGroup->{ATOMS}{$j}{$_} += $OFFSET->{$_};
	    }
	}
	($atoms->{STRUCTURE}, $bonds->{STRUCTURE}) = 
	    CombineMols($atoms->{STRUCTURE}, $currGroup->{ATOMS}, 
			$bonds->{STRUCTURE}, $currGroup->{BONDS});
	push @{ $bonds->{STRUCTURE}{$i} }, $atmCounter;
	push @{ $bonds->{STRUCTURE}{$atmCounter} }, $i;
	$atmCounter += scalar(@tmp);
    }
}

sub getOffset {
    my ($offset) = $_[0];
    my ($sign);

    $offset = 0 if (! defined($offset));
    if ($randomize) {
	$offset = 1 if (! $offset);
	$sign = abs($offset)/$offset;
	$offset = $sign * rand($offset);
    }

    return $offset;
}

sub getAltPos {
    my ($atoms, $bonds, $currAtom) = @_;
    my ($alternate, $i, %tmp);

    $alternate = 1;
    for $i (@{ $bonds->{STRUCTURE}{$currAtom} }) {
	$tmp{$i} = 1;
    }

    for $i (sort numerically keys %tmp) {
	if (exists($atoms->{STRUCTURE}{$i}{ALTPOS})) {
	    $alternate = -1 * $atoms->{STRUCTURE}{$i}{ALTPOS};
	    last;
	}
    }

    return $alternate;
}

sub init {
    my (%OPTS, $atomSelect, $usage, $dim);

    getopt('sgadwrf',\%OPTS);
    
    for ("s", "g", "a") {
	if (! exists($OPTS{$_})) {
	    $usage = &showUsage;
	    die "$usage\n";
	}
    }

    ($FILES->{STRUCTURE}, $FILES->{GROUP}, $atomSelect, $dim, $FILES->{SAVE}, $randomize, $isAlternate) = 
	($OPTS{s}, $OPTS{g}, $OPTS{a}, $OPTS{d}, $OPTS{w}, $OPTS{r}, $OPTS{f});

    print "Initializing...";
    for ("STRUCTURE", "GROUP") {
	FileTester($FILES->{$_});
    }

    if ($atomSelect =~ /\s+/) {
        @{ $selection } = split /\s+/, $atomSelect;
    } else {
        $selection->[0] = $atomSelect;
    }
    if (! defined($FILES->{SAVE})) {
	$FILES->{SAVE} = basename($FILES->{STRUCTURE});
	$FILES->{SAVE} =~ s/\.\w+$/_mod\.bgf/;
    }
    if (! defined($dim) or $dim !~ /(\w+):(\-?\d+)/) {
	$dim = "ZCOORD:0.5";
    }

    while ($dim =~ /(X|Y|Z|XCOORD|YCOORD|ZCOORD):(\-?\d+\.?\d*)/gi) {
	$DIM->{uc(substr($1, 0, 1)) . "COORD"} = $2;
    }
    die "ERROR: Expected (x|y|z):val for dimension. Got $dim\n"
	if (! defined($DIM));
    $randomize = 0 if (! defined($randomize));

    print "Done\n";
}

sub showUsage {
    return "usage: $0 -s structure file -g group file\n" . 
	"\t\t\t -a atom selection -d (dimesion) -w (save name) -r (randomize)\n" .
	"Options:\n\t-s structure file: location of bgf file\n" .
	"\t-g group file: location of group file to add to atoms in structure\n" .
	"\t-a (atom selection): [:][I|N|T][a|r]x[-y:z]\n" .
	"\t-d dimesion: Optional. Dimension to make attachment. Expected dim:offset\n" . 
	"\t\tDefault is +0.5 in z. Enclose multiple dimension in quotes\n" .
	"\t-w (save name): Optional. Save name of modified file. Will be structure_mod.bgf as default\n" .
	"\t-r (randomize): Optional. Specify whether to randomize the placement of group. Default is 0\n" .
	"\t-f (flip alternating atoms): Optional. Specifies whether to flip alternative (even/odd) atoms. Default 0\n";
}

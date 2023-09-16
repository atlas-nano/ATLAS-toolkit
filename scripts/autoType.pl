#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(ReadFFs FileTester LoadFFs LoadElements AddElementField AddRingField);
use FileFormats qw(addHeader createBGF ParseStructFile);
use ManipAtoms qw(SelectAtoms BuildAtomSelectionString CreateBondsByDistance);
use CERIUS2 qw(GenerateUFFParms);
use Getopt::Std qw(getopt);
use File::Basename qw(basename fileparse);
use BOX qw(GetBox);
use Storable qw(dclone);

sub numerically { ($a<=>$b) }

my ($inputFile, $ff, $selection, $saveFile, $itype, $typeBonds);
my ($ATOMS, $BONDS, $HEADER, $SELECT, $FFILES, $FF, $ffList, $ELEMENTS, $BOX, $CRYSTX);

$|++;
&init;
print "Getting atom data from $itype file $inputFile...";
($ATOMS, $BONDS, $HEADER) = ParseStructFile($inputFile, 1);
&GetMOL2Header($ATOMS, \@{ $HEADER }) if (exists($ATOMS->{HEADER}));
$SELECT = SelectAtoms($selection, $ATOMS);
die "ERROR: No valid atoms found!\n" if (! keys %{ $SELECT });
print "Done\n";
($FFILES, undef) = ReadFFs($ffList);
$FF = LoadFFs($FFILES);
print "Typing Atoms...";
&AddElementField($ATOMS, $ELEMENTS);
&CreateBondsByDistance($ATOMS,\%{ $BONDS },undef,$ATOMS,2,4) if ($itype =~ /(XYZ|CIF)/);
&setHybridzation($ATOMS, $BONDS);
#&AddRingField($ATOMS, $BONDS, $SELECT);
&typeAtoms($ATOMS, $BONDS, $FF->{TYPING}, $SELECT);
if($typeBonds) {
	print "Done\n...Typing Bonds...";
	&typeBonds($ATOMS, $BONDS, $FF->{TYPING}, $SELECT);
	$BOX = GetBox($ATOMS, undef, $HEADER);
	if (defined($BOX)) {
		if (exists($FF->{PARMS}{UFF})) {
			&updateUsedAtomTypes($FF, $ATOMS);
			&GenerateUFFParms($ATOMS, $BONDS, $FF) if (exists($FF->{PARMS}{UFF}));
		}
		&CreateBondsByDistance($ATOMS, \%{ $BONDS }, $BOX, $ATOMS, $FF);
	}
}
print "Done\nCreating ${saveFile}....";
&addHeader($ATOMS,$HEADER);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub updateUsedAtomTypes {
	my ($ff, $atoms) = @_;
	my ($i);

	for $i (values %{ $atoms }) {
		next if (!exists($ff->{ATOMTYPES}{$i->{FFTYPE}}));
		$ff->{ATOMTYPES}{$i->{FFTYPE}}{USED} = 1;
	}
}

sub typeBonds {
	my ($atoms, $bonds, $parms, $select) = @_;
	my ($i, $j, , $jI, $iI, $c);

	for $i (keys %{ $select }) {
		delete $atoms->{$i}{ORDER};
		$atoms->{$i}{DIFF} = $atoms->{$i}{ELEMENT}{NBONDS}-$atoms->{$i}{NUMBONDS};
	}

	for $i (keys %{ $select }) {
		next if($atoms->{$i}{DIFF} == 0);
		for $jI (0 .. $#{ $bonds->{$i} }) {
			$j = $bonds->{$i}[$jI];
			next if ($atoms->{$j}{DIFF} == 0);
			$atoms->{$i}{DIFF} > $atoms->{$j}{DIFF} ? $c = $atoms->{$j}{DIFF} : $c = $atoms->{$i}{DIFF};
			$atoms->{$i}{DIFF} -= $c;
			$atoms->{$i}{NUMBONDS} += $c;
			$atoms->{$j}{DIFF} -= $c;
			$atoms->{$j}{NUMBONDS} += $c;
			$iI = findBondIndex($bonds->{$j}, $i);
			$atoms->{$i}{ORDER}[$jI] = $c+1;
			$atoms->{$j}{ORDER}[$iI] = $c+1;
			last if ($atoms->{$i}{DIFF} == 0);
		}
	}

	for $i (keys %{ $select }) {
		next if (! exists($atoms->{$i}{ORDER}));
		for $j (0 .. $#{ $bonds->{$i} }) {
			$atoms->{$i}{ORDER}[$j] = 1 if (!defined($atoms->{$i}{ORDER}[$j]));
		}
	}

}

sub typeAtoms {
	my ($atoms, $bonds, $parms, $select) = @_;
	my ($i, $j, $plist, $fftype, $tmp, $maxmatch);

	for $i (keys %{ $select }) {
		$plist = ();
		undef $fftype;
		$maxmatch = -1;
		for $j (keys %{ $parms }) {
			&searchParm($atoms, $i, \@{ $parms->{$j} }, \@{ $plist });
		}
		die "ERROR: Cannot find any atom type for atom # $i ($atoms->{$i}{ATMNAME})\n"
			if ($#{ $plist } == -1);
		for $j (@{ $plist }) {
			next if ($j->{NBCOUNT} <= $maxmatch);
			if($j->{NBCOUNT} == 0) {
				$maxmatch = $j->{NBCOUNT};
				$fftype = $j->{FFTYPE};
				next;
			}
			#check that additional rules are staisfied
			$j->{VISITED} = ();	
			$tmp = checkParm($atoms, $bonds, $i, $j);
			if (defined($tmp) and $tmp ne "") {
				$maxmatch = $j->{NBCOUNT};
				$fftype = $tmp;
			}
		}
		die "ERROR: Cannot find any atom type for atom # $i ($atoms->{$i}{ATMNAME})\n"
			if ($maxmatch == -1);
		$atoms->{$i}{FFTYPE} = $fftype;	
	}
}

sub checkParm {
	my ($atoms, $bonds, $i, $plist) = @_;
	my ($j, $k, $l, $fftype, $tlist1, $tlist2, $nmatch);

	$plist->{VISITED}{$i} = 1;
	$nmatch = 0;
	for $j (@{ $plist->{BOND} }) {
		$tlist1 = ();
		push @{ $tlist1 }, dclone($j);
		undef $fftype;
		#check for match
		for $k (@{ $bonds->{$i} }) {
			next if (exists($plist->{VISITED}{$k}));
			$tlist2 = ();
			&searchParm($atoms, $k, $tlist1,\@{ $tlist2 });
			next if ($#{ $tlist2 } == -1);
			$plist->{VISITED}{$k} = 1;
			for $l (@{ $tlist2 }) {
				$fftype = checkParm($atoms, $bonds, $k, $l);
				last if(defined($fftype)); #match
			}
			last;
		}
		$nmatch++ if(defined($fftype) and $fftype ne "");
	}
	if($nmatch == $plist->{NBCOUNT}) {
		return $plist->{FFTYPE};
	}
}

sub searchParm {
	my ($atoms, $i, $plist, $slist) = @_;
	my ($j, $ele, $hyb);

	$ele = $atoms->{$i}{ELEMENT}{SYMBOL};
	$hyb = $atoms->{$i}{HYBRIDIZATION};
	for $j (@{ $plist }) {
		if ((($j->{ELE} eq "$ele") or ($j->{ELE} eq "any")) and ($#{ $plist } == 0 or (($j->{HYBRID} == 0) or ($j->{HYBRID} == $hyb)))) {
			$j->{NMATCH} = 0;
			push @{ $slist }, $j;
		}
	}
}

sub setHybridzation {
	my ($atoms, $bonds) = @_;
	my ($i, $j, $nb, $lp);

	for $i (keys %{ $atoms } ) {
		$nb = $#{ $bonds->{$i} } + 1; #number of bonds to atom
		$lp = $atoms->{$i}{ELEMENT}{LONEPAIRS};
		$lp = 0 if ($lp < 0);
		$atoms->{$i}{HYBRIDIZATION} = $nb+$lp-1;
		$atoms->{$i}{HYBRIDIZATION} = 3 if($nb > $atoms->{$i}{ELEMENT}{NBONDS}); #hypervalency
		$atoms->{$i}{LONEPAIRS} = $lp;
	}
}

sub findBondIndex {
	my ($bonds, $atom) = @_; 
	my ($i, $index);

	for $i (0 .. $#{ $bonds }) {
		if($bonds->[$i] == $atom) {
			$index = $i; 
			last;
		}
	}   

	return $index;
}

sub init {
    my (%OPTS, $atomSelect);

    getopt('itfabs',\%OPTS);
    for ("i", "f") {
	die "usage: $0 -i input_file -f force field -b (type bonds=no) -t (filetype=bgf|mol2|pdb|xyz) -a (atom selection = all) -s (bgf_file_savename)\n"
	    if (! exists($OPTS{$_}));
    }

    print "Initializing...";
    ($inputFile, $ffList, $atomSelect, $saveFile, $itype, $typeBonds) = ($OPTS{i}, $OPTS{f}, $OPTS{a}, $OPTS{s}, $OPTS{t}, $OPTS{b});
    FileTester($inputFile);

    $atomSelect =  "index > 0" if (! defined($atomSelect));
    $selection = BuildAtomSelectionString($atomSelect);
	$ELEMENTS = LoadElements();
	if (!defined($saveFile)) {
		$saveFile = basename($inputFile);
		$saveFile =~ s/\.\w+$//;
		$saveFile .= ".typed.bgf";
	}
	if (! defined($itype)) {
		(undef,undef,$itype) = fileparse($inputFile,qr"[^.]*$");
	}
	$typeBonds = 0 if (! defined($typeBonds) or $typeBonds =~ /0|no/i);
	$typeBonds = 1 if ($typeBonds =~ /1|yes/i);
    print "Done\n";

}

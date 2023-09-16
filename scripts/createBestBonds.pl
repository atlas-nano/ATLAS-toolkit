#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester IsInteger GetBondLength CoM Rotate);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);
use BOX qw(GetBox);
use Storable qw(dclone);
use Superimpose;
use Math::MatrixReal;

sub usage;
sub createBond;
sub doRotate;
sub numerically { ($a<=>$b); }

my ($bgfFile, $saveName, $atom1, $atom2, $rotate);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $SELECT1, $SELECT2);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&GetMols($ATOMS, $BONDS);
$BOX = GetBox($ATOMS, undef, $HEADERS);
print "Done\nSelecting relevant atoms...";
$SELECT1 = SelectAtoms($atom1, $ATOMS);
$SELECT2 = SelectAtoms($atom2, $ATOMS);
print "Done\nCreating Bonds...";
&createBond($ATOMS, $BONDS, $BOX, $SELECT1, $SELECT2, $rotate);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub createBond {
    my ($atoms, $bonds, $box, $sel1, $sel2, $rotopt) = @_;
	my ($i, $j, $k, $darray, $dist, $used, $tot, $order, $tmp, $tmp1, $counter, $nbond);

	$tot = scalar keys %{ $sel1 };

	$dist = 9999;
	for $i (keys %{ $sel1 }) {
		for $j (keys %{ $sel2 }) {
			$darray->{$i}{$j} = GetBondLength($atoms->{$i}, $atoms->{$j}, $box);
			if($darray->{$i}{$j} < $dist) {
				$dist = $darray->{$i}{$j};
				$tmp = ();
				$tmp->{i} = $i;
				$tmp->{j} = $j;
				$tmp->{dist} = $dist;
			}
		}
	}

	#now get order to progress. start at random and get nearest neighbors
	$i = $tmp->{i}; #start at shortest distance
	push @{ $order }, $i;
	print "shortest_dist: $tmp->{dist}\n";
	$counter = 1;
	while ($counter<$tot) {
		$used->{$i} = 1;
		$tmp = ();
		for $j (keys %{ $sel1 }) {
			next if ($i == $j or exists($used->{$j}));
			$dist = GetBondLength($atoms->{$i}, $atoms->{$j}, $box);
			$tmp->{$dist} = $j;
		}
		@{ $tmp1 } = sort numerically keys %{ $tmp };
		$j = $tmp1->[0];
		push @{ $order }, $tmp->{$j};
		$i = $tmp->{$j};
		$counter++;
	}

	print "order: @{ $order }\n";
	$used = ();
	for $i (@{ $order }) {
		for $j (sort {$darray->{$i}{$a} <=> $darray->{$i}{$b}} (keys %{ $darray->{$i} })) {
			next if (exists($nbond->{$j}));
			$nbond->{$j}{$i} = 1;
			last;
		}
	}

	&doRotate($atoms, $sel1, $sel2, $box, $nbond) if ($rotopt);
		
	for $i (keys %{ $nbond }) {
		for $j (keys %{ $nbond->{$i} }) {
			push @{ $bonds->{$i} }, $j;
			push @{ $bonds->{$j} }, $i;
		}
	}

}
sub doRotate {
	my ($atoms, $sel1, $sel2, $box, $blist) = @_;
	my ($ref, $mov, $i, $j, $rotMat, $CoM, $mol);

	for $i (keys %{ $sel1 }) {
		$ref->{$i} = dclone($atoms->{$i});
		$ref->{$i}{MASS} = 1;
	}
	for $i (keys %{ $sel2 }) {
		for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
			next if exists($mol->{$j});
			$mol->{$j} = \%{ $atoms->{$j} };
			$mol->{$j}{MASS} = 1;
		}
		$mov->{$i} = dclone($atoms->{$i});
	}

	$CoM = CoM($ref);
	&Kearsley::TranslateAllAtoms($ref, $CoM, -1);
	$CoM = CoM($mov);
	&Kearsley::TranslateAllAtoms($mov, $CoM, -1);
	$rotMat = Kearsley::GetRotMatrix($ref, $mov);

	$CoM = CoM($mol);
	&Kearsley::TranslateAllAtoms($mol, $CoM, -1);
	&Kearsley::RotateMolecule($mol, $rotMat);
	&Kearsley::TranslateAllAtoms($mol, $CoM, 1);
}

sub doRotateOld {
	my ($atoms, $sel1, $sel2, $box, $blist) = @_;
	my ($i, $j, $k, $mol, $CoM, $angles, $dim, $rot);
	my ($blen, $clen);

	for $i (keys %{ $sel2 }) {
		$mol->{$i} = dclone($atoms->{$i});
	}

	$CoM = CoM($mol);
	$clen = getBlen($atoms, $mol, $box, $blist);

	@{ $dim } = ("X","Y","Z");
	for $i (0 .. 2) {
		$j = 2;
		while($j<30) {
			@{ $angles } = (0, 0, 0);
			$angles->[0] = $j if ($dim->[$i] eq "X");
			$angles->[1] = $j if ($dim->[$i] eq "Y");
			$angles->[2] = $j if ($dim->[$i] eq "Z");
			&resetXYZ($atoms, $mol);
			&transmol($mol, $CoM, -1);
			&Rotate($mol, $angles, $i);
			&transmol($mol, $CoM, 1);
			$blen = getBlen($atoms, $mol, $box, $blist);
			if($blen < $clen) {
				$rot->{$dim->[$i]} = $j;
				$clen = $blen;
			}
			$j += 2;
		}
	}

	for $i (keys %{ $mol }) {
		for $j (keys %{ $mol->{$i}{MOLECULE}{MEMBERS} }) {
			$mol->{$j} = \%{ $atoms->{$j} };
		}
	}
	$CoM = CoM($mol);
	&transmol($mol, $CoM, -1);
	for $i (keys %{ $rot }) {
		@{ $angles } = (0,0,0);
		if ($i eq "X") {
			$angles->[0] = $rot->{$i};
			$j = 0;
		} elsif ($i eq "Y") {
			$angles->[1] = $rot->{$i};
			$j = 1;
		} else {
			$angles->[2] = $rot->{$i};
			$j = 2;
		}
		&Rotate($mol, $angles, $j);
	}
	&transmol($mol, $CoM, 1);
}

sub resetXYZ {
	my ($atoms, $mol) = @_;
	my ($i, $j);


	for $i (keys %{ $mol }) {
		for $j ("XCOORD","YCOORD","ZCOORD") {
			$mol->{$i}{$j} = $atoms->{$i}{$j};
		}
	}
}

sub getBlen {
	my ($atoms, $mol, $box, $blist) = @_;
	my ($dist, $i, $j);


	for $i (keys %{ $blist }) {
		$atom1 = $mol->{$i};
		for $j (keys %{ $blist->{$i} }) {
			$atom2 = $atoms->{$j};
			$dist += GetBondLength($atom1, $atom2, $box);
		}
	}

	return $dist;
}
		
sub transmol {
	my ($atoms, $com, $phase) = @_;
	my ($i, $j);

	for $j ("XCOORD","YCOORD","ZCOORD") {
		for $i (keys %{ $atoms }) {
			$atoms->{$i}{$j} = $atoms->{$i}{$j} + $phase*$com->{$j};
		}
	}
}

sub init {
    my (%OPTS, $sel1, $sel2);

    getopt('bijsr',\%OPTS);
    ($bgfFile, $sel1, $sel2, $rotate, $saveName) = ($OPTS{b},$OPTS{i},$OPTS{j},$OPTS{r},$OPTS{s});
    for ($bgfFile, $sel1, $sel2) {
        &usage if (! defined($_));
    }
    print "Initializing...";
    FileTester($bgfFile);
	$atom1 = BuildAtomSelectionString($sel1);
	$atom2 = BuildAtomSelectionString($sel2);
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$//;
		$saveName .= ".combined.bgf";
    }
	$rotate = 0 if (!defined($rotate) or $rotate !~ /1|yes/i);
	$rotate = 1 if ($rotate =~ /1|yes/i);
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -i 'atom1_selection' -j 'atom2_selection' -r (rotate_option=no) -s (save_name)
DATA

die "\n";

}

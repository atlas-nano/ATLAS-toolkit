#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms insertHeaderRemark);
use General qw(FileTester GetBondLength ShuffleArray);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);
use BOX qw(GetBox);
use Bulk qw(CalcDipole);
use Storable qw(dclone);

sub usage;
sub checkSystem;
sub getNeighbors;
sub generateIceStructure;
sub createAtom;

my ($bgfFile, $saveName, $selection, $max_iter, $cutoff, $oh_dist, $weight_f);
my ($SELECTIONS, $ATOMS, $BONDS, $tmp, $HEADERS, $BOX);
my ($PAIRS);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
$BOX = GetBox($ATOMS, undef, $HEADERS);
&fixBox($BOX);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
($ATOMS, $BONDS, $tmp) = GetBGFAtoms($SELECTIONS, $ATOMS, $BONDS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $BONDS });
&checkSystem($ATOMS, $BONDS);
print "Done\nGetting oxygen pair information...";
$PAIRS = getNeighbors($ATOMS, $cutoff, $BOX);
print "Done\n";
($ATOMS, $BONDS) = generateIceStructure($ATOMS, $BONDS, $BOX, $PAIRS, $oh_dist, $max_iter, $weight_f);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub generateIceStructure {
	my ($atoms, $bonds, $box, $pairs, $dist, $max, $kappa) = @_;
	my ($pairList, $totP, $i, $j, $totH, $lAtoms, $lBonds, $oxygens,$headers);
	my ($index, $currPair, $parent, $atomH, $violated, $totP_tmp, $icounter);
	my ($mols, $dM, $best, $prob, $count,$iter, $k, $sel, $pair_tmp, $prefix);

	$prefix = $saveName;
	$prefix =~ s/\.\w+$//;
	$pairList = ();
	$totP = 0;
	for $i (keys %{ $pairs }) {
		$oxygens->{$i} = 1;
		for $j (keys %{ $pairs->{$i} }) {
			next if ($j < $i);
			$oxygens->{$j} = 1;
			$pairList->[$totP][0] = $i;
			$pairList->[$totP][1] = $j;
			$totP++;
		}
	}

	$totH = scalar(keys %{ $pairs })* 2; #total number of hydrogens to place

	print "Generating Ice structures\n=========================\n";
	printf "%8s %10s\n","Iteration","Dipole Moment";

	$i = 1;
	$iter = 0;
	while ($i < $max+1) {
		$lAtoms = dclone($atoms);
		undef $lBonds;
		&ShuffleArray($pairList);
		#now place tot hydrogens by randomely selecting from pair array
		$count = scalar keys %{ $atoms };
		$totP_tmp = $totP;
		@{ $pair_tmp } = @{ $pairList };
		for $j (1 .. $totH) {
			$icounter = 0;
			while ($icounter < (1000*$totP)) {
				$icounter++;
				$index = int(rand() * $totP_tmp);
				$currPair = $pair_tmp->[$index];
				#figure out whether to place on atom a or b
				$parent = 0;
				$parent = 1 if (rand() > 0.5);
				if ($lAtoms->{ $currPair->[$parent] }{BONDS}{COUNT} < 2) {
					splice @{ $pair_tmp }, $index, 1;
					$totP_tmp--;
					last;
				}
			}
			last if ($icounter == (1000*$totP));
			last if ($totP_tmp==0);
			#create hydrogen atom and add to oxygen
			$count++;
			$atomH = createAtom($lAtoms, $box, $currPair, $parent, $dist, $count);
			$lAtoms->{$count} = $atomH;
		}
		#check that each oxygen has 2 hydrogens
		$violated = 0;
		$iter++;
		for $j (keys %{ $oxygens }) {
			if ($lAtoms->{$j}{BONDS}{COUNT} != 2) {
				$violated = 1;
				print "Iteration: $iter. Violation: Atom $j #bonds: $lAtoms->{$j}{BONDS}{COUNT} is not 2\r";
				last;
			}
		}
		next if ($violated);
		for $j (keys %{ $lAtoms }) {
			@{ $lBonds->{$j} } = keys %{ $lAtoms->{$j}{BONDS}{LIST} };
		}
		for $i (keys %{ $lAtoms }) {
			$sel->{$i} = 1;
		}
		$mols = GetMols($lAtoms, $lBonds);
		$dM = CalcDipole($lAtoms, $mols, undef, -1);
		#now see if dipolemoment is lower for new structure
		if (! defined($best) or $dM->{u} < $best->{dM}{u}) {
			$best->{dM} = \%{ $dM };
			$best->{atoms} = \%{ $lAtoms };
			$best->{bonds} = \%{ $lBonds };
		} elsif (defined($best)) {
			#accept with some probability
			$prob = rand() * 1/(($dM->{u} - $best->{dM}{u})*$kappa);
			next if ($prob < 0.5);
			$best->{dM} = \%{ $dM };
			$best->{atoms} = \%{ $lAtoms };
			$best->{bonds} = \%{ $lBonds };
		}
		printf "\nIteration: %d Sucess($i/$max). Dipole Moment %10.3f\n",$iter,$dM->{u};
		#write structure
		@{ $headers } = @{ $HEADERS };
		&insertHeaderRemark($headers, "REMARK DIPOLE MOMENT = $dM->{u}");
		&addHeader($lAtoms, $headers);
		&createBGF($lAtoms, $lBonds, "${prefix}.snap${i}.iter${iter}.bgf");
		$i++;
	}
	return($best->{atoms}, $best->{bonds});
}

sub createAtom {
	my ($atoms, $box, $pair, $parent_index, $del, $counter) = @_;
	my ($i, $parent_atom, $new_atom, $vec, $vlen, $atom_i, $atom_j);

	$parent_atom = \%{ $atoms->{ $pair->[$parent_index] } };
	$atom_i = $atoms->{ $pair->[0] };
	$atom_j = $atoms->{ $pair->[1] };
	if ($parent_index) {
		$atom_i = $atoms->{ $pair->[1] };
		$atom_j = $atoms->{ $pair->[0] };
	}
	for $i (keys %{ $parent_atom }) {
		next if ($i =~ /DISP|BOND|MOL/);
		$new_atom->{$i} = $parent_atom->{$i};
	}
	$new_atom->{ATMNAME} = "H" . $counter . "";
	$new_atom->{INDEX} = $counter;
	$new_atom->{FFTYPE} = "H_";
	$new_atom->{CHARGE} = -0.25*($atom_i->{CHARGE} + $atom_j->{CHARGE});
	$new_atom->{BONDS}{COUNT} = 1;
	$new_atom->{BONDS}{LIST}{ $pair->[$parent_index] } = 1;
	$parent_atom->{BONDS}{COUNT}++;
	$parent_atom->{BONDS}{LIST}{$counter} = 1;
	for $i ("XCOORD", "YCOORD", "ZCOORD") {
		$vec->{$i} = $atom_j->{$i} - $atom_i->{$i};
		$vec->{$i} -= $box->{$i}{len} if (abs($vec->{$i} - $box->{$i}{len}) < abs($vec->{$i}));
		$vec->{$i} += $box->{$i}{len} if (abs($vec->{$i} + $box->{$i}{len}) < abs($vec->{$i}));
		$vlen += $vec->{$i}**2;
	}
	$vlen = sqrt($vlen);
	for $i ("XCOORD", "YCOORD", "ZCOORD") {
		$new_atom->{$i} += ($del/$vlen)*$vec->{$i};
	}
	return $new_atom;
}

sub getNeighbors {
	my ($atoms, $max_dist, $box) = @_;
	my ($i, $j, $dist, $pairs);

	for $i (keys %{ $atoms }) {
		for $j (keys %{ $atoms }) {
			next if ($i >= $j);
			$dist = GetBondLength($atoms->{$i}, $atoms->{$j}, $box);
			next if ($dist > $max_dist);
			$pairs->{$i}{$j} = 1;
			$pairs->{$j}{$i} = 1;
		}
		#die "ERROR: No neighbours within $max_dist found for atom $i $atoms->{$i}{ATMNAME}\n"
			#if (! exists($pairs->{$i}) or ! keys %{ $pairs->{$i} });
	}
	return $pairs;
}

sub checkSystem {
	my ($atoms, $bonds) = @_;
	my ($i, $j);

	for $i (keys %{ $atoms }) {
		die "ERROR: Atom $i has bonds!"
			if ($#{ $bonds->{$i} } > -1);
		for $j (keys %{ $atoms->{$i} }) {
			delete $atoms->{$i}{$j} if ($j =~ /MOL/);
		}
	}
}

sub fixBox {
	my ($box) = $_[0];
	my ($i);

	for $i ("X","Y","Z") {
		%{ $box->{"${i}COORD"} }= %{ $box->{$i} };
	}
}

sub init {
	my (%OPTS, $atomSel);
	getopt('bsmocdp',\%OPTS);
	($bgfFile,$saveName,$max_iter,$atomSel,$cutoff,$oh_dist,$weight_f) = 
		($OPTS{b},$OPTS{s},$OPTS{m},$OPTS{o},$OPTS{c},$OPTS{d},$OPTS{p});
	for ($bgfFile, $atomSel) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	FileTester($bgfFile);
	$selection = BuildAtomSelectionString($atomSel);
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/_mod\.bgf/;
	}
	$max_iter = 100 if (!defined($max_iter) or $max_iter !~ /^\d+$/);
	$cutoff = 3.2 if (! defined($cutoff) or $cutoff !~ /^\d+\.?\d*$/);
	$oh_dist = 1.0 if (! defined($oh_dist) or $oh_dist !~ /^\d+\.?\d*$/);
	$weight_f = 10 if (! defined($weight_f) or $weight_f !~ /^\d+\.?\d*$/);
	srand(time() ^($$ + ($$ <<15)));
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -o oxygen_atoms -p (penalty_factor) -c (cutoff) -d (oh_distance) -m (max_num) -s (save_name)
Arguments:
  bgf_file: name of bgf_file
  oxygen_atoms: selection string for oxygen atoms
  penalty_factor: exponential penalty factor for rejecting strucutres = 10
  cutoff: oxygen - oxygen cutoff neighbor distance = 3.2 
  oh_distance: oxygen - hydrogen bond distance = 1.0
  save_name: name of file to save
  max_num: number of structures to generate before stopping
DATA
die "\n";

}

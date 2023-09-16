#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use MolData;
use Getopt::Std qw(getopt);
use General qw(GetSelections ShowSelectionInfo STDev GetAngle);

sub init;
sub showusage;
sub getAtoms;
sub calc;
sub calcAngle;

my ($molStructure, $fileName, $fileType, $verbose, $calcCom);
my ($SELECT1, $SELECT2, $atoms1, $atoms2, $disp, $dims, $atoms3, $SELECT3);

$|++;
&init;
print "Gettting data from $fileType file $fileName..." if ($verbose);
$molStructure->read($fileName, $fileType);
print "Done\nSelecting atoms/residues..." if ($verbose);
$atoms1 = getAtoms($molStructure, $SELECT1);
$atoms2 = getAtoms($molStructure, $SELECT2);
$atoms3 = getAtoms($molStructure, $SELECT3) if (defined($SELECT3));
print "Done\nCalculating COM displacement..." if ($verbose);
$disp = calc($atoms1, $atoms2, $dims, $calcCom, $atoms3);
print "Done\nCOM displacement $disp Angstoms\n" if ($verbose);
print "$disp\n" if (! $verbose);

sub calc {
	my ($mol1, $mol2, $dim, $calcCom, $mol3) = @_;
	my ($center, $atom, $atomMass, $molMass, $i, $dist);
	my ($j, $k, $distStr, $stdev);

	$dist = 0;
	if (! $calcCom) {
		for $i (keys %{ $mol1 }) {
			for $j (keys %{ $mol2 }) {
				next if ($mol1->{$i} == $mol2->{$j});
				if (! defined($mol3)) {
					$dist = calcDistance($mol1->{$i}, $mol2->{$j}, $dim);
					$distStr .= "$dist ";
				} else {
					for $k (keys %{ $mol3 }) {
						next if ($mol3->{$k} == $mol1->{$i} or $mol3->{$k} == $mol2->{$j});
						$dist = calcAngle($mol1->{$i}, $mol2->{$j}, $mol3->{$k}, $dim);
						$distStr .= "$dist ";
					}
				}
			}
		}
		chop $distStr;
		($dist, $stdev, undef) = STDev($distStr);
		#return sprintf("%.2f %.2f", $dist, $stdev);
		return sprintf("%.2f", $dist);
	} else {
		for ($mol1, $mol2) {
			$center = undef;
			$molMass = 0;
			for $atom (values %{ $_ }) {
				$atomMass = $atom->mass;
				$atomMass = 1 if (! defined($atomMass));
				$molMass += $atomMass;
				for $i ("x", "y", "z") {
					$center->{$i} += ($atom->$i * $atomMass);
				}
			}
			for $i ("x", "y", "z") {
				$center->{$i} /= $molMass;
			}
			$_->{CoM} = $center;
		}
		for $i (keys %{ $dim }) {
			$dist += ($mol1->{CoM}{$i} - $mol2->{CoM}{$i})**2;
		}
		return sprintf("%.2f",sqrt($dist));
	}

}

sub calcDistance {
	my ($mol1, $mol2, $dim) = @_;
	my ($dist, $k);
	
	for $k (keys %{ $dim }) {
		$dist += ($mol1->$k - $mol2->$k)**2;
	}
	$dist = sqrt($dist);
	
	return $dist;
}

sub calcAngle {
	my ($mol1, $mol2, $mol3, $dim) = @_;
	my ($dist, $atom1, $atom2, $atom3, $i);

	for ("x", "y", "z") {
		$i = uc($_);
		$atom1->{$i} = 0;
		$atom2->{$i} = 0;
		$atom3->{$i} = 0;
		if (exists($dim->{$_})) {
			$atom1->{$i} = $mol1->$_;
			$atom2->{$i} = $mol2->$_;
			$atom3->{$i} = $mol3->$_;
		}
	}
	$dist = GetAngle($atom1, $atom2, $atom3);
	
	return $dist;

}

sub getAtoms {
	my ($struct, $selection) = @_;
	my ($atomList, $i, $j, $k, $field, $fieldList, $atomId, $atoms);

	for $i (keys %{ $selection }) { # the outer loop is for or
		for $j (keys %{ $selection->{$i} }) { # this is the fields loop - assumed that different fields represent and
			$field = lc($j);
			$field = "resid" if ($field eq "resnum");
			next if (! exists($struct->shash->{$field}));
			for $k (keys %{ $selection->{$i}{$j} }) { # these have to be or
				next if (! exists($struct->shash->{$field}{$k}));
				for $atomId (keys %{ $struct->shash->{$field}{$k} }) {
					$atoms->{$atomId} = $struct->shash->{$field}{$k}{$atomId};
				}
			}
			if (! defined($fieldList)) {
				$fieldList = $atoms;
			} else {
				for $atomId (keys %{ $fieldList }) {
					delete $fieldList->{$atomId} if (! exists($atoms->{$atomId}));
				}
			}
			$atoms = ();
		}
		if (! defined($atomList)) {
			$atomList = $fieldList;
		} else {
			for $atomId (keys %{ $fieldList }) {
				$atomList->{$atomId} = $fieldList->{$atomId};
			}
		}
		undef($fieldList);
	}

	die "ERROR: No atom in file matched selection!\n" if (! $atomList);

	return $atomList;
}

sub init {
	my (%OPTS, $sel1, $sel2, $dimOpt, $sel3);

	getopt('ftmnvdca', \%OPTS);
	for ("f", "m", "n") {
		&showusage if (! exists($OPTS{$_}));
	}

   ($fileName, $fileType, $sel1, $sel2, $verbose, $dimOpt, $calcCom, $sel3) = 
	   ($OPTS{f}, $OPTS{t}, $OPTS{m}, $OPTS{n}, $OPTS{v}, $OPTS{d}, $OPTS{c}, $OPTS{a});
	$verbose = 0 if (! defined($verbose) or $verbose !~ /^\s*(1|yes)\s*$/i);
	$verbose = 1 if ($verbose =~ /^\s*(1|yes)\s*$/i);
	print "Initialzing..." if ($verbose);
	$dims = ({"x" => 1, "y" => 1, "z" => 1}) if (! defined($dimOpt) or $dimOpt !~ /(x|y|z)/i);
	$calcCom = 0 if (! defined($calcCom) or $calcCom !~ /^\s*(1|yes)\s*$/i);
	$calcCom = 1 if ($calcCom =~ /^\s*(1|yes)\s*$/i);
	
	if (defined($dimOpt)) {
		while ($dimOpt =~ /(x|y|z)/ig) {
			$dims->{$1} = 1;
		}
	}
	$molStructure =  MolData->new();
	$molStructure->testFile($fileName);
	if (! defined($fileType)) {
		$fileType = "bgf";
		if ($fileName =~ /\.(\w+)$/) {
			$fileType = lc $1;
		}
	}

	$SELECT1 = getSelOpts($sel1);
	$SELECT2 = getSelOpts($sel2);
	$SELECT3 = getSelOpts($sel3) if (defined($sel3));
	die "ERROR: Atom selections are the same!\n" if ($SELECT1 == $SELECT2);
	print "Done\n" if ($verbose);
}

sub getSelOpts {
	my ($selection) = $_[0];
	my ($SELECT, $group, $count, $tmpStr);

	$count = 1;
	$selection =~ s/^\s*//;
	$selection =~ s/\s*$//;

	if ($selection !~ /.*\(.+\)/) {
		while ($selection =~ /(\S+)/g) {
			$tmpStr .= "($1) ";
		}
		$selection = $tmpStr;
	}
	while ($selection =~/\((.+)\)/g) { #everything that is group in brackets represents a and group
		@{ $group } = split /\s/, $1;
		$SELECT->{$count} = GetSelections($group,0);
		$count++;
	}

	return $SELECT;
}

sub showusage {
	my ($usage) = "usage: $0 -f file name -m atom(s)1 selection -n atom(s)2 selection -a (angle atom(s)) -d (dimension) -c (calc com) -t (file type) -v (verbose = no)\n" .
	"options\n:\t-f filename: location of file\n\t-m/-n atom(s) selection: the atom(s) to be compared (see below for syntax)\n" .
	"\t-t file type: the file format. Will be guessed from the suffix if not selected\n\t-v verbose: print helpful line. default no\n" .
	"\t-d dimension: x|y|z for any combination of the three. Default to xyz\n\t-c calc com: group atoms into com or do per atom measurements. default is no\n\t-a angle atom(s): atom to select to measure atom1-atom2-angle_atom angle. default is none\n" . 
	&ShowSelectionInfo;
	die "$usage";
}

#!/usr/bin/perl -w

use strict;
use FindBin ();
use lib "$FindBin::Bin";
use Crossovers;

sub Initialize;
sub Rotate;
sub getDistance;
sub getCrossovers;
sub Numerically;
sub getMinDist;
sub rotateHelix;
sub writeData;

die "usage: $0 helix1 helix2 specsFile rotateOne\n"
    if (! @ARGV or $#ARGV < 3);

my ($helix1, $helix2, $specsFile, $rotateOne) = @ARGV;
my (%ATOMS, $opts, $cATOMS, $CROSS, %CENTER, $tolerance, $trans0); 
my ($helix, $searchAtom, $h1A, $h2A, $nF, $converged, $angle);
my ($radian) = cRadDegrees(1, 0);

Initialize;
$opts = readSpecsFile($specsFile);
readNamotPDB($helix1, \%ATOMS, 1, $trans0);
readNamotPDB($helix2, \%ATOMS, 2, $opts->{"Trans"});
%{ $CENTER{1} } = %{ $ATOMS{"HELIX1"} };
%{ $CENTER{2} } = %{ $ATOMS{"HELIX2"} };
		     
delete $ATOMS{"HELIX1"};
delete $ATOMS{"HELIX2"};

$ATOMS{1}{2} = reverseRes(\%{ $ATOMS{1}{2} });
$ATOMS{2}{2} = reverseRes(\%{ $ATOMS{2}{2} });

$searchAtom = "P";
($cATOMS, $CROSS, $nF) = getCrossovers(\%ATOMS, \@{ $opts->{"crossovers"} });
$tolerance = 2;
$converged = $angle = $h1A = $h2A = 0;

while (! $converged) {
    if (! $rotateOne) {
	$angle = getMinAngle(1);
	$converged = 1
	    if (($angle < $tolerance) or ((360 - $angle) < $tolerance));
	$h1A += $angle;
	$angle *= $radian;
	rotateHelix(1, $angle);
    }
    last
	if ($converged);

    $angle = getMinAngle(2);
    $converged = 1
	if (($angle < $tolerance) or ((360 - $angle) < $tolerance));
    $h2A += $angle;
    $angle *= $radian;
    rotateHelix(2, $angle);
    last
	if ($rotateOne);
}

$h1A = ($h1A % 360);
$h2A = ($h2A % 360);

print "Rotation Angles: Helix 1-$h1A Helix 2-$h2A...Done\n";

$ATOMS{1}{2} = reverseRes(\%{ $ATOMS{1}{2} });
$ATOMS{2}{2} = reverseRes(\%{ $ATOMS{2}{2} });

writeData(\%{ $ATOMS{1} }, $helix1, $trans0);
writeData(\%{ $ATOMS{2} }, $helix2, $opts->{"Trans"});

sub getMinAngle {
    my ($myHelix) = $_[0];
    my ($i, $currCross, $cPair, $a1, $a2, $atom1);
    my ($chain, $res, $base, $dist, $atom2, $sDist);
    $angle = 0;
    $sDist = 9999;

    for $i (1 .. 360) {
	$dist = 0;
	for $currCross (@{ $CROSS }) {
	    for $cPair (@{ $currCross }) {
		($a1, $a2) = split /\|/, $cPair;
		($chain, $res, $base) = split /,/, $cATOMS->{$a1};
		%{ $atom1 } = %{ $ATOMS{1}{$chain}{$res}{$base} };
		
		($chain, $res, $base) = split /,/, $cATOMS->{$a2};		
		%{ $atom2 } = %{ $ATOMS{2}{$chain}{$res}{$base} };
		
		if ($myHelix == 1) {
		    Rotate(\%{ $atom1 }, ($i * $radian), \%{ $CENTER{1} });
		} else {
		    Rotate(\%{ $atom2 }, ($i * $radian), \%{ $CENTER{2} });
		}
		
		$dist += getDistance($atom1, $atom2);
	    }
	}

	if ($dist < $sDist) {
	    $angle = $i;
	    $sDist = $dist;
	}

    }

    print "Helix $myHelix: Distance=" . ($sDist/$nF) . " Angle=$angle..";
    return ($sDist, $angle);
}

sub getCrossovers {
    my ($atomList, $crossList) = @_;
    my (%crossAtoms, @Crossovers, $counter, $c1, $c2);
    my ($atom1, $atom2, $res, $k1, $k2, $k3, $k4, $i);

    $i = 0;
    for $counter (@{ $crossList }) {
	if ($counter->{"BASE"} > 0 and $counter->{"BASE"} < $opts->{"Helix1"}{"TotalBases"} and 
	    $counter->{"BASE"} < $opts->{"Helix2"}{"TotalBases"}) {
	    $i++;
	    $atom1 = $atom2 = 0;
	    $c1 = getChain($counter->{"BASE"}, 1, $opts);
	    $c2 = getChain($counter->{"BASE"}, 2, $opts);
	    $res = $counter->{"BASE"};

	    $atom1 = findAtom(\%{ $atomList->{1}{$c1}{$res} }, $searchAtom);
	    $atom2 = findAtom(\%{ $atomList->{2}{$c2}{$res} }, $searchAtom);
	    next
		if (! $atom1 or ! $atom2);
	    $k1 = "100" . (10 * $res) . $atom1;
	    $k3 = "100" . (10 * ($res + 1)) . $atom1;
	    $k2 = "200" . (10 * $res) . $atom2;
	    $k4 = "200" . (10 * ($res + 1)) . $atom2;

	    print "Calculating distance between " . $ATOMS{1}{$c1}{$res}{$atom1}{"ATOMNUM"} .
		" and " . $ATOMS{2}{$c2}{$res}{$atom2}{"ATOMNUM"} . "\n";
	    $crossAtoms{$k1} = "$c1,$res,$atom1";
	    $crossAtoms{$k2} = "$c2,$res,$atom2";
	    $crossAtoms{$k3} = $c1 . "," . ($res + 1) . ",$atom1";
	    $crossAtoms{$k4} = $c2 . "," . ($res + 1) . ",$atom2";
	    if ($counter->{"TYPE"} == 1) {
		$Crossovers[($#Crossovers + 1)][0] = $k1 . "|" . $k4;
		$Crossovers[$#Crossovers][1] = $k2 . "|" . $k3;
	    } else {
		$Crossovers[($#Crossovers + 1)][0] = $k1 . "|" . $k2;
		$Crossovers[$#Crossovers][1] = $k3 . "|" . $k4;
	    }
	}
    }
    
    die "ERROR: No crossovers information recorded\n"
	if (! @Crossovers or ! %crossAtoms);

    return (\%crossAtoms, \@Crossovers, ($i * 2));
	    
}
	
sub rotateHelix {
    my ($myHelix, $angle) = @_;
    my ($atoms, $res, $chain, $cAtom);

    for $chain (keys %{ $ATOMS{$myHelix} }) {
	for $res (keys %{ $ATOMS{$myHelix}{$chain} }) {
	    for $atoms (keys %{ $ATOMS{$myHelix}{$chain}{$res} }) {
		$cAtom = \%{ $ATOMS{$myHelix}{$chain}{$res}{$atoms} };
		Rotate(\%{ $cAtom }, $angle, \%{ $CENTER{$myHelix} });
		print "";
	    }
	}
    }
}
    
sub Rotate {
    my ($atom, $angle, $CENTER) = @_;
    my ($returnVal, $counter, $i, $j);
    my (@rZ) = (
		[cos($angle),sin($angle),0],
		[-sin($angle),cos($angle),0],
		[0,0,1],
		);
    
# get internal coordinate molecule
    for ("XCOORD", "YCOORD", "ZCOORD") {
	$atom->{$_} -= $CENTER->{$_};
    }

# apply rotation
    $counter = 0;
    for ("XCOORD", "YCOORD", "ZCOORD") {
	$returnVal->{$_} = $atom->{"XCOORD"} * $rZ[$counter][0] + 
	    $atom->{"YCOORD"} * $rZ[$counter][1] +
	    $atom->{"ZCOORD"} * $rZ[$counter][2];
	$counter++;
    }

# map back to real space
    for ("XCOORD", "YCOORD", "ZCOORD") {
	$atom->{$_} = $returnVal->{$_} + $CENTER->{$_};
    }

}


sub Initialize {

    for ($helix1, $helix2, $specsFile) {
	die "ERROR: Cannot access $_: $!\n"
	    if (! -e $_ or ! -r $_ or ! -T $_);
    }

    $trans0 = (
	       {
		   "XCOORD" => 0,
		   "YCOORD" => 0,
		   "ZCOORD" => 0,
	       }
	       );
}

#!/usr/bin/perl -w
use strict;
sub readNamotPDB;
sub readSpecsFile;
sub Initialize;
sub Rotate;
sub getDistance;
sub getCrossovers;
sub findAtom;
sub getChain;
sub Trim;
sub cRadDegrees;
sub ReverseRes;
sub Numerically;

die "usage: $0 helix1 helix2 which_helix specsFile\n"
    if (! @ARGV or $#ARGV < 3);

my ($helix1, $helix2, $myHelix, $specsFile) = @ARGV;
my (%ATOMS, $opts, $cATOMS, $CROSS, $i, $atom1, $atom2, %CENTER, $angle); 
my ($cPair, $dist, $currCross, $helix, $searchAtom, $a1, $a2, $sDist);

Initialize;
$opts = readSpecsFile($specsFile);
readNamotPDB($helix1, \%ATOMS, 1, 0, 0, 0);
readNamotPDB($helix2, \%ATOMS, 2, $opts->{"Trans"}{"x"}, $opts->{"Trans"}{"y"}, $opts->{"Trans"}{"z"});
%{ $CENTER{1} } = %{ $ATOMS{"HELIX1"} };
%{ $CENTER{2} } = %{ $ATOMS{"HELIX2"} };
		     
delete $ATOMS{"HELIX1"};
delete $ATOMS{"HELIX2"};

$ATOMS{1}{2} = ReverseRes(\%{ $ATOMS{1}{2} });
$ATOMS{2}{2} = ReverseRes(\%{ $ATOMS{2}{2} });

$searchAtom = "P";
($cATOMS, $CROSS) = getCrossovers(\%ATOMS, \@{ $opts->{"crossovers"} });
my ($radian) = cRadDegrees(1, 0);

$sDist = 99999;

for $i (1 .. 360) {
    $dist = 0;
    for $currCross (@{ $CROSS }) {
	for $cPair (@{ $currCross }) {
	    ($a1, $a2) = split /\s+/, $cPair;
	    %{ $atom1 } = %{ $cATOMS->{$a1} };
	    %{ $atom2 } = %{ $cATOMS->{$a2} };
	    
	    if ($myHelix == 1) {
		Rotate(\%{ $atom1 }, ($i * $radian), \%{ $CENTER{1} });
	    } else {
		Rotate(\%{ $atom2 }, ($i * $radian), \%{ $CENTER{2} });
	    }

	    $dist = getDistance($atom1, $atom2);
	    print "REQUESTED:Distance is $dist\n";
	}
    }
#    if ($dist < $sDist) {
#	$angle = $i;
#	$sDist = $dist;
#    }


    print "REQUESTED:TER\n";
}

$sDist /= 4;

#print "Smallest distance: $sDist angle: $angle\n";

sub getCrossovers {
    my ($atomList, $crossList) = @_;
    my (%crossAtoms, @Crossovers, $counter, $c1, $c2);
    my ($atom1, $atom2, $res, $k1, $k2, $k3, $k4);

    for $counter (@{ $crossList }) {
	if ($counter->{"BASE"} > 0 and $counter->{"BASE"} < $opts->{"Helix1"}{"TotalBases"} and 
	    $counter->{"BASE"} < $opts->{"Helix2"}{"TotalBases"}) {
	    $atom1 = $atom2 = 0;
	    $c1 = getChain($counter->{"BASE"}, 1);
	    $c2 = getChain($counter->{"BASE"}, 2);
	    $res = $counter->{"BASE"};

	    $atom1 = findAtom(\%{ $atomList->{1}{$c1}{$res} }, $searchAtom);
	    $atom2 = findAtom(\%{ $atomList->{2}{$c2}{$res} }, $searchAtom);
	    next
		if (! $atom1 or ! $atom2);
	    $k1 = "100" . (10 * $res) . $atom1;
	    $k3 = "100" . (10 * ($res + 1)) . $atom1;
	    $k2 = "200" . (10 * $res) . $atom2;
	    $k4 = "200" . (10 * ($res + 1)) . $atom2;

	    %{ $crossAtoms{$k1} } = %{ $atomList->{1}{$c1}{$res}{$atom1} };
	    %{ $crossAtoms{$k2} } = %{ $atomList->{2}{$c2}{$res}{$atom2} };
	    %{ $crossAtoms{$k3} } = %{ $atomList->{1}{$c1}{($res + 1)}{$atom1} };
	    %{ $crossAtoms{$k4} } = %{ $atomList->{2}{$c2}{($res + 1)}{$atom2} };
	    if ($counter->{"TYPE"} == 1) {
		$Crossovers[($#Crossovers + 1)][0] = "$k1 $k4";
		$Crossovers[$#Crossovers][1] = "$k2 $k3";
	    } else {
		$Crossovers[($#Crossovers + 1)][0] = "$k1 $k2";
		$Crossovers[$#Crossovers][1] = "$k3 $k4";
	    }
	}
    }
    
    die "ERROR: No crossovers information recorded\n"
	if (! @Crossovers or ! %crossAtoms);

    return (\%crossAtoms, \@Crossovers);
	    
}
	
sub findAtom {
    my ($res, $atomName) = @_;
    my ($result, $aC, $atom, $foundAtom);

    $result = $foundAtom = 0;

    for $aC (keys %{ $res }) {
	$atom = \%{ $res->{$aC} };
	if (Trim(lc($atom->{"NAME"})) eq lc($atomName)) {
	    $foundAtom = $aC;
	    last;
	}
    }

    return $foundAtom;
}

sub getDistance {
    my ($atom1, $atom2) = @_;
    my ($dist);

    $dist = 0;
    for ("XCOORD","YCOORD","ZCOORD") {
	$dist += ($atom1->{$_} - $atom2->{$_})**2;
    }

    $dist = sqrt($dist);

    return $dist;
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

    if ($myHelix =~ /(1|2)/) {
	$myHelix = $1;
    } else {
	die "ERROR: Expected 1 or 2 for which_helix, got $myHelix\n";
    }
}

sub readSpecsFile {
    my ($cfile) = $_[0];
    my (%specs, $validCounter, $myKey, $isValid, $line_in, $rec);

    $validCounter = $isValid = 0;
    open CROSSFILE, $cfile or die "Error while reading $cfile: $!\n";
    while (<CROSSFILE>) {
	chomp;
	$line_in = $_;
	if ($line_in =~ /^Helix (\d+)\s+(\d+):(\d+)\s+(\d+)\s+bases/) {
	    $myKey = "Helix" . $1;
	    $specs{$myKey} = (
			      {
				  "MajorGroove"   => $2,
				  "MinorGroove"   => $3,
				  "Periodicity"   => ($2 + $3),
				  "TotalBases"    => $4,
			      }
			      );
	    $validCounter++;
	} elsif ($line_in =~ /^Transalation: (\S+)\s+(\S+)\s+(\S+)/) {
	    $specs{"Trans"} = (
			       {
				   "x" => $1,
				   "y" => $2,
				   "z" => $3,
			       }
			       );
	    $validCounter++;
	} elsif ($line_in =~ /^Parallel: ([1|0])/) {
	    $specs{"isParallel"} = $1;
	    $validCounter++;
	} elsif ($line_in =~ /^5PrimeToLeft: ([1|0])/) {
	    $specs{"5prime"} = $1;
	    $validCounter++;
	}elsif ($line_in =~ /CROSSOVERS/) {
            $isValid = 1; # file must have word "CROSSOVERS"
	    $validCounter++;
        }elsif ($isValid and $line_in =~ /^(\d+):(\d+)->(\d+):(\d+) TYPE:([1|2])$/) { # valid crossover spec
	    $rec = (
		    {
			"TYPE" => $5,
			"BASE" => $2,
		    }
		    );
	    push @{ $specs{"crossovers"} }, $rec;
	}
    }
    close CROSSFILE;

    die "ERROR: No crossvers found while reading file $cfile\n"
	if (! $isValid or $validCounter < 5);
    return \%specs;
}

sub readNamotPDB {
    my ($filenm, $AtomData, $helix, $xoff, $yoff, $zoff) = @_;
    my ($inStr, $curr_chain, @tmp, $curr_res, $rec, $counter, $atomCounter);

    $counter = $curr_chain = 1;
    $curr_res = -1;
    $rec = ();
    $atomCounter = 0;
    open NAMOTFILE, $filenm or die "ERROR: Cannot open $filenm: $!\n";
    while (<NAMOTFILE>) {
	chomp;
	$inStr = $_;
	if ($inStr =~ /^ATOM\s+\d+(.{5})(.{4})\s+(\w)\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
	    $atomCounter++;
	    if ($4 != $curr_res and $curr_res > -1) {
		$AtomData->{$helix}{$curr_chain}{$curr_res} = $rec;
		$rec = ();
		$counter = 1;
		$curr_res = $4;
	    } elsif ($curr_res == -1) {
		$curr_res = $4;
	    }
	    
	    $rec->{$counter} = (
				{
				    "NAME"    => $1,
				    "RESNAME" => $2,
				    "CHAIN"   => $3,
				    "RESNUM"  => $4,
				    "XCOORD"  => $5 + $xoff,
				    "YCOORD"  => $6 + $yoff,
				    "ZCOORD"  => $7 + $zoff,
				}
				);
	    for ("XCOORD","YCOORD","ZCOORD") {
		$AtomData->{("HELIX" . $helix)}{$_} += $rec->{$counter}{$_};
	    }

	    $counter++;
	} elsif ($inStr =~ /^TER/) {
	    $AtomData->{$helix}{$curr_chain}{$curr_res} = $rec;
	    $rec = ();
	    $curr_res = -1;
	    $curr_chain++;
	}
    }

    close NAMOTFILE;

    die "ERROR: No atom information found while reading Namot2 PDB file $filenm\n"
	if ($atomCounter == 0);
    for ("XCOORD","YCOORD","ZCOORD") {
	$AtomData->{("HELIX" . $helix)}{$_} /= $atomCounter;
    }
    
    
}

sub getChain {
# GetChain - gets the chain of a molecule
    my ($curr_group, $whichHelix) = @_;
    my ($returnval, $half_turn, $periodicity, $is5prime, $isParallel);

    $half_turn = $opts->{"Helix1"}{"MinorGroove"};
    $periodicity = $opts->{"Helix1"}{"MajorGroove"} + $half_turn;
    $is5prime = $opts->{"5prime"};
    $isParallel = $opts->{"isParallel"};

    $curr_group -= ($half_turn/2);

    while ($curr_group > $periodicity) {
	$curr_group -= $periodicity;
    }

    #print "5Prime " . $is5prime . "\n";

    if ($whichHelix == 1) {
	if ($is5prime and $isParallel) { #Parallel (PX)
	    if ($curr_group <= $half_turn) {
		$returnval = 1;
	    } else {
		$returnval = 2;
	    }
	} elsif (! $is5prime and $isParallel) { # Parallel (anti-PX)
	    if ($curr_group <= $half_turn) {
		$returnval = 2;
	    } else {
		$returnval = 1;
	    }
	} elsif ($is5prime and ! $isParallel) { # AntiParallel (DX)
	    if ($curr_group <= $half_turn) {
		$returnval = 1;
	    } else {
		$returnval = 2;
	    }
	} elsif (! $is5prime and ! $isParallel) { # AntiParallel (anti-DX)
	    if ($curr_group <= $half_turn) {
		$returnval = 2;
	    } else {
		$returnval = 1;
	    }
	}
    } else {
	if ($is5prime and $isParallel) { #Parallel (PX)
	    if ($curr_group <= $half_turn) {
		$returnval = 1;
	    } else {
		$returnval = 2;
	    }
	} elsif (! $is5prime and $isParallel) { # Parallel (anti-PX)
	    if ($curr_group <= $half_turn) {
		$returnval = 2;
	    } else {
		$returnval = 1;
	    }
	} elsif ($is5prime and ! $isParallel) { # AntiParallel (DX)
	    if ($curr_group <= $half_turn) {
		$returnval = 2;
	    } else {
		$returnval = 1;
	    }
	} elsif (! $is5prime and ! $isParallel) { # AntiParallel (anti-DX)
	    if ($curr_group <= $half_turn) {
		$returnval = 1;
	    } else {
		$returnval = 2;
	    }
	}
    }
	    
    return $returnval;
}

sub Trim {
    my ($inStr) = $_[0];
 
    $inStr =~ s/^\s+//;
    $inStr =~ s/\s+$//;

    return $inStr;
}

sub cRadDegrees {
                                                                                                                  
    my ($inital_angle, $convertToDegrees) = @_;

    my ($resulting_angle) = 0.0;                                                                      
    my ($pi) = atan2(1,1) *4;
                                                                                                                 
    if ($convertToDegrees) { $resulting_angle = $inital_angle * 180 / $pi; }
    else { $resulting_angle = $inital_angle * $pi/180; }
                                                                                                                  
    return $resulting_angle;
}

sub ReverseRes {
    my ($chainData) = $_[0];
    my ($resNum, $atom, $counter, @tmp, %tmpHash);

    @tmp = sort Numerically keys %{ $chainData };

    for $counter (0 .. $#tmp) {
	$resNum = ($#tmp - $counter) + 1;
	%{ $tmpHash{($resNum)} } = %{ $chainData->{$tmp[$counter]} };
    }

    return \%tmpHash;
}

sub Numerically {
    ($a<=>$b);
}

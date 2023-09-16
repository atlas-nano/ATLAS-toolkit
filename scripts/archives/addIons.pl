#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo createBGF addHeader GetResAtoms DeleteAtoms AddAtom PrintCoord UpdateBGF GetSystemCharge);
use Packages::General qw(FileTester IsInteger IsDecimal CoM);
use Packages::BOX qw(GetBox GetRadii CreateGrid GetSurface GetNeighbours);
use Packages::CERIUS2 qw(parseCerius2FF LoadFFs);
use Getopt::Std qw(getopt);

sub usage;
sub init;
sub getSurface;
sub addAtomRadii;
sub getIonOpts;
sub determineIonPlacement;
sub placeIons;
sub getIonParms;
sub findWater;
sub numerically;
sub replaceWat;
sub createIon;
sub getAllOverlaps;
sub placeRandom;
sub checkAtoms;

my ($bgfFile, $ionType, $ionConc, $saveName);
my ($GRID, $BOX, $atm_counter, $PARMS, $BBOX, %PLACED, $FFILES);
my ($IONS, $SURFACE, $tmp, $ATOMS, $BONDS, $HEADERS, $molCharge);

my($start) = time();
$|++;
($IONS, $FFILES) = &init;
$PARMS = LoadFFs($FFILES);
print "Getting atom data from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
checkAtoms($ATOMS, $PARMS);
getIonParms($IONS, $PARMS);
print "Done\nGetting System dimensions...";
$BOX = GetBox($ATOMS, $PARMS, $HEADERS);
print "Done\n";

print "Creating grid...\n";
($GRID, $BBOX, $atm_counter) = CreateGrid($ATOMS, 0, $BOX, 2.5, 1);
($SURFACE, $molCharge) = GetSurface($GRID);
print "\n";

getIonOpts($IONS, $BOX, $ionConc, $molCharge);
print "Determining Ion Placement...";
determineIonPlacement($GRID, $SURFACE, $IONS);
print "Done\n";

placeIons($ATOMS, $BONDS, $GRID, \%{ $IONS->{1} });
if (exists($IONS->{2}{"ATMNAME"})) {
    placeIons($ATOMS, $BONDS, $GRID, \%{ $IONS->{2} });
}
($ATOMS, $BONDS) = UpdateBGF($ATOMS, $BONDS);
print "Creating BGF file $saveName...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

$molCharge = GetSystemCharge($ATOMS);
printf "Total System Charge: %.3f\n", $molCharge;

my ($end) = time();

print "Elapsed: " . ($end - $start) . " seconds\n";

sub checkAtoms {
    my ($atoms, $parms) = @_;
    my ($i, $ffType);

    for $i (keys %{ $atoms }) {
        $ffType = $atoms->{$i}{FFTYPE};
        die "ERROR: Force field type $ffType not found in forcefield(s). Aborting\n"
            if (! exists($parms->{ATOMTYPES}{$ffType}));
        $atoms->{$i}{MASS} = $parms->{ATOMTYPES}{$ffType}{MASS};
        $atoms->{$i}{RADII} = $parms->{VDW}{$ffType}{$ffType}{VALS}[1];
    }
}

sub placeIons {
    my ($atoms, $bonds, $grid, $ions) = @_;
    my ($atmList, $currAtom, $solvAtoms, $atomC, @tmp);
    my ($CoM, $currIon, $resId, $overlap, $total, $resname);
    
    $total = $ions->{"total"};
    if ($#{ $ions->{"ATOM"} } <= $#tmp) {
        pop @{ $ions->{"ATOM"} };
    }
    for $currAtom (@{ $ions->{"ATOM"} }) {
        $solvAtoms = ();
        $atmList = GetResAtoms($currAtom, $bonds, ());
        for $atomC (keys %{ $atmList }) {
    	%{ $solvAtoms->{$atomC} } = %{ $atoms->{$atomC} };
        }
        $resId = $atoms->{$currAtom}{"RESNUM"};
        $resname = $atoms->{$currAtom}{"RESNAME"};
        
        $CoM = CoM($solvAtoms);
        $currIon = createIon($currAtom, $ions, $CoM, $resId);
        %{ $overlap } = %{ $atmList }; #getAllOverlaps($grid, $currIon, $bonds, \%{ $PLACED{$resId} });
        DeleteAtoms($overlap, $atoms, $bonds);
        AddAtom($currIon, $currAtom, $atoms, $bonds);
        print "Replacing $resname (residue# $resId) centered at " . PrintCoord($CoM) . " with " . $ions->{"ATMNAME"} . "\n";
    }
}

sub getAllOverlaps {
    my ($GRID, $atom, $bonds, $currCell) = @_;
    my ($x, $y, $z, $CLIST, $oLap, $i, %overlaps);
    my ($CELL, $counter, $dist, $index);
    $x = $currCell->{"XINDEX"};
    $y = $currCell->{"YINDEX"};
    $z = $currCell->{"ZINDEX"};

    $CLIST = GetNeighbours($GRID, $x, $y, $z);
    for $CELL (@{ $CLIST }) {
      for $counter (@{ $CELL->{"WATERS"} }) {
	  $dist = 0;
	  for $index ("XCOORD", "YCOORD", "ZCOORD") {
	      $dist += (($atom->{$index} - $counter->{$index}) ** 2);
	  }
	  $dist = sqrt($dist);
	  if ($dist <= ($atom->{"RADII"} + $counter->{"RADII"})) {
	      print "overlap with residue " . $counter->{"RESNUM"} . "..";
	      $oLap = GetResAtoms($counter->{"INDEX"}, $bonds, ());
	      for $i (keys %{ $oLap }) {
		  $overlaps{$i} = 1;
	      }
	  }
      }
    }
    
    return \%overlaps;

}

sub createIon {
    my ($ionIndex, $ionParms, $ionPos, $ionResID) = @_;
    my (%ION);
    
    %ION = %{ $ionPos };
    $ION{"INDEX"} = $ionIndex;
    $ION{"ATMNAME"} = $ionParms->{"ATMNAME"};
    $ION{"RESNAME"} = $ionParms->{"ATMNAME"};
    $ION{"RESNUM"} = $ionResID;
    $ION{"FFTYPE"} = $ionParms->{"FFTYPE"};
    $ION{"NUMBONDS"} = 0;
    $ION{"LONEPAIRS"} = 0;
    $ION{"CHARGE"} = $ionParms->{"Charge"};
    $ION{"RADII"} = $ionParms->{"RADII"};
    $ION{"LABEL"} = "HETATM";
    return \%ION;
}

sub determineIonPlacement {
    my ($grid, $surface, $ions) = @_;
    my ($counter, %CHARGE, $ionC, $zI, @tmp, $placed);
 
    for $counter (keys %{ $surface }) {
	if (! $surface->{$counter}{"BURRIED"}) {
	    $CHARGE{ $surface->{$counter}{"CHARGE"} } = $counter;
	}
    }

    $ions->{1}{"placed"} = 0;
    $ions->{2}{"placed"} = 0;
    if ($ions->{1}{"Charge"} > 0) {
	@tmp = sort numerically keys %CHARGE;
    } else {
	@tmp = reverse sort numerically keys %CHARGE;
    }
    findWater($grid, \%{ $ions->{1} }, \@tmp, \%CHARGE, $surface);
    $ions->{1}{"placed"} = $#{ $ions->{1}{"ATOM"}  } + 1;
    
    if (exists($ions->{2}{"ATMNAME"})) {
	@tmp = reverse sort numerically keys %CHARGE;
	#findWater($grid, \%{ $ions->{2} }, \@tmp, \%CHARGE, $surface);
        #$ions->{2}{"placed"} = $#{ $ions->{2}{"ATOM"} } + 1;
    } else {
	$ions->{2}{"total"} = 0;
	$ions->{2}{"placed"} = 0;
    }
    
    
    if ($ions->{1}{"total"} > $ions->{1}{"placed"}) {
        placeRandom(\%{ $ions->{1} }, $grid, $ions->{1}{"total"} - $ions->{1}{"placed"});
    }
    
    if ($ions->{2}{"total"} > $ions->{2}{"placed"}) {
        placeRandom(\%{ $ions->{2} }, $grid, $ions->{2}{"total"} - $ions->{2}{"placed"} - 1);
    }        
}

sub placeRandom {
    my ($ion, $grid, $amount) = @_;
    my ($counter, $xI, $yI, $zI);
    my ($xTotal, $yTotal, $zTotal, @tmp);
    
    
    @tmp = keys %{ $grid };
    $xTotal = $#tmp - 3;
    @tmp = keys %{ $grid->{1} };
    $yTotal = $#tmp - 3;
    @tmp = keys %{ $grid->{1}{1} };
    $zTotal = $#tmp - 3;

    $counter = 0;
    for (;;) {
	last unless ($amount > $counter);
        $xI = sprintf("%.0f", rand($xTotal) + 2);
	$yI = sprintf("%.0f", rand($yTotal) + 2);
	$zI = sprintf("%.0f", rand($zTotal) + 2);
	
	next if (exists($grid->{$xI}{$yI}{$zI}{"used"}));
		
	if (replaceWat($grid, $ion, $xI, $yI, $zI)) {
	    $counter++;
	}
	last if ($counter == $amount);
    }
}    
    
sub findWater {
    my ($grid, $ions, $tmp, $CHARGE, $surface) = @_;
    my ($i, $cellC, @tmp, $x, $y, $z, $placed);

    $placed = $ions->{"placed"};
    while ($ions->{"total"} > $placed and $#{ $tmp } > -1) {
	$cellC = pop @{ $tmp };
	$x = $surface->{ $CHARGE->{$cellC} }{"X"};
	$y = $surface->{ $CHARGE->{$cellC} }{"Y"};
	$z = $surface->{ $CHARGE->{$cellC} }{"Z"};

	next
	    if (exists($grid->{$x}{$y}{$z}{"used"}));
	
	$i = replaceWat($grid, $ions, $x, $y, $z);
	if ($i) {
	    $placed++;
	    delete $CHARGE->{$cellC};
	}
    }
}

sub replaceWat {
    my ($grid, $ions, $x, $y, $z) = @_;
    my ($currCell, $CELLS, $i, $watIndex, $resIndex, $isValid);
    
    $isValid = 0;
    $currCell = $grid->{$x}{$y}{$z};
    $currCell->{"used"} = 1;
    $CELLS = GetNeighbours($grid, $currCell);
    for $i (@{ $CELLS }) { # chose the first one with waters
	next
	    if (exists($i->{"used"}));
	if (exists( $i->{"WATERS"} )) {
	    if ($#{ $i->{"WATER"} } > 0) {
		$watIndex = sprintf("%.0f",rand($#{ $i->{"WATERS"} }));
	    } else {
		$watIndex = 0;
	    }
	    $resIndex = $i->{"WATERS"}[$watIndex]{"RESNUM"};
	    next
		if (exists($PLACED{$resIndex}));
	    $PLACED{$resIndex}{"XINDEX"} = $i->{"XINDEX"};
	    $PLACED{$resIndex}{"YINDEX"} = $i->{"YINDEX"};
	    $PLACED{$resIndex}{"ZINDEX"} = $i->{"ZINDEX"};
	    push @{ $ions->{"ATOM"} }, $i->{"WATERS"}[$watIndex]{"INDEX"};
	    $i->{"used"} = 1;
	    $isValid = 1;
	    last;
	}
    }
    
    return $isValid;
}

sub getIonParms {
    my ($ionInfo, $parms) = @_;
    my ($counter, $element, $radii);
    
    for $counter (keys %{ $parms->{"ATOMTYPES"} }) {
	$element = $parms->{"ATOMTYPES"}{$counter}{"ATOM"};
	if ($element eq $ionInfo->{1}{"ATMNAME"}) {
	    $ionInfo->{1}{"FFTYPE"} = $counter;
	    $ionInfo->{1}{"RADII"} = GetRadii(\%{ $ionInfo->{1} }, $parms);
	}
	if (exists($ionInfo->{2}) and $element eq $ionInfo->{2}{"ATMNAME"}) {
	    $ionInfo->{2}{"FFTYPE"} = $counter;
	    $ionInfo->{2}{"RADII"} = GetRadii(\%{ $ionInfo->{2} }, $parms);
	}
    }
    
    if (! exists($ionInfo->{1}{"FFTYPE"})) {
	die "ERROR: The " . $ionInfo->{1}{"ATMNAME"} .
	" ion does not have any parameters in the provided forcefield\n";
    } elsif (exists($ionInfo->{2}) and ! exists($ionInfo->{2}{"FFTYPE"})) {
	die "ERROR: The " . $ionInfo->{2}{"ATMNAME"} .
	" ion does not have any parameters in the provided forcefield\n";
    }
}


sub getIonOpts {
    my ($ionInfo, $BOX, $num_ions, $charge) = @_;
    my ($mole, $vol, $dim, $molar);
    
    $vol = 1;
    $mole = 6.0221415E23;
    for $dim ("X", "Y", "Z") {
	$vol = $vol * ($BOX->{$dim}{"hi"} - $BOX->{$dim}{"lo"});
    }
    printf "CELL VOLUME %.3G A^3", $vol;
    $vol *= 1E-30; #Angstroms^3 to meter^3
    $vol *= 1000; #meter^3 to liter
    printf "(%.3G liter)\n", $vol;
    $ionInfo->{1}{"Charge"} = $ionInfo->{1}{"total"} = 0;
    $ionInfo->{2}{"Charge"} = $ionInfo->{2}{"total"} = 0;

    if (! IsInteger(abs($charge))) {
	print "WARNING: System has net charge!\n";
	$charge = sprintf("%.0f", $charge);
    }

    if ($ionInfo->{1}{"ATMNAME"} =~ /Na|K/) { # +1 charge
	$ionInfo->{1}{"Charge"} = "+1";
    } elsif (lc($ionInfo->{1}{"ATMNAME"} eq "cl")) {
	$ionInfo->{1}{"Charge"} = -1;
    } else {
	$ionInfo->{1}{"Charge"} = "+2";
    }

    if (IsDecimal($num_ions) and $num_ions > 0) { #placing a concentration
	$ionInfo->{1}{"total"} = sprintf("%.0f",($mole * $vol * $num_ions));
	$molar = $num_ions;
    } elsif ($num_ions > 0 and ! IsDecimal($num_ions)) {
	$ionInfo->{1}{"total"} = $num_ions;
	$molar = $ionInfo->{1}{"total"} / ($mole * $vol);
    } elsif ($num_ions == 0 and ! IsDecimal($num_ions)) { # neutralize
	if ($ionInfo->{1}{"ATMNAME"} !~ /Na|K/) {
	    if (($charge % 2) == 1) {
		print "WARNING: Charge is not even but divalent cation chosen! System will have net charge\n";
	    }
	    $ionInfo->{1}{"total"} = abs(int($charge/2));
	} else {
	    $ionInfo->{1}{"total"} = abs($charge);
	}
	$molar = $ionInfo->{1}{"total"} / ($mole * $vol);
    } else {
	die "ERROR: Invalid option for ion/salt type and concentration\n";
    }

    if (exists($ionInfo->{2}{"ATMNAME"})) {
	$ionInfo->{2}{"Charge"} = -1;
	$ionInfo->{2}{"total"} = $ionInfo->{1}{"total"} * $ionInfo->{1}{"Charge"} + $charge + 1;
	if ($ionInfo->{2}{"total"} < 0) {
	    $ionInfo->{2}{"total"} = 0;
	}
    }
    
    if (! exists($ionInfo->{2}{"ATMNAME"})) {
	printf "Placing " . $ionInfo->{1}{"total"} . " molecule(s) of " . $ionInfo->{1}{"ATMNAME"} . 
	    $ionInfo->{1}{"Charge"} . " ion (%.3G molar)\n", $molar;
    } else {
	printf "Placing " . $ionInfo->{1}{"total"} . " molecule(s) of " . $ionInfo->{1}{"ATMNAME"} . 
	    $ionInfo->{2}{"ATMNAME"} . " salt (%.3G molar)\n", $molar;
	print "Placing " . $ionInfo->{2}{"total"} . " " . $ionInfo->{2}{"ATMNAME"} . "\n";
    }
}

sub init {
    my ($ION, %OPTS, $FF, @FFILES, @tmp, $i);

    getopt('bfins',\%OPTS);
    ($bgfFile, $FF, $ionType, $ionConc, $saveName) = 
	($OPTS{b}, $OPTS{f}, $OPTS{i}, $OPTS{n}, $OPTS{s});

    if (! defined($bgfFile) || ! defined($ionType)) {
	&usage;
	die "\n";
    }
    print "Initializing...";
    FileTester($bgfFile);

    die "ERROR: Available ion types are - Na, Cl, Mg, K, Ca " .
	"and the corresponding salts (e.g. KCl)\n"
	if ($ionType !~ /Na|Cl|Mg|K|Ca/i);
    if (! $saveName) {
	$saveName = $bgfFile;
	$saveName =~ s/\.w+$/_mod\.bgf/;
    }

    if (! defined($ionConc)) {
	$ionConc = 0;
    }

    if ($ionType !~ /^(Na|Mg|K|Ca|Cl)(Cl)?$/i) {
	print "ERROR: Invalid ion/salt type $ionType\n";
	usage;
	die "\n";
    } else {
	$ION->{1}{"ATMNAME"} = $1;
	if ($2 and $1 ne $2) {
	    $ION->{2}{"ATMNAME"} = $2;
	}
    }

    if ($ionConc eq "0" and length($ionType) > 2) {
	die "ERROR: Cannot use a salt to neutralize the system\n";
    }

    die "ERROR: ionConc has to be integer of decimal\n"
	if (! IsInteger($ionConc) and ! IsDecimal($ionConc));

   if ($FF =~ /\s+/) {
        @tmp = split /\s+/, $FF;
    } else {
        @tmp = ($FF);
    }

    for $i (@tmp) {
        if (-e $i &&  -r $i && -T $i) {
            push @FFILES, $i;
        }
    }

    die "ERROR: No valid CERIUS2 forcefield found!\n" if (! @FFILES);

    print "Done\n";
    return ($ION, \@FFILES);
}

sub usage {
    print <<USAGE
usage: $0 -b bgf file -f cerius2 force field(s) -i ion type -n [\# ions] -s [save_name]
options:
    bgf_file: the name of the BGF formatted file
    cerius2_ff: 1 (or more if in quotes) cerius2 formatted forcefile for the bgf file
    ion_type: Na, Cl, Mg, K, Ca and the combination salts (e.g. NaCl)
    [num_ions] : number of ions
                 0 - neutralize the system (default) - ions only
                 X - X number of ions/salt
                 X.x - X.x molar ion/salt
    [save_name]: name of output BGF file (optional)
USAGE
}

sub numerically {
    ($a<=>$b);
}

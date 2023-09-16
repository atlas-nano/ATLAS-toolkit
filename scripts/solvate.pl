#!/usr/bin/perl -w
#use diagnostics;

# solvate.pl: This script will take an input bgf_file, determine the surface atoms and add water molecules 
# of the desired radii around it in a random packing arrangement, in order to achieve the required density

use constant PI => 3.14159265358979;
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester GetBondLength Rotate CoM LoadFFs);
use FileFormats qw(GetBGFFileInfo createBGF addHeader addBoxToHeader createHeaders);
use BOX qw(GetBox CreateGrid CenterAtoms GetGridDims GetNeighbours);
use ManipAtoms qw(CenterSystem MapOrigin);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
sub init;
sub storeSurfaceAtm;
sub addSolv;
sub addRes;
sub numerically;
sub getLastRes;
sub getRadii;
sub getRandomNumber;
sub deflateGrid;
sub createWatBox;
sub placeOnGrid;
sub loadFF;
sub usage;
sub checkAtoms;
sub getWatOpts;
sub updateBox;

my ($bgfFile, $waterType, $density, $dimOpt, $saveName);
my ($dimMin, $dimMax, $start, $end, $GRIDOPTS, $loadWat);
my ($ATOMS, $BONDS, $SURFACE, $BBOX, $MOL_VOL, $vol, $WAT, $CENTER);
my ($atm_counter, $water_res, $GRID, $num_cells, $HEADERS, $num_wats);
my ($res_counter, $PARMS, $i, %WATS, $grid_len, $FFILES, $waterOpt);

$start = time();
$|++;
$FFILES = &init;
$PARMS = LoadFFs($FFILES);
getWatOpts if ($loadWat);
print "Parsing BGF File $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
&checkAtoms($ATOMS, $PARMS);
$CENTER = GetBox($ATOMS, undef, undef);
CenterSystem($ATOMS, $CENTER, ());
$MOL_VOL = GetBox($ATOMS, $PARMS, undef);
&MapOrigin($MOL_VOL, CoM($ATOMS));
print "Done\nCreating Grid...\n";
$GRIDOPTS = GetGridDims($MOL_VOL, $dimOpt, $dimMin, $dimMax);
($GRID, $BBOX, $atm_counter) = CreateGrid($ATOMS, $WAT->{"RADII"}, $GRIDOPTS, $grid_len, 1);
$res_counter = getLastRes($ATOMS);
for $i (keys % { $BBOX }) {
    $vol *= ($BBOX->{$i}{"hi"} - $BBOX->{$i}{"lo"});
}

$ATOMS->{"tot"} = $atm_counter;
$ATOMS->{"solute"} = $atm_counter;
$num_wats = int(($density * $vol * 6.023E-1)/18) if (! defined($num_wats));
print "Attempting to place $num_wats solvent molecules...";
$num_wats = createWatBox($ATOMS, $GRID, $BBOX, $grid_len, $WAT);
print "Done\n";
#$ATOMS->{"tot"} = $atm_counter;

#$num_wats = addSolv($ATOMS, $BONDS, $GRID, $WAT, $num_wats);
delete $ATOMS->{"tot"};
delete $ATOMS->{"solute"}; 

print "Total Simulation Box Size: $vol A^3\n";
print "Placed $num_wats solvent molecules\n";
$density = ((18 * $num_wats)/($vol * 6.023E-01));
print "Total solvent density achieved: $density (g/cm^3)\n";
print "Saving files $saveName...";
CenterAtoms($ATOMS, $BBOX, undef);
$HEADERS = createHeaders(updateBox($BBOX), $saveName);
addBoxToHeader($HEADERS, $BBOX);
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";
$end = time();
print "Elapsed: " . ($end - $start) . " seconds\n";

sub getWatOpts {
    my ($i, $CENTER, $dim);
    ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($waterOpt,0);
    checkAtoms($WAT->{"ATOMS"},$PARMS);
    $CENTER = CoM($WAT->{"ATOMS"});
    for $i (keys %{ $WAT->{ATOMS} }) {
        for $dim ("XCOORD", "YCOORD", "ZCOORD") {
            $WAT->{"ATOMS"}{$i}{$dim} -= $CENTER->{$dim};
        }
    }
    $WAT->{"RADII"} = getRadii($WAT->{"ATOMS"}, $WAT->{"BONDS"});
    printf "SOLVENT RADII: %.1f A\n", $WAT->{"RADII"};
    $grid_len = $WAT->{RADII} + 1.0;
}

sub checkAtoms {
    my ($atoms, $parms) = @_;
    my ($i, $ffType);

    for $i (keys %{ $atoms }) {
	$ffType = $atoms->{$i}{FFTYPE};
	die "ERROR: Force field type $ffType not found in forcefield(s). Aborting\n"
	    if (! exists($parms->{ATOMTYPES}{$ffType}));
	$atoms->{$i}{MASS} = $parms->{ATOMTYPES}{$ffType}{MASS};
	$atoms->{$i}{RADII} = $parms->{VDW}{$ffType}{$ffType}{1}{VALS}[1];
	die "ERROR: No VDW parameters found for atom type $ffType. Aborting\n"
	    if (! exists($atoms->{$i}{RADII}));
    }
}

sub usage {
    my ($printStr);
    
    $printStr = <<DATA;
usage: $0 -b bgfFile -f "forcefield1 forcefield2..." -w water type -d required density(g/cm3) -c dimension -h dimHi -l dimLo -s [saveName] -n [num wats]
OPTIONS (all are requied except saveName):
    -b bgfFile: Name of bgfFile, typed according to the forcefield (below)
    -f forcefield: 1 (or more enclose in quotes) Cerius2/MPSim force fields
    -w solvent type: Current options are
                   TIP3: the original Jorgenson TIP3 water model (rigid hydrogens)
                   TIP3_CHARMM: TIP3 water model as implemented in CHARMM
                   F3C: F3C water model (no rigid hydrogens)
		   F3C_OPT: F3C water model with optimized charges and vib. spectrum
                   SPC: SPC water model
		   MESO: Water model for Meso-scale simulations (Valeria Molinero)
                       -or- you can specify your own solvent: give the path
    -d density: the density of solvent required in g/cm3
    -c dimension: the dimension in which to extend the solvent box by dimLo, dimHi (see below). Options are:
                   xyz: extend in all three dimensions
		   xy|xz|yz: extend in two specified dimensions (other dimension is fixed)
                   x|y|z: extend in only one dimension (other two dimensions are fixed)
    -h dimHi: extend this far past the rectanglar box bounding the bgfFile. Options are:
		   -h integer: make the total bound box integer Angstroms. 
		               NOTE: adjusted to be at least the solute bounding box.
                   -h +integer: pad the box with integer Angstroms worth of solvent molecules of specified density
    -l dimLo: same as above, but in the negative dimension direction
    -s [saveName]: Optional. Name of final system. If not specified will be {bgfFile}_solv.bgf
    -n [num wats]: Optional. Specifies the maximum number of water molecules to add. Overrides the density number.
DATA
    return $printStr;
}

sub init {
    my (%OPT, $myPath, $FF, @tmp, $i);
    my (@FFILES);

    $myPath = "/ul/tpascal/scripts/dat/WAT/";
    getopt('bwdchlsfn',\%OPT);
    ($bgfFile, $FF, $waterOpt, $density, $dimOpt, $dimMin, $dimMax, $saveName, $num_wats) = 
	($OPT{b},$OPT{f},$OPT{w},$OPT{d},$OPT{c},$OPT{l},$OPT{h},$OPT{s},$OPT{n});

    for $i ($bgfFile, $FF, $waterOpt, $density, $dimOpt, $dimMin, $dimMax) {
	die &usage . "\n" if (! defined($i));
    }
    
    print "Initializing...";
    FileTester($bgfFile);    
    die "ERROR: Expected decimal/interger for water_density, got $density\n"
	if ($density !~ /^\d+[\.\d+]?/);
    for $i ($density, $dimMin, $dimMax) {
	die "ERROR: Expected decimal/interger value. Got \"$i\"!\n" if ($i !~ /^[\+|\-]?\d+[\.\d+]?/);
    }
    ($dimMin, $dimMax) = ($dimMax, $dimMin) if ($dimMin > $dimMax);
    if (! $saveName) {
	$saveName = basename($bgfFile);
	$saveName =~ s/\.w+$//;
	$saveName .= "_solv\.bgf";
    }

    $loadWat = 0;
    if (-e $waterOpt && -r $waterOpt && -T $waterOpt) {
	$loadWat = 1;
    } else {
	$waterOpt = uc($waterOpt);
	if($waterOpt eq "TIP3") {
	    ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($myPath . "tip3.bgf",0);
	    $WAT->{"RADII"} = 1.8;
	} elsif ($waterOpt eq "F3C") {
	    ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($myPath . "F3C.bgf",0);
	    $WAT->{"RADII"} = 1.8;
        } elsif ($waterOpt eq "F3C_OPT") {
            ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($myPath . "F3C_OPT.bgf",0);
            $WAT->{"RADII"} = 1.8;
	} elsif ($waterOpt eq "MESO") {
	    ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($myPath . "Meso.bgf",0);
	    $WAT->{"RADII"} = 1.6;
	} elsif ($waterOpt eq "TIP3_CHARMM") {
	    ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($myPath . "tip3_charmm.bgf",0);
	    $WAT->{"RADII"} = 1.8;
	}
    }

    die "ERROR: Invalid selection for solventType: $waterType\n" if (! keys %{ $WAT } && ! $loadWat);
	    
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
    $atm_counter = $water_res = $res_counter = 0;
    $vol = 1;
    $WAT->{RADII} = 1 if (! exists($WAT->{RADII}));
    $grid_len = $WAT->{RADII} + 1;
    if (defined($num_wats)) {
	undef($num_wats) if ($num_wats !~ /^\d+/);
    }
    print "Done\n";
    return \@FFILES;
}

sub addSolv {
    my ($atoms, $bonds, $grid, $myWat, $requiredWat) = @_;
    my (@tmp, $i, $wC, $water, $addWat);

    @tmp = keys %WATS;
    $addWat = $#tmp + 1;

    if ($addWat > $requiredWat) {
	print "randomely removing " . ($addWat - $requiredWat) . " waters...";
	while ($addWat > $requiredWat) {
	    $i = int(rand($#tmp));
	    if (exists($WATS{$tmp[$i]})) {
		delete $WATS{$tmp[$i]};
		$addWat--;
	    }
	}
    }

    $addWat = 0;
    for $wC (keys %WATS) {
	$water = \%{ $WATS{$wC} };
        addRes($water, $atoms, $myWat);
	$addWat++;
    }

    return $addWat;
}

sub addRes {
    my ($watPos, $All_Atoms, $watData, $grid, $cellIndex) = @_;
    my ($counter, $currAtm, $atom, $atomC, $bonds, $Angles, %WAT, $offset, $currCell);

    $water_res++;
    $res_counter++;
    $offset = $atomC = $All_Atoms->{"tot"};
    
    # do rotation
    $Angles = getRandomNumber(360);
    Rotate(\%{ $watData->{"ATOMS"} }, $Angles, 3);
    
    $currCell = $grid->{$cellIndex->{X}}{$cellIndex->{Y}}{$cellIndex->{Z}};

    for $atom (sort numerically keys %{ $watData->{"ATOMS"} }) {
	$atomC++;
	%{ $currAtm } = %{ $watData->{"ATOMS"}{$atom} };
	for $counter ("XCOORD", "YCOORD", "ZCOORD") {
	    $currAtm->{$counter} += $watPos->{$counter};
	}
	$currAtm->{"RESNUM"} = $res_counter;
	$currAtm->{"LABEL"} = "HETATM";
	$currAtm->{"CHAIN"} = "X";
	for $counter (keys %{ $currAtm }) {
	    $All_Atoms->{$atomC}{$counter} = $currAtm->{$counter};
	}
	if (exists($watData->{"BONDS"}{$atom})) {
	    @{ $bonds } = @{ $watData->{"BONDS"}{$atom} };
	    
	    for (@{ $bonds }) {
		push @{ $BONDS->{$atomC} }, ($_ + $offset);
	    }
	}
	for ("X", "Y", "Z") {
	    $currAtm->{CELL}{$_ . "INDEX"} = $cellIndex->{$_};
	}
	push @{ $currCell->{WATERS} }, $currAtm;
    }

    $All_Atoms->{"tot"} = $atomC;
}

sub numerically {
    ($a <=> $b);
}

sub getLastRes {
    my ($AtomData) = @_;
    my (@atoms) = sort numerically keys %{ $AtomData };

    return $AtomData->{$atoms[$#atoms]}{"RESNUM"};

}

sub getRadii {
    my ($ATOMS, $BONDS) = @_;
    my ($radii, $i, $dim, $aRad);

    $radii = 0;
    for $i (keys %{ $ATOMS }) {
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    $aRad = $ATOMS->{$i}{$dim} + $ATOMS->{$i}{RADII}/2;
	    $radii = $aRad if ($aRad > $radii);
	}
    }
    return $radii;
}

sub getRandomNumber(@) {
    my ($max) = $_[0];
    my (@output_array, $rand_angle, $counter);

    for $counter (0 .. 2) {
        $rand_angle = (rand $max) * PI/180;
        $output_array[$counter] = sprintf("%.2f", $rand_angle);
    }

    return \@output_array;

}

sub createWatBox {
    my ($Atoms, $GRID, $BBOX, $gridLen, $myWat) = @_;
    my ($i, $j, $k, $xIndex, $yIndex, $zIndex, $currCell, $dim, @tmp);
    my ($totWat, $radii, $variance, $WAT_POS, $BOX, $iterations);
    
    $radii = sqrt(2) * $myWat->{"RADII"};
    # get 1/2 of radii for deviation
    $variance = 0.5 * $radii;

    #deflate the box by the total variance
    for $dim ("X", "Y", "Z") {
        $BOX->{$dim}{"lo"} = $BBOX->{$dim}{"lo"} + $variance;
        $BOX->{$dim}{"hi"} = $BBOX->{$dim}{"hi"} - $variance;
    }
    
    $xIndex = $yIndex = $zIndex = 1;
    #$Atoms->{"tot"}++;
    $totWat = $iterations = 0;
    while ($totWat < $num_wats) {
        $i = $BOX->{"X"}{"lo"};
    while ($i < $BOX->{"X"}{"hi"}) {
        $yIndex = 1;
        $j = $BOX->{"Y"}{"lo"};
	if ($xIndex == -1) {
	    $j += $radii/2;
	}
        while ($j < $BOX->{"Y"}{"hi"}) {
            $zIndex = 1;
            $k = $BOX->{"Z"}{"lo"};
            if ($xIndex * $yIndex == -1) {
                $k += $radii/2;
            }
            while ($k < $BOX->{"Z"}{"hi"}) {
		$WAT_POS = ();
		$WAT_POS->{"XCOORD"} = $i + rand(2 * $variance) - $variance;
		$WAT_POS->{"YCOORD"} = $j + rand(2 * $variance) - $variance;
		$WAT_POS->{"ZCOORD"} = $k + rand(2 * $variance) - $variance;
		$WAT_POS->{"RADII"} = $myWat->{"RADII"};
		$WAT_POS->{"INDEX"} = $Atoms->{"tot"};
                if (placeOnGrid($GRID, $WAT_POS, $BBOX, $gridLen, $iterations)) {
		    $totWat++;
		}
		$k += $radii;
            }
            $yIndex *= -1;
            $j += $radii;
        }
        $xIndex *= -1;
	$yIndex = $xIndex;
        $i += $radii;
    }
    $iterations++;
    if ($iterations > 10) {
	print "giving up after 10 iterations! Max density reached?...";
	last;
    }
    }
    for $dim ("X", "Y", "Z") {
        $BOX->{$dim}{"lo"} -= (2.5*$variance);
        $BOX->{$dim}{"hi"} += (2.5*$variance);
    }
    return $totWat;
}

sub placeOnGrid {
    my ($GRID, $WATER, $BOX, $grid_len, $iteration) = @_;
    my ($dim, $Index, $currCell, $noCollide, $spacing);    
    my ($currAtom, $dist, $CLIST, $cListCell, $atom, $i);

    $spacing = 1.95 - (0.1 * $iteration);
    for $dim ("X", "Y", "Z") {
	$Index->{$dim} = int(($WATER->{$dim . "COORD"} - $BOX->{$dim}{"lo"})/$grid_len) + 1;
    }
	$noCollide = 1;
	$CLIST = GetNeighbours($GRID, \%{ $GRID->{$Index->{"X"}}{$Index->{"Y"}}{$Index->{"Z"} } });
	unshift @{ $CLIST }, \%{ $GRID->{$Index->{"X"}}{$Index->{"Y"}}{$Index->{"Z"}} };
      CLOOP: for $currCell (@{ $CLIST }) {
	  if (exists($currCell->{"ATOMS"}) or exists($currCell->{WATERS})) {
	      for $atom (@{ $currCell->{"ATOMS"}}, @{ $currCell->{WATERS} }) {
		  $dist = GetBondLength($WATER, $atom);
		  if ($dist < ($WATER->{"RADII"} + $atom->{"RADII"} + $spacing)) {
		      $noCollide = 0;
		      last CLOOP;
		  }
	      }
	  }
      }

    if ($noCollide) {
	#%{ $WATS{ $WATER->{"INDEX"} } } = %{ $WATER };
	&addRes($WATER, $ATOMS, $WAT, $GRID, $Index);
    }
    return $noCollide;
    
}

sub updateBox {
    my ($box) = $_[0];

    for ("X", "Y", "Z") {
	$box->{$_ . "COORD"}{hi} = $box->{$_}{hi};
	$box->{$_ . "COORD"}{lo} = $box->{$_}{lo};
    }
}


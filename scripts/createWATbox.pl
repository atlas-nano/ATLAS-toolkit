#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use constant PI => 3.14159265358979;
use General qw(GetBondLength Rotate CoM);
use FileFormats qw(createBGF addHeader createHeaders GetBGFFileInfo);
use BOX qw(CreateGrid GetGridDims GetNeighbours);
use Getopt::Std qw(getopt);
use strict;
use File::Basename qw(basename);
#use diagnostics;

# solvate.pl: This script will take an input bgf_file, determine the surface atoms and add water molecules 
# of the desired radii around it in a random packing arrangement, in order to achieve the required density

sub init;
sub addSolv;
sub numerically { ($a <=> $b); }
sub getRadii;
sub getRandomNumber;
sub createWatBox;
sub usage;
sub createBox;
sub getWatOpts;
my ($density, $dimMax, $num_wats, $saveName);
my ($start, $end, $GRIDOPTS, $loadWat);
my ($ATOMS, $BONDS, $BBOX, $MOL_VOL, $vol, $WAT);
my ($GRID, $num_cells, $HEADERS, $i);
my ($res_counter, $counter, %WATS, $grid_len, $FFILES, $waterOpt);

$start = time();
$|++;
$FFILES = &init;
getWatOpts if ($loadWat);
print "Creating bounding box...";
$MOL_VOL = &createBox($dimMax);
print "Done\nCreating Grid...\n";
$GRIDOPTS = GetGridDims($MOL_VOL, "xyz", 0, 0);
($GRID, $BBOX, $counter) = CreateGrid(undef, 0, $GRIDOPTS, $grid_len, 1);
$res_counter = 0;
for $i (keys % { $BBOX }) {
    $vol *= ($BBOX->{$i}{"hi"} - $BBOX->{$i}{"lo"});
    $BBOX->{"${i}COORD"}{"hi"} = $BBOX->{$i}{"hi"};
    $BBOX->{"${i}COORD"}{"lo"} = $BBOX->{$i}{"lo"};
}

$ATOMS->{"tot"} = 0;
$num_wats = int(($density * $vol * 6.023E-1)/18) if (! $num_wats);
&createWatBox($ATOMS, $GRID, $BBOX, $grid_len, $WAT);
$ATOMS->{"tot"} = $counter;
$num_wats = addSolv($ATOMS, $BONDS, $GRID, $WAT, $num_wats);
delete $ATOMS->{"tot"};

print "Total Simulation Box Size: $vol A^3\n";
print "Placed $num_wats solvent molecules\n";
$density = ((18 * $num_wats)/($vol * 6.023E-01));
print "Total solvent density achieved: $density (g/cm^3)\n";
print "Saving files $saveName...";
$HEADERS = createHeaders($BBOX, "water box");
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";
$end = time();
print "Elapsed: " . ($end - $start) . " seconds\n";

sub getWatOpts {
    my ($i, $CENTER, $dim);
    ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($waterOpt,0);
    $CENTER = CoM($WAT->{"ATOMS"});
    for $i (keys %{ $WAT->{ATOMS} }) {
        for $dim ("XCOORD", "YCOORD", "ZCOORD") {
            $WAT->{"ATOMS"}{$i}{$dim} -= $CENTER->{$dim};
        }
    }
    $WAT->{"RADII"} = 1.6;
    printf "SOLVENT RADII: %.1f A\n", $WAT->{"RADII"};
    $grid_len = $WAT->{RADII} + 1.0;
}

sub usage {
    my ($printStr);
    
    $printStr = <<DATA;
usage: $0 -w water type -l dim length -d [density] -n [num wats] -s [saveName] 
OPTIONS (all are requied except saveName):
    -w solvent type: Current options are
                   TIP3: the original Jorgenson TIP3 water model (rigid hydrogens)
                   TIP3_CHARMM: TIP3 water model as implemented in CHARMM
                   F3C: F3C water model (no rigid hydrogens)
		   F3C_OPT: F3C water model with optimized charges and vib. spectrum
                   SPC: SPC water model
		   MESO: Water model for Meso-scale simulations (Valeria Molinero)
                       -or- you can specify your own solvent: give the path
    -l box length: extend this in the dimension specified above. 1 number creates a cubic box,
		   multiple numbers in quotes creates an ortorhombic box
    -d [density]: the density of solvent required in g/cm3. Defaults to 1 gm/cm3.
    -n [num wats]: Optional. Specifies the maximum number of water molecules to add. Overrides the density number.
    -s [saveName]: Optional. Name of final system. If not specified will be {bgfFile}_solv.bgf
DATA
    return $printStr;
}

sub init {
    my (%OPT, $myPath, $FF, @tmp, $i);
    my (@FFILES);

    $myPath = "/ul/tpascal/scripts/dat/WAT/";
    getopt('wdlsn',\%OPT);
    ($waterOpt, $density, $dimMax, $saveName, $num_wats) = 
	($OPT{w},$OPT{d},$OPT{l},$OPT{s},$OPT{n});

    for $i ($waterOpt, $dimMax) {
	die &usage . "\n" if (! defined($i));
    }
    
    print "Initializing...";
    $density = 1 if (! defined($density) or $density !~ /^\d+[\.\d+]?/);
    
    if ($dimMax =~ /\s/) {
	@tmp = split /\s+/, $dimMax;
	if ($#tmp < 2) {
	    $tmp[1] = $tmp[0];
	    $tmp[2] = $tmp[0];
	}
	for $i (0 .. 2) {
	    die "ERROR: Expected decimal/interger value for dimension. Got \"$i\"!\n" if ($tmp[$i] !~ /^\+?\d+[\.\d+]?/);
	}
	$dimMax = ();
	$dimMax->{X} = $tmp[0];
	$dimMax->{Y} = $tmp[1];
	$dimMax->{Z} = $tmp[2];
    }

    $loadWat = 0;
    if (-e $waterOpt && -r $waterOpt && -T $waterOpt) {
	$loadWat = 1;
	$WAT->{TYPE} = "wat";
    } else {
	$waterOpt = uc($waterOpt);
	if($waterOpt eq "TIP3") {
	    ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($myPath . "tip3.bgf",0);
	    $WAT->{"RADII"} = 1.8;
	    $WAT->{TYPE} = "tip3";
	} elsif ($waterOpt eq "F3C") {
	    ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($myPath . "F3C.bgf",0);
	    $WAT->{"RADII"} = 1.8;
	    $WAT->{TYPE} = "f3c";
        } elsif ($waterOpt eq "F3C_OPT") {
            ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($myPath . "F3C_OPT.bgf",0);
            $WAT->{"RADII"} = 1.8;
	    $WAT->{TYPE} = "f3c";
	} elsif ($waterOpt eq "MESO") {
	    ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($myPath . "Meso.bgf",0);
	    $WAT->{"RADII"} = 1.6;
	    $WAT->{TYPE} = "meso";
	} elsif ($waterOpt eq "TIP3_CHARMM") {
	    ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo($myPath . "tip3_charmm.bgf",0);
	    $WAT->{"RADII"} = 1.8;
	    $WAT->{TYPE} = "tip3";
	} elsif ($waterOpt eq "DREIDING3") {
	    ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo("${myPath}dreiding3Wat.bgf",0);
            $WAT->{"RADII"} = 1.8;
            $WAT->{TYPE} = "dreiding";
        } elsif ($waterOpt eq "MW") {
            ($WAT->{"ATOMS"}, $WAT->{"BONDS"}) = GetBGFFileInfo("${myPath}mW.bgf",0);
            $WAT->{"RADII"} = 1.8;
            $WAT->{TYPE} = "mW";
	}
    }

    die "ERROR: Invalid selection for solventType: $waterOpt\n" if (! keys %{ $WAT } && ! $loadWat);
	    
    if (! $saveName) {
	$saveName = $WAT->{TYPE};
	$saveName .= "_" . $dimMax->{X} . "x";
	$saveName .= "_" . $dimMax->{Y} . "y";
	$saveName .= "_" . $dimMax->{Z} . "z";
	$saveName .= ".bgf";
    }

    $res_counter = 0;
    $vol = 1;
    $WAT->{RADII} = 1 if (! exists($WAT->{RADII}));
    $grid_len = $WAT->{RADII} + 1;
    if (defined($num_wats)) {
	$num_wats = 0 if ($num_wats !~ /^\d+/);
    } else {
	$num_wats = 0;
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
	$wC = 1;
	while ($wC <= ($addWat - $requiredWat)) { #randomely remove waters
	    $i = int(rand($#tmp));
	    if (exists($WATS{$tmp[$i]})) {
		delete $WATS{$tmp[$i]};
		$wC++;
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
    my ($watPos, $All_Atoms, $watData) = @_;
    my ($counter, $currAtm, $atom, $atomC, $bonds, $Angles, %WAT, $offset);

    $res_counter++;
    $offset = $atomC = $All_Atoms->{"tot"};
    
    # do rotation
    $Angles = getRandomNumber(360);
    Rotate(\%{ $watData->{"ATOMS"} }, $Angles, 3);
    
    for $atom (sort numerically keys %{ $watData->{"ATOMS"} }) {
	$atomC++;
	%{ $currAtm } = %{ $watData->{"ATOMS"}{$atom} };
	for $counter ("XCOORD", "YCOORD", "ZCOORD") {
	    $currAtm->{$counter} += $watPos->{$counter};
	}
	$currAtm->{"RESNUM"} = $res_counter;
	$currAtm->{"LABEL"} = "HETATM";
	for $counter (keys %{ $currAtm }) {
	    $All_Atoms->{$atomC}{$counter} = $currAtm->{$counter};
	}
	if ($watData->{"ATOMS"}{$atom}{NUMBONDS} > 0) {
	    @{ $bonds } = @{ $watData->{"BONDS"}{$atom} };
	    
	    for (@{ $bonds }) {
		push @{ $BONDS->{$atomC} }, ($_ + $offset);
	    }
	}
    }

    $All_Atoms->{"tot"} = $atomC;
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
    my ($totWat, $radii, $variance, $WAT_POS, $BOX);
    
    $radii = sqrt(2) * $myWat->{"RADII"};
    for $dim ("X", "Y", "Z") {
        $BOX->{$dim}{"lo"} = $BBOX->{$dim}{"lo"} + $radii/2;
        $BOX->{$dim}{"hi"} = $BBOX->{$dim}{"hi"} - $radii/2;
    }
    
    # get 10% of radii for deviation
    $variance = .15 * $myWat->{"RADII"};

    $i = $BOX->{"X"}{"lo"};
    
    $xIndex = $yIndex = $zIndex = 1;
    $Atoms->{"tot"}++;
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
                if (placeOnGrid($GRID, $WAT_POS, $BBOX, $gridLen)) {
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
    
    for $dim ("X", "Y", "Z") {
        $BOX->{$dim}{"lo"} -= $radii/2;
        $BOX->{$dim}{"hi"} += $radii/2;
    }
}

sub placeOnGrid {
    my ($GRID, $WATER, $BOX, $grid_len) = @_;
    my ($dim, $Index, $currCell, $noCollide, $waterEff);    
    my ($currAtom, $dist, $CLIST, $cListCell, $atom);

    $noCollide = 0;
    for $dim ("X", "Y", "Z") {
	$Index->{$dim} = int(($WATER->{$dim . "COORD"} - $BOX->{$dim}{"lo"})/$grid_len) + 1;
    }
    if (! $GRID->{$Index->{"X"}}{$Index->{"Y"}}{$Index->{"Z"}}{"Surface"}) {
	$noCollide = 1;
    } else {
	$CLIST = GetNeighbours($GRID, \%{ $GRID->{$Index->{"X"}}{$Index->{"Y"}}{$Index->{"Z"} } });
	unshift @{ $CLIST }, \%{ $GRID->{$Index->{"X"}}{$Index->{"Y"}}{$Index->{"Z"}} };
      CLOOP: for $currCell (@{ $CLIST }) {
	  if (exists($currCell->{"ATOMS"})) {
	      for $atom (@{ $currCell->{"ATOMS"}}) {
		  $dist = GetBondLength($WATER, $atom);
		  if ($dist < ($WATER->{"RADII"} + $atom->{"RADII"})) {
		      $noCollide = 0;
		      last CLOOP;
		  }
	      }
	  }
      }
    }

    if ($noCollide) {
	%{ $WATS{ $WATER->{"INDEX"} } } = %{ $WATER };
	$ATOMS->{"tot"}++;
    }
    return $noCollide;
    
}

sub createBox {
    my ($dim) = $_[0];
    my (%BOX);

    %BOX = (
                "X" => {
                    "hi"    => $dim->{X},
                    "lo"    => 0,
                    "len"   => 0,
                    "angle" => 90,
                },
                "Y" => {
                    "hi"    => $dim->{Y},
                    "lo"    => 0,
                    "len"   => 0,
                    "angle" => 90,
                },
                "Z" => {
                    "hi"    => $dim->{Z},
                    "lo"    => 0,
                    "len"   => 0,
                    "angle" => 90,
                },
            );
    
    return \%BOX;
}

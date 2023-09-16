package Packages::BOX;

require Exporter;
use strict;
use Cwd;

our (@ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = qw(CreateGrid GetNeighbours GetBox MapGrid PrintBox GetRadii CenterAtoms GetSurface 
	     GetGridDims ConvertBox MakeBox  MoveAtomsToOrigin PlaceAtomsOnGrid);
$VERSION = "1.00";

sub getBoxDims {
    my ($BBOX, $grid_len, $printGrid) = @_;
    my ($counter, $bLen);

    if ($printGrid) {
	PrintBox("Total bounding box for atom centers:", $BBOX);
    }

    for $counter (keys %{ $BBOX }) {
	$bLen = $BBOX->{$counter}{"hi"} - $BBOX->{$counter}{"lo"};
	if ($bLen =~ /\d+\.\d+/ || $bLen % $grid_len > 0) {
	    $BBOX->{$counter}{"hi"} = ((int($bLen/$grid_len) + 1) * $grid_len) + $BBOX->{$counter}{"lo"};
	}
    }

}

sub CreateGrid {
    my ($ATOMS, $radii, $BBOX, $grid_len, $printGrid) = @_;
    my ($counter, $CURR_BOX, $Index, %GRID, $AtomC, $total_cells, $bLen);
    my ($currCell, $burried, $exclude);

    if (! defined($printGrid)) {
	$printGrid = 0;
    }

    $AtomC = $total_cells = $exclude = $burried = 0;
    if ($radii) {
	getBoxDims($BBOX, $grid_len, $printGrid);
    }

    $CURR_BOX->{x1} = $BBOX->{X}{lo};
    $CURR_BOX->{y1} = $BBOX->{Y}{lo};
    $CURR_BOX->{z1} = $BBOX->{Z}{lo};
    $Index->{x} = 0;
    $Index->{y} = 0;
    $Index->{z} = 0;

    while ($CURR_BOX->{x1} < $BBOX->{X}{hi}) {
	$Index->{x} += 1;
	$Index->{y} = 0;
	$CURR_BOX->{x2} = $CURR_BOX->{x1} + $grid_len;
	$CURR_BOX->{y1} = $BBOX->{Y}{lo};
	while ($CURR_BOX->{y1} < $BBOX->{Y}{hi}) {
	    $Index->{y} += 1;
	    $Index->{z} = 0;
	    $CURR_BOX->{y2} = $CURR_BOX->{y1} + $grid_len; 
	    $CURR_BOX->{z1} = $BBOX->{Z}{lo};
	    while ($CURR_BOX->{z1} < $BBOX->{Z}{hi}) {
		$total_cells++;
		$Index->{z} += 1;
		$CURR_BOX->{z2} = $CURR_BOX->{z1} + $grid_len; 
		$currCell = \%{ $GRID{$Index->{x}}{$Index->{y}} };
		$currCell->{$Index->{z}} = (
					      {
						  "XINDEX" => $Index->{x},
						  "YINDEX" => $Index->{y},
						  "ZINDEX" => $Index->{z},
						  "X"      => (
							       {
								   "lo" => $CURR_BOX->{x1},
								   "hi" =>  $CURR_BOX->{x2},
							       },
							       ),
						  "Y"      => (
							       {
								   "lo" => $CURR_BOX->{y1},
								   "hi" =>  $CURR_BOX->{y2},
							       },
							       ),
						  "Z"      => (
							       {
								   "lo" => $CURR_BOX->{z1},
								   "hi" =>  $CURR_BOX->{z2},
							       },
							       ),
					      }
					      );
		$CURR_BOX->{z1} += $grid_len;
	    }
	    $CURR_BOX->{y1} += $grid_len;
	}
	$CURR_BOX->{x1} += $grid_len;
    }
    ($GRID{X}{tot}, $GRID{Y}{tot}, $GRID{Z}{tot}) = ($Index->{x}, $Index->{y}, $Index->{z});
    ($AtomC, $exclude, $burried) = PlaceAtomsOnGrid($ATOMS, \%GRID, $BBOX, $grid_len) if (defined($ATOMS));

    if ($printGrid) {
	PrintBox("Total vdw box size:", $BBOX);
	print "Total Atoms: $AtomC\n";
	print "Total Cells: $total_cells\n";
	print "Cells with Atoms: $exclude (" . sprintf("%.2f", (100 * ($exclude/$total_cells))) . ")%\n";
	print "Inital burried Cells:  $burried (" . sprintf("%.2f", (100 * ($burried/$total_cells))) . "%)\n";
    }

#    $burried = setAllBurriedCells(\%GRID);
    if ($printGrid) {
	print "Total Burried Cells: $burried (" . sprintf("%.2f", (100 * ($burried/$total_cells))) . "%)\n";
	$exclude -= $burried;
	print "Cells defining molecular surface: $exclude (" . sprintf("%.2f", (100 * ($exclude/$total_cells))) . "%)\n";
    }

    return (\%GRID, $BBOX, $AtomC);
}

sub setAllBurriedCells {
    my ($GRID) = $_[0];
    my ($x, $y, $z, $cell, $burried, $CLIST, $nCells, $newBurried, $j, $newCell);
    my ($xIndex, $yIndex, $zIndex);

    $burried = 0;
    $newBurried = 1;
    while ($newBurried) {
	$newBurried = 0;
	for $x (keys %{ $GRID }) {
	    for $y (keys %{ $GRID->{$x} }) {
	      CELLLOOP: for $z (keys %{ $GRID->{$x}{$y} }) {
		  $cell = \%{ $GRID->{$x}{$y}{$z} };
		  if ($cell->{Burried}) {
		      $cell->{Surface} = 0;
		      $burried++;
		      next CELLLOOP;
		  }
		  $CLIST = GetNeighbours($GRID, $cell);
		  for $nCells (@{ $CLIST }) {
		      $xIndex = $nCells->{XINDEX};
		      $yIndex = $nCells->{YINDEX};
		      $zIndex = $nCells->{ZINDEX};
		      $newCell = \%{ $GRID->{$nCells->{XINDEX}}{$nCells->{YINDEX}}{$nCells->{ZINDEX}} };
		      next CELLLOOP if (! exists($newCell->{VOL}) || $newCell->{VOL} == 0 || ! $newCell->{Burried});
		  }
		  $cell->{Burried} = 1; #if all of my neighbours are burried then i am as well
		  $cell->{Surface} = 0;
		  $burried++;
		  $newBurried = 1;
	      }
	  }
	}
    }
    return $burried;
}
		
sub PlaceAtomsOnGrid {
    my ($ATOMS, $GRID, $BOX, $grid_len) = @_;
    my ($counter, $AtomC, $Index, $dim, @deleted_keys, $i);  
    my ($atom, $cell, $maxVol, $exclude, $atomVol, $burried);
    my ($pi) = atan2(1,1) *4;
    
    $maxVol = $grid_len**3;

    $AtomC = $counter = $exclude = $burried = 0;
    @deleted_keys = keys %{ $ATOMS };
    while ($#deleted_keys > -1) {
	$counter = pop @deleted_keys;
	for $dim ("X", "Y", "Z") {
	    $i = int(($ATOMS->{$counter}{$dim . "COORD"} - $BOX->{$dim}{lo})/$grid_len) + 1;
	    $Index->{$dim} = $i;
	}
	$AtomC++;
	$atom = \%{ $ATOMS->{$counter} };
	next if (! exists($GRID->{$Index->{X}}) or
		 ! exists($GRID->{$Index->{X}}{$Index->{Y}}) or
		 ! exists($GRID->{$Index->{X}}{$Index->{Y}}{$Index->{Z}}));
	$cell = \%{ $GRID->{$Index->{X}}{$Index->{Y}}{$Index->{Z}} };
	$atom->{CELL}{XINDEX} = $Index->{X};
	$atom->{CELL}{YINDEX} = $Index->{Y};
	$atom->{CELL}{ZINDEX} = $Index->{Z};
	
	$atomVol  = 0;
	$cell->{VOL} = 0 if (! exists($cell->{VOL}));
	if ($atom->{IS_SOLVENT}) {
	    push @{ $cell->{WATERS} }, $ATOMS->{$counter};
	} elsif ($atom->{NUMBONDS} == 0) {
	    push @{ $cell->{IONS} }, $ATOMS->{$counter};
	} elsif ($atom->{RESNAME} !~ /(WAT|Na|Mg|Cl)/i) {
	    push @{ $cell->{ATOMS} }, $ATOMS->{$counter};
	    $atomVol  = (4 * $atom->{"RADII"}**3 * $pi)/3;
	    $exclude++ if ($#{ $cell->{ATOMS} } == 0);
	} elsif ($atom->{RESNAME} eq "WAT" or $atom->{RESNAME} eq "HOH") {
	    push @{ $cell->{WATERS} }, $ATOMS->{$counter};
	} else {
	    push @{ $cell->{IONS} }, $ATOMS->{$counter};
	}
				   
	$cell->{VOL} += $atomVol;
	
	if ($cell->{VOL} > 0) {
	    setExcludedVol($cell);
	    if ($cell->{VOL} >= $maxVol) {
		if (! $cell->{Burried}) {
		    $cell->{Burried} = 1;
		    $cell->{Surface} = 0;
		    $burried++;
		}
	    } else {
		$cell->{Burried} = 0;
		$cell->{Surface} = 1;
	    }
	}
    }
    return ($AtomC, $exclude, $burried);
}

sub setExcludedVol {
    my ($cell) = $_[0];
    my ($atom, $radii, $dim);

    for $atom (@{ $cell->{ATOMS} }, @{ $cell->{IONS} }) {
	$radii = $atom->{RADII};
	for $dim ("X", "Y", "Z") {
	    $cell->{EXCLUDE}{$dim}{hi} = ($atom->{$dim . "COORD"} + $radii)
		if (! exists($cell->{EXCLUDE}{$dim}{hi}) or 
		    ($atom->{$dim . "COORD"} + $radii) > $cell->{EXCLUDE}{$dim}{hi});
	    $cell->{EXCLUDE}{$dim}{"lo"} = ($atom->{$dim . "COORD"} - $radii)
		if (! exists($cell->{EXCLUDE}{$dim}{lo}) or 
		    ($atom->{$dim . "COORD"} - $radii) < $cell->{EXCLUDE}{$dim}{lo});
	}
    }
}

sub GetNeighbours {
    my ($GRID, $CELL, $layers) = @_;
    my (@CLIST, @range, $i);
    my ($xR, $yR, $zR, $xI, $yI, $zI);
    
    if (! defined($layers)) {
	@range = (-1,0,1);
    } else {
	for $i ((-1 * $layers) .. $layers) {
	    push @range, $i;
	}
    }
    if (! defined($CELL->{XINDEX}) ||
	! defined($CELL->{YINDEX}) ||
	! defined($CELL->{ZINDEX})) {
	    return ();
    }
    for $xR (@range) {
	$xI = $CELL->{XINDEX} + $xR;
	if (! exists($GRID->{$xI})) {
	    #print "\nwrapping... $xI";
	    if ($xI == 0) { # if at the end of box, wrap back to beginning
		$xI =  $GRID->{X}{tot};
	    }elsif ($xI > 0) { # if at the end of box, wrap back to beginning
		$xI -= $GRID->{X}{tot};
	    } else { # else at the beginning to wrap to end
	        $xI += $GRID->{X}{tot};
	    }
	    #print "-> $xI (x tot:  $GRID->{X}{tot}\n";
	}
	for $yR (@range) {
	    $yI = $CELL->{YINDEX} + $yR;
	    if (! exists($GRID->{$xI}{$yI})) {
		#print "\nwrapping... $yI";
		if ($yI == 0) { # if at the end of box, wrap back to beginning
		    $yI =  $GRID->{Y}{tot};
		}elsif ($yI > 0) { # if at the end of box, wrap back to beginning
		    $yI -= $GRID->{Y}{tot};
		} else { # else at the beginning to wrap to end
		    $yI += $GRID->{Y}{tot};
		}
		#print "-> $yI (y tot:  $GRID->{Y}{tot}\n";
	    }
	    for $zR (@range) {
		$zI = $CELL->{ZINDEX} + $zR;
		delete $GRID->{$xI}{$yI}{$zI} if (! keys %{ $GRID->{$xI}{$yI}{$zI} });
		if (! exists($GRID->{$xI}{$yI}{$zI})) {
		    #print "\nwrapping... $zI";
		    if ($zI == 0) { # if at the end of box, wrap back to beginning
			$zI =  $GRID->{Z}{tot};
		    }elsif ($zI > 0) { # if at the end of box, wrap back to beginning
			$zI -= $GRID->{Z}{tot};
		    } else { # else at the beginning to wrap to end
			$zI += $GRID->{Z}{tot};
		    }
		    #print "-> $zI (z tot:  $GRID->{Z}{tot}\n";
		}
		push @CLIST, \%{ $GRID->{$xI}{$yI}{$zI} };
	    }
	}
    }
    
    return \@CLIST;
}

sub IsInBox(@) {
    my ($Atom, $Box) = @_;
    my ($index, $returnval);

    $returnval = 0;
    if ($Atom->{"XCOORD"} >= $Box->{"x1"} and $Atom->{"XCOORD"} <= $Box->{"x2"}) {
	if ($Atom->{"YCOORD"} >= $Box->{"y1"} and $Atom->{"YCOORD"} <= $Box->{"y2"}) {
	    if ($Atom->{"ZCOORD"} >= $Box->{"z1"} and $Atom->{"ZCOORD"} <= $Box->{"z2"}) {
		$returnval = 1;
	    }
	}
    }

    return $returnval;

}

sub GetBox {
    my ($ATOMS, $PARMS, $HEADERS) = @_;
    my (%BOX, $isValid, $i, $j, $coord); 
    my ($cVal, @tmp, $radii, @vals);
    
    $isValid = 0;
    @tmp = (-9999,99999,90,-99999,99999,90,-99999,99999,90);
    if (defined($HEADERS)) {
	for $i (@{ $HEADERS }) {
	    if ($i =~ /^CRYSTX\s+(.+)/) {
		@vals = split /\s+/, $1;
		for $j (0 .. 2) {
		    $tmp[($j * 3)] = $vals[$j];
		    $tmp[($j * 3) + 1] = 0,
		    $tmp[($j * 3) + 2] = $vals[3+$j];
		}
		$isValid = 1;
	    }
	}
    }
    
    if (! $isValid) {
	for $i (keys %{ $ATOMS }) {
	    $radii = GetRadii($ATOMS->{$i}, $PARMS);
	    if (! defined($radii)) {
		$radii = 2;
	    }
	    
	    $j = 0;
	    for $coord ("X", "Y", "Z") {
		$cVal = $ATOMS->{$i}{$coord . "COORD"};
		$tmp[$j] = ($cVal + $radii)
		    if (($cVal + $radii) > $tmp[$j]);
		$tmp[$j + 1] = ($cVal - $radii)
		    if (($cVal - $radii) < $tmp[$j + 1]);
		$j += 3;
	    }
	}
    }

    %BOX = (
		"X" => {
		    "hi"    => $tmp[0],
		    "lo"    => $tmp[1],
		    "len"   => $tmp[0] - $tmp[1],
		    "angle" => $tmp[2],
		},
		"Y" => {
		    "hi"    => $tmp[3],
		    "lo"    => $tmp[4],
                    "len"   => $tmp[3] - $tmp[4],
		    "angle" => $tmp[5],
		},
		"Z" => {
		    "hi"    => $tmp[6],
		    "lo"    => $tmp[7],
                    "len"   => $tmp[6] - $tmp[7],
		    "angle" => $tmp[8],
		},
	    );
    
    return \%BOX;
}

sub GetRadii {
    my ($atom, $PAR) = @_;
    my ($returnVal, $atmName);

    $atmName = $atom->{"FFTYPE"};
    if (! $PAR or ! exists($PAR->{"VDW"}{$atmName}{$atmName})) {
	$atmName = $atom->{"ATMNAME"};
    }

    if (exists($PAR->{"VDW"}{$atmName}{$atmName})) {
	$returnVal = $PAR->{"VDW"}{$atmName}{$atmName}{1}{"VALS"}[1]/2;
    }

    return $returnVal;
}

sub MapGrid {
    my ($CELLS, $tot_cells, $water_rad, $BOX) = @_;
    my ($xindex, $yindex, $zindex, $counter, $holder, $index, $exclude);
    my ($xc, $yc, $zc, $radii, $cell_vol, $atoms_vol, $burried);


    $burried = $exclude = 0;
    for $xindex (keys %{ $CELLS }) {
	for $yindex (keys %{ $CELLS->{$xindex} }) {
	    for $zindex (keys %{ $CELLS->{$xindex}{$yindex} }) {
		$counter = \%{ $CELLS->{$xindex}{$yindex}{$zindex} };
		
		for $index ("X", "Y", "Z") {
		    if ($BOX->{$index}{"hi"} < $counter->{$index}{"hi"}) {
			$BOX->{$index}{"hi"} = $counter->{$index}{"hi"};
		    }
		    if ($BOX->{$index}{"lo"} > $counter->{$index}{"lo"}) {
			$BOX->{$index}{"lo"} = $counter->{$index}{"lo"};
		    }
		    
		}

		$cell_vol = ($counter->{"X"}{"hi"} - $counter->{"X"}{"lo"}) *
		     ($counter->{"Y"}{"hi"} - $counter->{"Y"}{"lo"}) *
		     ($counter->{"Z"}{"hi"} - $counter->{"Z"}{"lo"});
		$atoms_vol = 0;
		for $holder (0 .. $#{ $counter->{"ATOMS"} }) {
		    $xc = $counter->{"ATOMS"}[$holder]{"XCOORD"};
		    $yc = $counter->{"ATOMS"}[$holder]{"YCOORD"};
		    $zc = $counter->{"ATOMS"}[$holder]{"ZCOORD"};
		    $radii = $counter->{"ATOMS"}[$holder]{"RADII"};
		    
		    if (! defined($counter->{"EXCLUDE"})) {
			$counter->{"EXCLUDE"}{"xhi"} = $xc + $radii;
			$counter->{"EXCLUDE"}{"xlo"} = $xc - $radii;
			$counter->{"EXCLUDE"}{"ylo"} = $yc - $radii;
			$counter->{"EXCLUDE"}{"yhi"} = $yc + $radii;
			$counter->{"EXCLUDE"}{"zhi"} = $zc + $radii;
			$counter->{"EXCLUDE"}{"zlo"} = $zc - $radii;
		    } else {
			if (($xc + $radii) > $counter->{"EXCLUDE"}{"xhi"}) {
			    $counter->{"EXCLUDE"}{"xhi"} = $xc + $radii;
			}elsif (($xc - $radii) < $counter->{"EXCLUDE"}{"xlo"}) {
			    $counter->{"EXCLUDE"}{"xlo"} = $xc - $radii;
			}elsif (($yc + $radii) > $counter->{"EXCLUDE"}{"yhi"}) {
			    $counter->{"EXCLUDE"}{"yhi"} = $yc + $radii;
			}elsif (($yc - $radii) < $counter->{"EXCLUDE"}{"ylo"}) {
			    $counter->{"EXCLUDE"}{"ylo"} = $yc - $radii;
			}elsif (($zc + $radii) > $counter->{"EXCLUDE"}{"zhi"}) {
			    $counter->{"EXCLUDE"}{"zhi"} = $zc + $radii;
			}elsif (($zc - $radii) < $counter->{"EXCLUDE"}{"zlo"}) {
			    $counter->{"EXCLUDE"}{"zlo"} = $zc - $radii;
			}
		    }

		    $atoms_vol += ($radii * 2)^3;
		}
		if (defined($counter->{"EXCLUDE"})) {
		    $atoms_vol = ($counter->{"EXCLUDE"}{"xhi"} - $counter->{"EXCLUDE"}{"xlo"}) *
			($counter->{"EXCLUDE"}{"yhi"} - $counter->{"EXCLUDE"}{"ylo"}) *
			($counter->{"EXCLUDE"}{"zhi"} - $counter->{"EXCLUDE"}{"zlo"});
		    $atoms_vol += ($water_rad * 2)^3;
		    
		    if ($atoms_vol >= $cell_vol) {
			$burried++;
			$counter->{Surface} = 0;
		    } else {
			$counter->{Surface} = 1;
		    }
		    $exclude++;
		} else {
		    $counter->{Surface} = 1;
		}
	    }
	}
    }

    PrintBox("Total vdw box size:", $BOX);

    print "Total Cells: $tot_cells\n";
    print "Burried Cells: $burried (" . sprintf ("%.2f", (100 * ($burried/$tot_cells))) . " %)\n";
    print "Cells with Atoms: $exclude (" . sprintf ("%.2f", (100 * ($exclude/$tot_cells))) . " %)\n";

    return ($CELLS, $BOX);
}

sub PrintBox {
    my ($intext, $inBox) = @_;
    my ($counter, $out_string);

    $out_string = sprintf("%-40s", $intext);
    for $counter ("X", "Y", "Z") {
	$out_string .= sprintf("%8.3f ", ($inBox->{$counter}{"hi"} - $inBox->{$counter}{"lo"}));
    }

    print "$out_string\n";

}

sub CenterAtoms {
    my ($ATOMS, $BOX, $radii) = @_;
    my (%Offset, $dim, $atomC);

    for $dim ("X", "Y", "Z") {
	$Offset{$dim} = $BOX->{$dim}{"lo"};
    }

    for $atomC (keys %{ $ATOMS }) {
	for $dim ("X", "Y", "Z") {
	    $ATOMS->{$atomC}{$dim . "COORD"} -= $Offset{$dim};
	}
    }
   
}

sub GetSurface {
    my ($GRID) = $_[0];
    my ($i, $j, $k, $currCell, %SURFACE, $sCell, $atom, $charge, $totCharge, $dim, $ions);

    $sCell = $charge = $totCharge = $ions = 0;
    for $i (keys %{ $GRID }) {
	for $j (keys %{ $GRID->{$i} }) {
	    for $k (keys %{ $GRID->{$i}{$j} }) {
		$currCell = \%{ $GRID->{$i}{$j}{$k} };
		#next if (! $currCell->{Surface});
		if (defined($currCell->{ATOMS})) {
		    $sCell++;
		    $charge = 0;
                    for $dim ("X", "Y", "Z") {
                        $SURFACE{$sCell}{$dim} = $currCell->{$dim . "INDEX"};
                    }
		    for $atom (@{ $currCell->{ATOMS} }, @{ $currCell->{IONS} }) {
			$charge += $atom->{CHARGE};
		    }
		    $SURFACE{$sCell}{CHARGE} = $charge;
		    $SURFACE{$sCell}{BURRIED} = 0;
		    $totCharge += $charge;
		}elsif (defined($currCell->{IONS})) {
		    $ions++;
		    for $atom (@{ $currCell->{IONS} }) {
			$charge += $atom->{CHARGE};
		    }
		}
	    }
	}
    }

    print "Cells defining molecular surface: $sCell\n";
    print "Surface Cells containing ions: $ions\n";
    printf "Total charge of surface atoms: %.3f\n", $totCharge;

    return (\%SURFACE, $totCharge);

}

sub getTriclinicParms {
    my ($box, $cell) = @_;
    my ($hyp, $angle);
    my ($PI) = atan2(1,1) * 4;

    $angle = $box->{$cell->{3}}{angle};
    $angle = $angle % 180;
    $angle = 180 - $angle if ($angle > 90);
    $angle *= ($PI/180);
    $hyp = $box->{$cell->{2}}{len};

    return ($angle, $hyp);
}

sub GetGridDims {
    my ($molBox, $dimOpt, $min, $max) = @_;
    my ($i, %GRID, $j, $angle, $c1, $c2, $hyp, @TRICLINIC);

    @TRICLINIC = (
		  (
		   {
		       1 => "X",
		       2 => "Y",
		       3 => "Z",
		   },
		   ),
		  (
		   {
		       1 => "X",
		       2 => "Z",
		       3 => "Y",
		   },
		   ),
		  (
		   {
		       1 => "Y",
		       2 => "Z",
		       3 => "X",
		   },
		   )
		  );

#corrections for triclinic cell
    for $i (@TRICLINIC) {
	($angle, $hyp) = getTriclinicParms($molBox,$i);
	next if (cos($angle) < 0.01);
	$c1 = $hyp * cos($angle)/2;
	$c2 = $hyp * (1 - sin($angle)); 
	$molBox->{$i->{2}}{lo} += $c1;
	$molBox->{$i->{2}}{hi} -= $c1;
	$molBox->{$i->{2}}{len} -= ($c1 * 2);

	$molBox->{$i->{1}}{hi} -= $c2;
	$molBox->{$i->{1}}{lo} += $c2;
	$molBox->{$i->{1}}{len} -= ($c2 * 2);
    }

    for $i (keys %{ $molBox }) {
	for $j ("hi", "lo") {
	    $GRID{$i}{$j} = $molBox->{$i}{$j};
	}
    }

    for $i ("x","y","z") {
	$j = uc($i);
	if (lc($dimOpt) =~ /$i/) {
	    if ($min =~ /\+(\d+\.*\d*)/) {
		$GRID{$j}{lo} -= $1;
	    }elsif ($min < $molBox->{$j}{lo}) {
	        $GRID{$j}{lo} = $min;
	    }

            if ($max =~ /\+(\d+\.*\d*)/) {
                $GRID{$j}{hi} += $1;
            }elsif ($max > $molBox->{$j}{hi}) {
                $GRID{$j}{hi} = $max;
            }
	}
    }

    return \%GRID;
}

sub MakeBox { 
    my ($box) = $_[0];
    my (%BOX, $i, @dim);

    @dim = ("X", "Y", "Z");
    for $i (0 .. 2) {
        $BOX{$dim[$i]}{lo} = 0;
        $BOX{$dim[$i]}{hi} = $box->{$i + 2}{DATA};
        $BOX{$dim[$i] . "COORD"}{lo} = 0;
        $BOX{$dim[$i] . "COORD"}{hi} = $box->{$i + 2}{DATA};
        $BOX{$dim[$i] . "COORD"}{len} = $box->{$i + 2}{DATA};
    }

    return \%BOX;
}

sub ConvertBox {
    my ($lammpsBox) = $_[0];
    my (%amberBox, $i);

    $amberBox{1}{DATA} = 90;
    for $i (0 .. 2) {
        $amberBox{$i + 2}{DATA} = $lammpsBox->[$i]{hi} - $lammpsBox->[$i]{lo};
    }

    return \%amberBox;
}

sub MoveAtomsToOrigin {
    my ($atoms) = $_[0];
    my ($MIN, $i, $j);

    for my $i (keys %{ $atoms }) {
	for $j ("XCOORD", "YCOORD", "ZCOORD") {
	    $MIN->{$j} = $atoms->{$i}{$j} if (! exists($MIN->{$j}) or $atoms->{$i}{$j} < $MIN->{$j});
	}
    }

    for my $i (keys %{ $atoms }) {
        for $j ("XCOORD", "YCOORD", "ZCOORD") {
            $atoms->{$i}{$j} -= $MIN->{$j};
        }
    }

}
1;

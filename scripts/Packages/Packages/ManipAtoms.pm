package ManipAtoms;

require Exporter;
use strict;
use Storable qw(dclone);
use Cwd;
use General qw(CoM GetBondLength GetMinDist LoadElements FindElement FindElementByMass SearchAroundAtom);
use REPLICATE qw (GetBoxDisplacementTensor TransAtom);
use BOX qw(GetBox CreateGrid GetNeighbours);
no warnings 'recursion';

our (@ISA, @EXPORT, @EXPORT_OK, $VERSION);

@ISA = qw(Exporter);
$VERSION = "1.00";
@EXPORT = ();
@EXPORT_OK = qw(ImageAtoms ScaleAtoms UnwrapAtoms CenterSystem GetAtmList 
				RotateAbout GetSolvent SplitAtomsByMol GetAtmData GetMols MapOrigin 
				GetAbituraryRotationMatrix GroupAtomsByField MakeFieldSequential
				BuildAtomSelectionString SelectAtoms ReimageAtoms CreateBondsByDistance
				RemoveMolLongBonds ReimageMols SelectAtomsByField AddMolsToSelection);

sub numerically { ($a<=>$b); }

sub MapOrigin {
	my ($box, $com) = @_;
	my ($boxCenter, $offset);

	for ("X","Y","Z") {
		$boxCenter->{$_} = ($box->{$_}{hi} - $box->{$_}{lo})/2;
		$offset->{$_} = $boxCenter->{$_} - $com->{$_ . "COORD"};
		$box->{$_}{hi} -= $offset->{$_};
		$box->{$_}{lo} -= $offset->{$_};
	}
}

sub ReimageMols {
	my ($atoms, $mols, $box) = @_;
	my ($i, $j, $k, $cMol, $com, $tmp, $bCenter, $offset, @m, $t);

	&GetBoxDisplacementTensor($box);
	%{ $bCenter } = (XCOORD =>$box->{X}{center}, YCOORD=>$box->{Y}{center}, ZCOORD=>$box->{Z}{center});
	@m = ([-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1]);
	for $i (keys %{ $mols }) {
		$cMol = GetAtmData($atoms, $mols->{$i}{MEMBERS});
		$com = General::CoM($cMol, $box);
		%{ $tmp } = %{ $com };
		%{ $offset } = ();
		$t = 0;
		for $j ("X", "Y", "Z") {
			&findMinDist($bCenter,$tmp,$box,$j,\@{ $m[$t++] },\@{ $m[$t++] });
			$offset->{"${j}COORD"} = $tmp->{"${j}COORD"} - $com->{"${j}COORD"};
		}
		next if ($offset->{XCOORD}==0 and $offset->{YCOORD}==0 and $offset->{ZCOORD}==0);
		for $k (keys %{ $mols->{$i}{MEMBERS} }) {
			for $j ("XCOORD", "YCOORD", "ZCOORD") {
				$atoms->{$k}{$j} += $offset->{$j};
			}
		}
	}
}

sub RemoveMolLongBonds {
	my ($atoms, $bonds, $mols, $box) = @_;
	my ($i, $tmp);
	&GetBoxDisplacementTensor($box);
	for $i (keys %{ $mols }) {
		@{ $tmp } = %{ $mols->{$i}{MEMBERS} };
		&removeLongBonds($atoms, $bonds, $box, $tmp->[0]);
	}
}

sub removeLongBonds {
	my ($atoms, $bonds, $box, $currAtm) = @_;
	my (@m, $i, $j, $a1, $a2, $a3, $d1, $d2, $t, $dir);

	@m = ([-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1]);
	$a1 = $atoms->{$currAtm};
	$a1->{FIXED} = 1;
	for $i (@{ $bonds->{$currAtm} }) {
		$a2 = $atoms->{$i};
		next if (exists($a2->{FIXED}));
		$t=0;
		for $j ("X", "Y", "Z") {
			&findMinDist($a1,$a2,$box,$j,\@{ $m[$t++] },\@{ $m[$t++] });
		}
	}
	for $i (@{ $bonds->{$currAtm} }) {
		next if (exists($atoms->{$i}{FIXED}));
		&removeLongBonds($atoms, $bonds, $box, $i);
	}
}

sub findMinDist {
	my ($a1, $a2, $box, $d, $c1, $c2) = @_;
	my ($d1, $d2, $a3, $valid);

	$valid = 1;
	while ($valid) {
		$valid = 0;
		$d1 = $a1->{"${d}COORD"} - $a2->{"${d}COORD"}; #original cell distance
		$a3 = TransAtom($a2,$box,$c1);
		$d2 = $a1->{"${d}COORD"} - $a3->{"${d}COORD"}; #distance +cell
		if($d2*$d2>$d1*$d1) {
			$a3 = TransAtom($a2,$box,$c2);
			$d2 = $a1->{"${d}COORD"} - $a3->{"${d}COORD"}; #distance -cell
			$valid = 1 if ($d2*$d2<$d1*$d1);
		} else {
			$valid = 1;
		}
		$a2->{"${d}COORD"} = $a3->{"${d}COORD"} if ($valid); #update new atom position and loop
	}
}

sub CreateBondsByDistance {
	my ($atoms, $bonds, $box, $atomSel, $ff, $max_dist, $max_bond) = @_;
	my ($c, $d, $i, $j, $k, $l, $dist, $a1, $a2, $grid, $tol, $rec, $r0, $score, $v); 
	my ($i1, $i2, $t1, $t2, $aCells, $tmp, $dCells, $candidates, $ele, $eNum, $min_dist);

	srand (time ^ ($$ + ($$ << 11))); #initiate new random number generator
	$tol = 0.1; #10% tolerance
	$max_dist = 1.6 if (! defined($max_dist)); 

	$ele = &General::LoadElements();

	#first delete all bonds to selected atoms
	for $i (keys %{ $atoms }) {
		for $j (@{ $bonds->{$i} }) {
			for $k (0 .. $#{ $bonds->{$j} }) {
				splice @{ $bonds->{$j} }, $k, 1 if(exists($atomSel->{ $bonds->{$j}[$k] }));
			}
		}
	}
	for $i (keys %{ $atomSel }) {
		for $j ("DISPX", "DISPY", "DISPZ", "ORDER") {
			delete $atoms->{$i}{$j};
		}
		$atoms->{$i}{NUMBONDS} = 0;
		delete $atomSel->{$i};
		$j =0;
		if(!defined($max_bond)) {
			if(!defined($atoms->{$i}{ELEMENT})) {
				if(exists($atoms->{$i}{MASS})) {
					$eNum = General::FindElementByMass($atoms->{$i}{MASS}, $ele);
				} else {
					$eNum = General::FindElement($ele, $atoms->{$i}{ATMNAME});
				}
			} else {
				$eNum = $atoms->{$i}{ELEMENT}{NUM};
			}
			$atoms->{$i}{MAXBOND} = $ele->{$eNum}{NVALENCE};
		} else {
			$atoms->{$i}{MAXBOND} = $max_bond;
		}
		while($j <= $#{ $bonds->{$i} }) {
			$k = $bonds->{$i}[$j];
			if(exists($atomSel->{$k})) {
				splice @{ $bonds->{$i} }, $j, 1;
			} else {
				$j++;
			}
		}
		$atomSel->{$i} = $atoms->{$i};
	}

	#initialize grid for searching
	if(defined($box)) {
		#get box displacement tensor
		&GetBoxDisplacementTensor($box);
	} else {
		$box = GetBox($atomSel,undef,undef,10.0);
	}
	($grid, undef, undef) = CreateGrid($atomSel, 3, $box, 2.5, 0);

	#now find new bonds
	for $l (sort keys %{ $atomSel }) {
		$a1 = $atomSel->{$l};
		$t1 = $a1->{FFTYPE};
		$i1 = $a1->{INDEX};
		($aCells, $dCells) = GetNeighbours($grid, $a1->{CELL}, 1, undef, 1); #atoms and cell displacement
		$candidates = (); #hash for candidates to bond too
		for $c (0 .. $#{ $aCells }) {
			$i = $aCells->[$c];
			$d = $dCells->[$c];
			$tmp = (); #temporary array storing neighbors
			for ("ATOMS", "WATERS", "IONS") {
				push @{ $tmp }, @{ $i->{$_} } if (exists($i->{$_}));
			}
			for $a2 (@{ $tmp}) {
				$i2 = $a2->{INDEX};
				next if (exists($v->{$i1}{$i2}) or exists($v->{$i2}{$i1})); #skip if i,j pair already considered
				$v->{$i1}{$i2} = $v->{$i2}{$i1} = 1; #set considered flag
				$t2 = $a2->{FFTYPE};
				if(defined($ff) and exists($ff->{BONDS}{$t1}) and exists($ff->{BONDS}{$t1}{$t2})) {
					$r0 = $ff->{BONDS}{$t1}{$t2}{1}{VALS}[1];
				} elsif(defined($ff) and exists($ff->{BONDS}{$t2}) and exists($ff->{BONDS}{$t2}{$t1})) {
					$r0 = $ff->{BONDS}{$t2}{$t1}{1}{VALS}[1];
				} elsif(defined($ff)) {
					next;
				} else {
					$r0 = $max_dist
				}
				($dist, $k) = General::GetBondLength($a1,$a2,$box,1);
				$rec = ();
				if($dist < $r0*(1+$tol) and $dist > $r0*(1-$tol)) {
					$rec->{i2} = $i2;
					@{ $rec->{d} } = @{ $d };
					$rec->{dist} = $dist;
					$rec->{r0} = $r0;
					#score bond, assuming harmonic functions about equilibrium + random number to prevent overwriting
					$score = 100000*($r0 - $dist)**2+rand(); 
					%{ $candidates->{$score} } = %{ $rec };
				}
			}
		}
		for $c (sort numerically keys %{ $candidates } ) {
			$max_bond = $atoms->{$i1}{MAXBOND};
			last if (exists($a1->{NUMBONDS}) and $a1->{NUMBONDS} == $atoms->{$i1}{MAXBOND});
			$i2 = $candidates->{$c}{i2};
			next if (exists($atoms->{$i2}{NUMBONDS}) and $atoms->{$i2}{NUMBONDS} == $atoms->{$i2}{MAXBOND});
			$d = $candidates->{$c}{d};
			push @{ $bonds->{$i1} }, $i2;
			push @{ $bonds->{$i2} }, $i1;
			if($d->[0] != 0) {
				$atoms->{$i1}{DISPX}[ $#{ $bonds->{$i1} } ] =  $d->[0];
				$atoms->{$i2}{DISPX}[ $#{ $bonds->{$i2} } ] = -$d->[0];
			}
			if($d->[1] != 0) {
				$atoms->{$i1}{DISPY}[ $#{ $bonds->{$i1} } ] =  $d->[1];
				$atoms->{$i2}{DISPY}[ $#{ $bonds->{$i2} } ] = -$d->[1];
			}
			if($d->[2] != 0) {
				$atoms->{$i1}{DISPZ}[ $#{ $bonds->{$i1} } ] =  $d->[2];
				$atoms->{$i2}{DISPZ}[ $#{ $bonds->{$i2} } ] = -$d->[2];
			}
			$a1->{NUMBONDS}++;
			$atoms->{$i2}{NUMBONDS}++;
		}
	}
}

sub CreateBondsByDistanceNoGrid {
	my ($atoms, $bonds, $box, $atomSel, $ff, $max_bond) = @_;
	my ($i, $j, $k, $nbond, $dist, $disp, @tmp, $m, $n); 
	my ($a1, $a2, @tmat, $t1, $t2, $max_dist, $tol);

	$tol = 1.15;
	$max_bond = 9999 if (!defined($max_bond));
	
	#get box displacement tensor
	&GetBoxDisplacementTensor($box);

	#first delete all bonds to selected atoms
	for $i (keys %{ $atoms }) {
		for $j (@{ $bonds->{$i} }) {
			for $k (0 .. $#{ $bonds->{$j} }) {
				splice @{ $bonds->{$j} }, $k, 1 if(exists($atomSel->{ $bonds->{$j}[$k] }));
			}
		}
	}
	for $i (keys %{ $atomSel }) {
		$bonds->{$i} = ();
		for $j ("DISPX", "DISPY", "DISPZ", "ORDER") {
			delete $atoms->{$i}{$j};
		}
	}
	#initialize search vectors
	if (defined($box)) {
		for $i (-1, 0, 1) {
			for $j (-1,0,1) {
				for $k (-1,0,1) {
					push @tmat, [$i,$j,$k];
				}
			}
		}
	} else {
		@tmat = ([0,0,0]);
	}

	#now find new bonds
	@tmp = sort numerically keys %{ $atoms };
	for $m (0 .. $#tmp) {
		$i = $tmp[$m];
		next if (!exists($atomSel->{$i}));
		$a1 = $atoms->{$i};
		$nbond = 0;
		$t1 = $a1->{FFTYPE};
		for $n ($m+1 .. $#tmp) {
			last if ($nbond == $max_bond);
			$j = $tmp[$n];
			next if ($atoms->{$j}{NUMBONDS} == $max_bond);
			$t2 = $atoms->{$j}{FFTYPE};
			if(exists($ff->{BONDS}{$t1}) and exists($ff->{BONDS}{$t1}{$t2})) {
				$max_dist = $ff->{BONDS}{$t1}{$t2}{1}{VALS}[1]*$tol;
			} else {
				$max_dist = $ff->{BONDS}{$t2}{$t1}{1}{VALS}[1]*$tol;
			}
			for $k (@tmat) { 
				$a2 = TransAtom($atoms->{$j}, $box, $k);
				$dist = 0;
				for ("XCOORD", "YCOORD", "ZCOORD") {
					$dist += ($a1->{$_} - $a2->{$_})**2;
				}
				$dist = sqrt($dist);
				if($dist < $max_dist) {
					push @{ $bonds->{$i} }, $j;
					push @{ $bonds->{$j} }, $i;
					if($k->[0] != 0) {
						$atoms->{$i}{DISPX}[ $#{ $bonds->{$i} } ] = $k->[0];
						$atoms->{$j}{DISPX}[ $#{ $bonds->{$j} } ] = -$k->[0];
					}
					if($k->[1] != 0) {
						$atoms->{$i}{DISPY}[ $#{ $bonds->{$i} } ] = $k->[1];
						$atoms->{$j}{DISPY}[ $#{ $bonds->{$j} } ] = -$k->[1];
					}
					if($k->[2] != 0) {
						$atoms->{$i}{DISPZ}[ $#{ $bonds->{$i} } ] = $k->[2];
						$atoms->{$j}{DISPZ}[ $#{ $bonds->{$j} } ] = -$k->[2];
					}
				}
				$nbond++;
				$atoms->{$j}{NUMBONDS}++;
			}
		}
	}
}

sub CreateBondsByDistance_old {
	my ($atoms, $bonds, $box, $atomSel, $max_dist, $max_bond) = @_;
	my ($i, $j, $k, $nbond, $dist, $disp, @tmp, $m, $n);

	#first delete all bonds to selected atoms
	for $i (keys %{ $atoms }) {
		for $j (@{ $bonds->{$i} }) {
			for $k (0 .. $#{ $bonds->{$j} }) {
				splice @{ $bonds->{$j} }, $k, 1 if(exists($atomSel->{ $bonds->{$j}[$k] }));
			}
		}
	}
	for $i (keys %{ $atomSel }) {
		$bonds->{$i} = ();
		for $j ("DISPX", "DISPY", "DISPZ", "ORDER") {
			delete $atoms->{$i}{$j};
		}
	}
	#now find new bonds
	@tmp = sort numerically keys %{ $atomSel };
	for $m (0 .. $#tmp) {
		$i = $tmp[$m];
		$nbond = 0;
		for $n ($m+1 .. $#tmp) {
			#last if ($nbond == $max_bond or (exists($atoms->{$i}{ELEMENT}) and $nbond == $atoms->{$i}{ELEMENT}{NBONDS}));
			last if ($nbond == $max_bond);
			$j = $tmp[$n];
			($dist, $disp) = GetBondLength($atoms->{$i}, $atoms->{$j}, $box, 1);
			if($dist < $max_dist) {
				push @{ $bonds->{$i} }, $j;
				push @{ $bonds->{$j} }, $i;
				for $k ("DISPX", "DISPY", "DISPZ") {
					next if (!exists($disp->{$k}));
					$atoms->{$i}{$k}[ $#{ $bonds->{$i} } ] = $disp->{$k};
					$atoms->{$j}{$k}[ $#{ $bonds->{$j} } ] = -$disp->{$k};
				}
				$nbond++;
			}
		}
		for $j ("DISPX", "DISPY", "DISPZ") {
			next if (! defined($atoms->{$i}{$j}));
			$k = 0;
			for (@{ $atoms->{$i}{$j} }) {
				$k = 1 if($atoms->{$i}{$j}[$_] != 0);
			}
			delete $atoms->{$i}{$j} if (! $k);
		}
	}
}

sub CenterSystem {
	my ($atoms, $offset, $atomList) = @_;
	my ($i, $dim);

	for $i (keys %{ $atoms }) {
		next if (keys %{ $atomList } && ! exists($atomList->{$i}));
		for $dim ("XCOORD", "YCOORD", "ZCOORD") {
			$offset->{$dim} = $offset->{substr($dim,0,1)}{lo} if (! exists($offset->{$dim}));
			next if (! $offset->{$dim});
			$atoms->{$i}{$dim} -= $offset->{$dim};
		}
	}
}

sub FindElement_old {
	my ($fftype, $parms, $elements) = @_;
	my ($elementName, $elementNum, $i);

	#die "ERROR: $fftype is not a valid force field type!\n" if (! exists($parms->{$fftype}));
	$elementName = $parms->{$fftype}{ATOM};
	$elementNum = 0;

	for $i (keys %{ $elements }) {
		if (uc($elements->{$i}{SYMBOL}) eq uc($elementName)) {
			$elementNum = $i;
			last;
		}
	}

	#die "ERROR: Element $elementName is not a valid element!\n" if (! $elementNum);

	return ($elementName, $elementNum);
}

sub ImageAtoms_new {
	my ($atoms, $CoM, $box) = @_;
	my ($i, $j, $CENTER, $boxAtm, $OFFSET);

	for $i ("XCOORD", "YCOORD", "ZCOORD") {
		$boxAtm->{$i} = $box->{$i}{CENTER};
	}

	&GetMinDist($CoM, $boxAtm, $box, \%{ $CENTER});
	for $i ("X", "Y", "Z") {
		$OFFSET->{"${i}COORD"} = ($CoM->{"${i}COORD"} - $CENTER->{$i});
	}

	for $i (keys %{ $atoms }) {
		for $j ("XCOORD", "YCOORD", "ZCOORD") {
			$atoms->{$i}{$j} -= $CoM->{$j};
		}
		$atoms->{$i}{IMAGED} = 1;
	}

	return $CENTER;
}

sub ImageAtoms {
	my ($atoms, $CoM, $box) = @_;
	my ($MOL, $i, $OFFSET, $dim, $index, $CENTER);

	%{ $CENTER } = %{ $CoM };
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
		$OFFSET->{$dim} = 0;
		while ($CENTER->{$dim} > $box->{$dim}{hi}) {
			$OFFSET->{$dim}++;
			$CENTER->{$dim} -= $box->{$dim}{len};
		}
		while ($CENTER->{$dim} < $box->{$dim}{lo}) {
			$OFFSET->{$dim}--;
			$CENTER->{$dim} += $box->{$dim}{len};
		}
	}
	
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
		$index = $dim;
		$index =~ s/COORD/INDEX/;
		for $i (keys %{ $atoms }) {
			$atoms->{$i}{$dim} -= ($OFFSET->{$dim} * $box->{$dim}{len});
			$atoms->{$i}{$index} = $OFFSET->{$dim};
			if ($atoms->{$i}{$dim} < $box->{$dim}{lo}) {
				$atoms->{$i}{$index}--;
			} elsif ($atoms->{$i}{$dim} > $box->{$dim}{hi}) {
				$atoms->{$i}{$index}++;
			}
			$atoms->{$i}{IMAGED} = 1;
		}
	}

	return $CENTER;
}

sub ScaleAtoms {
	my ($atoms, $box) = @_;
	my ($i, $coord, $dim, $index, $pos);

	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
		$index = $dim;
		$index =~ s/COORD/INDEX/;
		for $i (keys %{ $atoms }) {
			$pos = $atoms->{$i}{$dim};
			$pos /= $box->{$dim}{len};
			if ($pos =~ /(\d+)(\.\d+)/) {
				if ($1 > 0) {
					$atoms->{$i}{$index} = $1 + 1;
				} else {
					$atoms->{$i}{$index} = $1;
				}
				$atoms->{$i}{$dim} = sprintf("%.6f ", $pos);
			}
		}
	}
}

sub UnwrapAtoms {
	my ($ATOMS, $BOX, $scaled) = @_;
	my ($atomC, $atom, $i, $index, $pos);

	for $atomC (keys %{ $ATOMS }) {
		$atom = \%{ $ATOMS->{$atomC} };
		for $i ("X", "Y", "Z") {
			$index = "${i}INDEX";
			next if(!exists($atom->{$index}));
			$index = $atom->{$index};
			if (! defined($scaled) or $scaled == 1) {
				$atom->{"${i}COORD"} *= $BOX->{"${i}COORD"}{len};
				$atom->{"${i}COORD"} += $BOX->{"${i}COORD"}{lo};
			}
			$atom->{"${i}COORD"} += ($index * $BOX->{"${i}COORD"}{len});
		}
	}
}

sub GetAtmList {
	my ($select, $ATOMS) = @_;
	my ($i, %LIST, $field, $val, $rec, $operator, $excluded, $tmp);
	my ($count);

	for $field (keys %{ $select }) {
		$tmp = "";
		for $val (keys %{ $select->{$field} }) {
			$excluded = 0;
			$operator = "";
			if ($val eq "*") {
				$operator = ">";
				$val = "0 and \$ATOMS->{\$i}{INDEX} < 999999";
			} elsif ($val =~ /^\d+/) { #integer
				$operator = "==";
			} elsif ($val =~ /^\!(\d+)/) {
				$operator = "!=";
				$val = $1;
				$excluded = 1;
			} elsif ($val =~ /^(>|<|=|\+|\-)(\d+\.*\d*)/) {
				$operator = $1;
				$val = $2;
			} elsif ($val =~ /^\!(>|<|=)(\d+\.*\d*)/) {
				$excluded = 1;
				if ($1 eq ">") {
					$operator = "<";
				} elsif ($1 eq "<") {
					$operator = ">";
				} else {
					$operator = "!=";
				}
				$val = $2;
			} elsif ($val =~ /^(\w+)$/) {
				$operator = "eq";
				$val = "\"${1}\"";
			} elsif ($val =~ /^\!(\w+)$/) {
				$excluded = 1;
				$operator = "ne";
				$val = "\"${1}\"";
			} elsif ($val =~ /^(\S+)/) {
				$operator = "=~";
				$val = "/$1/";
			}
			next if (! $operator);
			if ($tmp) {
				if (! $excluded) {
					$tmp .= "or \$ATOMS->{\$i}{${field}} $operator $val ";
				} else {
					$tmp .= "and \$ATOMS->{\$i}{${field}} $operator $val ";
				}
			} else {
				$tmp = "{${field}} $operator $val ";
			}
		}
		if ($rec) {
			$rec .= "and (\$ATOMS->{\$i}$tmp) ";
		} else {
			$rec = "$tmp";
		}
	}

	$count = 0;
	for $i (keys %{ $ATOMS }) {
		$excluded = 0;
		if (! eval('$ATOMS->{$i}' . $rec)) {
			$excluded = 1;
		}
		if (! $excluded) {
			$count++;
			$LIST{$i} = $count;
		}
	}
	
	return \%LIST;
}

sub GetSolvent {
	my ($atoms, $solvType) = @_;
	my (%SOLVOPTS, %SOLVATMS, $i);

	%SOLVOPTS = (
						"WATER" => {
										"MOLSIZE" => 3,
										"RESNAME" => "WAT|HOH|RES",
										"FFTYPE"  => "OW|HW|OT|HT|HO|OH|H_|O_3"
									}
				);

	return () if (! defined($SOLVOPTS{uc($solvType)}));

	for $i (keys %{ $atoms }) {
		if ($atoms->{$i}{MOLSIZE} == $SOLVOPTS{$solvType}{MOLSIZE}) {
			if ($atoms->{$i}{RESNAME} =~ /$SOLVOPTS{$solvType}{RESNAME}/ or
				$atoms->{$i}{FFTYPE} =~ /$SOLVOPTS{$solvType}{FFTYPE}/) {
					$atoms->{$i}{IS_SOLVENT} = 1;
					$SOLVATMS{$i} = 1;
			}
		}
	}

	return \%SOLVATMS;
}

sub getMolList {
	my ($atoms, $bonds, $atomID, $mol, $used) = @_;
	my ($i, $j, $k);

	$mol->{MEMBERS}{$atomID} = 1;
	$mol->{MOLSIZE}++;
	for $i ("MOLECULE", "MOLECULEID", "MOLSIZE") {
		delete $atoms->{$atomID}{$i};
	}
	$atoms->{$atomID}{MOLECULE} = $mol;
	$atoms->{$atomID}{MOLECULEID} = \$atoms->{$atomID}{MOLECULE}{INDEX};
	$atoms->{$atomID}{MOLSIZE} = \$atoms->{$atomID}{MOLECULE}{MOLSIZE};
	$used->{$atomID} = 1;
	for $i (@{ $bonds->{$atomID} }) {
		next if (exists($used->{$i}));
		$used->{$i} = 1;
		#add bond data
		for $j (keys %{ $atoms->{$atomID} }) {
			next if ($j =~ /(MOL|DISP|ORD|BONDATOM)/);
			push @{ $atoms->{$atomID}{BONDATOM}{uc $j}{$atoms->{$i}{$j}} }, \%{ $atoms->{$i} } if(exists($atoms->{$i}{$j}));
			push @{ $atoms->{$i}{BONDATOM}{uc $j}{$atoms->{$atomID}{$j}} }, \%{ $atoms->{$atomID} };
		}
		&getMolList($atoms, $bonds, $i, $mol, $used);
	 }
}

sub GetMols {
	my ($atoms, $bonds, $select) = @_;
	my ($i, $j, $MOLS, $counter, $msize, $USED);

	$select = $atoms if (! defined($select));
	$i = 1;
	while (! exists($atoms->{$i})) {
		$i++;
	}

	$counter = 0;
	for $i (sort numerically keys %{ $atoms }) {
		if (! exists($USED->{$i}) and exists($select->{$i})) {
			$counter++;
			$MOLS->{$counter} = ();
			$MOLS->{$counter}{INDEX} = $counter;
			$MOLS->{$counter}{MOLSIZE} = 0;
			&getMolList($atoms, $bonds, $i, $MOLS->{$counter}, $USED);
		}
	}

	return $MOLS;
}

sub GetMols_old {
	my ($ATOMS, $BONDS) = @_;
	my ($atomC, @tmp, $rec, $i, $min, $molNum);

	for $atomC (keys %{ $ATOMS }) {
		delete $ATOMS->{$atomC}{MOLECULE};
		delete $ATOMS->{$atomC}{MOLECULEID};
	}
	
	@tmp = sort numerically keys %{ $ATOMS };
	for $atomC (@tmp) {
		$rec = ();
		$rec->{MEMBERS}{$atomC} = 1;
		if ($ATOMS->{$atomC}{NUMBONDS} == 0 || ! $BONDS->{$atomC}) { #ions
			$rec->{INDEX} = $molNum;
			$ATOMS->{$atomC}{MOLECULE} = \%{ $rec };
			$molNum++;
		} else {
			$min = $atomC;
			for $i (@{ $BONDS->{$atomC} })  {
				$rec->{MEMBERS}{$i} = 1;
				$min = $i if ($i < $min);
			}
			if ($min < $atomC) { # found a member which has a lower index, so merge this rec with it, and update
				for $i (keys %{ $rec->{MEMBERS} }) {
					$ATOMS->{$min}{MOLECULE}{MEMBERS}{$i} = 1;
					next if ($min == $i or $i == $atomC);
					if (exists($ATOMS->{$i}{MOLECULE}{MEMBERS})) {
						for (keys %{ $ATOMS->{$i}{MOLECULE}{MEMBERS} }) {
							$ATOMS->{$min}{MOLECULE}{MEMBERS}{$_} = 1;
						}
					}
				}
				$rec = \%{ $ATOMS->{$min}{MOLECULE} };
			} else {
				$rec->{INDEX} = $molNum;
				$molNum++;
			}

			for $i (keys %{ $rec->{MEMBERS} })  {
				next if ($i == $min);
				$ATOMS->{$i}{MOLECULE} = \%{ $rec };
			}
		}
	}

	$molNum = 1;
	for $atomC (@tmp) {
		if (! exists($ATOMS->{$atomC}{MOLECULE}{INDEX}) || $ATOMS->{$atomC}{NUMBONDS} == 0 || ! $BONDS->{$atomC}) {
			$ATOMS->{$atomC}{MOLECULE}{INDEX} = $molNum;
			$molNum++;
		}
		$ATOMS->{$atomC}{MOLECULEID} = $ATOMS->{$atomC}{MOLECULE}{INDEX};
		$ATOMS->{$atomC}{MOLECULE}{SIZE} = scalar(keys %{ $ATOMS->{$atomC}{MOLECULE}{MEMBERS} });
		$ATOMS->{$atomC}{MOLSIZE} = $ATOMS->{$atomC}{MOLECULE}{SIZE};
	}
}

sub SplitAtomsByMol {
	my ($atoms, $selList) = @_;
	my ($i, %molList, $counter, $j, @tmp);

	$selList = $atoms if (! defined($selList));
	@tmp = sort { ($a<=>$b) } keys %{ $selList };

	for $i (@tmp) {
		next if (exists($atoms->{$i}{IS_SPLIT}));
		for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
			next if (! exists($selList->{$j}));
			$atoms->{$j}{IS_SPLIT} = 1;
			$counter = $atoms->{$j}{MOLECULEID};
			$molList{$counter}{$j} = $atoms->{$i}{MOLECULE}{MEMBERS}{$j};
		}
	}

	for $i (@tmp) {
		delete $atoms->{$i}{IS_SPLIT}
	}

	return \%molList;
}

sub GetAtmData {
	my ($allAtoms, $atomList) = @_;
	my (%ATOMS, $i);
																																	  
	for $i (keys %{ $atomList }) {
		$ATOMS{$i} = \%{ $allAtoms->{$i} };
	}
	return \%ATOMS;
}

sub vecLen {
	my ($v) = $_[0];
	my ($vLen);

	for ("XCOORD", "YCOORD", "ZCOORD") {
		$vLen += $v->{$_}*$v->{$_};
	}
	return sqrt($vLen);
}

sub GetAbituraryRotationMatrix {
	my ($v, $s, $c) = @_;
	my ($rM,$vx,$vy,$vz,$scale,$tmp);
	$tmp = vecLen($v);
	$scale = 1/$tmp;
	($vx,$vy,$vz) = ($v->{XCOORD}*$scale,$v->{YCOORD}*$scale,$v->{ZCOORD}*$scale);
	$rM = [
		   [ 1+(1-$c)*($vx*$vx-1), (1-$c)*$vx*$vy-$vz*$s, (1-$c)*$vx*$vz+$vy*$s],
		   [(1-$c)*$vx*$vy-$vz*$s,  1+(1-$c)*($vy*$vy-1), (1-$c)*$vy*$vz-$vx*$s],
		   [(1-$c)*$vx*$vz-$vy*$s, (1-$c)*$vy*$vz+$vx*$s,  1+(1-$c)*($vz*$vz-1)]
		  ];
	return $rM;
}

sub RotateAbout {
	my ($atoms, $rotM) = @_;
	my ($orig, $i, $j, $count);

	for $i (keys %{ $atoms }) {
		for $j ("XCOORD","YCOORD","ZCOORD") {
			$orig->{$i}{$j} = $atoms->{$i}{$j};
		}
	}

	for $i (keys %{ $atoms }) {
		$atoms->{$i}{XCOORD} =  $orig->{$i}{XCOORD}*$rotM->[0][0]+
								$orig->{$i}{YCOORD}*$rotM->[0][1]+
								$orig->{$i}{ZCOORD}*$rotM->[0][2];
		$atoms->{$i}{YCOORD} =  $orig->{$i}{XCOORD}*$rotM->[1][0]+
								$orig->{$i}{YCOORD}*$rotM->[1][1]+
								$orig->{$i}{ZCOORD}*$rotM->[1][2];
		$atoms->{$i}{ZCOORD} =  $orig->{$i}{XCOORD}*$rotM->[2][0]+
								$orig->{$i}{YCOORD}*$rotM->[2][1]+
								$orig->{$i}{ZCOORD}*$rotM->[2][2];
	}
}

sub TransAtoms {
	my ($atoms, $com, $d) = @_;
	my ($i, $j);

	$d = 1 if (! defined($d));
	for $i ("XCOORD", "YCOORD", "ZCOORD") {
		for $j (keys %{ $atoms }) {
			$atoms->{$j}{$i} += $d*$com->{$i};
		}
	}
}

sub updateAtomIndex {
	my ($atoms, $bonds, $select) = @_;
	my ($i, $j, $index, $newAtoms, $newBonds, @fields, @list);

	$newBonds = ();
	@list = keys %{ $atoms };
	@fields = grep {!/MOL/} keys %{ $atoms->{ $list[0] } };
	for $i (@list) {
		next if (! $select->{$i});
		if (! keys %{ $atoms->{$i} } ) {
			delete $atoms->{$i};
			next;
		}
		$index = $atoms->{$i}{INDEX};
		for $j (@fields) {
			$newAtoms->{$index}{$j} = $atoms->{$i}{$j};
		}
		$newBonds->{$index} = ();
		next if (! exists($newAtoms->{$index}{BONDS}) || ! $newAtoms->{$index}{BONDS});
		for $j (0 .. $#{ $newAtoms->{$index}{BONDS} }) {
			$newBonds->{$index}[$j] = $newAtoms->{$index}{BONDS}[$j]{INDEX};
		}
		delete $newAtoms->{$index}{BONDS};
		delete $atoms->{$i};
		delete $bonds->{$i};
	}
	@list = keys %{ $atoms };
	for $i (@list) {
		next if (! $select->{$i});
		delete $atoms->{$i};
		delete $bonds->{$i};
	}
	for $i (keys %{ $newAtoms }) {
		%{ $atoms->{$i} } = %{ $newAtoms->{$i} };
		$bonds->{$i} = ();
		@{ $bonds->{$i}  } = @{ $newBonds->{$i} } if($newBonds->{$i});
	}
}

sub SelectAtomsByField {
	my ($atoms, $bonds, $mode, $select) = @_;
	my ($i, $FIELD, $MOLS, $hasSelect);

	$hasSelect = 1;
	if(! defined($select)) {
		$hasSelect = 0;
		$select = $atoms;
	}
	
	$MOLS = &GetMols($atoms, $bonds, $select) if ($mode =~ /MOL/);

	for $i (keys %{ $atoms }) {
		next if (! exists($select->{$i}));
		if ($mode =~ /MOL/) {
			$FIELD->{ ${ $atoms->{$i}{$mode} } }{$i} = $i;
		} else {
			$FIELD->{ $atoms->{$i}{$mode} }{$i} = $i
		}
	}

	return $FIELD;
}

sub GroupAtomsByField {
	my ($atoms, $bonds, $mode, $select, $mode2, $reverse_opt) = @_;
	my ($i, $j, $FIELD, $index, @sorted, $sort_val, $sort_numeric); 
	my ($resid, $resNum, $MOLS, $tmp, $hasSelect);

	$reverse_opt = 0 if (! defined($reverse_opt));
	$hasSelect = 1;
	if(! defined($select)) {
		$hasSelect = 0;
		$select = $atoms;
	}
	$mode2 = "INDEX" if (! defined($mode2));

	$MOLS = &GetMols($atoms, $bonds, $select) 
		if ($mode =~ /MOL/ or $mode2 =~ /MOL/);

	for $i (keys %{ $atoms }) {
		next if (! exists($select->{$i}));
		next if (! keys %{ $atoms->{$i}});
		#store the hash key values based on {mode}{mode2}
		$sort_val = $atoms->{$i}{$mode2};
		$sort_val = ${ $atoms->{$i}{$mode2} } if ($mode2 =~ /MOL/);
		if ($mode =~ /MOL/) {
			$FIELD->{ ${ $atoms->{$i}{$mode} } }{$i} = $sort_val;
		} else {
			$FIELD->{ $atoms->{$i}{$mode} }{$i} = $sort_val;
		}
		#now store the original bond indices
		next if (! exists($bonds->{$i}) || ! $bonds->{$i});
		$index = 0;
		for $j (@{ $bonds->{$i} }) {
			$atoms->{$i}{BONDS}[$index] = \%{ $atoms->{$j} };
			$index++;
		}
	}

	if ($mode =~ /MOL|RESNUM|NUMBONDS|COORD/) {
		@sorted = sort {$a <=> $b} keys %{ $FIELD };
	} else {
		@sorted = sort {$a cmp $b } keys %{ $FIELD };
	}
	$sort_numeric = 0;
	$sort_numeric = 1 if ($mode2 =~ /MOL|RESNUM|NUMBONDS|COORD/);

	$index = $resid = $resNum = 0;
	if($hasSelect) {
		@{ $tmp } = sort {$a<=>$b} keys %{ $FIELD->{$sorted[0]} };
		$index = $atoms->{$tmp->[0]}{INDEX}-1;
		$resid = $resNum = $atoms->{$tmp->[0]}{RESNUM};
	}
	for $i (0 .. $#sorted) {
		if($sort_numeric) {
			@{ $tmp } = sort { $FIELD->{$sorted[$i]}{$a} <=> $FIELD->{$sorted[$i]}{$b} } keys %{ $FIELD->{ $sorted[$i]} };
		} else {
			@{ $tmp } = sort { $FIELD->{$sorted[$i]}{$a} cmp $FIELD->{$sorted[$i]}{$b} } keys %{ $FIELD->{ $sorted[$i]} };
		}
		@{ $tmp } = reverse @{ $tmp } if ($reverse_opt);
		for $j (@{ $tmp }) {
			$index++;
			$atoms->{$j}{INDEX} = $index;
			if (! $resNum || $atoms->{$j}{RESNUM} != $resNum) {
				$resid++;
				$resNum = $atoms->{$j}{RESNUM};
			}
			$atoms->{$j}{RESNUM} = $resid;
			if ($mode =~ /MOLECULEID/) {
				$atoms->{$j}{CHAIN} = "X";
				$atoms->{$j}{CHAIN} = chr(64+ ${ $atoms->{$j}{MOLECULEID} }) 
					if($atoms->{$j}{MOLECULEID} < 10);
			} elsif ($mode =~ /MOLSIZE/) {
				$atoms->{$j}{CHAIN} = chr(64+ $i+1);
			}
		}
	}
	$FIELD = ();
	&updateAtomIndex(\%{ $atoms },\%{ $bonds }, $select);
}

sub MakeFieldSequential {
	my ($atoms, $field) = @_;
	my ($i, $atomsByField, @list, $isnumeric, $index, $j);

	@list = keys %{ $atoms };
	$isnumeric = 1;
	die "ERROR: field $field does not exists in MakeFieldSequential\n"
		if (! exists($atoms->{$list[0]}{$field}));
	for $i (@list ) {
		push @{ $atomsByField->{ $atoms->{$i}{$field} } }, $atoms->{$i};
		$isnumeric=0 if ($atoms->{$i}{$field} !~ /^\-?\d+\.?\d*$/);
	}
	return 0 if (scalar(keys %{ $atomsByField} ) == 1);
	@list = sort numerically keys %{ $atomsByField } if($isnumeric);
	@list = sort { ($a cmp $b); } keys %{ $atomsByField } if(!$isnumeric);
	$index = shift @list;
	$index++;
	for $i (@list) {
		for $j (@{ $atomsByField->{$i} }) {
			$j->{$field} = $index;
		}
		$index++;
	}
}

sub BuildAtomSelectionString {
	my ($atomSel) = $_[0];
	my ($fields, $fl, $ud, $sel, $atm_no);
	my ($dist, $x, $y, $z, $oper, $sstr);

	$fields = (
				{
				"INDEX"		=> 0,
				"ATMNAME"	=> 1,
				"RESNAME"	=> 1,
				"CHAIN"		=> 1,
				"RESNUM"	=> 0,
				"XCOORD"	=> 0,
				"YCOORD"	=> 0,
				"ZCOORD"	=> 0,
				"FFTYPE"	=> 1,
				"NUMBONDS"	=> 0,
				"LONEPAIRS"	=> 0,
				"CHARGE"	=> 0,
				"MOLECULEID"=> 0,
				"MOLSIZE"	=> 0,
				"BONDATOM"	=> 0,
				"fa"		=> 0,
				"fb"		=> 0,
				"fc"		=> 0,
				}
			);
	$fl = "(\-?\-?\>?)(" . join("|",keys %{ $fields }) . ")";
	$sel = $atomSel;
	#replace all the field strings with $atoms->{$i}{FIELD}
	while ($atomSel =~ /$fl/gi) {
		if ($1 eq "") {
			$ud = '$atoms->{$i}{' . uc  $2 . '}';
			$ud = '${ $atoms->{$i}{' . uc  $2 . '} }' if (uc($2) eq "MOLSIZE" or uc($2) eq "MOLECULEID");
			$sel =~ s/$2/$ud/;
		} else {
			$ud = "{".  uc $2 . "}";
			$sel =~ s/${1}${2}/$ud/;
		}
	}
	#special case for if the dist keyword is used
	while ($sel =~ /dist\((\d+)\)/i) { 
		$atm_no = $1;
		$sel =~ s/dist\(\w+\)/GetBondLength\(\$atoms->\{\$i\}\,\$atoms->\{$atm_no\},\$box\)/i;
	}
	#special case for if the count keyword is used
	while ($sel =~ /count\((\S+)\)/i) {
		$atm_no = $1;
		$sel =~ s/count\((\S+)\)/scalar \@{ $atm_no }/i;
	}
	#special case if the with(in/out) keyword is used. syntax is with(in/out)(distance of atom_index) or with(in/out)(distance of x,y,z)
	while ($sel =~ /with(in|out)\s*\((\d+\.?\d*)\s+of\s+(\d+\.?\d*),(\d+\.?\d*),(\d+\.?\d*)\)/i) {
		$dist = $2;
		$x = $3; $y = $4; $z= $5;
		$oper = "<";
		$oper = ">" if($1 =~ /out/);
		$sel =~ s/with(in|out)\s*\((\d+\.?\d*)\s+of\s+(\d+\.?\d*),(\d+\.?\d*),(\d+\.?\d*)\)/GetBondLength\(\$atoms->\{\$i\}\,\"$x,$y,$z\"\, \$box\) ${oper}= $dist/i;
	}
	while ($sel =~ /with(in|out)\s*\((\d+\.?\d*)\s+of\s+(\d+)\)/i) {
		$dist = $2;
		$atm_no = $3;
		$oper = "<";
		$oper = ">" if($1 =~ /out/);
		$sel =~ s/with(in|out)\s*\((\d+\.?\d*)\s+of\s+(\d+)\)/GetBondLength\(\$atoms->\{\$i\}\,\$atoms->\{$atm_no\},\$box\) ${oper}= $dist/i;
	}
	while ($sel =~ /with(in|out)\s*\((\d+\.?\d*)\s+of\s+([^\)]+)\)/i) {
		$dist = $2;
		$sstr = $3;
		$oper = "<";
		$oper = ">" if($1 =~ /out/);
		$sel =~ s/with(in|out)\s*\((\d+\.?\d*)\s+of\s+([^\)]+)\)/SearchAroundAtom\('$sstr',\$box\,\$atoms,\"$oper\",$dist,\\\%{\$atomList})/i;
	}
	return $sel;
}

sub SelectAtoms {
	my ($selectionStr, $atoms, $box, $atomList) = @_;
	my ($i, $err, $sList, @tmp);

	#special case for SearchAroundAtom
	while($selectionStr =~ /(SearchAroundAtom\([^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+\))/gi) {
		eval($1);
	}

	$selectionStr =~ s/or SearchAroundAtom\([^\)]+\)//g;
	#now try to see if remaining expresion is valid
	$i = 1;
	#try
	$SIG{__WARN__} = sub {  };
	eval($selectionStr);
	if($@) {
		$err = $@;
		$err =~ s/^.* line \d+, near .../near \"/;
		die "ERROR: Invalid atom selection ${err}";
	}
	#now select atoms
	for $i (keys %{ $atoms }) {
		if (eval($selectionStr)) {
			$atomList->{$i} = 1;
		}
	}
	# make sure at least 1 atom record selection
	die "ERROR: No atoms matched selection!\n" if (! $atomList);
	return $atomList;
}

sub ReimageAtoms {
	my ($atoms, $bonds, $mols, $box, $selection) = @_;
	my ($i, $j, $k, $l, $dist, $factor, $altered, $com, $curr);

	$altered = ();

	for $i (keys %{ $mols }) {
		for $j (keys %{ $mols->{$i}{MEMBERS} }) {
			next if(!exists($selection->{$j}));
			for $k (@{ $bonds->{$j} }) {
				next if($j > $k);
				for $l ("X", "Y", "Z") {
					$dist = getdist($atoms->{$j}, $atoms->{$k}, $l);
					$factor = 1;
					$factor = -1 if ($atoms->{$k}{"${l}COORD"} < $atoms->{$j}{"${l}COORD"});
					while($dist>4) {
						$dist -= $box->{$l}{len};
						$atoms->{$k}{"${l}COORD"} -= $factor * $box->{$l}{len};
					}
				}
			}
		}
		$curr = GetAtmData($atoms, $mols->{$i}{MEMBERS});
		$com = General::CoM($curr);
		for $l ("X", "Y", "Z") {
			next if ($com->{"${l}COORD"} > 0 and $com->{"${l}COORD"} < $box->{$l});
			$factor = 0;
			if ($com->{"${l}COORD"} < 0) {
				while($com->{"${l}COORD"}<0) {
					$factor += $box->{$l}{len};
					$com->{"${l}COORD"} += $box->{$l}{len};
				}
			} else {
				while($com->{"${l}COORD"}>$box->{$l}{len}) {
					$factor -= $box->{$l}{len};
					$com->{"${l}COORD"} -= $box->{$l}{len};
				}
			}
			for $j (keys %{ $mols->{$i}{MEMBERS} }) {
				$atoms->{$j}{"${l}COORD"} += $factor;
			}
		}
	}
}

sub getdist {
	my ($atom1, $atom2, $dim) = @_;

	return (sqrt(($atom1->{"${dim}COORD"}-$atom2->{"${dim}COORD"})**2));
}

sub AddMolsToSelection {
	my ($select, $atoms) = @_;
	my ($i, $j, @tmp);

	@tmp = keys %{ $select };
	for $i (@tmp) {
		for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
			$select->{$j} = 1;
		}
	}
}

1;

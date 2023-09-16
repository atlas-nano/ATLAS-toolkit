package REPLICATE;

require Exporter;
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use Cwd;
use Storable qw(dclone);
use FileFormats qw(ParseStructFile);

#use General qw(PrintProgress CombineMols GetTime);
#use Superimpose qw(SuperimposeAtoms);

use constant PI => atan2(1,1) * 4;

our (@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(ReplicateCell GetBoxVol GetBoxDisplacementTensor TransAtom TransCellAtoms);
$VERSION = "1.00";

sub TransAtom {
	my ($atom, $box, $vec) = @_;
	my ($newAtm, @mat, $i);

	@mat = ([$box->{X}{DISP_V}{X},$box->{Y}{DISP_V}{X},$box->{Z}{DISP_V}{X}],
			[$box->{X}{DISP_V}{Y},$box->{Y}{DISP_V}{Y},$box->{Z}{DISP_V}{Y}],
			[$box->{X}{DISP_V}{Z},$box->{Y}{DISP_V}{Z},$box->{Z}{DISP_V}{Z}]);
	$newAtm->{XCOORD} = $atom->{XCOORD} + $mat[0][0]*$vec->[0]+$mat[0][1]*$vec->[1]+$mat[0][2]*$vec->[2];
	$newAtm->{YCOORD} = $atom->{YCOORD} + $mat[1][0]*$vec->[0]+$mat[1][1]*$vec->[1]+$mat[1][2]*$vec->[2];
	$newAtm->{ZCOORD} = $atom->{ZCOORD} + $mat[2][0]*$vec->[0]+$mat[2][1]*$vec->[1]+$mat[2][2]*$vec->[2];
	return $newAtm;
}

sub TransCellAtomsOld {
	my ($unit, $box, $vec) = @_;
	my ($i,%ATOMS, $j, $tmp);
	
	%ATOMS = %{ $unit };
	for $i (keys %{ $unit }) {
		$tmp = TransAtom($unit->{$i},$box, $vec);
		for $j ("XCOORD","YCOORD","ZCOORD") {
			$ATOMS{$i}{$j} = $tmp->{$j};
		}
	}
	return \%ATOMS;
}

sub TransCellAtoms {
	my ($unit, $box, $vec) = @_;
	my ($i, $ATOMS, $j, $tmp);
	
	$ATOMS = dclone($unit);
	for $i (keys %{ $ATOMS }) {
		$tmp = TransAtom($ATOMS->{$i}, $box, $vec);
		for $j ("XCOORD","YCOORD","ZCOORD") {
			$ATOMS->{$i}{$j} = $tmp->{$j};
		}
	}
	return $ATOMS;
}

sub invertMol {
	my ($atoms, $dim) = @_;
	my ($i);

	for $i (keys %{ $atoms }) {
		$atoms->{$i}{"${dim}COORD"} *= -1;
	}
}

sub GetBoxVol {
		my ($box) = $_[0];

		&GetBoxDisplacementTensor($box);

}

sub GetBoxDisplacementTensor {
	my ($box) = $_[0];
	return if (exists($box->{X}{DISP_V}));

	my ($lx,$ly,$lz,$a,$b,$c,$cos_alpha,$alpha,$beta,$gamma);
	my ($ax,$bx,$cx,$by,$cy,$cz);

	$lx = $box->{X}{len};
	$ly = $box->{Y}{len};
	$lz = $box->{Z}{len};
	$lx = 1 if ($lx == 0);
	$ly = 1 if ($ly == 0);
	$lz = 1 if ($lz == 0);
	$alpha = $box->{X}{angle}*PI/180;
	$beta  = $box->{Y}{angle}*PI/180;
	$gamma = $box->{Z}{angle}*PI/180;

	# 3 x 3 displacement tensor
	#(a b c) =  (ax bx cx)
	#				(0  by cy)
	#				(0  0  cz)
	#ax = lx
	#bx = ly cos(gamma)
	#cx = lz cos(beta)
	#cy = ly*lz*cos(alpha) - bx*cx/by
	#cz = sqrt(lz*lz - cx*cx - cy*cy)
	$ax = $lx; $bx = $ly*cos($gamma); $cx = $lz*cos($beta);
	$by = $ly*sin($gamma); $cy = ($ly*$lz*cos($alpha)-$bx*$cx)/$by;
	$cz = sqrt($lz*$lz - $cx*$cx - $cy*$cy);
	$box->{X}{DISP_V}{X} = $ax; $box->{X}{DISP_V}{Y} = 0;   $box->{X}{DISP_V}{Z} = 0;
	$box->{Y}{DISP_V}{X} = $bx; $box->{Y}{DISP_V}{Y} = $by; $box->{Y}{DISP_V}{Z} = 0;
	$box->{Z}{DISP_V}{X} = $cx; $box->{Z}{DISP_V}{Y} = $cy; $box->{Z}{DISP_V}{Z} = $cz;
	#update box center
	$box->{X}{center} = ($ax+$bx+$cx)/2;
	$box->{Y}{center} = ($by+$cy)/2;
	$box->{Z}{center} = $cz/2;
}

sub ReplicateCell {
	my ($atoms, $bonds, $box, $dims, $centerMol, $str, $updateResNum, $offset1D, $qAdjust, $setPBC) = @_;
	my ($i, $j, $k, $l, $m, $cellAtoms, $cellBonds, $unitAtoms, $pbcbonds, $nAtoms); 
	my ($strLen, $offset, $tot, $count, $total, $start, $vec, $outdir);
	my ($outAtomFile, $outBondFile, $f, $chain);
	
	$qAdjust->{X} = $qAdjust->{Y} = $qAdjust->{Z} = 0 
		if (! defined($qAdjust));
	$setPBC = 1 if (! defined($setPBC));	
	&GetBoxDisplacementTensor($box);
	&setPBCbonds($atoms,$bonds,$box) if ($setPBC);
	print "${str}calculating time remaining\r";
	$start = time();
	if ($centerMol) {
		@{ $vec } = (int(($dims->{X}-1)/2),int(($dims->{Y}-1)/2),int(($dims->{Z}-1)/2));
		$atoms = TransCellAtoms($atoms, $box, $vec);
	}

	$total = 1;
	for $i ("X", "Y", "Z") {
		$total *= $dims->{$i};
	}

	$unitAtoms = dclone($atoms);
	$cellBonds = dclone($bonds);

	$outdir = "_tmpdir_$$";
	mkdir $outdir;
	$outAtomFile = "${outdir}/__atoms.dat";
	$outBondFile = "${outdir}/__bonds.dat";
	open  OUTATOM,"> $outAtomFile" or die "ERROR: Cannot create $outAtomFile: $!\n";
	print OUTATOM "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,f10.5)\n";
	open OUTBOND, "> $outBondFile" or die "ERROR: Cannot create $outBondFile: $!\n";
	print OUTBOND "FORMAT CONECT (a6,12i6)\n";

	$nAtoms = scalar(keys %{ $atoms });

	$count = 0;
	$offset = 0;
	for $i (0 .. ($dims->{X}-1)) {
		$vec->[0] = $i;
		for $j (0 .. ($dims->{Y}-1)) {
			$vec->[1] = $j;
			for $k (0 .. ($dims->{Z}-1)) {
				$vec->[2] = $k;
				$count++;
				$cellAtoms = TransCellAtoms($unitAtoms, $box, $vec);
				for $l (1 .. $nAtoms) {
					$chain = "X";
					$chain = $atoms->{$l}{CHAIN} if (exists($atoms->{$l}{CHAIN}));
					printf OUTATOM "%-6s %5d %-5s %3s %1s %5d%10.5f%10.5f%10.5f %-5s %2d %1d %9.6f\n",
					$atoms->{$l}{LABEL}, $l + $offset, $atoms->{$l}{ATMNAME}, $atoms->{$l}{RESNAME}, $chain, $atoms->{$l}{RESNUM},
					$cellAtoms->{$l}{XCOORD}, $cellAtoms->{$l}{YCOORD}, $cellAtoms->{$l}{ZCOORD}, 
					$atoms->{$l}{FFTYPE}, $atoms->{$l}{NUMBONDS}, $atoms->{$l}{LONEPAIRS}, $atoms->{$l}{CHARGE};
					printf OUTBOND "%-6s%6d","CONECT", $l + $offset;
					for $m (@{ $bonds->{$l} }) {
						 printf OUTBOND "%6d", $m + $offset;
					}
					printf OUTBOND "\n"; 
					for $f ("ORDER", "DISPX", "DISPY", "DISPZ") {
						next if(! exists($atoms->{$l}{$f}));
						printf OUTBOND "%-6s%6d",$f, $l + $offset;
						for $m (@{ $atoms->{$l}{$f} }) {
							printf OUTBOND "%6d", $m;
						}
						printf OUTBOND "\n";
					}
				}
				$strLen = printProgress($count, $total, $start, $str);
				$offset += $nAtoms;
			}
		}
	}
	close OUTATOM;
	close OUTBOND;
	for $i ("X", "Y", "Z") {
		$box->{$i}{hi} = ($box->{$i}{len} * $dims->{$i});
		$box->{$i}{lo} = 0;
		$box->{$i}{len} = ($box->{$i}{len} * $dims->{$i});
	}
	$total = getTime(time() - $start);
	printf "$str%-${strLen}s\n", "${total}s elapsed...Done";

	#now consolidate all the temporary atoms and bonds
	system("cat ${outdir}/__atoms.dat ${outdir}/__bonds.dat > ${outdir}/tmp.bgf");
	($atoms, $bonds) = FileFormats::ParseStructFile("${outdir}/tmp.bgf",0);
	system("rm -fr ${outdir}");
	return ($atoms, $bonds, $box);
}

sub ReplicateCellOld {
	my ($atoms, $bonds, $box, $dims, $centerMol, $str, $updateResNum, $offset1D, $qAdjust, $setPBC) = @_;
	my ($i, $j, $k, $cellAtoms, $cellBonds, $unitAtoms, $pbcbonds); 
	my ($strLen, $offset, $tot, $count, $total, $start, $vec);
	
	$qAdjust->{X} = $qAdjust->{Y} = $qAdjust->{Z} = 0 
		if (! defined($qAdjust));
	$setPBC = 1 if (! defined($setPBC));	
	&GetBoxDisplacementTensor($box);
	&setPBCbonds($atoms,$bonds,$box) if ($setPBC);
	print "${str}calculating time remaining\r";
	$start = time();
	if ($centerMol) {
		@{ $vec } = (int(($dims->{X}-1)/2),int(($dims->{Y}-1)/2),int(($dims->{Z}-1)/2));
		$atoms = TransCellAtomsOld($atoms, $box, $vec);
	}

	$total = 1;
	for $i ("X", "Y", "Z") {
		$total *= $dims->{$i};
	}

	for $i ("X", "Y", "Z") {
		$unitAtoms = dclone($atoms);
		$cellBonds = dclone($bonds);
		$tot = scalar(keys %{ $atoms });
		$offset = 0;
		for $j (1 .. ($dims->{$i} - 1)) {
			$count++;
			$strLen = printProgress($count, $total, $start, $str);
			$offset += $tot;
			@{ $vec } = (0,0,0);
			$vec->[0] = 1 if ($i eq "X");
			$vec->[1] = 1 if ($i eq "Y");
			$vec->[2] = 1 if ($i eq "Z");
			$cellAtoms = TransCellAtomsOld($unitAtoms, $box, $vec);
			if($offset1D) {
				for $k ("X", "Y", "Z") {
					next if ($k eq $i);
					$cellAtoms = TransCellAtomsOld($cellAtoms, $box, $vec);
				}
			}
			($atoms, $bonds) = combineMols($atoms, $cellAtoms, $bonds, $cellBonds, $updateResNum);
			&updatePBCbonds($atoms, $bonds, $offset, $tot, $i, $qAdjust);
		}
		$box->{$i}{hi} = ($box->{$i}{len} * $dims->{$i});
		$box->{$i}{lo} = 0;
		$box->{$i}{len} = ($box->{$i}{len} * $dims->{$i});
	}
	$tot = getTime(time() - $start);
	printf "$str%-${strLen}s\n", "${tot}s elapsed...Done";
	return ($atoms, $bonds, $box);
}

sub setPBCbonds {
	my ($atoms, $bonds, $box) = @_;
	my ($i,$j,$d,$bl,$pbl,$disp,$counter);

	for $i (keys %{ $atoms }) {
		for $d ("X","Y","Z") {
			next if (exists($atoms->{$i}{"DISP${d}"}));
			next if (!$bonds->{$i});
			$atoms->{$i}{"DISP${d}"} = ();
			$counter = 0;
			for $j (0 .. $#{ $bonds->{$i} }) {
				$disp = 0;
				$bl = abs($atoms->{$i}{"${d}COORD"}-$atoms->{ $bonds->{$i}[$j] }{"${d}COORD"});
				#look in +dim direction for bond
				$pbl = abs($atoms->{$i}{"${d}COORD"}-$atoms->{ $bonds->{$i}[$j] }{"${d}COORD"}+$box->{$d}{DISP_V}{$d});
				$disp = -1 if($pbl<$bl);
				#loos in -dim direction for bond
				$pbl = abs($atoms->{$i}{"${d}COORD"}-$atoms->{ $bonds->{$i}[$j] }{"${d}COORD"}-$box->{$d}{DISP_V}{$d});
				$disp = 1 if($pbl<$bl);
				push @{ $atoms->{$i}{"DISP${d}"} }, $disp;
				$counter += abs($disp);
			}
			delete $atoms->{$i}{"DISP${d}"} if (!$counter);
		}
	}
}

sub updatePBCbonds {
	my ($atoms, $bonds, $aOffset, $cellatoms, $dim, $dQ) = @_;
	my ($i, $j, $disp, $counter);

	for $i (keys %{ $atoms }) {
		next if (! exists($atoms->{$i}{"DISP${dim}"}));
		$counter = 0;
		for $j (0 .. $#{ $bonds->{$i} }) {
			$disp = $atoms->{$i}{"DISP${dim}"}[$j];
			next if (!$disp);
			if ($i <= $aOffset) { #atoms in original cell
				if($disp<0) {
					$bonds->{$i}[$j] -= $disp*$cellatoms;
					##$atoms->{$i}{"DISP${dim}"}[$j]++;
				} else {
					$bonds->{$i}[$j] += $disp*$aOffset;
					$atoms->{$i}{"DISP${dim}"}[$j]-- if ($disp>0);
					##$atoms->{$i}{"DISP${dim}"}[$j]--;
				}
			} else {
				if($disp<0) {
					$bonds->{$i}[$j] += $disp*$cellatoms;
					$atoms->{$i}{"DISP${dim}"}[$j]++;
				} else {
					$bonds->{$i}[$j] -= $disp*$aOffset;
					##$atoms->{$i}{"DISP${dim}"}[$j]--;
				}
				$atoms->{$i}{CHARGE} += $dQ->{$dim}/2;
				$atoms->{ $bonds->{$i}[$j] }{CHARGE} += $dQ->{$dim}/2;
			}
			$counter += abs($atoms->{$i}{"DISP${dim}"}[$j]);
		}
		delete $atoms->{$i}{"DISP${dim}"} if ($counter == 0);
	}
}

sub updatePBCbonds_old {
	my ($atoms, $bonds, $bondlist, $atomOffset, $tot, $dim) = @_;
	my ($i, $j, $atom1, $atom2, $atom3, $atom4, @list, $pos);

	for $i (keys %{ $bondlist }) {
		@list = keys %{ $bondlist->{$i} };
		for $j (0 .. $#list) {
			$pos = $bondlist->{$i}{ $list[$j] };
			$atom1 = $i;
			$atom2 = $list[$j];
			$atom3 = $atom1 + $atomOffset;
			$atom4 = $atom2 + $tot;
			#for atom1 , delete bond to atom2 and form bond to atom4
			$bonds->{$atom1}[$pos->[0]] = $atom4;
			#for atom2, delete bond to atom1 and form bond to atom3
			$bonds->{$atom2}[$pos->[1]] = $atom3;
			delete $atoms->{$atom2}{"DISP${dim}"};
			#$atoms->{$atom2}{"DISP${dim}"}[$pos->[1]] = 0;
			#for atom3, delete bond to atom4 and form bond to atom2
			$bonds->{$atom3}[$pos->[0]] = $atom2;
			delete $atoms->{$atom3}{"DISP${dim}"};
			#$atoms->{$atom3}{"DISP${dim}"}[$pos->[0]] = 0;
			#updateBond($atoms->{$atom3}, $bonds->{$atom3}, $atom4, $atom2);
			#delete $atoms->{$atom3}{"DISP${dim}"};
			#for atom4, delete bond to atom3 and form bond to atom1
			$bonds->{$atom4}[$pos->[1]] = $atom1;
			#updateBond($atoms->{$atom4}, $bonds->{$atom4}, $atom3, $atom1);
			#now update bondlist
			delete $bondlist->{$i}{$atom2};
			$bondlist->{$i}{$atom4}[0] = $pos->[0];
			$bondlist->{$i}{$atom4}[1] = $pos->[1];
		}
	}
	print "";
}

sub getPBCbonds {
	my ($atoms, $bonds, $box, $k) = @_;
	my ($i, $j, $l, $dist, $atom1, $atom2, $bondlist, $sign,$flag2);

	#search for multiple bond entries to same atom and image flags
	for $i (keys %{ $bonds }) {
		$atom1 = $i;
		for $j (0 .. $#{ $bonds->{$i} }) {
			$sign = 0;
			$sign = $atoms->{$i}{"DISP${k}"}[$j] if (exists($atoms->{$i}{"DISP${k}"}));
			next if($sign>-1);
			#find symmetric entry
			undef($flag2);
			$atom2 =  $bonds->{$i}[$j];
			for $l (0 .. $#{ $bonds->{$atom2} }) {
				if ($bonds->{$atom2}[$l] == $atom1) {
					if (exists($atoms->{$atom2}{"DISP${k}"}) and $atoms->{$atom2}{"DISP${k}"}[$l] == -1*$sign) {
						$flag2 = $l;
						last;
					}
				}
			}
			next if (! defined($flag2));
			$bondlist->{$i}{$bonds->{$i}[$j]} = ([$j, $flag2]);
		}
	}

	return $bondlist if(defined($bondlist));
	#search for bonds by distance
	for $i (keys %{ $bonds }) {
		$atom1 = \%{ $atoms->{$i} };
		for $j (0 .. $#{ $bonds->{$i} }) {
			$atom2 = \%{ $atoms->{$bonds->{$i}[$j]} };
			for $k ("X","Y","Z") {
				$dist = $atom1->{"${k}COORD"} - $atom2->{"${k}COORD"};
				next if ($dist == 0);
				$sign = $dist/abs($dist);
				if ($sign < 0) {
					$dist *= -1;
				}
				if ($dist > 5 and $dist > $box->{$k}{len}/3) { #search for long bonds
					$bondlist->{$bonds->{$i}[$j]}{$i} = ([$j,$sign]);
				}
			}
		}
	}

	return $bondlist;
}
#helper functions
sub numerically { ($a<=>$b); }
sub combineMols {
	my ($mol1, $mol2, $CONN, $connections, $updateRes) = @_;
	my ($atom, $tot_atoms, $tot_res, @tmp, $ATOMS, @BONDS, %CONS, $bond);

	$updateRes = 1 if (! defined($updateRes));
	@tmp = sort numerically keys %{ $mol1 };
	$tot_atoms = $tmp[$#tmp];
	$tot_res = $mol1->{$tot_atoms}{"RESNUM"};
	@tmp = sort numerically keys %{ $mol2 };

	$ATOMS = $mol1;
	%CONS = %{ $CONN };
	for $atom (@tmp) {
		$ATOMS->{($atom + $tot_atoms)} = dclone($mol2->{$atom});
		$ATOMS->{($atom + $tot_atoms)}{"RESNUM"} += $tot_res if($updateRes);
		$ATOMS->{($atom + $tot_atoms)}{"INDEX"} = $atom + $tot_atoms;
		@BONDS = ();
		@BONDS = @{ $connections->{$atom} } if (defined($connections->{$atom}));
		for $bond (@BONDS) {
			push @{ $CONS{($atom + $tot_atoms)} }, ($bond + $tot_atoms);
		}
	}

	return ($ATOMS, \%CONS);
}

sub printProgress {
	my ($currPos, $total, $start, $pStr) = @_;
	my ($progress, $str, $end);
	
	$end = time();

	$progress = $currPos/$total;
	
	$str = sprintf("%.2f%% complete %s\r", 
				   100*$progress, getEta(($end - $start), $progress));
	
	print "${pStr}${str}" if (defined($pStr));
	return length($str);
}

sub getTime {
	my ($timeLeft) = $_[0];
	my ($returnStr);
 
	if ($timeLeft > 60) {
		if ($timeLeft > 3600) {
			if ($timeLeft > 86400) {
				$returnStr = int($timeLeft/86400) . "d ";
				$timeLeft = $timeLeft % 86400;
			}
			$returnStr .= int($timeLeft/3600) . "h ";
			$timeLeft = $timeLeft % 3600;
		}
		$returnStr .= int($timeLeft/60) . "m ";
		$timeLeft = $timeLeft % 60;
	}
	$returnStr .= sprintf("%.0f", $timeLeft);
	return $returnStr;
}

sub getEta {
	my ($elapsed, $percentage) = @_;
	my ($totalTime) = $elapsed/$percentage;
	my ($timeLeft) = $totalTime - $elapsed;
	my ($returnStr) = "(";
	
	$returnStr .= getTime($timeLeft) . "s remaining)		  ";
	
	return $returnStr;
}
1;

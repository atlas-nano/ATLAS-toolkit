#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester CoM CombineMols LoadFFs ReadFFs GetBondLength CoM GetFFTypeStr GetFileTypeStr);
use FileFormats qw(ParseStructFile createBGF createHeaders addHeader
					addBoxToHeader GetBGFAtoms AddMass);
use BOX qw(GetBox CreateGrid GetNeighbours Cart2Frac);
use ManipAtoms qw(GetMols MapOrigin SplitAtomsByMol);
use Getopt::Std qw(getopt);

use File::Basename qw(basename);

my ($DATA, $FFILES, $ATOMS, $BONDS, $HEADERS, $BOX, $PARMS, $reversePlace);
my ($i, $soluStruct, $centerType, $soluAtms, $solvStruct, $saveName, $overlap);

$|++;
&init;
$PARMS = LoadFFs($FFILES) if (defined($FFILES));
print "Parsing solvent/membrane bgf $solvStruct...";
($DATA->{SOLVENT}{ATOMS}, $DATA->{SOLVENT}{BONDS}, $DATA->{SOLVENT}{HEADERS}) = ParseStructFile($solvStruct, 1);
print "Done\nParsing solute bgf $soluStruct...";
($DATA->{SOLUTE}{ATOMS}, $DATA->{SOLUTE}{BONDS}, $DATA->{SOLUTE}{HEADERS}) = ParseStructFile($soluStruct, 1);
$soluAtms = scalar(keys %{ $DATA->{SOLUTE}{ATOMS} });
print "Done\n";
&checkAtomTypes($DATA, $PARMS) if (defined($FFILES));
&AddMass($DATA->{SOLVENT}{ATOMS}, $PARMS) if (defined($FFILES));
&AddMass($DATA->{SOLUTE}{ATOMS}, $PARMS) if (defined($FFILES));
for $i ("SOLUTE", "SOLVENT") {
	print "Computing ${centerType} for $i...\r";
	$DATA->{$i}{BOX} = GetBox($DATA->{$i}{ATOMS}, undef, $DATA->{$i}{HEADERS});
}
print "Computing ${centerType}s...Completed\n";
print "Embedding systems...";
&centerMols($DATA) if ($centerType ne "NONE");
if (! $reversePlace) {
	($ATOMS, $BONDS) = CombineMols($DATA->{SOLUTE}{ATOMS}, $DATA->{SOLVENT}{ATOMS},
		$DATA->{SOLUTE}{BONDS}, $DATA->{SOLVENT}{BONDS});
} else {
	($ATOMS, $BONDS) = CombineMols($DATA->{SOLVENT}{ATOMS}, $DATA->{SOLUTE}{ATOMS},
		$DATA->{SOLVENT}{BONDS}, $DATA->{SOLUTE}{BONDS});
}
&GetMols($ATOMS, $BONDS);
$BOX = GetBox($ATOMS, $PARMS, $DATA->{SOLVENT}{HEADERS});
&MapOrigin($BOX, CoM($ATOMS));
&addRadii($ATOMS, $PARMS);
print "Done\n";
($ATOMS, $BONDS) = checkOverlaps($ATOMS, $BONDS, $soluAtms, $BOX, $reversePlace) if($overlap);
print "Creating BGF file $saveName...";
$HEADERS = createHeaders(undef, $saveName);
&addBoxToHeader($HEADERS, $BOX);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub checkOverlaps {
	my ($atoms, $bonds, $nsolu, $box, $rev) = @_;
	my ($i, $solvAtms, $solvMols, $TMP, $totAtms, $nsolv);

	print "Removing overlaps...";
	$totAtms = scalar(keys %{ $atoms } );
	$nsolv = 0;
	if ($rev) {
		for $i (1 .. ($totAtms - $nsolu)) {
			$solvAtms->{$i} = 1;
			$nsolv++;
		}
	} else {
		for $i (($nsolu + 1) .. $totAtms) {
			$solvAtms->{$i} = 1;
			$nsolv++;
		}
	}
	$solvMols = GetMols($atoms, $bonds, $solvAtms);
	$TMP = &removeOverlaps($atoms, $nsolu, $nsolv, $solvMols, $box, $rev);
	($atoms, $bonds) = GetBGFAtoms($TMP, $atoms, $bonds);
	print "Done\n";
	return ($atoms, $bonds);
}

sub removeOverlaps {
	my ($atoms, $nsolu, $nsolv, $molSolv, $box, $rev) = @_;
	my ($i, $j, $k, $l, $catoms, $mgrid, $solvCoM, $gspacing, $pos, $molData, $solvGIndex, $solvGNeigh, $nremove);

	for $i ("X","Y","Z") {
		$gspacing->{$i} = 2.5/$BOX->{$i}{len};
	}
	#place solute atoms on grid
	for $i (keys %{ $atoms }) {
		$catoms->{$i} = 1;
		next if($i>$nsolu and ! $rev);
		next if($i<=$nsolv and $rev);
		$pos = ();
		%{ $pos } = (XCOORD=>$atoms->{$i}{XCOORD},YCOORD=>$atoms->{$i}{YCOORD},ZCOORD=>$atoms->{$i}{ZCOORD});
		$mgrid->{int($atoms->{$i}{FA}/$gspacing->{X})}{int($atoms->{$i}{FB}/$gspacing->{Y})}{int($atoms->{$i}{FC}/$gspacing->{Z})}{SOLUTE}{$i} = \%{ $pos };
		$atoms->{$i}{GRID}{X} = int($atoms->{$i}{FA}/$gspacing->{X});
		$atoms->{$i}{GRID}{Y} = int($atoms->{$i}{FB}/$gspacing->{Y});
		$atoms->{$i}{GRID}{Z} = int($atoms->{$i}{FC}/$gspacing->{Z});
	}
	
	$nremove = 0;
	for $i (keys %{ $molSolv }) {
		($molData,undef) = GetBGFAtoms($molSolv->{$i}{MEMBERS},$atoms,$BONDS);
		$solvCoM->{1} = CoM($molData);
		&Cart2Frac($solvCoM, $box);
		@{ $solvGIndex } = (int($solvCoM->{1}{FA}/$gspacing->{X}),int($solvCoM->{1}{FB}/$gspacing->{Y}),int($solvCoM->{1}{FC}/$gspacing->{Z}));
		$solvGNeigh = getNeighbors($mgrid,@{ $solvGIndex });
		for $j (@{ $solvGNeigh }) {
			foreach $k (keys %{ $j->{SOLUTE} }) {
				next if (!exists($catoms->{$k}));
				if(GetBondLength($solvCoM->{1},$atoms->{$k},$box) < 3) {
					$nremove++;
					for $l (keys %{ $molSolv->{$i}{MEMBERS} }) {
						delete $catoms->{$l};
					}
				}
			}
		}
	}

	print "removed $nremove solvent molecules...";
	return $catoms;
}

sub getNeighbors {
	my ($mgrid, $ix, $iy, $iz) = @_;
	my ($cgrid, $i, $j, $k);

	for $i (-1 .. 1) {
		next if !(exists($mgrid->{$ix+$i}));
		for $j (-1 .. 1) {
			next if !(exists($mgrid->{$ix+$i}{$iy+$j}));
			for $k (-1 .. 1) {
				next if !(exists($mgrid->{$ix+$i}{$iy+$j}{$iz+$k}));
				push @{ $cgrid }, $mgrid->{$ix+$i}{$iy+$j}{$iz+$k};
			}
		}
	}
	
	return $cgrid;
}

sub removeOverlapsOld {
	my ($bgfAtoms, $grid, $offset, $molList, $box) = @_;
	my ($i, $x, $y, $z, $CELL, $solvent, $CLIST, $molID); 
	my ($dist, $index, $currCell, $atoms, $j, $solAtom, $GRID);

	%{ $atoms } = %{ $bgfAtoms };
	for $i (keys %{ $atoms }) {
		 next if (! exists($atoms->{$i}) or $atoms->{$i}{IS_SOLVENT});
		$x = $atoms->{$i}{CELL}{"XINDEX"};
		$y = $atoms->{$i}{CELL}{"YINDEX"};
		$z = $atoms->{$i}{CELL}{"ZINDEX"};
		next if (! $x or ! $y or ! $z);
		$currCell = $grid->{$x}{$y}{$z};
		next if (exists($currCell->{VISITED}));
		$CLIST = GetNeighbours($GRID, $currCell);
		$currCell->{VISITED} = 1;
		for $solAtom (@{ $currCell->{ATOMS} }) {
			for $CELL (@{ $CLIST }) {
				for $solvent (@{ $CELL->{WATERS} }) {
					next if (! $solvent->{IS_SOLVENT});
					$dist = GetBondLength($solAtom,$solvent,$box);
					if ($dist < 3) { #remove entire solvent molecule
						$molID = $solvent->{MOLECULEID};
						for $j (keys %{ $molList->{$molID} }) {
							next if (! exists($atoms->{($j + $offset)}));
							delete $atoms->{($j + $offset)};
						}
					}
				}
			}
		}
	}
	
	return $atoms;
}

sub removeOverlapsFinal {
	my ($bgfAtoms, $grid, $offset, $molList, $box) = @_;
	my ($i, $x, $y, $z, $CELL, $solvent, $CLIST, $molID); 
	my ($dist, $index, $currCell, $atoms, $j, $solAtom, $GRID);

	my ($solvStart, $tot, $totMol, $solvStartMol);
	$solvStart = scalar(keys %{ $DATA->{SOLUTE}{ATOMS} });
	$solvStartMol = $bgfAtoms->{$solvStart}{MOLECULEID};
	$tot = scalar(keys %{ $bgfAtoms });
	$totMol = $bgfAtoms->{$tot}{MOLECULEID};

	for $i ($solvStartMol .. $tot) {
	}
	
	%{ $atoms } = %{ $bgfAtoms };
	for $i (keys %{ $atoms }) {
		 next if (! exists($atoms->{$i}) or $atoms->{$i}{IS_SOLVENT});
		$x = $atoms->{$i}{CELL}{"XINDEX"};
		$y = $atoms->{$i}{CELL}{"YINDEX"};
		$z = $atoms->{$i}{CELL}{"ZINDEX"};
		next if (! $x or ! $y or ! $z);
		$currCell = $grid->{$x}{$y}{$z};
		next if (exists($currCell->{VISITED}));
		$CLIST = GetNeighbours($GRID, $currCell);
		$currCell->{VISITED} = 1;
		for $solAtom (@{ $currCell->{ATOMS} }) {
			for $CELL (@{ $CLIST }) {
				for $solvent (@{ $CELL->{WATERS} }) {
					next if (! $solvent->{IS_SOLVENT});
					$dist = GetBondLength($solAtom,$solvent,$box);
					if ($dist < 3) { #remove entire solvent molecule
						$molID = $solvent->{MOLECULEID};
						for $j (keys %{ $molList->{$molID} }) {
							next if (! exists($atoms->{($j + $offset)}));
							delete $atoms->{($j + $offset)};
						}
					}
				}
			}
		}
	}
	
	return $atoms;
}
sub init {
	my (%OPTS, @tmp, @FFILES, $i, $FF);
	
	getopt('msfcwor',\%OPTS);
	($solvStruct, $soluStruct,$FF,     $centerType,$saveName,$overlap,$reversePlace) = (
	  $OPTS{m},   $OPTS{s},   $OPTS{f},$OPTS{c},   $OPTS{w}, $OPTS{o},$OPTS{r});
	
	die &usage if (! defined($solvStruct) || ! defined($soluStruct));
	
	print "Initializing...";

	$centerType = "COM" if (! defined($centerType) or $centerType !~ /^com|cog|cop|none/);
	$centerType = uc($centerType);
	$centerType = "COP" if ($centerType !~/COM|NONE/);

	if (! defined($saveName)) {
		$soluStruct =~ /^\s*(\S+)/;
		$saveName = basename($1);
		$saveName =~ s/\.\w+$/_embed\.bgf/;
	}

	#if ($FF =~ /\s+/) {
	#	@tmp = split /\s+/, $FF;
	#} else {
	#	@tmp = ($FF);
	#}
  
	#for $i (@tmp) {
	#	if (-e $i &&  -r $i && -T $i) {
	#		push @FFILES, $i;
	#	}
	#}

    ($FFILES, undef) = ReadFFs($FF, 0);

	$overlap = 1 if (! defined($overlap));
	$overlap = 1 if ($overlap !~ /^0|no/i);
	$overlap = 0 if ($overlap =~ /^0|no/i);

	print "Done\n";
}

sub addRadii {
	my ($atoms, $parms) = @_;
	my ($i, $ffType);

	for $i (keys %{ $atoms }) {
		$ffType = $atoms->{$i}{FFTYPE};
		if (! defined($PARMS)) {
			$atoms->{$i}{RADII} = 1;
		} else {
			$atoms->{$i}{RADII} = $parms->{VDW}{$ffType}{$ffType}{1}{VALS}[1];
		}
	}
}

sub checkAtomTypes {
	my ($mols, $parms) = @_;
	my ($i, $j, $ffType);

	for $i ("SOLUTE", "SOLVENT") {
		for $j (keys %{ $mols->{$i}{ATOMS} }) {
			$ffType = $mols->{$i}{ATOMS}{$j}{FFTYPE};
			$ffType = $parms->{OVERWRITE_FFTYPES}{$ffType} 
				if (exists($parms->{OVERWRITE_FFTYPES}) and exists($parms->{OVERWRITE_FFTYPES}{$ffType}));
			die "ERROR: Force field type $ffType not found in forcefield(s). Aborting\n" 
				if (! exists($parms->{ATOMTYPES}{$ffType}));
			if ($centerType eq "COM") {
				$mols->{$i}{ATOMS}{$j}{RADII} = $parms->{VDW}{$ffType}{$ffType}{1}{VALS}[1];
				$mols->{$i}{ATOMS}{$j}{MASS} = $parms->{ATOMTYPES}{$ffType}{MASS};
			}
			if ($i eq "SOLVENT") {
				$mols->{$i}{ATOMS}{$j}{IS_SOLVENT} = 1;
			} else {
				$mols->{$i}{ATOMS}{$j}{IS_SOLUTE}  = 1;
			}
		}
	}
}

sub embeddMols_old {
	my ($mols) = $_[0];
	my ($i, $j, @dims);

	@dims = keys %{ $mols->{SOLVENT}{BOX} };
	for $i (@dims) {
		$mols->{SOLUTE}{OFFSET}{$i} = (($mols->{SOLUTE}{BOX}{$i}{hi} + $mols->{SOLUTE}{BOX}{$i}{lo})/2 - 
									  ($mols->{SOLVENT}{BOX}{$i}{hi} + $mols->{SOLVENT}{BOX}{$i}{lo})/2);
	}

	for $i (keys %{ $mols->{SOLUTE}{ATOMS} }) {
		for $j (@dims) {
			$mols->{SOLUTE}{ATOMS}{$i}{"${j}COORD"} -= $mols->{SOLUTE}{OFFSET}{$j};
		}
	}
}

sub centerMols {
	my ($mols) = $_[0];
	my ($i, $j, @dims, $coms);

	$coms->{SOLUTE} = CoM($mols->{SOLUTE}{ATOMS});
	$coms->{SOLVENT} = CoM($mols->{SOLVENT}{ATOMS});

	@dims = ("X","Y","Z");
	for $i (@dims) {
		$mols->{SOLUTE}{OFFSET}{$i} = ($coms->{SOLVENT}{"${i}COORD"} - $coms->{SOLUTE}{"${i}COORD"});
	}

	for $i (keys %{ $mols->{SOLUTE}{ATOMS} }) {
		for $j (@dims) {
			$mols->{SOLUTE}{ATOMS}{$i}{"${j}COORD"} += $mols->{SOLUTE}{OFFSET}{$j};
		}
	}
}

sub usage {
	my $fTypeStr = GetFileTypeStr;
	my $ffTypeStr = GetFFTypeStr;
	my ($usage) = <<DATA;	
usage: $0 -s solu_struct -m solv_struct_file -f "forcefield(s)" -c (centering) -o (remove_overlap) -r (reverse_placement) -w (saveName) 
Required arguments:
	-s solu_struct: Location of the solute structure file.
$fTypeStr	-m solv_struct: Location of solvent/membrane structure file (see above)
	-f \"forcefield(s)": 1 or more Cerius2|Polygraf|ITP|CHARMM_PRM formatted forcefields
$ffTypeStr
Optional arrguments:
	-c centering: Ways of centering the solute. Expected:
		com (center of mass - default)
		cog (center of gravity)
		cop (center of position)
		none
	-o overlap: Flag specifying whether to remove overlapping solvent molecules. Expected yes|1 (default) or no|0
	-r reverse_placement: Flag for whether to reverse the order of the final structure, to be solvent:solute.
		Expected no|0 (default) or yes|1
	-w saveName: Name of embedded structure file (in BGF format)	
DATA

}

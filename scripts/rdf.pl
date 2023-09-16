#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use warnings;
no warnings "recursion";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo AddMass createBGF createHeaders addHeader);
use General qw(STDev FileTester CoM GetSelections IsInteger LoadFFs
						 IsDecimal TrjSelections CenterOnMol GetBondLength GetStats);
use AMBER qw(ParseAmberTrj GetAmberByteOffset);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use BOX qw(GetBox ConvertBox MakeBox CreateGrid GetNeighbours);
use ManipAtoms qw(UnwrapAtoms ImageAtoms GetAtmList SplitAtomsByMol 
							GetMols SelectAtoms BuildAtomSelectionString);
use FSearch qw(StoreCoords GetGridATOMS SetResiduals);
use constant PI => atan2(1,1) * 4;

sub numerically { ($a<=>$b); }
sub saveFile;
sub calcRDF;
sub getResiduals;
sub setResiduals;
sub storeCoords;
sub getGridATOMS;
sub init;
sub boxConvert;
sub showUsage;

my ($SOLUTE, $SOLVENT, $trjFile, $BGF, $BONDS, $rMax, $printStr, $gridLen);
my ($getSnapshot, $getByteOffset, $trjType, $LAMMPSOPTS, $delR, $nBins, $rounding);
my ($BBOX, $reImage, $SELECT, $saveFile, $field, $DATA, $PARMS, $allMOLS);

$|++;
&init;
$gridLen = 5;

if (! $trjType) {
	print "Calculating RDF...";
	&calcRDF($BGF, $BBOX, 1, undef);
	print "Done\n";
} else {
	$field = scalar keys %{ $BGF };
	$getByteOffset->($SELECT, $trjFile, $field);
	if ($trjType == 2) {
		&GetLammpsTrjType($SELECT, $trjFile, "atom", \%{ $LAMMPSOPTS });
		$field = "atom";
	}
	$printStr = "Calculating RDF from $trjFile...";
	$getSnapshot->($BGF, $trjFile, $SELECT, $field, \&calcRDF, $printStr, undef);
}

print "Saving RDF data to $saveFile...";
saveFile($DATA, $saveFile);
print "Done\n";

sub saveFile {
	my ($data, $saveName) = @_;
	my ($frame, $bin, $rlower, $rcenter, $rupper); 
	my ($nSolu, $nSolv, $ncoord_stats, $nideal);
	my (@ncoord_ave, @gr_ave, $gr, $ncoord); 
	my ($gr_stats, $avgFile);

	$nSolu = scalar(keys %{ $SOLUTE });
	$nSolv = scalar(keys %{ $SOLVENT });
	open OUTFILE, "> $saveName" or die "ERROR: Cannot create file $saveName: $!\n";
	for $frame (keys %{ $data }) {
		print OUTFILE "\#Snapshot \# $frame\n";
		$ncoord = 0;
		for $bin (0 .. $nBins) {
			$rlower = $bin*$delR;
			$rcenter = $rlower + 0.5*$delR;
			$rupper = $rlower + $delR;
			$nideal = $data->{$frame}{CONST} * $nSolv *
						($rupper*$rupper*$rupper - $rlower*$rlower*$rlower);
			$gr = $data->{$frame}{VALS}[$bin]/($nSolu*$nideal);
			$ncoord += $gr*$nideal;
			print OUTFILE "$rcenter $gr $ncoord\n";
			$gr_ave[$bin] += $gr;
			$ncoord_ave[$bin] += $ncoord;
		}
		print OUTFILE "\n";
	}
	close OUTFILE;

	if ($saveName =~ /(.*)\.(.+)$/) {
		$avgFile = "${1}_avg.${2}";
	} else {
		$avgFile = "${saveName}_avg.rdf";
	}

	open OUTFILE, "> $avgFile" or die "ERROR: Cannot create $avgFile:$!\n";
	print OUTFILE "\#AVERAGE\n";
	for $bin (0 .. $nBins) {
		$gr_ave[$bin] /= $nBins;
		$ncoord_ave[$bin] /= $nBins;
		$rlower = $bin*$delR;
		$rcenter = $rlower + 0.5*$delR;
		$rupper = $rlower + $delR;
		print OUTFILE "$rcenter $gr_ave[$bin] $ncoord_ave[$bin]\n";
	}
	close OUTFILE;
}

sub calcRDF {
	my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
	my ($sCoords, $BBOX, $i, @tmp, $j, $x, $y, $z);
	my ($gridATOMS, $soluAtoms, $solvAtoms, $count);
	my ($residualList, $vRes, $dist, $tot, $bin);

	$count = 0;
	if ($trjType == 2) { #LAMMPS
		$frameNum = $ATOMS->{TIMESTEP}[0];
		$BOX = ConvertLammpsBox($ATOMS->{"BOX BOUNDS"});
		$BBOX = $BOX;
		$ATOMS = $ATOMS->{ATOMS};
		if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
			UnwrapAtoms($DATA->{ATOMS}, $BOX, $LAMMPSOPTS->{scaled});
		}
		for $i (keys %{ $ATOMS }) {
			@tmp = grep {!/COORD/i} keys %{ $BGF->{$i} };
			for $j (@tmp) {
				next if ($j =~ /COORD/i);
				$ATOMS->{$i}{$j} = $BGF->{$i}{$j};
			}
		}		
	} elsif ($trjType == 1) {
		$BBOX = MakeBox($BOX);
	} else {
		$BBOX = $BOX;
	}
	
	($gridATOMS, $soluAtoms, $solvAtoms) = GetGridATOMS($ATOMS, $SOLUTE, $SOLVENT);
	$sCoords = StoreCoords($gridATOMS, $BOX);
	$tot = scalar(keys %{ $gridATOMS });
	for $i (0 .. $nBins) {
		$DATA->{$frameNum}{VALS}[$i] = 0; #initialize all elements of array
	}

	for $j (@{ $soluAtoms }) {
		next if (! $j);
		$x = $gridATOMS->{$j}{SORT}{X};
		$residualList = ();
		$vRes = ();
		&SetResiduals($solvAtoms, $rMax, $j, $tot, \@{ $residualList }, \@{ $vRes });
		next if (! getResiduals($x, $sCoords->{X}, $residualList, $vRes, $tot, $BBOX->{XCOORD}{len}));
		$y = $gridATOMS->{$j}{SORT}{Y};
		next if (! getResiduals($y, $sCoords->{Y}, $residualList, $vRes, $tot, $BBOX->{YCOORD}{len}));
		$z = $gridATOMS->{$j}{SORT}{Z};
		next if (! getResiduals($z, $sCoords->{Z}, $residualList, $vRes, $tot, $BBOX->{ZCOORD}{len}));
		for $i (1 .. $#{ $vRes }) {
			next if (! $vRes->[$i]);
			$dist = sqrt(($rMax - $residualList->[$i]));
			$bin = sprintf("%${rounding}f",($dist/$delR));
			$DATA->{$frameNum}{VALS}[$bin]++;
			$count++;
		}
	}
	$DATA->{$frameNum}{TOTAL} = $count;
	$DATA->{$frameNum}{CONST} =  4.0*PI/
		(3 * $BBOX->{XCOORD}{len}*$BBOX->{YCOORD}{len}*$BBOX->{ZCOORD}{len});
}

sub getResiduals {
	my ($arrayIndex, $distArray, $residuals, $validRes, $tot, $boxLen) = @_;
	my ($i, $baseVal, $offset, $j, $atomC, $tmp, $count);

	$baseVal = $distArray->[$arrayIndex][0];
	$count = 0;
	
	for $i (1 .. $tot) {
		$tmp->[$i] = $validRes->[$i];
		$validRes->[$i] = 0;
	}
 
	for ($j = ($arrayIndex + 1); $j < $tot; $j++) {
		$atomC = $distArray->[$j][1];
		next if (! $tmp->[$atomC]);
		$offset = $distArray->[$j][0] - $baseVal;
		$offset -= $boxLen while ($offset > $boxLen);
		$residuals->[$atomC] -= ($offset * $offset);
		if ($residuals->[$atomC] > 0) {
		$validRes->[$atomC] = $atomC;
		$count++; }
	}

	for ($j = ($arrayIndex - 1); $j > -1; $j--) {
		$atomC = $distArray->[$j][1];
		next if (! $tmp->[$atomC]);
		$offset = ($baseVal - $distArray->[$j][0]);
		$offset -= $boxLen while ($offset > $boxLen);
		$residuals->[$atomC] -= ($offset * $offset);
		if ($residuals->[$atomC] > 0) {
		$validRes->[$atomC] = $atomC;
		$count++; }
	}

	return $count;
}

sub setResiduals {
	my ($solvAtoms, $distMax, $skip, $tot, $residuals, $vRes) = @_;
	my ($i);
	
	for $i (1 .. $tot) {
		$residuals->[$i] = $distMax;
		$vRes->[$i] = $i;
		$vRes->[$i] = 0 if (! $solvAtoms->[$i]);
	}
	$vRes->[$skip] = 0;
}

sub storeCoords {
	my ($atoms, $box) = @_;
	my (%cSORTED, $DATA, $i, $d, %atomMap, $index, @arryIndex, $tot);

	@arryIndex =keys %{ $atoms };
	$tot = $#arryIndex;
	for $i (@arryIndex) {
		$d = $atoms->{$i}{XCOORD};
		$d += 0.00001 while (exists($DATA->{X}{$d}));
		$DATA->{X}{$d} = $i;

		$d = $atoms->{$i}{YCOORD};
		$d += 0.00001 while (exists($DATA->{Y}{$d}));
		$DATA->{Y}{$d} = $i;

		$d = $atoms->{$i}{ZCOORD};
		$d += 0.00001 while (exists($DATA->{Z}{$d}));
		$DATA->{Z}{$d} = $i;
	}

	for $i ("X", "Y", "Z") {
		$index = 0;
		$#{ $cSORTED{$i} } = $tot;
		for $d (sort numerically keys %{ $DATA->{$i} }) {
			$cSORTED{$i}[$index][0] = $d;
			$cSORTED{$i}[$index][1] = $DATA->{$i}{$d};
			$atoms->{ $DATA->{$i}{$d} }{SORT}{$i} = $index;
			$index++;
		}
	}

	return \%cSORTED;
}

sub getGridATOMS {
	my ($allAtoms, $soluMols, $solvMols) = @_;
	my ($i, $count, %atomCOMS, @solu, @solv, $MOL, %SOLVSOLU, @tmp);
	
	$count = 0;
	$solu[0] = $solv[0] = 0;
	for $i (keys %{ $soluMols }) {
		$count++;
		$MOL = GetAtmData($allAtoms, $soluMols->{$i});
		$atomCOMS{$count} = CoM($MOL);
		$atomCOMS{$count}{IS_SOLUTE} = 1;
		$atomCOMS{$count}{MOLECULEID} = $i;
		$atomCOMS{$count}{RESNAME} = "SOL";
		$solu[$count] = $count;
		@tmp  = keys %{ $soluMols->{$i} };
		if ($allAtoms->{ $tmp[0] }{IS_SOLVENT}) {
			$SOLVSOLU{ $allAtoms->{ $tmp[0] }{SOLVENTID} } = 1;
			$atomCOMS{$count}{RESNAME} = "WAT";
			$solv[$count] = $count;
		} else {
			$solv[$count] = 0;
		}
	}
		
	for $i (keys %{ $solvMols }) {
		next if (exists($SOLVSOLU{$i}));
		$count++;
		$MOL = GetAtmData($allAtoms, $solvMols->{$i});
		$atomCOMS{$count} = CoM($MOL);
		$atomCOMS{$count}{IS_SOLUTE} = 1;
		$atomCOMS{$count}{MOLECULEID} = $i;
		$atomCOMS{$count}{RESNAME} = "WAT";
		$solv[$count] = $count;
		$solu[$count] = 0;
	}

	return (\%atomCOMS, \@solu, \@solv);
}
		
sub init {
	my (%OPTS, @tmp, $i, $j, $list, $solvTmp, $HEADERS);
	my ($soluSel, $solvSel, $bgfFile, $tSel, $usage, $FF, $prec);
	getopt('mvbtdrlswif',\%OPTS);

	$usage = &showUsage;

	for ($OPTS{m},$OPTS{v},$OPTS{b}) {
		die "$usage" if (! defined($_));
	}

	print "Initializing...";
	($soluSel,$solvSel,$trjFile,$bgfFile,$rMax,$delR,$trjType,$tSel,$saveFile,$reImage,$FF) = 
		($OPTS{m},$OPTS{v},$OPTS{t},$OPTS{b},$OPTS{r},$OPTS{d},$OPTS{l},$OPTS{s},$OPTS{w},$OPTS{i},$OPTS{f});

	if (defined($trjFile) and -e $trjFile and -r $trjFile and -T $trjFile) {
		if (! defined($trjType) or lc($trjType) ne "lammps") {
			$trjType = 1;
			$getSnapshot = \&ParseAmberTrj;
			$getByteOffset = \&GetAmberByteOffset;
		} else {
			$trjType = 2;
			$getSnapshot = \&ParseLAMMPSTrj;
			$getByteOffset = \&GetLammpsByteOffset;
		}
		if (! defined($tSel)) {
			$tSel = "*";
		}
		$list = TrjSelections($tSel);
		for $j (keys %{ $list }) {
			$SELECT->{$j} = $list->{$j};
		}
		
		die "ERROR: No valid frames selected with selection $tSel!\n"
			if (! keys %{ $SELECT } and $tSel ne "*");
	} else {
		$trjType = 0;
	}
		
	FileTester($bgfFile);

	$soluSel = "*" if (!defined($soluSel));
	$solvSel = "*" if (!defined($solvSel));
	$soluSel = BuildAtomSelectionString($soluSel);
	$solvSel = BuildAtomSelectionString($solvSel);

	$delR = 0.1 if (! defined($delR) or (! IsInteger($delR) and ! IsDecimal($delR)));

	while ($delR =~ /(\d|.)/g) {
		if ($1 ne ".") {
			if ($1 > 0) {
				$prec .= "1";
			} else {
				$prec .= "0";
			}
		} else {
			$prec .= ".";
		}
	}

	$rounding = $prec;
	chop $rounding;
	$rounding .= "0";

	if (! defined($saveFile)) {
		$saveFile = basename($bgfFile);
		$saveFile =~ s/\.\w+$/\.rdf/;
	}
	$reImage = 1 if (! defined($reImage) or $reImage !~ /^0/);
	print "Done\nParsing BGF file $bgfFile...";
	($BGF, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
	&GetMols($BGF, $BONDS);

	$SOLUTE = SelectAtoms($soluSel, $BGF);
	for (keys %{ $SOLUTE }) {
		$BGF->{$_}{IS_SOLUTE} = 1;
	}
	$SOLUTE = SplitAtomsByMol($BGF, $SOLUTE);
	die "ERROR: Solute atoms does not correspond to any in BGF file!\n"
		if (! keys %{ $SOLUTE });
	for $i (keys %{ $SOLUTE }) {
		for $j (keys %{ $SOLUTE->{$i} }) {
			$BGF->{$j}{SOLUTEID} = $i;
		}
	}
									  
	$SOLVENT = SelectAtoms($solvSel, $BGF);
	for (keys %{ $SOLVENT }) {
		$BGF->{$_}{IS_SOLVENT} = 1;
	}
	$SOLVENT = SplitAtomsByMol($BGF, $SOLVENT);
	die "ERROR: Solvent atoms does not correspond to any in BGF file!\n"
		if (! keys %{ $SOLVENT });
	for $i (keys %{ $SOLVENT }) {
		for $j (keys %{ $SOLVENT->{$i} }) {
			$BGF->{$j}{SOLVENTID} = $i;
		}
	}
	
	$allMOLS = SplitAtomsByMol($BGF, $BGF);

	$BBOX =GetBox($BGF, undef, $HEADERS);
	&boxConvert($BBOX);
	if (! defined($rMax) or (! IsInteger($rMax) and ! IsDecimal($rMax))) {
		$rMax = 12;
	}
	$nBins = $rMax/$delR; # total number of bins
	$rMax *= $rMax;
	print "Done\n";
	if (exists($OPTS{f})) {
		if ($FF =~ /\s/) {
			@tmp = split /\s+/, $FF;
		} else {
			$tmp[0] = $FF;
		}
		
		$PARMS = LoadFFs(\@tmp);
		AddMass($BGF, $PARMS);
	}
}

sub boxConvert {
	my ($box) = $_[0];
	
	for ("X", "Y", "Z") {
		%{ $box->{$_ . "COORD"} } = %{ $box->{$_} };
	}
}

sub showUsage {
	my ($usage) = "usage: $0 -m solute atoms -v solvent atoms -b bgf file " . 
		"\nOptional parameters:" . 
		"\n\t-w save name -r max distance -d del r -f cerius2 forcefield" . 
		"\n\t-t trajectory -s traj selection -l traj type (amber(default)|lammps)\n";
	return $usage;
}

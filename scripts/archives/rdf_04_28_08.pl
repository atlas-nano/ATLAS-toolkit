#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use warnings;
no warnings "recursion";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetBGFFileInfo AddMass createBGF createHeaders addHeader);
use Packages::General qw(STDev FileTester CoM GetSelections IsInteger LoadFFs
			 IsDecimal TrjSelections CenterOnMol GetBondLength);
use Packages::AMBER qw(ParseAmberTrj GetAmberByteOffset);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use Packages::BOX qw(GetBox ConvertBox MakeBox CreateGrid GetNeighbours);
use Packages::ManipAtoms qw(UnwrapAtoms ImageAtoms GetAtmList SplitAtomsByMol 
			    GetAtmData);
use constant PI => atan2(1,1) * 4;

sub init;
sub calcRDF;
sub saveFile;
sub getAtoms;
sub boxConvert;
sub saveFile;
sub numerically { ($a<=>$b); }
sub showUsage;
sub getSoluSolv;
sub getGridAtoms;
sub setResiduals;
sub getResiduals;
sub storeCoords;

my ($SOLUTE, $SOLVENT, $trjFile, $BGF, $BONDS, $rMax, $printStr, $gridLen);
my ($getSnapshot, $getByteOffset, $trjType, $LAMMPSOPTS, $delR, $nBins);
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
    my ($i, $j, $totSolu, $totSolv, $STATS, $norm); 
    my ($avg, $stdev, $rl, $ru, $avgFile);

    $totSolv = scalar keys %{ $SOLVENT };
    $totSolu = scalar keys %{ $SOLUTE };
    open OUTFILE, "> $saveName" or die "ERROR: Cannot create file $saveName: $!\n";
    for $i (sort numerically keys %{ $data }) {
	print OUTFILE "\#Snapshot \# $i\n";
	for $j (0 .. $#{ $data->{$i}{VALS} }) {
	    $rl = $j*$delR;
	    $ru = $rl + $delR;
	    $norm = $data->{$i}{CONST} * $data->{$i}{TOTAL} * ($ru*$ru*$ru - $rl*$rl*$rl);
	    $data->{$i}{VALS}[$j] /= ($norm * $totSolu);
	    print OUTFILE "$rl $data->{$i}{VALS}[$j]\n";
	    $STATS->[$j] .= "$data->{$i}{VALS}[$j] ";
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
    for $j (0 .. $#{ $STATS }) {
	chomp $STATS->[$j];
	($avg, $stdev, $i) = STDev($STATS->[$j]);
	print OUTFILE ($j*$delR) . " $avg $stdev\n";
    }
    close OUTFILE;
}

sub getAtoms {
    my ($allAtoms, $atomList) = @_;
    my (%ATOMS, $i);

    for $i (keys %{ $atomList }) {
        $ATOMS{$i} = $allAtoms->{$i};
    }

    return \%ATOMS;
}

sub calcRDF {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($MOL, $CENTER, @tmp, $i, $j, $k, $SOLCENTER, $solvMols, %USED); 
    my ($bin, $dist, $count, $prec, $rounding, $soluMols);
    my ($GRID, $currCell, $numNeigh, $x, $y, $z, $nCells);
    my ($comATOMS, $solvAtoms, $soluAtoms, $density, $vRes);
    my ($residualList, $sCoords, $tot, $atomMap, $d, $gridATOMS);

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

    if ($trjType == 2) { #LAMMPS
	$frameNum = $ATOMS->{TIMESTEP}[0];
        $BOX = ConvertLammpsBox($ATOMS->{"BOX BOUNDS"});
        $BBOX = $BOX;
        $ATOMS = $ATOMS->{ATOMS};
	if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
	    UnwrapAtoms($DATA->{ATOMS}, $BOX, $LAMMPSOPTS->{scaled});
	}
        @tmp = grep {!/COORD/i} keys %{ $BGF->{1} };
        for $i (keys %{ $ATOMS }) {
            for $j (@tmp) {
                next if ($j =~ /COORD/i or ! exists($BGF->{$i}{$j}));
                $ATOMS->{$i}{$j} = $BGF->{$i}{$j};
            }
        }	
    } elsif ($trjType == 1) {
	$BBOX = MakeBox($BOX);
    } else {
	$BBOX = $BOX;
    }
    
    for $i (0 .. $nBins) {
        $DATA->{$frameNum}{VALS}[$i] = 0; #initialize all elements of array
    }

    $MOL = GetAtmData($ATOMS, $SOLUTE->{1});
    $SOLCENTER = CoM($MOL);

#    if ($reImage) {
#        @tmp = ("X", "Y", "Z");
#        for $j (@tmp) {
#            $BBOX->{$j}{CENTER} = $BBOX->{$j . "COORD"}{len}/2;
#	    $BBOX->{$j . "COORD"}{CENTER} = $BBOX->{$j}{CENTER};
#            $BBOX->{$j}{hi} = $BBOX->{$j . "COORD"}{len};
#	    $BBOX->{$j . "COORD"}{hi} = $BBOX->{$j}{hi};
#            $BBOX->{$j}{lo} = 0;
#	    $BBOX->{$j . "COORD"}{lo} = $BBOX->{$j}{lo};
#            for $i (keys %{ $ATOMS }) {
#                $ATOMS->{$i}{$j . "COORD"} += ($BBOX->{$j}{CENTER} - $SOLCENTER->{$j . "COORD"});
#            }
#        }
        #CenterOnMol($ATOMS, $SOLCENTER);
#        for $i (keys %{ $allMOLS }) {
#            $MOL = GetAtmData($ATOMS, $allMOLS->{$i});
#            $CENTER = CoM($MOL);
#            ImageAtoms($MOL, $CENTER, $BOX);
#         }

#    }

#    $tot = createHeaders(undef, "test.bgf");
#    &addHeader($ATOMS, $tot);
#    createBGF($ATOMS, $BONDS, "test.bgf");
    ($gridATOMS, $soluAtoms, $solvAtoms) = getGridATOMS($ATOMS, $SOLUTE, $SOLVENT);
    $sCoords = storeCoords($gridATOMS);
    $count = 0;
    $tot = scalar(keys %{ $solvAtoms }) - 1;
    for $j (keys %{ $soluAtoms }) {
	$x = $soluAtoms->{$j}{SORT}{X};
	($residualList, $vRes) = setResiduals($solvAtoms, $rMax, $j);
	next if (! getResiduals($x, $sCoords->{X}, $residualList, $vRes, $tot, $BBOX->{XCOORD}{len}));
	$y = $soluAtoms->{$j}{SORT}{Y};
	next if (! getResiduals($y, $sCoords->{Y}, $residualList, $vRes, $tot, $BBOX->{YCOORD}{len}));
	$z = $soluAtoms->{$j}{SORT}{Z};
	next if (! getResiduals($z, $sCoords->{Z}, $residualList, $vRes, $tot, $BBOX->{ZCOORD}{len}));
	for $i (keys %{ $vRes }) {
	    #next if (! exists($comATOMS->{$i}{IS_SOLVENT}));
	    $dist = sqrt(($rMax - $residualList->{$i}));
	    #next if ($dist < $delR);
	    $bin = sprintf("%${rounding}f",($dist/$delR));
	    $DATA->{$frameNum}{VALS}[$bin]++;
	    $count++;
	}
    }

    $DATA->{$frameNum}{TOTAL} = $tot + 1;
    $DATA->{$frameNum}{CONST} =  4.0*PI/
	(3 * $BBOX->{XCOORD}{len}*$BBOX->{YCOORD}{len}*$BBOX->{ZCOORD}{len});
}

sub getResiduals {
    my ($arrayIndex, $distArray, $residuals, $validRes, $tot, $boxLen) = @_;
    my ($i, $baseVal, $offset, $j, $atomC, $tmp);

    $baseVal = $distArray->[$arrayIndex][0];
    $j = $arrayIndex + 1;
    $tot++;
    while ($j < $tot) {
	$atomC = $distArray->[$j][1];
	if (! exists($validRes->{$atomC})) {
	    $j++;
	    next;
	}
	$offset = ($distArray->[$j][0] - $baseVal);
	if (abs($offset) > $boxLen) {
	    $tmp = $offset/$boxLen;
	    if ($tmp =~ /(\..*)/) {
		$offset = $1 * $boxLen;
	    } else {
		$offset = 0;
	    }
	}

	$offset *= $offset;
	$residuals->{$atomC} -= $offset;
	if ($residuals->{$atomC} < 0) {
	    delete $validRes->{$atomC};
	    $i = $j + 1;
	    while ($i < $tot) {
		$atomC = $distArray->[$i][1];
		delete $validRes->{$atomC} if (exists($validRes->{$atomC}));
		$i++;
	    }
	    last;
	} 
	$j++;
    }

    $j = $arrayIndex - 1;
    while ($j > -1) {
        $atomC = $distArray->[$j][1];
        if (! exists($validRes->{$atomC})) {
	    $j--;
	    next;
	}
        $offset = $distArray->[$j][0] - $baseVal;
	if (abs($offset) > $boxLen) {
	    $tmp = $offset/$boxLen;
	    if ($tmp =~ /(\..*)/) {
		$offset = $1 * $boxLen;
	    } else {
		$offset = 0;
	    }
	}
        $offset *= $offset;
        $residuals->{$atomC} -= $offset;
        if ($residuals->{$atomC} < 0) {
            delete $validRes->{$atomC};
            $i = $j - 1;
            while ($i > -1) {
                $atomC = $distArray->[$i][1];
                delete $validRes->{$atomC} if (exists($validRes->{$atomC}));
		$i--;
            }
            last;
	}
	$j--;
    }

    return 0 if (! %{ $validRes });
    return 1;
}

sub setResiduals {
    my ($atomsData, $distMax, $skip) = @_;
    my (%residuals, $i, %vRes);

    for $i (keys %{ $atomsData }) {
	next if ($i == $skip);
	$residuals{$i} = $distMax;
	$vRes{$i} = $i;
    }

    return (\%residuals, \%vRes);
}

sub storeCoords {
    my ($atoms) = @_;
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
	for $d (sort keys %{ $DATA->{$i} }) {
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
    my ($i, $count, %atomCOMS, %solu, %solv, $MOL, %SOLVSOLU, @tmp);
    
    $count = 0;
    for $i (keys %{ $soluMols }) {
	$count++;
	$MOL = GetAtmData($allAtoms, $soluMols->{$i});
	$atomCOMS{$count} = CoM($MOL);
	$atomCOMS{$count}{IS_SOLUTE} = 1;
	$atomCOMS{$count}{MOLECULEID} = $i;
	$atomCOMS{$count}{RESNAME} = "SOL";
	$solu{$count} = $atomCOMS{$count};
	@tmp  = keys %{ $soluMols->{$i} };
	if ($allAtoms->{ $tmp[0] }{IS_SOLVENT}) {
	    $SOLVSOLU{ $allAtoms->{ $tmp[0] }{SOLVENTID} } = 1;
	    $atomCOMS{$count}{RESNAME} = "WAT";
	    $solv{$count} = $atomCOMS{$count};
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
	$solv{$count} = $atomCOMS{$count}
    }

    return (\%atomCOMS, \%solu, \%solv);
}
	
sub getSoluSolv {
    my ($atomData, $currCell, $nCells) = @_;
    my (@SOLV, $i, $cell, @SOLU);

    for $i (@{ $currCell->{IONS} }, @{ $currCell->{WATERS} }, @{ $currCell->{ATOMS} }) {
	push @SOLU, $i if ($i->{IS_SOLUTE});
    }

    for $cell (@{ $nCells }) {
	for $i (@{ $cell->{WATERS} }) {
	    push @SOLV, $i if ($i->{IS_SOLVENT});
	}	
    }
    
    return (\@SOLU, \@SOLV);
}

sub init {
    my (%OPTS, @tmp, $i, $j, $list, $solvTmp, $HEADERS);
    my ($solSel, $solvSel, $bgfFile, $tSel, $usage, $FF);
    getopt('mvbtdrlswif',\%OPTS);

    $usage = &showUsage;

    for ($OPTS{m},$OPTS{v},$OPTS{b}) {
	die "$usage" if (! defined($_));
    }

    print "Initializing...";
    ($solSel,$solvSel,$trjFile,$bgfFile,$rMax,$delR,$trjType,$tSel,$saveFile,$reImage,$FF) = 
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

    if ($solSel =~ /\s+/) {
	@tmp = split /\s+/, $solSel;
    } else {
	$tmp[0] = $solSel;
    }
    $list = GetSelections(\@tmp, 0);
    for $j (keys %{ $list }) {
	$SOLUTE->{$j} = $list->{$j};
    }
    die "ERROR: No valid atoms selected for SOLUTE with selection $solSel!\n"
	if (! keys %{ $SOLUTE });

    
    if ($solvSel =~ /\s+/) {
	@tmp = split /\s+/, $solvSel;
    } else {
	$tmp[0] = $solvSel;
    }
    $list = GetSelections(\@tmp, 0);
    for $j (keys %{ $list }) {
	$SOLVENT->{$j} = $list->{$j};
    }
    die "ERROR: No valid atoms selected for SOLVENT with selection $solvSel!\n"
	if (! keys %{ $SOLVENT });

    $delR = 0.1 if (! defined($delR) or (! IsInteger($delR) and ! IsDecimal($delR)));

    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$/\.rdf/;
    }
    $reImage = 1 if (! defined($reImage) or $reImage !~ /^0/);
    print "Done\nParsing BGF file $bgfFile...";
    ($BGF, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);

    $SOLUTE = GetAtmList($SOLUTE, $BGF);
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
                                      
    $SOLVENT = GetAtmList($SOLVENT, $BGF);
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

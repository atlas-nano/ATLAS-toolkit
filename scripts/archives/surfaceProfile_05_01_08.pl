#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use warnings;
no warnings "recursion";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetBGFFileInfo AddMass);
use Packages::General qw(STDev FileTester CoM GetSelections TrjSelections GetBondLength
			 CenterOnMol IsInteger IsDecimal ShowSelectionInfo STDev LoadFFs);
use Packages::AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use Packages::BOX qw(GetBox CreateGrid CenterAtoms GetGridDims GetNeighbours MoveAtomsToOrigin);
use Packages::ManipAtoms qw(UnwrapAtoms ImageAtoms GetAtmList GetSolvent SplitAtomsByMol GetAtmData);

sub init;
sub showUsage;
sub boxConvert;
sub calcProfile;
sub getAtoms;
sub getAtmsInCell;
sub getGridList;
sub numerically { ($a<=>$b); }
sub getResRange;
sub writeDensity;
sub writeShellData;

my ($SOLUTE, $SOLVENT, $trjFile, $BGF, $BONDS, $rMax, $printStr, $axis);
my ($getSnapshot, $getByteOffset, $trjType, $LAMMPSOPTS, $delR);
my ($SELECT, $saveFile, $field, $DATA, $PARMS, $prec, $rounding);
my ($MOL_VOL, $HEADERS, $GRIDOPTS, $BBOX, $DIM, $calcType, $shellData);

$|++;
$calcType = 1;
&init;

if (! $trjType) {
    print "Calculating solvent <-> surface profile...";
    &calcProfile($BGF, $BBOX, 1, undef);
    print "Done\n";
} else {
    $field = scalar keys %{ $BGF };
    $getByteOffset->($SELECT, $trjFile, $field);
    if ($trjType == 2) {
        &GetLammpsTrjType($SELECT, $trjFile, "atom", \%{ $LAMMPSOPTS });
        $field = "atom";
    }
    $printStr = "Calculating solvent <-> surface profile from $trjFile...";
    $getSnapshot->($BGF, $trjFile, $SELECT, $field, \&calcProfile, $printStr, undef);
}

&writeDensity($DATA, $saveFile) if ($calcType == 1);
&writeShellData($DATA, $saveFile) if ($calcType != 1);

sub writeShellData {
    my ($data, $saveName) = @_;
    my (%SHELLMOLS, $shellName, $shellID, $i, $prefix);
    my (%SHELLTRANS, $j, $resID, $resData, $k, $tot, %SHELLCOUNT);

    $prefix = $_[1];
    $prefix =~ s/\.\w+$//;
    for $i (keys %{ $data }) {
	if ($#{ $data->{$i}{SHELLS} } == 0) { # stayed in one shell entire time
	    ($shellName, $shellID) = ($data->{$i}{SHELLS}[0]{NAME}, $data->{$i}{SHELLS}[0]{ID});
	    $SHELLMOLS{$shellName}{$i} = $data->{$i}{TRAJ};
	} else { # moved in between shells
	    for $j (1 .. ($#{ $data->{$i}{SHELLS} } - 1)) {
		$shellName = $data->{$i}{SHELLS}[$j]{NAME} . "_" . $data->{$i}{SHELLS}[($j + 1)]{NAME};
		$k = $data->{$i}{SHELLS}[$j]{STOP} - $data->{$i}{SHELLS}[$j]{START};
		$SHELLTRANS{$shellName} .= "$k ";
		$SHELLCOUNT{$shellName}++;
	    }
	}
    }

    for $shellName (keys %SHELLMOLS) {
	print "Writing molecules occupying $shellName...";
	system("mkdir -p ${prefix}/${shellName}");
	$resData = ();
	for $i (keys %{ $SHELLMOLS{$shellName} }) {
	    for $j (keys %{ $SOLVENT->{$i} }) {
		$resID = $BGF->{$j}{RESNUM};
		last;
	    }
	    $resData->{$resID} = 1;
	    $saveName = "${prefix}/${shellName}/${resID}.dat";
	    open PROFDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
	    for $j (sort numerically keys %{ $SHELLMOLS{$shellName}{$i} }) {
		print PROFDATA $SHELLMOLS{$shellName}{$i}{$j}{DIST} .  "\n";
	    }
	    close PROFDATA;
	}
	$j = getResRange($resData);
	$saveName = "${prefix}_${shellName}.dat";
	open SHELLDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
	for $k (@{ $j }) {
	    print SHELLDATA " $k";
	}
	close SHELLDATA;
	print "Done\n";
    }

    $saveName = "${prefix}_shelltrans.dat";
    print "Writing shell transition statistics to $saveName...\n";
    open SHELLTRANS, "> $saveName" or die "ERROR: Cannot write to $saveName: $!\n";
    for $i (keys %SHELLTRANS) {
	$shellName = $i;
	$shellName =~ s/_/ \-> /;
	chomp $SHELLTRANS{$i};
	($j, $k, $tot) = STDev($SHELLTRANS{$i});
	printf SHELLTRANS "%-20s %12.3f %8.3f %8d %8d\n", $shellName, $j, $k, $tot, $SHELLCOUNT{$i};
	printf "%-20s %12.3f %8.3f %8d %8d\n", $shellName, $j, $k, $tot, $SHELLCOUNT{$i};
    }
    close SHELLTRANS;
}

sub writeDensity {
    my ($data, $saveName) = @_;
    my ($prefix, $rT, $frame, $count, $bin, $i1, $i2);
    my (%binData, $avg, $stdev, $tot);

    $prefix = $_[1];
    $prefix =~ s/\.\w+$//;
    $rT = 1;
    print "Writing data to $saveFile...";
    for $frame (keys %{ $data }) {
	$count++;
	for $bin (keys %{ $data->{$frame} }) {
	    for $i1 (keys %{ $data->{$frame}{$bin} }) {
		for $i2 (keys %{ $data->{$frame}{$bin}{$i1} }) {
		    next if (! exists($data->{$frame}{$bin}{$i1}{$i2}));
		    $binData{$bin} .= $data->{$frame}{$bin}{$i1}{$i2} . " ";
		}
	    }
	}
    }
    
    open OUTDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    printf OUTDATA "%-${prec}s %8s %8s %5s\n","#", "TOTAL", "AVG", "STDEV";
    for $bin (sort numerically keys %binData) {
	chop $binData{$bin};
	($avg, $stdev, $tot) = STDev($binData{$bin});
	printf OUTDATA "%${prec}f %8d %8.3f %5.3f %8d\n", ($bin * $delR), $tot/$avg, $avg, $stdev, $tot;
    }
    close OUTDATA;
    print "Done\n";
}

sub calcProfile {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($MOL, $CENTER, @tmp, $i, $j, $SOLCENTER, $SOLATMS, $SOLVATMS, $index);
    my ($bin, $dist, $count, $GRID, $shellName, $shellID, $solvGRIDLIST);
    my ($x, $y, $z, $k, $solGRIDLIST, $SOLVCENTER, $minDist, $i1, $i2);

    if ($trjType == 2) { #LAMMPS
        $frameNum = $ATOMS->{TIMESTEP}[0];
        $BOX = ConvertLammpsBox($ATOMS->{"BOX BOUNDS"});
        $ATOMS = $ATOMS->{ATOMS};
	if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
	    UnwrapAtoms($ATOMS, $BBOX, $LAMMPSOPTS->{scaled});
	}
	@tmp = grep {!/COORD|VISITED/i} keys %{ $BGF->{1} };
	for $i (keys %{ $BGF }) {
	    for $j (@tmp) {
		$ATOMS->{$i}{$j} = $BGF->{$i}{$j};
	    }
	    $ATOMS->{$i}{IS_SOLVENT} = 1 if (exists($BGF->{$i}{IS_SOLVENT}));
	}
    } elsif ($trjType == 1) {
        $BOX = ConvertAmberBox($BOX);
    } 

    $MOL = GetAtmData($ATOMS, $SOLUTE);
    $CENTER = CoM($MOL);
    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    for $j (@tmp) {
	$BOX->{$j}{CENTER} = $BOX->{$j}{len}/2;
	$BOX->{$j}{hi} = $BOX->{$j}{len};
	$BOX->{$j}{lo} = 0;
	for $i (keys %{ $ATOMS }) {
	    $ATOMS->{$i}{$j} += ($BOX->{$j}{CENTER} - $CENTER->{$j});
	}
    }
     
    for $i (keys %{ $SOLVENT }) {
	$MOL = GetAtmData($ATOMS, $SOLVENT->{$i});
	$CENTER = CoM($MOL);
	ImageAtoms($MOL, $CENTER, $BOX);
    }

   for $i ("X", "Y", "Z") {
	$BOX->{$i} = $BOX->{"${i}COORD"};
   }
   ($GRID, $i, $j) = CreateGrid($ATOMS, 0, $BOX, 2.5, 0);

    $count = 0;
    for $i (keys %{ $SOLUTE }) {
	($x, $y, $z) = ($ATOMS->{$i}{CELL}{XINDEX}, $ATOMS->{$i}{CELL}{YINDEX},$ATOMS->{$i}{CELL}{ZINDEX});
	next if (exists($GRID->{$x}{$y}{$z}{VISITED}) and $GRID->{$x}{$y}{$z}{VISITED});
	$GRID->{$x}{$y}{$z}{VISITED} = 1;
	$i1 = eval('$' . $DIM->[0]);
	$i2 = eval('$' . $DIM->[1]);

	($solGRIDLIST, $solvGRIDLIST) = getGridList($GRID, $axis, $x, $y, $z);
	($SOLATMS, $SOLVATMS) = getAtmsInCell($ATOMS, $SOLUTE, $solGRIDLIST, $solvGRIDLIST);

	for $j (keys %{ $SOLVATMS }) {
	    $SOLVCENTER = CoM($SOLVATMS->{$j});
	    undef($minDist);
	    for $k (keys %{ $SOLATMS }) {
		$dist = GetBondLength($SOLATMS->{$k}, $SOLVCENTER);
		$minDist = $dist if (! defined($minDist) or ($minDist > $dist));
	    }
	    next if (! defined($minDist));
	    $bin = sprintf("%${rounding}f",($minDist/$delR) + 1);
	    if ($calcType == 1)  {
		$DATA->{$frameNum}{$bin}{$i1}{$i2}++;
	    } else {
		undef($shellID);;
		for $k (@{ $shellData }) {
		    next if ($minDist > $k->{END});
		    $shellName = $k->{NAME};
		    $shellID = $k->{ID};
		    last;
		}
		($shellName, $shellID) = ("bulk", 0) if (! defined($shellID));
		    
		$DATA->{$j}{TRAJ}{$frameNum} = (
						{
						    "NAME" => $shellName,
						    "DIST" => $minDist,
						}
						);
		if (! exists($DATA->{$j}{SHELLS})) {
		    $DATA->{$j}{SHELLS}[0]{NAME} = $shellName;
		    $DATA->{$j}{SHELLS}[0]{START} = $frameNum;
		    $DATA->{$j}{SHELLS}[0]{ID} = $shellID;
		    $index = 0;
		} else { 
		    $index = $#{ $DATA->{$j}{SHELLS} };
		    if ($DATA->{$j}{SHELLS}[$index]{NAME} ne $shellName) {
			$index++;
			$DATA->{$j}{SHELLS}[$index]{ID} = $shellID;
			$DATA->{$j}{SHELLS}[$index]{NAME} = $shellName;
			$DATA->{$j}{SHELLS}[$index]{START} = $frameNum;
		    }
		}
		$DATA->{$j}{SHELLS}[$index]{STOP} = $frameNum;
	    }
	    $count++;
	}		 
    }
}

sub getGridList {
    my ($grid, $dim, $x, $y, $z) = @_;
    my (@solvGRIDS, @solGRIDS, $i, $j);

    if ($dim eq "x") {
	for $i (keys %{ $grid }) {
	    push @solvGRIDS, $grid->{$i}{$y}{$z};
	}
       
	for $i (($y - 1) .. ($y + 1)) {
	    for $j (($z - 1) .. ($z + 1)) {
		next if (! exists($grid->{$x}{$i}) or ! exists($grid->{$x}{$i}{$j}));
		push @solGRIDS, $grid->{$x}{$i}{$j};
	    }
	}
    } elsif ($dim eq "y") {
	for $i (keys %{ $grid->{$x} }) {
	    push @solvGRIDS, $grid->{$x}{$i}{$z};
	}

	for $i (($x - 1) .. ($x + 1)) {
	    for $j (($z - 1) .. ($z + 1)) {
		next if (! exists($grid->{$i}) or ! exists($grid->{$i}{$y}) or ! exists($grid->{$i}{$y}{$j}));
		push @solGRIDS, $grid->{$i}{$y}{$j};
	    }
	}
    } else {
	for $i (keys %{ $grid->{$x}{$y} }) {
	    push @solvGRIDS, $grid->{$x}{$y}{$i};
	}

	for $i (($x - 1) .. ($x + 1)) {
	    for $j (($y - 1) .. ($y + 1)) {
		next if (! exists($grid->{$i}{$j}{$z}));
		push @solGRIDS, $grid->{$i}{$j}{$z};
	    }
	}

    }
    return (\@solGRIDS, \@solvGRIDS);
}

sub getAtmsInCell {
    my ($atoms, $solute, $solCells, $solvCells) = @_;
    my (%SOLUTEATMS, %SOLVATMS, $i, $counter, $j, $cell, $k);

    $counter = 1;
    for $cell (@{ $solCells }) {
	for $i (@{ $cell->{ATOMS} }) {
	    next if (! exists($solute->{ $i->{INDEX} }) or $atoms->{$i->{INDEX}}{IS_SOLVENT});
	    $SOLUTEATMS{$counter} = $atoms->{$i->{INDEX}};
	    $counter++;
	}
    }

    for $cell (@{ $solvCells }) {
	for $i (@{ $cell->{ATOMS} }, @{ $cell->{WATERS} }, @{ $cell->{IONS} }) {
	    next if (exists($atoms->{$i->{INDEX}}{VISITED}));
	    if (exists($i->{IS_SOLVENT}) and $i->{IS_SOLVENT}) {
		for $j (keys %{ $i->{MOLECULE} }) {
		    $SOLVATMS{$i->{MOLECULEID}}{$j} = $atoms->{$j};
		    $atoms->{$j}{VISITED} = 1;
		}
	    }
	}
    }

    return (\%SOLUTEATMS, \%SOLVATMS)
    
}

sub getAtoms {
    my ($allAtoms, $atomList) = @_;
    my (%ATOMS, $i);

    for $i (keys %{ $atomList }) {
        $ATOMS{$i} = $allAtoms->{$i};
    }

    return \%ATOMS;
}

sub init {
    my (%OPTS, @tmp, $i, $j, $list, $solvTmp, $BOX, $nBins, $count);
    my ($solSel, $solvSel, $bgfFile, $tSel, $usage, $FF, $solvShells);

    getopt('mvbtdrlswafc',\%OPTS);

    $usage = &showUsage;

    for ($OPTS{a},$OPTS{b}) {
        die "$usage" if (! defined($_));
    }

    print "Initializing...";
    ($solSel,$solvSel,$trjFile,$bgfFile,$rMax,$delR,$trjType,$tSel,$saveFile,$axis,$FF, $solvShells) =
        ($OPTS{m},$OPTS{v},$OPTS{t},$OPTS{b},$OPTS{r},$OPTS{d},$OPTS{l},$OPTS{s},$OPTS{w},$OPTS{a},$OPTS{f},$OPTS{c});

    if ($axis =~ /(X|Y|Z)/i) {
	$axis = lc($1);
    } else {
	die "ERROR: Expected X|Y|Z for axis. Got \"$axis\"!\n";
    }

    @{ $DIM } = grep {!/$axis/} ("x","y","z");

    if (defined($trjFile) and -e $trjFile and -r $trjFile and -T $trjFile) {
	if (! defined($trjType)) {
	    if ($trjFile =~ /\.lammpstrj/) {
		$trjType = "lammps";
	    } else {
		$trjType = "amber";
	    }
	}

        if (lc($trjType) ne "lammps") {
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

    if (defined($solSel)) {
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
    }

    if (defined($solvSel)) {
	if ($solvSel =~ /\s+/) {
	    @tmp = split /\s+/, $solvSel;
	} else {
	    $tmp[0] = $solvSel;
	}
	$list = GetSelections(\@tmp, 0);
	for $j (keys %{ $list }) {
	    $solvTmp->{$j} = $list->{$j};
	}
	
	die "ERROR: No valid atoms selected for SOLVENT with selection $solvSel!\n"
	    if (! keys %{ $solvTmp });
    }

    $delR = 0.1 if (! defined($delR) or (! IsInteger($delR) and ! IsDecimal($delR)));

    if (! defined($saveFile)) {
        $saveFile = basename($bgfFile);
        $saveFile =~ s/\.\w+$//;
	$saveFile .= ".rdf";
    }

    print "Done\nParsing BGF file $bgfFile...";
    ($BGF, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
    if (defined($SOLUTE)) {
	$SOLUTE = GetAtmList($SOLUTE, $BGF);
    } else {
	$SOLUTE = \%{ $BGF->{1}{MOLECULE} };
    }
    die "ERROR: Solute atoms does not correspond to any in BGF file!\n"
        if (! keys %{ $SOLUTE });
    if (defined($solvTmp)) {
	$SOLVENT = GetAtmList($solvTmp, $BGF);
    } else {
	$SOLVENT = GetSolvent($BGF,"WATER");
    }
    die "ERROR: Solvent atoms does not correspond to any in BGF file!\n"
        if (! keys %{ $SOLVENT });
    $SOLVENT = SplitAtomsByMol($BGF, $SOLVENT);
    for $i (keys %{ $SOLVENT }) {
	for $j (keys %{ $SOLVENT->{$i} }) {
	    $BGF->{$j}{IS_SOLVENT} = 1;
	}
    }
    
    $BOX = GetBox($BGF, undef, $HEADERS);
    %{ $BBOX } = %{ $BOX };
    &boxConvert($BBOX);
    if (! defined($rMax) or (! IsInteger($rMax) and ! IsDecimal($rMax))) {
        $rMax = 0;
        for $i ("X", "Y", "Z") {
            $rMax = $BBOX->{$i}{len} if ($BBOX->{$i}{len} > $rMax);
        }
    }

    $nBins = $rMax/$delR; # total number of bins
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

    $count = 0;
    if (defined($solvShells)) {
	@tmp = split /\s+/, $solvShells;
	for $i (@tmp) {
	    if ($i =~ /^\d+\.?\d*$/) {
		$count++;
		push @{ $shellData }, (
				       {
					   "END" => $i,
					   "ID"  => $count,
					   "NAME" => "shell${count}",
				       }
				       );
		$calcType = 2;
	    }
	}
	if (defined($shellData)) {
	    push @{ $shellData }, (
				   {
				       "END" => $rMax,
				       "ID"  => 0,
				       "NAME" => "bulk",
				   }
				   );
	}
    }
    print "Done\n";

    if (exists($OPTS{f})) {
        if ($FF =~ /\s/) {
            @tmp = split /\s+/, $FF;
        } else {
            $tmp[0] = $FF;
        }
	$list = ();
	for $i (@tmp) {
	    push @{ $list }, $i if (-e $i);
	}
 
        $PARMS = LoadFFs($list);
        AddMass($BGF, $PARMS);
    }
}

sub boxConvert {
    my ($box) = $_[0];

    for ("X", "Y", "Z") {
        %{ $box->{$_ . "COORD"} } = %{ $box->{$_} };
    }
}

sub getResRange {
    my ($resData) = $_[0];
    my ($i, @tmp, @RANGE, $start, $prev);

    @tmp = sort numerically keys %{ $resData };
    $i = 0;
    $start = $prev = -1;
    for $i (0 .. $#tmp) {
	if ($start == -1) {
	    $start = $tmp[$i];
	} elsif (($tmp[$i] - $prev) > 1) {
	    if (($start - $prev) > 0) {
		push @RANGE, ":${start}-${prev}";
	    } else {
		push @RANGE, ":${start}";
	    }
	    $start = $tmp[$i];
	}
	$prev = $tmp[$i];
    }
    if (($start - $prev) > 0) {
	push @RANGE, ":${start}-${prev}";
    } else {
	push @RANGE, ":${start}";
    }
    return \@RANGE;
}

sub showUsage {
    my ($usage) = "usage: $0 -b bgf file -a axis (non periodic surface dimension)" .
        "\nOptional parameters:" .
	"\n\t-m solute atoms (see selection info) -v solvent atoms (see selection info) -w save name " .
        "\n\t-r max distance from surface -d del r for binning -f cerius2 forcefield -t trajectory file " .
        "\n\t-s trajectory selection (see selection info) -l traj type (amber(default)|lammps) -c \"solv shell peaks\"\n" . 
	&ShowSelectionInfo;
    return $usage;
}

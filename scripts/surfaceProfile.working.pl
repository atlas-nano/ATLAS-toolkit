#!/usr/bin/perl

use strict;
use warnings;
no warnings "recursion";
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo AddMass createHeaders addHeader createBGF);
use General qw(STDev FileTester CoM GetSelections TrjSelections ShowSelectionInfo LoadFFs Round GetSigDigits);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox CreateLAMMPSTrj);
use BOX qw(GetBox CreateGrid);
use ManipAtoms qw(UnwrapAtoms ImageAtoms GetAtmList GetSolvent SplitAtomsByMol GetAtmData);

sub numerically { ($a<=>$b); }
sub writeShellCount;
sub writeShellData;
sub writeIFT;
sub writeDensity;
sub calcProfile;
sub findSolute;
sub storeSoluCoords;
sub calcIFT;
sub getCorrectShell;
sub saveShellData;
sub saveBGF;
sub saveLammpsShellTraj;
sub init;
sub boxConvert;
sub getResRange;
sub findSurface;
sub writeIFTFrameData;
sub openIFTFile;

my ($SOLUTE, $SOLVENT, $trjFile, $BGF, $BONDS, $rMax, $printStr, $axis);
my ($getSnapshot, $getByteOffset, $trjType, $LAMMPSOPTS, $delR, $offset);
my ($SELECT, $saveFile, $field, $DATA, $PARMS, $prec, $rounding, $i, $iftFile);
my ($HEADERS, $BBOX, $DIM, $calcType, $shellData, $dumpBGFs, $numProcs);

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
        &GetLammpsTrjType($SELECT, $trjFile, "ift", \%{ $LAMMPSOPTS });
        $field = "ift";
    }
    &openIFTFile($saveFile) if ($LAMMPSOPTS->{ISIFT});
    $printStr = "Calculating solvent <-> surface profile from $trjFile...";
    $getSnapshot->($BGF, $trjFile, $SELECT, $field, \&calcProfile, $printStr, $iftFile);
    close $iftFile if ($LAMMPSOPTS->{ISIFT});
}

if ($calcType == 1 or ! defined($SOLUTE)) {
    &writeDensity($DATA->{DENSITY}, $saveFile);
    &writeIFT($DATA->{IFT}, $saveFile) if ($LAMMPSOPTS->{ISIFT});
} elsif ($calcType == 2 and defined($SOLUTE)) {
    &writeShellData($DATA, $saveFile);
} elsif (defined($SOLUTE)) {
    &writeShellCount($DATA, $saveFile);
}

die "\n";

sub writeShellCount {
    my ($data, $prefix) = @_;
    my ($i, $saveName, @frames, $count, @shells, %shellData, $j);

    while ($data->{FRAMES} =~ /(\d+)/g) {
	push @frames, $1;
    }
    delete $data->{FRAMES};
    @shells = sort {($a cmp $b ) } keys %{ $data };

    print "Writing # molecules in shells...";
    $prefix =~ s/\.\w+$//;
    $saveName = "${prefix}_shells_count.dat";
    open SHELLCOUNT ,"> $saveName" or die "ERROR: Cannot create $saveName: $!\n";

    printf SHELLCOUNT "%-12s ","#TIMESTEP";
    for $i (@shells) {
	printf SHELLCOUNT "%8s ", $i;
    }
    print SHELLCOUNT "\n";

    for $i (@shells) {
	$count = 0;
	while ($data->{$i} =~ /(\d+)/g) {
	    $shellData{$i}{$frames[$count]} = $1;
	    $count++;
	}
    }

    for $i (0 .. $count -1) {
	printf SHELLCOUNT "%12d ", $frames[$i];
	for $j (@shells) {
	    printf SHELLCOUNT "%8.3f ", $shellData{$j}{$frames[$i]};
	}
	print SHELLCOUNT "\n";
   }
   close SHELLCOUNT;
   print "Done\n";
}

sub writeShellData {
    my ($data, $saveName) = @_;
    my (%SHELLMOLS, $shellName, $shellID, $i, $prefix, $count, $vacData);
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
		next if ($k == 0); # next if lasted less than 1ps
		$SHELLTRANS{$shellName} .= "$k ";
		$SHELLCOUNT{$shellName}++;
	    }
	}
    }

    $vacData .= "Group 1 Atoms " . scalar(keys %{ $SOLUTE }) . "\n";
    $vacData .= "1-" . scalar(keys %{ $SOLUTE }) . "\n";
    $count = 2;
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
	($k, $j) = getResRange($resData);

	$vacData .= "Group $count Atoms " . $k->{TOT} . "\n";
        $vacData .= $k->{ATOMS} . "\n";
	$count++;

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
	printf SHELLTRANS "%-20s %12.2f %8.2f %8d\n", $shellName, $j, $k, $SHELLCOUNT{$i};
	printf "%-20s %12.2f %8.2f %8d\n", $shellName, $j, $k, $SHELLCOUNT{$i};
    }
    close SHELLTRANS;

    $vacData = "Total groups $count\n$vacData";
    open VACGRPFILE, "> ${prefix}_grps.vac" or die "ERROR: Cannot create ${prefix}_grps.vac: $!\n";
    print VACGRPFILE $vacData;
    close VACGRPFILE;
}

sub writeIFT {
    my ($data, $prefix) = @_;
    my ($saveName, %sigma, $bin);
    my $atomic2gpcc = 1.66030;
    my $atm2dynepscm=1013250;
    my $A2cm=1e-8;
    my $atomic2bulkVol = 0.6023;
    my $nBins = $rMax/$delR;

    return if (! exists($data->{Pt}));

    $prefix =~ s/\.\w+$//;
    $saveName = "${prefix}_ift_profile.dat";
    print "Writing IFT Profile data to $saveName...";
    open IFTPROFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    $prec += 8;
    printf IFTPROFILE "%-${prec}s %12s %12s %12s %12s %12s\n","#BIN", "IFT", "Pn", "Pt", "DENSITY", "PRESSURE";
    for $bin (sort numerically keys %{ $data->{Pt} }) {
	$data->{Pn}{$bin} /= $data->{count};
	$data->{Pt}{$bin} /= $data->{count};
	$data->{Press}{$bin} /= $data->{count};
	$data->{density}{$bin} /= $data->{count};
	# now compute ift
	$sigma{$bin} = ($data->{Pn}{$bin}-$data->{Pt}{$bin})*$delR*$atm2dynepscm*$A2cm/2;
        printf IFTPROFILE "%-${prec}f %12.7f %12.5f %12.5f %12.5f %12.3f\n",($bin*$delR),$sigma{$bin},
			  $data->{Pn}{$bin},$data->{Pt}{$bin},$data->{density}{$bin},$data->{Press}{$bin};
    }
    close IFTPROFILE;
    print "Done\n";
}

sub writeDensity {
    my ($data, $saveName) = @_;
    my ($prefix, $rT, $frame, $count, $bin, $i1, $i2);
    my (%binData, $avg, $stdev, $tot, $rec);

    $prefix = $_[1];
    $prefix =~ s/\.\w+$//;
    $rT = 1;
    print "Writing data to $saveName...";
    for $frame (keys %{ $data }) {
	$count++;
	for $bin (keys %{ $data->{$frame} }) {
	    $binData{$bin} .= "$data->{$frame}{$bin} ";
	}
    }
    
    $count = 0;
    open OUTDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    printf OUTDATA "%-${prec}s %-8s %-5s %-8s %-8s %-8s\n","#", "AVG", "STDEV", "TOTAL", "PTS", "INT";
    for $bin (sort numerically keys %binData) {
	chop $binData{$bin};
	($avg, $stdev, $tot) = STDev($binData{$bin});
	$count += $avg;
	printf OUTDATA "%${prec}f %8.3f %5.3f %8d %8d %8.3f\n", ($bin * $delR), $avg, $stdev, $tot, $tot/$avg, $count;
    }
    close OUTDATA;
    print "Done\n";
}

sub calcProfile {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($count, $i, $j, $k, $tot,$offset, $dM, $minDist, $bin, $dist);
    my ($SOLATMS, $chain, $shellName, $shellID, $index, $soluArray);
    my ($SHELLTRJDATA, %SHELLCOUNT, %IFT, $CENTER, $MOL, $iftDATA);

    $dM = uc($axis) . "COORD";
    if ($trjType == 2) { #LAMMPS
        $frameNum = $ATOMS->{TIMESTEP}[0];
	$tot = $ATOMS->{"NUMBER OF ATOMS"}[0];
        $BOX = ConvertLammpsBox(\%{ $ATOMS->{"BOX BOUNDS"} });
        $ATOMS = $ATOMS->{ATOMS};
	if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
	    UnwrapAtoms($ATOMS,  $BOX, $LAMMPSOPTS->{scaled});
	}
    } elsif ($trjType == 1) {
        $BOX = ConvertAmberBox(\%{ $BOX });
	$tot = scalar(keys %{ $ATOMS });
    } 

    $soluArray = findSurface($ATOMS, $SOLUTE, $BOX, $DIM, $dM) if (defined($SOLUTE));

    for $i (keys %{ $SOLVENT }) {
	$MOL = ();
	$CENTER = ();
	$index = 0;
	for $j (keys %{ $SOLVENT->{$i} }) {
	    $ATOMS->{$j}{MASS} = 1;
	    $ATOMS->{$j}{MASS} = $BGF->{$j}{MASS} if (exists($BGF->{$j}{MASS}));
	    $MOL->{$j} = \%{ $ATOMS->{$j} };
	    for $k ("XCOORD", "YCOORD", "ZCOORD") {
		$CENTER->{$k} += $ATOMS->{$j}{$k} * $ATOMS->{$j}{MASS};
	    }
	    $index += $ATOMS->{$j}{MASS}; #total mass
        }
	for $k ("XCOORD", "YCOORD", "ZCOORD") {
	    $CENTER->{$k} /= $index;
	}
	$CENTER = ImageAtoms($MOL, $CENTER, $BOX);
	$SOLATMS = ();
	$SOLATMS = findSolute($soluArray, $CENTER, $BOX, $DIM, 2) if (defined($SOLUTE));
	$j = $index = 0;
	$minDist = $CENTER->{$dM};
	for $k (keys %{ $SOLATMS }) {
	    $dist = ($ATOMS->{$k}{XCOORD} - $CENTER->{XCOORD})**2 + 
		($ATOMS->{$k}{YCOORD} - $CENTER->{YCOORD})**2 +
		($ATOMS->{$k}{ZCOORD} - $CENTER->{ZCOORD})**2;
	    if (! $j or $minDist > $dist) {
		$minDist = $dist;
		$j = $k;
		$index = $ATOMS->{$j}{$dM};
	    }
	}
	$bin = sprintf("%${rounding}f",((sqrt($minDist)/$delR) + 1) * ($CENTER->{$dM} <=> $index));

	if ($calcType == 1)  {
	    $DATA->{DENSITY}{$frameNum}{$bin}++;
	    $IFT{$bin}{$i} = 1;
	} elsif ($calcType == 2 and defined($SOLUTE)) {
	    ($index, $shellName) = &saveShellData($shellData, $minDist, $frameNum, \%{ $DATA->{$i} });
	    $DATA->{$i}{SHELLS}[$index]{STOP} = $frameNum;
	    next if (! $dumpBGFs);
	    ($shellName, $index) = &getCorrectShell($shellData, $minDist);
	    $chain = chr(66 + $index);
	    for $k (keys %{ $MOL }) {
		$SHELLTRJDATA->{ATOMS}{$k} = $ATOMS->{$k};
		$SHELLTRJDATA->{BONDS}{$k} = $BONDS->{$k};
		$SHELLTRJDATA->{ATOMS}{$k}{CHAIN} = $chain;
	    }
	} elsif ($calcType == 3 and defined($SOLUTE)) {
	    ($shellName, $index) = &getCorrectShell($shellData, $minDist);
	    $SHELLCOUNT{$shellName}++;
	    for $k (keys %{ $MOL }) {
		$SHELLTRJDATA->{$shellName}{$k} = $ATOMS->{$k};
	    }
	}		
	$count++;
    }		 

    if ($calcType == 1) {
	if ($LAMMPSOPTS->{ISIFT}) {
	    for $bin (keys %IFT) {
		($iftDATA->{IFT}{$bin}, $iftDATA->{CELLVOL}, $iftDATA->{SLABVOL}) = 
		&calcIFT($ATOMS, $SOLVENT, \%{ $IFT{$bin} }, $BOX, $axis);
	    }
	    &writeIFTFrameData(\%{ $DATA->{IFT} }, $iftDATA, $frameNum, $iftFile);
	}
    } elsif ($calcType == 2 and $dumpBGFs and defined($SOLUTE)) {
        for $i (keys %{ $SOLUTE })  {
	    $SHELLTRJDATA->{ATOMS}{$i} = $ATOMS->{$i};
	    $SHELLTRJDATA->{BONDS}{$i} = $BONDS->{$i};
	}
	&saveBGF($SHELLTRJDATA, $BOX, $frameNum);
    } elsif ($calcType == 3 and defined($SOLUTE)) {
	&saveLammpsShellTraj(\%{ $SHELLTRJDATA }, \%{ $shellData }, $frameNum, $tot, $BOX, $axis);
	$DATA->{FRAMES} .= "$frameNum ";
	for $i (keys %SHELLCOUNT) {
	    $DATA->{$i} .= "$SHELLCOUNT{$i} ";
	}
    }
    
    %SHELLCOUNT = ();
    $SHELLTRJDATA = ();
    $ATOMS = ();
    $SOLATMS = ();
    $soluArray = ();
    $iftDATA = ();
}

sub writeIFTFrameData {
    my ($ift, $data, $frame, $iftFHandle) = @_;
    my ($invol, $cellVol, $bin);
    my $atomic2gpcc = 1.66030;
    my $atm2dynepscm=1013250;
    my $atomic2bulkVol = 0.6023;
    my $A2cm=1e-8;
    my $nBins = $rMax/$delR;

    $ift->{count}++;
    $invol = 1/$data->{SLABVOL};
    $cellVol = $data->{CELLVOL};
    for $bin (keys %{ $data->{IFT} }) {
        $ift->{Pn}{$bin} += -$invol*$data->{IFT}{$bin}{Pn};
        $ift->{tot}{Pn} += -$invol*$data->{IFT}{$bin}{Pn};
        $ift->{Pt}{$bin} += -$invol*$data->{IFT}{$bin}{Pt}/2;
        $ift->{tot}{Pt} += -$invol*$data->{IFT}{$bin}{Pt}/2;
        $ift->{density}{$bin} += $invol*$atomic2bulkVol*$data->{IFT}{$bin}{mass};
        $ift->{tot}{density} += $invol*$atomic2bulkVol*$data->{IFT}{$bin}{mass};
        $ift->{Press}{$bin} += $data->{IFT}{$bin}{P}/(-3*$cellVol);
        $ift->{tot}{Press} += $data->{IFT}{$bin}{P}/(-3*$cellVol);
    }
    $ift->{tot}{sigma} = ($ift->{tot}{Pn}-$ift->{tot}{Pt})*$delR*$atm2dynepscm*$A2cm/2;
    printf $iftFHandle "%-12d %12.7f %12.5f %12.5f %12.5f %12.3f\n",
    $frame,($ift->{tot}{sigma}/$ift->{count}),($ift->{tot}{Pn}/$ift->{count}),
    ($ift->{tot}{Pt}/$ift->{count}),($ift->{tot}{density}/$ift->{count}),
    ($ift->{tot}{Press}/$ift->{count});
}

sub findSurface {
    my ($atoms, $solute, $box, $store, $dM) = @_;
    my ($MOL, $CENTER, @dims, $i, $j, $k, $i1, $i2, $dist, $count);
    my ($atomC, $surface, $coordGrid, $soluCOORDS, $soluArray, $tol);

    $MOL = GetAtmData($atoms, $solute);
    $CENTER = CoM($MOL);
    $tol = 0;
    @dims = ("XCOORD", "YCOORD", "ZCOORD");
    for $i (@dims) {
        $box->{$i}{CENTER} = $box->{$i}{len}/2;
        $box->{$i}{hi} = $box->{$i}{len};
        $box->{$i}{lo} = 0;
        for $j (keys %{ $atoms }) {
            $atoms->{$j}{$i} += ($box->{$i}{CENTER} - $CENTER->{$i});
        }
    }
    for $i (keys %{ $solute }) {
        $j = ($atoms->{$i}{$store->[0]});
        $k = ($atoms->{$i}{$store->[1]});
        $i1 = Round($j, $tol);
        $i2 = Round($k, $tol);
        $dist = $atoms->{$i}{$dM};
        if (! exists($coordGrid->{$i1}{$i2})) {
            $coordGrid->{$i1}{$i2}{MAX}{dist} = $dist;
            $coordGrid->{$i1}{$i2}{MAX}{atom} = $i;
            $coordGrid->{$i1}{$i2}{MIN}{dist} = $dist;
            $coordGrid->{$i1}{$i2}{MIN}{atom} = $i;
        } elsif ($dist > $coordGrid->{$i1}{$i2}{MAX}{dist}) {
            $atomC = $coordGrid->{$i1}{$i2}{MAX}{atom};
            delete $surface->{$atomC} if ($coordGrid->{$i1}{$i2}{MIN}{atom} != $atomC);
            $coordGrid->{$i1}{$i2}{MAX}{dist} = $dist;
            $coordGrid->{$i1}{$i2}{MAX}{atom} = $i;
        } elsif ($dist < $coordGrid->{$i1}{$i2}{MIN}{dist}) {
            $atomC = $coordGrid->{$i1}{$i2}{MIN}{atom};
            delete $surface->{$atomC} if ( $coordGrid->{$i1}{$i2}{MAX}{atom} != $atomC);
            $coordGrid->{$i1}{$i2}{MIN}{dist} = $dist;
            $coordGrid->{$i1}{$i2}{MIN}{atom} = $i;
        } else {
            next;
        }
        $surface->{$i}{1} = $j;
        $surface->{$i}{2} = $k;
    }
    for $i (keys %{ $surface }) {
        $soluCOORDS->{ $surface->{$i}{1} }{VALS}{ $surface->{$i}{2} }{VALS} = $i;
        $count++;
    }
    $soluArray = storeSoluCoords($soluCOORDS, $box, $store, 2);

    return $soluArray;
}


sub findSolute {
    my ($soluArray, $solvent, $box, $axes, $axisLen) = @_;
    my ($i1, $i2, $s1, $s2, $i, $j, $k, $soluList);
    
    $i1 = int(($solvent->{ $axes->[0] } - $box->{ $axes->[0] }{lo})/$axisLen);
    $i2 = int(($solvent->{ $axes->[1] } - $box->{ $axes->[1] }{lo})/$axisLen);

    if ($i1 > 0 and $i1 < $#{ $soluArray }) {
	@{ $s1 } = (-1 + $i1, $i1, 1 + $i1);
    } elsif ($i1 == 0) {
	@{ $s1 } = ($i1, 1 + $i1, 2 + $i1);
    } elsif ($i1 > 2) {
	@{ $s1 } = (-2 + $i1, -1 + $i1, $i1);
    } else {
	@{ $s1 } = ($i1);
    }

    if ($i2 > 0 and $i2 < $#{ $soluArray }) {
	@{ $s2 } = (-1 + $i2, $i2, 1 + $i2);
    } elsif ($i2 == 0) {
	@{ $s2 } = ($i2, 1 + $i2, 2 + $i2);
    } elsif ($i2 > 2) {
	@{ $s2 } = (-2 + $i2, -1 + $i2, $i2);
    } else {
	@{ $s2 } = ($i2);
    }

    for $i (@{ $s1 }) {
	for $j (@{ $s2 }) {
	    for $k (@{ $soluArray->[$i][$j] }) {
		$soluList->{$k} = 1;
	    }
	}
    }
    
    return $soluList;
}

sub storeSoluCoords {
    my ($soluCoords, $box, $axes, $axisLen) = @_;
    my ($i, $j, $soluArray, $i1, $i2);

    for $i (keys %{ $soluCoords }) {
	$i1 = int(($i - $box->{ $axes->[0] }{lo})/$axisLen);
	next if ($i1 < 0);
	for $j (keys %{ $soluCoords->{$i}{VALS} } ) {
	    $i2 = int(($j - $box->{ $axes->[1] }{lo})/$axisLen);
	    next if ($i2 < 0);
	    push @{ $soluArray->[$i1][$i2] }, $soluCoords->{$i}{VALS}{$j}{VALS};
	    $soluCoords->{$i}{INDEX} = $i1;
	    $soluCoords->{$i}{VALS}{$j}{INDEX} = $i2;
	}
    }
    
    return $soluArray;
}

sub calcIFT {
    my ($atoms, $mols, $iftData, $BOX, $axis) = @_;
    my ($i, @atmList, $atomC, $vol, @tmp, %map);
    my (%IFT, $slab_vol, $tot_vol, $stress, $j);

    $axis = uc($axis);
    @tmp = ("X", "Y", "Z");
    %map = ("X" => 0, "Y" => 1, "Z" => 2);

    for $i (@tmp) {
        $vol->{$i} = $BOX->{"${i}COORD"}{hi} - $BOX->{"${i}COORD"}{lo};
    }
    $tot_vol = $vol->{X}*$vol->{Y}*$vol->{Z}; #cell volume
    $vol->{$axis} = $delR;
    $slab_vol = $vol->{X}*$vol->{Y}*$vol->{Z}; #slab volume

    for $i (keys %{ $iftData }) {
        @atmList = keys %{ $mols->{$i} };
        for $atomC (@atmList) {
	    for $j (@tmp) {
		$stress = $atoms->{$atomC}{"S" . $j . $j};
		$IFT{P} += $stress; # slab scalar total pressure*vol (sum of ensembles)
		if ($j eq $axis) {
		    $IFT{Pn} += $stress; # slab normal pressure*vol
		} else {
		    $IFT{Pt} += $stress; # slab tangential pressure*vol
		}
	    }
	    $IFT{mass} += $atoms->{$atomC}{MASS}; # slab mass
	}
    }

    return (\%IFT, $tot_vol, $slab_vol);
}

sub getCorrectShell {
    my ($shellData, $minDist) = @_;
    my ($k, $shellName, $shellID);

    for $k (@{ $shellData }) {
        next if ($minDist > $k->{END});
        $shellName = $k->{NAME};
        $shellID = $k->{ID};
        last;
    }
    ($shellName, $shellID) = ("bulk", 0) if (! defined($shellID));

    return ($shellName, $shellID);
}

sub saveShellData {
    my ($shellData, $minDist, $frameNum, $DATA) = @_;
    my ($shellID, $k, $shellName, $index);

    ($shellName, $shellID) = getCorrectShell($shellData, $minDist);

    $DATA->{TRAJ}{$frameNum} = (
				{
				    "NAME" => $shellName,
				    "DIST" => $minDist,
				}
			       );
    if (! exists($DATA->{SHELLS})) {
        $DATA->{SHELLS}[0]{NAME} = $shellName;
        $DATA->{SHELLS}[0]{START} = $frameNum;
        $DATA->{SHELLS}[0]{ID} = $shellID;
        $index = 0;
    } else {
        $index = $#{ $DATA->{SHELLS} };
        if ($DATA->{SHELLS}[$index]{NAME} ne $shellName) {
	    $index++;
	    $DATA->{SHELLS}[$index]{ID} = $shellID;
	    $DATA->{SHELLS}[$index]{NAME} = $shellName;
	    $DATA->{SHELLS}[$index]{START} = $frameNum;
	}
    }

    return ($index, $shellName);
}

sub saveBGF {
    my ($shellData, $BOX, $frame) = @_;
    my ($saveName, $HEADERS, $ATOMS, $BONDS, $uBGF, $i, $j);

    $saveName = $saveFile;
    $saveName =~ s/\.\w+$//;
    $saveName .= "_ts${frame}.bgf";
    for $i (keys %{ $BGF }) {
	for $j (keys %{ $shellData->{ATOMS}{$i} }) {
	    $BGF->{$i}{$j} = $shellData->{ATOMS}{$i}{$j};
	}
    }
    %{ $uBGF } = %{ $BGF };
    $BONDS = $shellData->{BONDS};
    $HEADERS = createHeaders($BOX, $saveName);
    &addHeader($uBGF, $HEADERS);
    &createBGF($uBGF, $BONDS, $saveName);
}

sub saveLammpsShellTraj {
    my ($trjData, $shellData, $timeStep, $solvAtms, $BOX, $axis) = @_;
    my ($i, $trjName, $prefix, $headers, $shellName, $prevShell, @tmp, %map);

    $prefix = $saveFile;
    $prefix =~ s/\.\w+$//;
    $prevShell = $offset;
    $axis = uc($axis);
    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    %map = ("X" => 0, "Y" => 1, "Z" => 2);

    $shellData->[ $#{ $shellData } ]{END} = $BOX->{"${axis}COORD"}{hi}/2;
    $headers->{TIMESTEP}[0] = $timeStep;
    $headers->{"NUMBER OF ATOMS"}[0] = $solvAtms;
    for $i (0 .. $#tmp) {
        $headers->{"BOX BOUNDS"}[$i]{lo} = $BOX->{$tmp[$i]}{lo};
        $headers->{"BOX BOUNDS"}[$i]{hi} = $BOX->{$tmp[$i]}{hi};
    }

    for $i (0 .. $#{ $shellData }) {
	$shellName = $shellData->[$i]{NAME};
        next if (! exists($trjData->{$shellName}));
	$prevShell = $shellData->[$i-1]{END} if ($i > 0);
	$headers->{"BOX BOUNDS"}[$map{$axis}]{lo} = $prevShell;
	$headers->{"BOX BOUNDS"}[$map{$axis}]{hi} = $shellData->[$i]{END};
	$trjName = "${prefix}_${shellName}.lammpstrj";
	open LAMMPSTRJ, ">> $trjName" or die "ERROR: Cannot write to $trjName: $!\n";
	&CreateLAMMPSTrj($trjData->{$shellName}, $headers, \*LAMMPSTRJ);
	close LAMMPSTRJ;
    }

}

sub init {
    my (%OPTS, @tmp, $i, $j, $list, $solvTmp, $BOX, $nBins, $count);
    my ($solSel, $solvSel, $bgfFile, $tSel, $usage, $FF, $solvShells);
    my ($prefix, $fileName);

    getopt('mvbtdrlswafcogp',\%OPTS);

    $usage = &showUsage;

    for ($OPTS{a},$OPTS{b}) {
        die "$usage" if (! defined($_));
    }

    print "Initializing...";
    ($solSel,$solvSel,$trjFile,$bgfFile,$rMax,$delR,$trjType,$tSel,$saveFile,$axis,$FF, $solvShells, $offset, $dumpBGFs,$numProcs) =
        ($OPTS{m},$OPTS{v},$OPTS{t},$OPTS{b},$OPTS{r},$OPTS{d},$OPTS{l},$OPTS{s},$OPTS{w},$OPTS{a},$OPTS{f},$OPTS{c},$OPTS{o},$OPTS{g},$OPTS{p});

    $numProcs = 1 if (! defined($numProcs) or $numProcs !~ /^\d+$/);
    $offset = 0 if (! defined($offset) or $offset !~ /^\d+\.?\d*$/);
    if ($axis =~ /(X|Y|Z)/i) {
	$axis = lc($1);
    } else {
	die "ERROR: Expected X|Y|Z for axis. Got \"$axis\"!\n";
    }

    @{ $DIM } = grep {!/^$axis/i} ("XCOORD","YCOORD","ZCOORD");

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

    $delR = 0.1 if (! defined($delR) or $delR !~ /^\d+\.?\d*/);

    if (! defined($saveFile)) {
        $saveFile = basename($bgfFile);
        $saveFile =~ s/\.\w+$//;
	$saveFile .= "_density.dat";
    }

    print "Done\nParsing BGF file $bgfFile...";
    ($BGF, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
    if (defined($SOLUTE)) {
	$SOLUTE = GetAtmList($SOLUTE, $BGF);
	die "ERROR: Solute atoms does not correspond to any in BGF file!\n"
	    if (! keys %{ $SOLUTE });
	for $i (keys %{ $SOLUTE}) {
	    $BGF->{$i}{IS_SOLUTE} = 1;
	}
    }
    if (defined($solvTmp)) {
	$SOLVENT = GetAtmList($solvTmp, $BGF);
    } else {
	$SOLVENT = $BGF;
    }
    die "ERROR: Solvent atoms does not correspond to any in BGF file!\n"
        if (! keys %{ $SOLVENT });
    $SOLVENT = SplitAtomsByMol($BGF, $SOLVENT);
    for $i (keys %{ $SOLVENT }) {
	for $j (keys %{ $SOLVENT->{$i} }) {
	    $BGF->{$j}{IS_SOLVENT} = 1;
	}
    }
    $dumpBGFs = 0 if (! defined($dumpBGFs));
    $dumpBGFs = 1 if ($dumpBGFs =~ /(1|yes)/i);

    $BOX = GetBox($BGF, undef, $HEADERS);
    %{ $BBOX } = %{ $BOX };
    &boxConvert($BBOX);
    if (! defined($rMax) or ($rMax !~ /^\d+\.?\d*/)) {
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
    $calcType = 3 if ($calcType == 2 and $offset > 0);
    if ($calcType == 3) {
        $prefix = basename($saveFile);
        $prefix =~ s/\.\w+$//;

        for $i (@{ $shellData }) {
            $fileName = "${prefix}_" . $i->{NAME} . ".lammpstrj";
            #eval('my $' . $i->{NAME});
            open TESTTRJ, "> $fileName" or die "ERROR: Cannot create $fileName: $!\n";
            close TESTTRJ;
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

sub openIFTFile {
   my ($prefix) = $_[0];
   $prefix =~ s/\.\w+$//;
   my $saveName = "${prefix}_ift.dat";
   open $iftFile, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
   printf $iftFile "#%-11s %12s %12s %12s %12s %12s\n","TIMESTEP", "IFT", "Pn", "Pt", "DENSITY", "PRESSURE";
}

sub boxConvert {
    my ($box) = $_[0];

    for ("X", "Y", "Z") {
        %{ $box->{$_ . "COORD"} } = %{ $box->{$_} };
    }
}

sub getResRange {
    my ($resData) = $_[0];
    my ($i, $j, @tmp, @tmp1, @RANGE, $start, $end, $prev, %ATMRANGE, %RESDATA);

    @tmp = sort numerically keys %{ $resData };
    $i = 0;
    for $i (keys %{ $BGF }) {
	next if (exists($RESDATA{$BGF->{$i}{RESNUM}}));
        $RESDATA{$BGF->{$i}{RESNUM}} = $BGF->{$i}{MOLECULE};
    }
    $start = $prev = -1;
    $ATMRANGE{count} = $ATMRANGE{TOT} = 0;
    for $i (0 .. $#tmp) {
	if ($start == -1) {
	    $start = $tmp[0];
	} elsif (($tmp[$i] - $prev) > 1) {
	    if (($prev - $start) > 0) {
		push @RANGE, ":${start}-${prev}";
	        @tmp1 = sort numerically keys %{ $RESDATA{$start} };
		$ATMRANGE{start} = $tmp1[0];
		@tmp1 = sort numerically keys %{ $RESDATA{$prev} };
		$ATMRANGE{end} = $tmp1[$#tmp1];
		$ATMRANGE{ATOMS} .= $ATMRANGE{start} . " - " . $ATMRANGE{end} . ", ";
		$ATMRANGE{TOT} += ($ATMRANGE{end}-$ATMRANGE{start}) + 1;
		$ATMRANGE{count}++;
		if ($ATMRANGE{count}  % 5 == 0) {
		    $ATMRANGE{ATOMS} = substr($ATMRANGE{ATOMS}, 0, -2);
		    $ATMRANGE{ATOMS} .= "\n";
		}
	    } else {
		push @RANGE, ":${start}";
                @tmp1 = sort numerically keys %{ $RESDATA{$start} };
                $ATMRANGE{start} = $tmp1[0];
                $ATMRANGE{end} = $tmp1[$#tmp1];
                $ATMRANGE{ATOMS} .= $ATMRANGE{start} . " - " . $ATMRANGE{end} . ", ";
		$ATMRANGE{TOT} += ($ATMRANGE{end}-$ATMRANGE{start}) + 1;
		$ATMRANGE{count}++;
                if ($ATMRANGE{count}  % 5 == 0) {
	            $ATMRANGE{ATOMS} = substr($ATMRANGE{ATOMS}, 0, -2);
                    $ATMRANGE{ATOMS} .= "\n";
                }		
	    }
	    $start = $tmp[$i];
	}
	$prev = $tmp[$i];
    }
    if (($start - $prev) > 0) {
	push @RANGE, ":${start}-${prev}";
        @tmp = sort numerically keys %{ $RESDATA{$start} };
        $ATMRANGE{start} = $tmp[0];
        @tmp = sort numerically keys %{ $RESDATA{$prev} };
        $ATMRANGE{end} = $tmp[$#tmp];
        $ATMRANGE{ATOMS} .= $ATMRANGE{start} . " - " . $ATMRANGE{end} . ", ";
	$ATMRANGE{TOT} += ($ATMRANGE{end}-$ATMRANGE{start}) + 1;
	$ATMRANGE{count} += ($prev - $start) + 1;
    } else {
	push @RANGE, ":${start}";
	@tmp = sort numerically keys %{ $RESDATA{$start} };
	$ATMRANGE{start} = $tmp[0];
        $ATMRANGE{end} = $tmp[$#tmp];
        $ATMRANGE{ATOMS} .= $ATMRANGE{start} . " - " . $ATMRANGE{end} . ", ";
	$ATMRANGE{TOT} += ($ATMRANGE{end}-$ATMRANGE{start}) + 1;
	$ATMRANGE{count}++;
    }

    $ATMRANGE{ATOMS} = substr($ATMRANGE{ATOMS}, 0, -2);
    return (\%ATMRANGE, \@RANGE);
}

sub showUsage {
    my ($usage) = "usage: $0 -b bgf file -a axis (non periodic surface dimension)" .
        "\nOptional parameters:" .
	"\n\t-m solute atoms (see selection info) -v solvent atoms (see selection info) -w save name " .
        "\n\t-r max distance from surface -d del r for binning -f cerius2 forcefield -t trajectory file " .
        "\n\t-s trajectory selection (see selection info) -l traj type (amber(default)|lammps) -c \"solv shell peaks\"" . 
	"\n\t-o offset: integer. The distance from the solute surface-solvent distance" .
	"\n\t-g dump bgf: Flag to dump bgf files every timestep based on shell. Default 0\n" .
	&ShowSelectionInfo;
    return $usage;
}

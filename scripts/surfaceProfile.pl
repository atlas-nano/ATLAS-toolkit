#!/usr/bin/perl

use strict;
use warnings;
no warnings "recursion";
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo AddMass createHeaders addHeader createBGF);
use General qw(STDev FileTester CoM GetSelections TrjSelections ShowSelectionInfo LoadFFs Round GetSigDigits DotProduct);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox CreateLAMMPSTrj);
use BOX qw(GetBox CreateGrid);
use ManipAtoms qw(UnwrapAtoms ImageAtoms GetAtmList GetSolvent SplitAtomsByMol GetAtmData);
use Bulk qw(CalcMolOrientCorrelation CalcDipole CalcDielectricConstant CalcDiffusionConstant GetOrientOpt
		      WriteData WriteStats SetOpts);
use Math::MatrixReal;
use constant q2C => 1.602176487E-19; # q = 1.602176487 x 10^-19 C
use constant a2m => 10E-10; # 1 A = 10^-10 m

sub numerically { ($a<=>$b); }
sub acos;
sub writeShellCount;
sub writeShellData;
sub writeIFT;
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
sub centerCell;
sub calcDensity;
sub calcMolOrient;
sub calcChargeDensity;
sub writeBinData;
sub calcElecPot;

my ($SOLUTE, $SOLVENT, $trjFile, $BGF, $BONDS, $rMax, $printStr, $axis, $ORIENT, $tStep, $sysTemp);
my ($SELECT, $getSnapshot, $getByteOffset, $trjType, $LAMMPSOPTS, $delR, $offset, $tStart);
my ($HEADERS, $savePrefix, $field, $DATA, $PARMS, $prec, $rounding, $i, $BULK, $oPts);
my ($BBOX, $DIM, $calcType, $shellData, $dumpBGFs, $numProcs, $iftFile, $COMstart);

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
    &openIFTFile($savePrefix) if ($LAMMPSOPTS->{ISIFT} and ($calcType == 1 or ! defined($SOLUTE)));
    $printStr = "Calculating atom slab density/ift profile from $trjFile...";
    $printStr = "Calculating solvent <-> surface profile from $trjFile..." if (defined($SOLUTE));
    $getSnapshot->($BGF, $trjFile, $SELECT, $field, \&calcProfile, $printStr, $iftFile);
    close $iftFile if ($LAMMPSOPTS->{ISIFT} and ($calcType == 1 or ! defined($SOLUTE)));
}

if ($calcType == 1 and $LAMMPSOPTS->{ISIFT}) {
    &writeIFT($DATA->{IFT}, $savePrefix);
} elsif ($calcType == 1 or ! defined($SOLUTE)) {
    &writeBinData($DATA->{DENSITY}, "atomic_density", $savePrefix);
    &writeBinData($DATA->{ORIENT}{a}, "mol_alpha_orientation", $savePrefix);
    &writeBinData($DATA->{ORIENT}{t}, "mol_theta_orientation", $savePrefix);
    &writeBinData($DATA->{CHARGE_DENSITY}, "charge_density", $savePrefix);
    print "Calculating electrostatic potential...";
    $DATA->{ELEC_POT} = calcElecPot($DATA->{CHARGE_DENSITY});
    print "Done\n";
    &writeBinData($DATA->{ELEC_POT}, "electrostatic_potential", $savePrefix);
} elsif ($calcType == 2 and defined($SOLUTE)) {
    &writeShellData($DATA, $savePrefix);
} elsif (defined($SOLUTE)) {
    &writeShellCount($DATA, $savePrefix);
    print "Writing data...";
    for $i (grep { !/VOLUME|TEMPERATURE/i } keys %{ $BULK }) {
	&SetOpts(\%{ $BULK->{$i} });
    	&WriteData($BULK->{$i}, "${savePrefix}_${i}", 1);
    }
    $i = "MOL_CORRELATION";
    &SetOpts(\%{ $BULK->{$i} });
    &WriteData($BULK->{$i}, "${savePrefix}_${i}", 0);
    print "Done\n";
    &WriteStats($BULK, "${savePrefix}_stats.dat", 1);
}

sub writeShellCount {
    my ($data, $prefix) = @_;
    my ($i, $saveName, @frames, $count, @shells, %shellData, $j);

    while ($data->{FRAMES} =~ /(\d+)/g) {
	push @frames, $1;
    }
    delete $data->{FRAMES};
    @shells = sort {($a cmp $b ) } keys %{ $data };

    print "Writing # molecules in shells...";
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
    my ($data, $prefix) = @_;
    my (%SHELLMOLS, $shellName, $shellID, $i, $saveName, $count, $vacData);
    my (%SHELLTRANS, $j, $resID, $resData, $k, $tot, %SHELLCOUNT, $resRange);

    for $i (keys %{ $data }) {
	if ($#{ $data->{$i}{SHELLS} } == 0) { # stayed in one shell entire time
	    ($shellName, $shellID) = ($data->{$i}{SHELLS}[0]{NAME}, $data->{$i}{SHELLS}[0]{ID});
	    $SHELLMOLS{$shellName}{$i} = $data->{$i}{TRAJ};
	} else { # moved in between shells
	    $SHELLMOLS{trans}{$i} = $data->{$i}{TRAJ};
	    for $j (1 .. $#{ $data->{$i}{SHELLS} }) {
		$shellName = $data->{$i}{SHELLS}[$j]{NAME} . "_" . $data->{$i}{SHELLS}[($j - 1)]{NAME};
		$k = ($data->{$i}{SHELLS}[$j]{STOP} - $data->{$i}{SHELLS}[$j]{START}) * $delR/1000;
		next if ($k == 0); # next if lasted less than 1ps
		$SHELLTRANS{$shellName}{x1} += $k;
		$SHELLTRANS{$shellName}{x2} += $k**2;
		$SHELLCOUNT{$shellName}++;
	    }
	}
    }

    $saveName = "${prefix}_shelltrans.dat";
    print "Writing shell transition statistics to $saveName\n";
    open SHELLTRANSFILE, "> $saveName" or die "ERROR: Cannot write to $saveName: $!\n";
    for $i (keys %SHELLTRANS) {
        $j = $SHELLTRANS{$i}{x1}/$SHELLCOUNT{$i};
        $k = ($SHELLTRANS{$i}{x2}/$SHELLCOUNT{$i}) - $j**2;
        $tot = $SHELLCOUNT{$i};
        $shellName = $i;
        $shellName =~ s/_/ \-> /;
        printf SHELLTRANSFILE "%-20s %12.2f %8.2f %8d\n", $shellName, $j, $k, $SHELLCOUNT{$i};
        printf "%-20s %12.2f %8.2f %8d\n", $shellName, $j, $k, $SHELLCOUNT{$i};
    }
    close SHELLTRANSFILE;

    $vacData = "Group 1 Atoms " . scalar(keys %{ $SOLUTE }) . "\n";
    $vacData .= "1-" . scalar(keys %{ $SOLUTE }) . "\n";
    $count = 1;
    for $shellName (sort {$a cmp $b } keys %SHELLMOLS) {
	print "Writing molecules occupying $shellName...";
	system("mkdir -p ${prefix}/${shellName}");
	$resData = ();
	for $i (keys %{ $SHELLMOLS{$shellName} }) {
	    for $j (keys %{ $SOLVENT->{$i} }) {
		$resID = $BGF->{$j}{RESNUM};
		last;
	    }
	    $resData->{$i} = $resID;
	    $saveName = "${prefix}/${shellName}/${resID}.dat";
	    open PROFDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
	    for $j (sort numerically keys %{ $SHELLMOLS{$shellName}{$i} }) {
		print PROFDATA $SHELLMOLS{$shellName}{$i}{$j}{DIST} .  "\n";
	    }
	    close PROFDATA;
	}
	($k, $j) = getResRange($resData, $SOLVENT);

	$count++;
	$vacData .= "Group $count Atoms $k\n";
	$saveName = "${prefix}_${shellName}_mols.dat";
	open SHELLDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
	$k = 0;
	for $resRange (@{ $j }) { 
	    $vacData .= "$resRange, ";
	    print SHELLDATA " $resRange";
	    $k++;
	    if ($k == 10) {
		$vacData =~ s/\, $//;
		$vacData .= "\n";
		$k = 0;
	    }
	}
	$vacData .= "\n" if ($k > 0);
	close SHELLDATA;
	print "Done\n";
    }
    print "Writing VAC group file to ${prefix}_grps.vac...";
    $vacData = "Total groups $count\n$vacData";
    open VACGRPFILE, "> ${prefix}_grps.vac" or die "ERROR: Cannot create ${prefix}_grps.vac: $!\n";
    print VACGRPFILE $vacData;
    print "Done\n";
    undef($data);
}

sub writeIFT {
    my ($data, $prefix) = @_;
    my ($saveName, %sigma, $bin, $myPrec);
    my $atomic2gpcc = 1.66030;
    my $atm2dynepscm=1013250;
    my $A2cm=1e-8;
    my $atomic2bulkVol = 0.6023;
    my $nBins = $rMax/$delR;

    return if (! exists($data->{Pt}));

    $saveName = "${prefix}_ift_profile.dat";
    print "Writing IFT Profile data to $saveName...";
    open IFTPROFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    $myPrec = $prec;
    $myPrec += 7;
    printf IFTPROFILE "#%-${myPrec}s %12s %12s %12s %12s %12s\n","BIN", "IFT", "Pn", "Pt", "DENSITY", "PRESSURE";
    for $bin (sort numerically keys %{ $data->{Pt} }) {
	$data->{Pn}{$bin} /= $data->{count};
	$data->{Pt}{$bin} /= $data->{count};
	$data->{Press}{$bin} /= $data->{count};
	$data->{density}{$bin} /= $data->{count};
	# now compute ift
	$sigma{$bin} = ($data->{Pn}{$bin}-$data->{Pt}{$bin})*$delR*$atm2dynepscm*$A2cm/2;
        printf IFTPROFILE "%-${myPrec}f %12.7f %12.5f %12.5f %12.5f %12.3f\n",($bin*$delR),$sigma{$bin},
			  $data->{Pn}{$bin},$data->{Pt}{$bin},$data->{density}{$bin},$data->{Press}{$bin};
    }
    close IFTPROFILE;
    print "Done\n";
}

sub writeBinData {
    my ($data, $suffix, $prefix) = @_;
    my ($saveName, $frame, $count, $bin, $i1, $i2);
    my (%binData, $sigma, $x1, $x2, $tot);

    $x1 = $x2 = ();
    $sigma = $count = 0;
    $saveName = "${prefix}_${suffix}.dat";
    print "Writing data to $saveName...";
    for $frame (keys %{ $data }) {
	$count++;
	for $bin (keys %{ $data->{$frame} }) {
	    $x1->{$bin} += $data->{$frame}{$bin};
	    $x2->{$bin} += $data->{$frame}{$bin}**2;
	}
    }
    for $bin (keys %{ $x1 } ) {
	$x1->{$bin} /= $count;
	$x2->{$bin} /= $count;
	$sigma = ($x2->{$bin} - $x1->{$bin}**2);
        $sigma = 0 if ($sigma < 0.0000001);
	$sigma = sqrt($sigma);
	$binData{$bin} = {"AVG" => $x1->{$bin}, "STDEV" => $sigma, "TOT" => ($x1->{$bin} * $count)};
    }
    
    open OUTDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    printf OUTDATA "%-8s %-8s %-5s %-8s %-8s %-8s\n","#BIN", "AVG", "STDEV", "TOTAL", "PTS", "INT";
    for $bin (sort numerically keys %binData) {
	$tot += $binData{$bin}{AVG};
	printf OUTDATA "%${prec}f %8.3f %5.3f %8d %8d %8.3f\n", ($bin * $delR), $binData{$bin}{AVG}, 
	               $binData{$bin}{STDEV}, $count, $binData{$bin}{TOT}, $tot;
    }
    close OUTDATA;
    print "Done\n";
}

sub calcProfile {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($count, $i, $j, $k, $tot,$offset, $dM, $minDist, $bin, $dist);
    my ($SOLATMS, $chain, $shellName, $shellID, $index, $soluArray);
    my ($SHELLTRJDATA, %SHELLCOUNT, %IFT, $CENTER, $MOL, $iftDATA);
    my ($slab_vol, $cell_vol, $delT, $shell_vol);

    $dM = uc($axis) . "COORD";
    $cell_vol = $slab_vol = 1;
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

    $tStart = $frameNum if (! defined($tStart));
    $delT = $tStep * ($frameNum - $tStart); #time in ps
    $BULK->{VOLUME}{tStep}{$delT} = $BOX->{XCOORD}{len}*$BOX->{YCOORD}{len}*$BOX->{ZCOORD}{len};
    $BULK->{TEMPERATURE}{tStep}{$delT} = $sysTemp;

    &centerCell($ATOMS, $SOLUTE, $SOLVENT, $BOX, $frameNum);
    $soluArray = findSurface($BGF, $SOLUTE, $BOX, $DIM, $dM, $LAMMPSOPTS) if (defined($SOLUTE));
    for $i ("XCOORD", "YCOORD", "ZCOORD") {
        $cell_vol *= $BOX->{$i}{len};
        if ($i =~ /^$axis/i) {
            $slab_vol *= $delR;
        } else {
            $slab_vol *= $BOX->{$i}{len};
        }
    }

    for $i (keys %{ $SOLVENT }) {
	$MOL = GetAtmData($BGF, $SOLVENT->{$i});
	$CENTER = CoM($MOL);
	$CENTER = ImageAtoms($MOL, $CENTER, $BOX);
	$SOLATMS = ();
	$SOLATMS = findSolute($soluArray, $CENTER, $BOX, $DIM, 2) if (defined($SOLUTE));
	$j = $index = 0;
	$minDist = $CENTER->{$dM}**2;
	for $k (keys %{ $SOLATMS }) {
	    $dist = ($BGF->{$k}{XCOORD} - $CENTER->{XCOORD})**2 + 
		($BGF->{$k}{YCOORD} - $CENTER->{YCOORD})**2 +
		($BGF->{$k}{ZCOORD} - $CENTER->{ZCOORD})**2;
	    if (! $j or $minDist > $dist) {
		$minDist = $dist;
		$j = $k;
		$index = $BGF->{$j}{$dM};
	    }
	}
        $minDist = sqrt($minDist);
	$bin = sprintf("%${rounding}f",(($minDist/$delR) + 1) * ($CENTER->{$dM} <=> $index));

	if ($calcType == 1)  {
	    $DATA->{DENSITY}{$frameNum}{$bin}++;
	    $IFT{$bin}{$i} = $CENTER;
	} elsif ($calcType == 2 and defined($SOLUTE)) {
	    ($index, $shellName) = &saveShellData($shellData, $minDist, $frameNum, \%{ $DATA->{$i} });
	    $DATA->{$i}{SHELLS}[$index]{STOP} = $frameNum;
	    next if (! $dumpBGFs);
	    ($shellName, $index) = &getCorrectShell($shellData, $minDist);
	    $chain = chr(66 + $index);
	    for $k (keys %{ $MOL }) {
		$SHELLTRJDATA->{ATOMS}{$k} = $BGF->{$k};
		$SHELLTRJDATA->{BONDS}{$k} = $BONDS->{$k};
		$SHELLTRJDATA->{ATOMS}{$k}{CHAIN} = $chain;
	    }
	} elsif ($calcType == 3 and defined($SOLUTE)) {
	    ($shellName, $index) = &getCorrectShell($shellData, $minDist);
	    $SHELLCOUNT{$shellName}{COUNT}++;
	    $SHELLCOUNT{$shellName}{BIN}{$bin}++;
	    for $k (keys %{ $MOL }) {
		$SHELLTRJDATA->{$shellName}{$i}{$k} = $BGF->{$k};
		$SHELLCOUNT{$shellName}{MASS} += $BGF->{$k}{MASS};
	    }
	}		
	$count++;
    }		 

    if ($calcType == 1) {
	($iftDATA->{CELLVOL}, $iftDATA->{SLABVOL}) = ($cell_vol, $slab_vol);
	if ($LAMMPSOPTS->{ISIFT}) {
	    delete $DATA->{DENSITY};
	    for $bin (keys %IFT) {
		$iftDATA->{IFT}{$bin} = calcIFT($BGF, $SOLVENT, \%{ $IFT{$bin} }, $BOX, $axis);
	    }
	    &writeIFTFrameData(\%{ $DATA->{IFT} }, $iftDATA, $frameNum, $iftFile);
	} else {
	    for $bin (keys %IFT) {
		($i, $j) = calcMolOrient($BGF, $SOLVENT, \%{ $IFT{$bin} }, $oPts, $slab_vol);
		$DATA->{ORIENT}{a}{$frameNum}{$bin} = $i;
		$DATA->{ORIENT}{t}{$frameNum}{$bin} = $j;
		&calcChargeDensity(\%{ $DATA->{CHARGE_DENSITY}{$frameNum} }, \%{ $IFT{$bin} }, $bin, $SOLVENT, $BGF, $delR, $dM, $slab_vol);
	    }
	}
    } elsif ($calcType == 2 and $dumpBGFs and defined($SOLUTE)) {
        for $i (keys %{ $SOLUTE })  {
	    $SHELLTRJDATA->{ATOMS}{$i} = $BGF->{$i};
	    $SHELLTRJDATA->{BONDS}{$i} = $BONDS->{$i};
	}
	&saveBGF($SHELLTRJDATA, $BOX, $frameNum);
    } elsif ($calcType == 3 and defined($SOLUTE)) {
        $DATA->{FRAMES} .= "$delT ";
        for $k (keys %SHELLCOUNT) {
	    #first save the shell count info
            $DATA->{$k} .= "$SHELLCOUNT{$k}{COUNT} ";

	    #density
	    ($BULK->{$k}{DENSITY}{tStep}{$delT}{T}, $shell_vol) = calcDensity(\%{ $SHELLCOUNT{$k} }, $k, $BOX, $shellData, $dM);

	    # diffusion constant
	    $BULK->{$k}{DIFFUSION_CONSTANT}{tStep}{$delT} =
	        &CalcDiffusionConstant($ATOMS, $SHELLTRJDATA->{$k}, \%{ $BULK->{$k}{DIFFUSION_CONSTANT}{DAT} }, $delT, 1);

	    # dielectric constant
	    $BULK->{$k}{DIELECTRIC_CONSTANT}{VOLUME}{tStep}{$delT} = $shell_vol;
	    $BULK->{$k}{DIELECTRIC_CONSTANT}{tStep}{$delT} =
		&CalcDielectricConstant($BGF, $SHELLTRJDATA->{$k}, \%{ $BULK->{$k}{DIELECTRIC_CONSTANT}{DAT} }, 
					$BULK->{$k}{DIELECTRIC_CONSTANT}{VOLUME}, $BULK->{TEMPERATURE});
	    $BULK->{$k}{DIELECTRIC_CONSTANT}{tStep}{$delT}{u} =
		&CalcDipole($BGF, $SHELLTRJDATA->{$k}, \%{ $BULK->{$k}{DIELECTRIC_CONSTANT}{DAT} }, 1);

	    # molecular vector correlation
	    next if (! defined($ORIENT));
	    $BULK->{$k}{MOL_CORRELATION}{tStep}{$delT} =
		&CalcMolOrientCorrelation($BGF, $SHELLTRJDATA->{$k}, \%{ $BULK->{$k}{MOL_CORRELATION}{DAT} }, $ORIENT);
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
        $ift->{density}{$bin} += $invol*$atomic2bulkVol*$data->{IFT}{$bin}{mass}/$nBins;
        $ift->{tot}{density} += $invol*$atomic2bulkVol*$data->{IFT}{$bin}{mass}/$nBins;
        $ift->{Press}{$bin} += $data->{IFT}{$bin}{P}/(-3*$cellVol);
        $ift->{tot}{Press} += $data->{IFT}{$bin}{P}/(-3*$cellVol);
    }
    $ift->{tot}{sigma} = ($ift->{tot}{Pn}-$ift->{tot}{Pt})*$delR*$atm2dynepscm*$A2cm/2;
    printf $iftFHandle "%-12d %12.7f %12.5f %12.5f %12.5f %12.3f\n",
    $frame,($ift->{tot}{sigma}/$ift->{count}),($ift->{tot}{Pn}/$ift->{count}),
    ($ift->{tot}{Pt}/$ift->{count}),($ift->{tot}{density}/$ift->{count}),
    ($ift->{tot}{Press}/$ift->{count});
}

sub centerCell {
    my ($atoms, $solute, $solvent, $box, $frame) = @_;
    my ($MOL, $CENTER, @dims, $i, $j, $k, $offset, $solvID);
    
    $CENTER = $offset = {"XCOORD" => 0, "YCOORD" => 0, "ZCOORD" => 0};
    if (defined($solute)) {
	$MOL = GetAtmData($atoms, $solute);
	$CENTER = CoM($MOL);
	#if (! defined($COMstart)) {
	    #$COMstart = $CENTER;
	#} else {
	    #for $j (keys %{ $CENTER }) {
		#$CENTER->{$j} -= $COMstart->{$j};
	    #}
	#}
    }
    @dims = ("XCOORD", "YCOORD", "ZCOORD");
    for $i (@dims) {
        $box->{$i}{CENTER} = $box->{$i}{len}/2;
        $box->{$i}{hi} = $box->{$i}{len};
        $box->{$i}{lo} = 0;
    }

    for $i (keys %{ $atoms }) {
	for $j (keys %{ $atoms->{$i} }) {
	    $BGF->{$i}{$j} = $atoms->{$i}{$j};
	}
        for $j (@dims) {
            $BGF->{$i}{$j} += ($box->{$j}{CENTER} - $CENTER->{$j});
        }
    }
}

sub findSurface {
    my ($atoms, $solute, $box, $store, $dM, $lammpsopts) = @_;
    my ($i, $j, $k, $i1, $i2, $dist, $count, $imageFlag, $index, $mass);
    my ($atomC, $surface, $coordGrid, $soluCOORDS, $soluArray, $tol);

    $tol = 0;
    
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

sub calcChargeDensity {
    my ($data, $iftBinData, $bin, $solvInfo, $atoms, $binSize, $coord, $vol) = @_;
    my ($i, $j, $currBin, $solvID, $atomID, $coordDist, $factor);

    $factor = 10 * q2C * 10E-9/($vol * a2m**3); # charge density in C/m^3 x 10^9
    for $solvID (keys %{ $iftBinData }) {
	for $atomID (keys %{ $solvInfo->{$solvID} }) {
	    $coordDist = $iftBinData->{$solvID}{$coord} - $atoms->{$atomID}{$coord};
	    $currBin = $bin + int($coordDist/$binSize);
	    $data->{$currBin} += $factor * $atoms->{$atomID}{CHARGE}; # charge density in C/A^3 x 10^9
	}
    }
}

sub calcIFT {
    my ($atoms, $mols, $iftData, $BOX, $axis) = @_;
    my ($i, @atmList, $stress, %IFT, $j, $atomC);

    for $i (keys %{ $iftData }) {
        @atmList = grep { /\d+/ } keys %{ $mols->{$i} };
        for $atomC (@atmList) {
	    for $j ("X", "Y", "Z") {
		$stress = $atoms->{$atomC}{"S" . $j . $j};
		$IFT{P} += $stress; # slab scalar total pressure*vol (sum of ensembles)
		if ($j =~  /$axis/i) {
		    $IFT{Pn} += $stress; # slab normal pressure*vol
		} else {
		    $IFT{Pt} += $stress; # slab tangential pressure*vol
		}
	    }
	    $IFT{mass} += $atoms->{$atomC}{MASS}; # slab mass
	}
    }

    return \%IFT;
}

sub calcMolOrient {
    my ($atoms, $mols, $iftData, $orientPts, $vol) = @_;
    my ($i, @atmList, $j, $dM, $angle, $alpha, $k, $tot, $totAngle, $rec, $density);

    $alpha = $tot = 0;
    for $i (keys %{ $iftData }) {
	$tot++;
        @atmList = grep { /\d+/ } keys %{ $mols->{$i} };
	$dM = (); #molecular dipole moment
	for $j (@atmList) {
	    for $k ("X", "Y", "Z") {
		$dM->{$k} += $atoms->{$j}{"${k}COORD"}*$atoms->{$j}{CHARGE};
	    }
	}
	$angle = DotProduct($dM, $orientPts); # cosine of angle between vectors
	$alpha += 2.5*(3*$angle**2 - 1); #5/2*(3 (cos angle)^2 -1)
	$totAngle += $angle;
    }
    return ($alpha/$tot, $totAngle/$tot);
}

sub calcElecPot {
    my ($data) = $_[0];
    my ($frame, $chargeDensity, %binData, $bin, $i, $j, $e0, $A, $LR, $f); 
    my ($zeroPos, $elecPot, @tmp, $maxBin, $minBin, $avg, $del, $elecPotArray);
    my ($lower, $diag, $upper, $extrema);
    $e0 = 8.854187817; #e0 = 1/(u0c^2) = 8.854 x 10^-12 A^2 s^4/(kg m^3)
    for $frame (keys %{ $data }) {
        for $bin (keys %{ $data->{$frame} }) {
            $binData{$bin} .= "$data->{$frame}{$bin} ";
        }
    }

    @tmp = sort numerically keys %binData;
    ($maxBin, $minBin) = ($tmp[$#tmp], $tmp[0]);

    $zeroPos = 0;
    for $i (1 .. $#tmp) {
	if ($tmp[$i] == 0 or ($tmp[$i - 1] < 0 and $tmp[$i] > 0)) {
	    $zeroPos = $i;
	    last;
	} 
    }

    $extrema->{hi} = int(($#tmp - $zeroPos) * 0.1); #remove 10% from end of array to remove end effects
    $extrema->{lo} = int($zeroPos * 0.1); # remove 10% from start of array to remove end effects
    for $i (1 .. $extrema->{hi}) {
	pop @tmp;
    }
    for $i (1 .. $extrema->{lo}) {
	shift @tmp;
    }
    for $bin (0 .. $#tmp) {
        chop $binData{$tmp[$bin]};
        ($avg, undef, undef) = STDev($binData{$tmp[$bin]});
	$chargeDensity->{$tmp[$bin]} = $avg;
    }
    $zeroPos -= $extrema->{lo};

    $del = 1/(($#tmp - $zeroPos) * $delR); # grid spacing on the [0,1] interval in the + dim
    $del *= $del; # grid spacing squared
    $del = -1 * $del/$e0; # multiplication factor for charge density in equation Au = f 
			  #f = -h^2*charge_density(x)/e0, A is tri diagonal matrix for finite difference: 1 2 1
			  #u is electrostatic potential, so u = A^-1 x f
    $f = new Math::MatrixReal(($#tmp - $zeroPos + 1), 1);
    for $i ($zeroPos .. $#tmp) {
	$f->assign(($i - $zeroPos +1),1,($chargeDensity->{ $tmp[$i] } * $del));
	push @{ $lower }, 1;
	push @{ $upper }, 1;
	push @{ $diag }, -2;
    }

    shift @{ $upper };
    shift @{ $lower };
    $A = Math::MatrixReal->new_tridiag($lower, $diag, $upper);
    $LR = $A->decompose_LR();
    ($j, $elecPot, $i) = $LR->solve_LR($f);
    for $i (1 .. $elecPot->[1]) {
	$elecPotArray->{1}{ $tmp[$i + $zeroPos - 1] } = $elecPot->element($i,1) * (($#tmp - $zeroPos) * $delR)**2;
    }

    $zeroPos--;
    $extrema = $#tmp - $extrema;
    $upper = $lower = $diag = ();
    undef($f);
    $f = new Math::MatrixReal(($zeroPos +1), 1);
    for $i (reverse 0 .. $zeroPos) {
	$f->assign(($i+1),1,($chargeDensity->{ $tmp[$i] } * $del));
	push @{ $lower }, 1;
	push @{ $upper }, 1;
	push @{ $diag }, -2;
    }
 
    shift @{ $upper };
    shift @{ $lower };
    $A = Math::MatrixReal->new_tridiag($lower, $diag, $upper);
    $LR = $A->decompose_LR();
    ($j, $elecPot, $i) = $LR->solve_LR($f);
    for $i (1 .. $elecPot->[1]) {
        $elecPotArray->{2}{ $tmp[$i - 1] } = $elecPot->element($i,1) * (($#tmp - $zeroPos) * $delR)**2;
    }

    return $elecPotArray;
}

sub calcDensity {
    my ($shellData, $shellName, $box, $shellInfo, $axis) = @_;
    my ($vol, $i, $mass, $density, $start, $stop, $dist, $dim);

    $vol = 1;
    $mass = $shellData->{MASS};
    for $i (grep {!/$axis/i} keys %{ $box }) {
	$vol *= $box->{$i}{len};
    }

    for $i (keys %{ $shellData->{BIN} }) {
	$dist = $i*$delR;
	$dim = "pos";
	$dim = "neg" if ($dist < 0);
	$start->{$dim} = $dist if (! exists($start->{$dim}) or $start->{$dim} > $dist);
	$stop->{$dim} = $dist if (! exists($stop->{$dim}) or $stop->{$dim} < $dist);
    }
    $dim = 0;
    for $i ("pos", "neg") {
	next if (! exists($start->{$i}));
	$dim += sqrt(($start->{$i} - $stop->{$i})**2);
    }
    $vol *= $dim;
    $density = $mass/($vol * 0.6023);
    return ($density, $vol);

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

    $saveName = "${savePrefix}_ts${frame}.bgf";
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

    $prefix = $savePrefix;
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
    my ($fileName, $orientOpt);

    getopt('mvbtdrlswafcogpik',\%OPTS);

    $usage = &showUsage;

    for ($OPTS{a},$OPTS{b}) {
        die "$usage" if (! defined($_));
    }

    print "Initializing...";
    ($solSel,$solvSel,$trjFile,$bgfFile,$rMax,$delR,$trjType,$tSel,$savePrefix,$axis,$FF, $solvShells, $offset, $dumpBGFs,$orientOpt,$tStep,$sysTemp) =
    ($OPTS{m},$OPTS{v},$OPTS{t},$OPTS{b},$OPTS{r},$OPTS{d},$OPTS{l},$OPTS{s},$OPTS{w},$OPTS{a},$OPTS{f},$OPTS{c},$OPTS{o},$OPTS{g},$OPTS{p},$OPTS{i}, $OPTS{k});

    #$numProcs = 1 if (! defined($numProcs) or $numProcs !~ /^\d+$/);
    $tStep = 0.001 if (! defined($tStep) or ($tStep !~ /^(?:\d+|\d*\.\d+)$/));
    $offset = 0 if (! defined($offset) or $offset !~ /^\d+\.?\d*$/);
    $sysTemp = 300 if (! defined($sysTemp) or $sysTemp !~ /^\d+\.?\d*$/);
    if ($axis =~ /(X|Y|Z)/i) {
	$axis = lc($1);
    } else {
	die "ERROR: Expected X|Y|Z for axis. Got \"$axis\"!\n";
    }

    if ($axis eq "x") {
	$oPts = {"XCOORD" => 1, "YCOORD" => 0, "ZCOORD" => 0};
    } elsif ($axis eq "y") {
	$oPts = {"XCOORD" => 0, "YCOORD" => 1, "ZCOORD" => 0};
    } else {
	$oPts = {"XCOORD" => 0, "YCOORD" => 0, "ZCOORD" => 1};
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

    $savePrefix = basename($bgfFile) if (! defined($savePrefix));
    $savePrefix =~ s/\.\w+$//;

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
	    $BGF->{$j}{MOLECULEID} = $i;
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
    } else {
	for $i (keys %{ $BGF }) {
	    $BGF->{$i}{MASS} = 1;
	}
    }

    $ORIENT = GetOrientOpt($BGF, $SOLVENT, $orientOpt) if (defined($orientOpt));
}

sub openIFTFile {
   my ($prefix) = $_[0];
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
    my ($resData, $solvInfo) = @_;
    my ($i, $j, @atomIDs, $start, $prev, $count, $data, @range, $tmp);

    for $i (keys %{ $resData }) {
	for $j (keys %{ $solvInfo->{$i} }) {
	    $tmp->{$j} = $j;
	}
    }

    @atomIDs = sort numerically keys %{ $tmp };
    $start = $prev = shift @atomIDs;
    $i = $count = 0;
    while ($i < $#atomIDs) {
	if (($atomIDs[$i] - $prev) > 1) {
	    if (($prev - $start) > 1) {
		$data = "${start}-${prev}";
	    } else {
		$data = "${start}";
	    }
	    $start = $atomIDs[$i];
	    $count += length($data);
	    if ($count > 100) {
		#$data .= "\n";
		$count = 0;
	    }
	    push @range, $data;
	}
	$prev = $atomIDs[$i];
	$i++;
    }
    #add the last entry info
    if (($atomIDs[$#atomIDs] - $prev) > 1) {
	$data = "${start}-" . $atomIDs[$#atomIDs];
	push @range, $data;
    } else {
	$data = "${start}-${prev}";
	push @range, $data;
	push @range, $atomIDs[$#atomIDs];
    }
    return ((scalar(@atomIDs) + 1), \@range);
}

sub showUsage {
    my ($usage) = "usage: $0 -b bgf file -a axis (non periodic surface dimension)" .
        "\nOptional parameters:" .
	"\n\t-m solute atoms (see selection info) -v solvent atoms (see selection info) -w save name " .
        "\n\t-r max distance from surface -d del r for binning -f cerius2 forcefield -t trajectory file " .
        "\n\t-s trajectory selection (see selection info) -l traj type (amber(default)|lammps) -c \"solv shell peaks\"" . 
	"\n\t-o offset: integer. The distance from the solute surface-solvent distance" .
        "\n\t-p mol orientation options (none default): For the molecule, expected 1_2:uM for correlation of 1->2 vector" .
        "\n\t\t labeled uM or {name1}_x_uM:p for the vector cross product of predefined name1 (e.g. 1->2) and uM called P" .
        "\n\t\t Use dm for dipole moment. Enclose multiple in quotes" .
	"\n\t-i simulation timestep (in ps). Default is 0.001 (1 fs) -k system temperature (default = 300K)" .
	"\n\t-g dump bgf: Flag to dump bgf files every timestep based on shell. Default 0\n" .
	&ShowSelectionInfo;
    return $usage;
}

sub acos {
    my($x) = $_[0];
    if (abs($x) > 1.0) {
        return 0;
    } else {
        return atan2(sqrt(1 - $x * $x), $x);
    }
}


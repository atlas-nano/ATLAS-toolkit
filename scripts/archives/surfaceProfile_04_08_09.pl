#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use warnings;
no warnings "recursion";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetBGFFileInfo AddMass createHeaders addHeader createBGF);
use Packages::General qw(STDev FileTester CoM GetSelections TrjSelections ShowSelectionInfo LoadFFs);
use Packages::AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox CreateLAMMPSTrj);
use Packages::BOX qw(GetBox CreateGrid);
use Packages::ManipAtoms qw(UnwrapAtoms ImageAtoms GetAtmList GetSolvent SplitAtomsByMol GetAtmData);

sub init;
sub showUsage;
sub boxConvert;
sub calcProfile;
sub getGridAtms;
sub numerically { ($a<=>$b); }
sub getResRange;
sub writeDensity;
sub writeShellData;
sub saveLammpsShellTraj;
sub getSoluNeighCells;
sub saveShellData;
sub writeShellCount;
sub calcIFT;
sub writeIFT;
sub getCorrectShell;
sub saveBGF;
sub splitTrjSelect;

my ($SOLUTE, $SOLVENT, $trjFile, $BGF, $BONDS, $rMax, $printStr, $axis);
my ($getSnapshot, $getByteOffset, $trjType, $LAMMPSOPTS, $delR, $offset);
my ($SELECT, $saveFile, $field, $DATA, $PARMS, $prec, $rounding, $i);
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
        &GetLammpsTrjType($SELECT, $trjFile, "atom", \%{ $LAMMPSOPTS });
        $field = "ift";
    }
    $printStr = "Calculating solvent <-> surface profile from $trjFile...";
    $getSnapshot->($BGF, $trjFile, $SELECT, $field, \&calcProfile, $printStr, undef);
}

if ($calcType == 1) {
    &writeDensity($DATA, $saveFile);
    &writeIFT($DATA, $saveFile);
} elsif ($calcType == 2) {
    &writeShellData($DATA, $saveFile);
} else {
    &writeShellCount($DATA, $saveFile);
}

die "\n";
sub splitTrjSelect {
    my ($selection, $procs) = @_;
    my ($i, $SELECT, @tmp, $count, $j, $key, $k);

    @tmp = sort numerically keys %{ $selection };
    $count = scalar @tmp;
    $j = int($count/$procs);
    for $i (1 .. $procs) {
	for $k (1 .. $j) {
	    last if (! @tmp);
	    $key = shift @tmp;
	    $SELECT->{$i}{$key} = $selection->{$key};
	}
    }

    if (@tmp) {
	for $key (@tmp) {
	    $SELECT->{1}{$key} = $selection->{$key};
	}
    }
    return $SELECT;
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
    my ($saveName, %Pn, %Pt, %density, %sigma, %tot, %Press);
    my ($bin, $count, $invol, $frame, $cellVol, @frames, @tmp1);
    my $atomic2gpcc = 1.66030;
    my $atm2dynepscm=1013250;
    my $A2cm=1e-8;
    my $nBins = $rMax/$delR;

    @frames = sort numerically keys %{ $data };
    @tmp1 = keys %{ $data->{$frames[0]}{IFT} };
    return if (! exists($data->{$frames[0]}{IFT}{$tmp1[0]}));

    $prefix =~ s/\.\w+$//;
    $saveName = "${prefix}_ift.dat";
    print "Writing IFT data to $saveName...";
    $count = 0;

    open IFTDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    printf IFTDATA "#%-11s %12s %12s %12s %12s %12s\n","TIMESTEP", "IFT", "Pn", "Pt", "DENSITY", "PRESSURE";
    for $frame (@frames) {
	$count++;
	$invol = 1/$data->{$frame}{SLABVOL};
	$cellVol = $data->{$frame}{CELLVOL};
	for $bin (keys %{ $data->{$frame}{IFT} }) {
	    $Pn{$bin} += -$invol*$data->{$frame}{IFT}{$bin}{Pn};
	    $tot{Pn} += -$invol*$data->{$frame}{IFT}{$bin}{Pn};
	    $Pt{$bin} += -$invol*$data->{$frame}{IFT}{$bin}{Pt};
	    $tot{Pt} += -$invol*$data->{$frame}{IFT}{$bin}{Pt};
	    $density{$bin} += $nBins*$atomic2gpcc*$data->{$frame}{IFT}{$bin}{mass}/$cellVol;
	    $tot{density} += $nBins*$atomic2gpcc*$data->{$frame}{IFT}{$bin}{mass}/$cellVol;
	    $Press{$bin} += $data->{$frame}{IFT}{$bin}{P}/(-3*$cellVol);
	    $tot{Press} += $data->{$frame}{IFT}{$bin}{P}/(-3*$cellVol);
	}
        $tot{sigma} = ($tot{Pn}-$tot{Pt})*$delR*$atm2dynepscm*$A2cm/2;
	printf IFTDATA "%-12d %12.7f %12.5f %12.5f %12.5f %12.3f\n",$frame,($tot{sigma}/$count),
	($tot{Pn}/$count),($tot{Pt}/$count),($tot{density}/$count),($tot{Press}/$count);
    }
    close IFTDATA;

    $saveName = "${prefix}_ift_profile.dat";
    print "Done\nWriting IFT Profile data to $saveName...";
    open IFTPROFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    $prec += 8;
    printf IFTPROFILE "%-${prec}s %12s %12s %12s %12s %12s\n","#BIN", "IFT", "Pn", "Pt", "DENSITY", "PRESSURE";
    for $bin (sort numerically keys %Pt) {
	$Pn{$bin} /= $count;
	$Pt{$bin} /= $count;
	$Press{$bin} /= $count;
	$density{$bin} /= $count;
	# now compute ift
	$sigma{$bin} = ($Pn{$bin}-$Pt{$bin})*$delR*$atm2dynepscm*$A2cm/2;
        printf IFTPROFILE "%-${prec}f %12.7f %12.5f %12.5f %12.5f %12.3f\n",($bin*$delR),
			  $sigma{$bin},$Pn{$bin},$Pt{$bin},$density{$bin},$Press{$bin};
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
	for $bin (keys %{ $data->{$frame}{DENSITY} }) {
	    $binData{$bin} .= "$data->{$frame}{DENSITY}{$bin} ";
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
    my ($MOL, $CENTER, @tmp, $i, $j, $SOLCENTER, $SOLATMS, $SOLVATMS, $index, $soluCells);
    my ($bin, $dist, $count, $GRID, $shellName, $shellID, $solvGRIDLIST, $tot, $gridAtms);
    my ($x, $y, $z, $k, $solGRIDLIST, $SOLVCENTER, $minDist, $i1, $i2, $SHELLTRJDATA);
    my ($soluGridAtms, $solvGridAtms, $atmOffset, $totMass, %SHELLCOUNT, %IFT, $isIFT);
    my ($chain);

    $isIFT = 0;
    if ($trjType == 2) { #LAMMPS
        $frameNum = $ATOMS->{TIMESTEP}[0];
	$tot = $ATOMS->{"NUMBER OF ATOMS"}[0];
        $BOX = ConvertLammpsBox(\%{ $ATOMS->{"BOX BOUNDS"} });
        $ATOMS = $ATOMS->{ATOMS};
	if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
	    UnwrapAtoms(\%{ $ATOMS }, \%{ $BBOX }, $LAMMPSOPTS->{scaled});
	}
	@tmp = grep {!/COORD|VISITED/i} keys %{ $BGF->{1} };
	for $i (keys %{ $BGF }) {
	    for $j (@tmp) {
		$ATOMS->{$i}{$j} = $BGF->{$i}{$j};
	    }
	    $ATOMS->{$i}{IS_SOLVENT} = 1 if (exists($BGF->{$i}{IS_SOLVENT}));
	    $isIFT = 1 if (exists($ATOMS->{$i}{SXX}));
	}
    } elsif ($trjType == 1) {
        $BOX = ConvertAmberBox(\%{ $BOX });
	$tot = scalar(keys %{ $ATOMS });
    } 

    $atmOffset = -1;
    $totMass = 0;
    for $i (sort numerically keys %{ $SOLUTE }) {
        $gridAtms->{$i} = \%{ $ATOMS->{$i} };
        for $j ("XCOORD", "YCOORD", "ZCOORD") {
	    $CENTER->{$j} += $ATOMS->{$i}{$j} * $ATOMS->{$i}{MASS};
        }
        $soluGridAtms->{$i} = \%{ $gridAtms->{$i} };
	$gridAtms->{$i}{NUMBONDS} = 3;
        $atmOffset = $i;
	$totMass += $ATOMS->{$i}{MASS};
    }
    $atmOffset++;

    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    for $j (@tmp) {
	$BOX->{$j}{CENTER} = $BOX->{$j}{len}/2;
	$BOX->{$j}{hi} = $BOX->{$j}{len};
	$BOX->{$j}{lo} = 0;
	for $i (keys %{ $ATOMS }) {
	    $ATOMS->{$i}{$j} += ($BOX->{$j}{CENTER} - $CENTER->{$j}/$totMass);
	}
    }
     
    for $i (keys %{ $SOLVENT }) {
	$MOL = GetAtmData(\%{ $ATOMS }, \%{ $SOLVENT->{$i} });
	$CENTER = CoM(\%{ $MOL });
	$k = $i + $atmOffset;
	$gridAtms->{$k} = ImageAtoms(\%{ $MOL }, \%{ $CENTER }, \%{ $BOX });
        $gridAtms->{$k}{IS_SOLVENT} = 1;
        $gridAtms->{$k}{RESNAME} = "WAT";
        $gridAtms->{$k}{RESNUM} = $i;
        $gridAtms->{$k}{CHARGE} = 0;
        $gridAtms->{$k}{RADII} = 1.8;
        $gridAtms->{$k}{NUMBONDS} = 3;
        $gridAtms->{$k}{INDEX} = $k;
	$gridAtms->{$k}{MOLECULE} = \%{ $SOLVENT->{$i} };
        $solvGridAtms->{$k} = \%{ $gridAtms->{$k} };
    }

   for $i ("X", "Y", "Z") {
	$BOX->{$i} = $BOX->{"${i}COORD"};
   }
   ($GRID, $i, $j) = CreateGrid(\%{ $gridAtms }, 5, \%{ $BOX }, 2, 0);
   $soluCells = ();

   for $i (keys %{ $soluGridAtms }) { 
	next if (! defined($soluGridAtms->{$i}{CELL}) or ! keys %{ $soluGridAtms->{$i}{CELL} });
	($x, $y, $z) = ($soluGridAtms->{$i}{CELL}{XINDEX}, $soluGridAtms->{$i}{CELL}{YINDEX}, $soluGridAtms->{$i}{CELL}{ZINDEX});
	$soluCells->{$x}{$y}{$z} = \%{ $GRID->{$x}{$y}{$z} };
   }
   &getSoluNeighCells($GRID, $soluCells, $axis);

   $count = 0;
   for $i (keys %{ $solvGridAtms }) {
	next if (exists($gridAtms->{$i}{PLACED}));
	($x, $y, $z) = ($gridAtms->{$i}{CELL}{XINDEX}, $gridAtms->{$i}{CELL}{YINDEX}, $gridAtms->{$i}{CELL}{ZINDEX});
	($SOLATMS, $SOLVATMS) = getGridAtms(\%{ $GRID }, $axis, $x, $y, $z, \%{ $soluCells }, $solvGridAtms);
	next if (! %{ $SOLVATMS } or ! %{ $SOLATMS });
	for $j (keys %{ $SOLVATMS }) {
	    $dist = ();
	    for $k (keys %{ $SOLATMS }) {
		$dist->{$k} = ($SOLATMS->{$k}{XCOORD} - $SOLVATMS->{$j}{XCOORD})**2 + 
			      ($SOLATMS->{$k}{YCOORD} - $SOLVATMS->{$j}{YCOORD})**2 +
			      ($SOLATMS->{$k}{ZCOORD} - $SOLVATMS->{$j}{ZCOORD})**2;
	    }
	    @tmp = sort numerically values %{ $dist };
	    $minDist = sqrt(shift @tmp);
	    $bin = sprintf("%${rounding}f",($minDist/$delR) + 1);
	    if ($calcType == 1)  {
		$DATA->{$frameNum}{DENSITY}{$bin}++;
		$IFT{$bin}{$j} = 1;
	    } elsif ($calcType == 2) {
		($index, $shellName) = &saveShellData($shellData, $minDist, $frameNum, \%{ $DATA->{$j} });
		$DATA->{$j}{SHELLS}[$index]{STOP} = $frameNum;
                next if (! $dumpBGFs);
		($shellName, $index) = &getCorrectShell($shellData, $minDist);
                $chain = chr(66 + $index);
                for $k (keys %{ $SOLVATMS->{$j}{MOLECULE} }) {
                    $SHELLTRJDATA->{ATOMS}{$k} = $ATOMS->{$k};
                    $SHELLTRJDATA->{BONDS}{$k} = $BONDS->{$k};
                    $SHELLTRJDATA->{ATOMS}{$k}{CHAIN} = $chain;
                }
	    } elsif ($calcType == 3) {
		($shellName, $index) = &getCorrectShell($shellData, $minDist);
		$SHELLCOUNT{$shellName}++;
                for $k (keys %{ $SOLVATMS->{$j}{MOLECULE} }) {
                    $SHELLTRJDATA->{$shellName}{$k} = $ATOMS->{$k};
                }
	    }		
	    $count++;
	}		 
    }
    if ($calcType == 1) {
	if ($isIFT) {
	    for $bin (keys %IFT) {
		($DATA->{$frameNum}{IFT}{$bin}, $DATA->{$frameNum}{CELLVOL}, $DATA->{$frameNum}{SLABVOL}) = 
		&calcIFT($ATOMS, $SOLVENT, \%{ $IFT{$bin} }, $BOX, $axis);
	    }
	}
    } elsif ($calcType == 2 and $dumpBGFs) {
        for $i (keys %{ $SOLUTE })  {
	    $SHELLTRJDATA->{ATOMS}{$i} = $ATOMS->{$i};
	    $SHELLTRJDATA->{BONDS}{$i} = $BONDS->{$i};
	}
	&saveBGF($SHELLTRJDATA, $BOX, $frameNum);
    } elsif ($calcType == 3) {
	&saveLammpsShellTraj(\%{ $SHELLTRJDATA }, \%{ $shellData }, $frameNum, $tot, $BOX, $axis);
	$DATA->{FRAMES} .= "$frameNum ";
	for $i (keys %SHELLCOUNT) {
	    $DATA->{$i} .= "$SHELLCOUNT{$i} ";
	}
    }
    
    %SHELLCOUNT = ();
    $SHELLTRJDATA = ();
    $ATOMS = ();
    $gridAtms = ();
    $soluGridAtms = ();
    $solvGridAtms = ();
    $GRID = ();
    $soluCells = ();
    $SOLATMS = ();
    $SOLVATMS = ();
}

sub calcIFT {
    my ($atoms, $mols, $iftData, $BOX, $axis) = @_;
    my ($i, @atmList, $atomC, $vol, @tmp, %map);
    my (%IFT, $slab_vol, $tot_vol, $stress, $j);

    $axis = uc($axis);
    @tmp = ("X", "Y", "Z");
    %map = ("X" => 0, "Y" => 1, "Z" => 2);

    for $i (@tmp) {
        $vol->{$i} = $BOX->{$i}{hi} - $BOX->{$i}{lo};
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

sub getSoluNeighCells {
    my ($grid, $soluCells, $caxis) = @_;
    my ($x, $y, $z, $cells, $k, $l, $m, $map, $currCell, $rec);

    $caxis = uc($caxis);

    @{ $cells->{X} } = sort numerically keys %{ $soluCells };
    for $x (0 .. $#{ $cells->{X} }) {
	$map->{X} = ();
	if ($caxis eq "X") {
	    $map->{X}[0] = $cells->{X}[$x];
	} else {
            for $k (-1, 0, 1) {
	        $l = $x+$k;
	        $l = $#{ $cells->{X} } if ($l < 0);
	        $l = 0 if ($l > $#{ $cells->{X} });
	        push @{ $map->{X} }, $cells->{X}[$l];
	    }
	}
	@{ $cells->{Y} } = sort numerically keys %{ $soluCells->{$cells->{X}[$x]} };
	for $y (0 .. $#{ $cells->{Y} }) {
	    $map->{Y} = ();
	    if ($caxis eq "Y") {
		$map->{Y}[0] = $cells->{Y}[$y];
	    } else {
                for $k (-1, 0, 1) {
                    $l = $y+$k;
                    $l = $#{ $cells->{Y} } if ($l < 0);
                    $l = 0 if ($l > $#{ $cells->{Y} });
                    push @{ $map->{Y} }, $cells->{Y}[$l];
		}
            }
            @{ $cells->{Z} } = sort numerically keys %{ $soluCells->{$cells->{X}[$x]}{$cells->{Y}[$y]} };
            for $z (0 .. $#{ $cells->{Z} }) {
	        $map->{Z} = ();
		if ($caxis eq "Z") {
		    $map->{Z}[0] = $cells->{Z}[$z];
		} else {
	            for $k (-1, 0, 1) {
        	        $l = $z+$k;
	                $l = $#{ $cells->{Z} } if ($l < 0);
        	        $l = 0 if ($l > $#{ $cells->{Z} });
	                push @{ $map->{Z} }, $cells->{Z}[$l];
	            }
		}
		$soluCells->{$cells->{X}[$x]}{$cells->{Y}[$y]}{$cells->{Z}[$z]}{NEIGHS} = ();
		$currCell = \%{ $soluCells->{$cells->{X}[$x]}{$cells->{Y}[$y]}{$cells->{Z}[$z]} };
		for $k (@{ $map->{X} }) {
		    for $l (@{ $map->{Y} }) {
			next if (! $soluCells->{$k}{$l});
			for $m (@{ $map->{Z} }) {
                            next if (! $soluCells->{$k}{$l}{$m} or ! $soluCells->{$k}{$l}{$m}{ATOMS});
			    $rec->{ATOMS} = ();
                            @{ $rec->{ATOMS} } = @{ $soluCells->{$k}{$l}{$m}{ATOMS} };
                            push @{ $currCell->{NEIGHS} }, $rec;
			}
		    }
		}
	    }
	}
    }
    for $k ("X", "Y", "Z") {
        delete $map->{$k};
        delete $cells->{$k};
    }
    $map = ();
    $cells = ();
    $currCell = ();

}

sub saveBGF {
    my ($shellData, $BOX, $frame) = @_;
    my ($saveName, $HEADERS, $ATOMS, $BONDS);

    $saveName = $saveFile;
    $saveName =~ s/\.\w+$//;
    $saveName .= "_ts${frame}.bgf";
    $ATOMS = $shellData->{ATOMS};
    $BONDS = $shellData->{BONDS};
    $HEADERS = createHeaders($BOX, $saveName);
    &addHeader($ATOMS, $HEADERS);
    &createBGF($ATOMS, $BONDS, $saveName);

}

sub saveLammpsShellTraj {
    my ($trjData, $shellData, $timeStep, $solvAtms, $BOX, $axis) = @_;
    my ($i, $trjName, $prefix, $headers, $shellName, $prevShell, @tmp, %map);

    $prefix = $saveFile;
    $prefix =~ s/\.\w+$//;
    $prevShell = $offset;
    $axis = uc($axis);
    @tmp = ("X", "Y", "Z");
    %map = ("X" => 0, "Y" => 1, "Z" => 2);

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
	$headers->{"BOX BOUNDS"}[$map{$axis}]{hi} = 2*$shellData->[$i]{END};
	$trjName = "${prefix}_${shellName}.lammpstrj";
	open LAMMPSTRJ, ">> $trjName" or die "ERROR: Cannot write to $trjName: $!\n";
	&CreateLAMMPSTrj($trjData->{$shellName}, $headers, \*LAMMPSTRJ);
	close LAMMPSTRJ;
    }

}

sub getGridAtms {
    my ($grid, $dim, $a, $b, $c, $soluCells, $solvAtmList) = @_;
    my (%solvATMS, %soluATMS, $i, $j, $k, $currCell, $isFound, $l);
    my ($map);

    $currCell = $grid->{$a}{$b}{$c};
    for $i (@{ $currCell->{WATERS} }) {
	$solvATMS{$i->{RESNUM}} = \%{ $i };
    }

    $isFound = 0;
    $map = ();
    if ($dim eq "x") {
	@{ $map } = ([1 .. $grid->{X}{tot}],[$b, ($b + 1), ($b - 1)],[$c, ($c + 1), ($c - 1)]);
    } elsif ($dim eq "y") {
	@{ $map } = ([$a, ($a + 1), ($a - 1)],[1 .. $grid->{Y}{tot}],[$c, ($c + 1), ($c - 1)]);
    } else {
	@{ $map } = ([$a, ($a + 1), ($a - 1)],[$b, ($b + 1), ($b - 1)],[1 .. $grid->{Z}{tot}]);
    }

    for $i (@{ $map->[0] }) {
	for $j (@{ $map->[1] }) {
	    for $k (@{ $map->[2] }) {
		next if (! exists($soluCells->{$i}) or ! exists($soluCells->{$i}{$j}) or ! exists($soluCells->{$i}{$j}{$k}));
		$currCell = $grid->{$i}{$j}{$k};
		for $l (@{ $currCell->{ATOMS} }) {
		    $soluATMS{ $l->{INDEX} } = \%{ $l };
		}
		#$isFound = 1;
		#last;
	    }
	    last if ($isFound);
	}
	last if ($isFound);
    }

    return (\%soluATMS, \%solvATMS);
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

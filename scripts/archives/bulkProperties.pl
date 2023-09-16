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
use Packages::General qw(FileTester CoM GetSelections TrjSelections ShowSelectionInfo DotProduct STDev CrossProduct LoadFFs);
use Packages::AMBER qw(ParseAmberTrj GetAmberByteOffset);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType);
use Packages::BOX qw(ConvertBox MakeBox);
use Packages::ManipAtoms qw(UnwrapAtoms GetAtmList GetSolvent SplitAtomsByMol);
use constant PI => atan2(1,1) * 4;

sub init;
sub showUsage;
sub calcProperties;
sub calcDipoleMoment;
sub calcDisplacement;
sub calcVolume;
sub getAtoms;
sub getTotMass;
sub saveProperties;
sub determineNumFrames;
sub numerically { ($a<=>$b) }
sub calcDielectricConstant;
sub calcDiffusionConstant;
sub calcMolOrientCorrelation;
sub writeData;
sub getOrientOpt;
sub getVecs;

my ($MOLS, $trjFile, $BGF, $BONDS, $printStr, $PARMS, $list, $i, $prevTime);
my ($getSnapshot, $getByteOffset, $trjType, $LAMMPSOPTS, $tStep, $ORIENT, $totAtms);
my ($SELECT, $saveFile, $field, $DATA, $HEADERS, $totMass, $totFrames, $totMols);

$|++;
&init;
$PARMS = LoadFFs($list);
AddMass($BGF, $PARMS);
$totMass = getTotMass($BGF, $MOLS);
$field = scalar keys %{ $BGF };
$getByteOffset->($SELECT, $trjFile, $field);
$totFrames = &determineNumFrames($SELECT);
if ($trjType == 2) {
    &GetLammpsTrjType($SELECT, $trjFile, "atom", \%{ $LAMMPSOPTS });
    $field = "atom";
}
$printStr = "Calculating bulk properties from $trjFile...";
$prevTime = 0;
$getSnapshot->($BGF, $trjFile, $SELECT, $field, \&saveProperties, $printStr, $MOLS);

delete $DATA->{DIFFUSION_CONSTANT}{START_COM};
delete $DATA->{MOL_CORRELATION}{START};

print "Writing data...";
writeData($DATA, $saveFile);
print "Done\n";

for $i ("DIELECTRIC_CONSTANT", "DIFFUSION_CONSTANT", "DENSITY") {
    printf "%-20s: %8.3g $DATA->{$i}{UNITS} +/- %8.5g\n", $i, $DATA->{$i}{STATS}{T}{AVG}, $DATA->{$i}{STATS}{T}{STDEV};
}

sub writeData {
    my ($data, $savePrefix) = @_;
    my ($i, $avg, $stdev, $tmp, $fileName, $type); 
    my ($avg_str, $stdev_str, @tList, $j, @index, $delT);

    $data->{DENSITY}{UNITS} = "g/cm3";
    $data->{DIELECTRIC_CONSTANT}{UNITS} = "";
    $data->{DIFFUSION_CONSTANT}{UNITS} = "cm2/s";
    @index = sort numerically keys %{ $data->{DIELECTRIC_CONSTANT}{DATA} };

    for $type ("DIFFUSION_CONSTANT", "DENSITY", "DIELECTRIC_CONSTANT", "MOL_CORRELATION") {
	@tList = ("T");
	for $i (keys %{ $data->{$type}{DATA}{$index[0]} }) {
	    next if ($i eq "T");
	    push @tList, $i;
	}
	$fileName = $savePrefix . "_" . lc($type) . ".dat";
	open DATA, "> $fileName" or die "ERROR: Cannot create $fileName: $!\n";
	printf DATA "#%7s ", "TIME";
	for $j (@tList) {
	    printf DATA "%12s ", $j;
	}
	print DATA "\n";
	$tmp = ();
	for $i (@index) {
	    printf DATA "%8.3f ", $i; 
	    for $j (@tList) {
		printf DATA "%12.5f ", $data->{$type}{DATA}{$i}{$j};
		if ($type ne "DIFFUSION_CONSTANT") {
		    $tmp->{$j} .= "$data->{$type}{DATA}{$i}{$j} ";
		} elsif ($i > $index[0]) {
		    $delT = $i - $index[0];
		    $tmp->{$j} .= ($data->{$type}{DATA}{$i}{$j} * 1e-03/(6 * $delT)) . " ";
		}
	    }
	    print DATA "\n";
	}
	
	if ($type ne "MOL_CORRELATION") {
	    $avg_str = sprintf("#%7s ", "AVG");
	    $stdev_str = sprintf("#%7s ", "STDEV");
	    for $j (@tList) {
		next if (! exists($data->{$type}{DATA}{$index[0]}{$j}));
		chop $tmp->{$j};
		($avg, $stdev, undef) = STDev($tmp->{$j});
		$data->{$type}{STATS}{$j}{AVG} = $avg;
		$data->{$type}{STATS}{$j}{STDEV} = $stdev;
		$avg_str .= sprintf("%12.5G ", $avg);
		$stdev_str .= sprintf("%12.8G ", $stdev);
	    }
	    print DATA "${avg_str}\n${stdev_str}\n";
	}
	close DATA;
    }

}
sub saveProperties {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($CENTER, @tmp, $i, $j, $volume, $MOLDATA);
    my ($DIPOLE, $DISPLACEMENT, $index);
    
    if ($trjType == 2) { #LAMMPS
        $frameNum = $ATOMS->{TIMESTEP}[0]/1000; #timestep in ps
        $BOX = $ATOMS->{"BOX BOUNDS"};
	$volume = ($BOX->[0]{hi} - $BOX->[0]{lo}) * ($BOX->[1]{hi} - $BOX->[1]{lo}) *
		  ($BOX->[2]{hi} - $BOX->[2]{lo});
        $ATOMS = $ATOMS->{ATOMS};
	$BOX = MakeBox(ConvertBox($BOX));
        if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
            UnwrapAtoms($ATOMS, $BOX, $LAMMPSOPTS->{scaled});
        }

        for $i (keys %{ $ATOMS }) {
            for $j ("MASS", "CHARGE") {
                $ATOMS->{$i}{$j} = $_[1]->{$i}{$j};
            }
        }
    } else {
        $volume = $BOX->{2}{DATA} * $BOX->{3}{DATA} * $BOX->{4}{DATA};
    }

    @tmp = ("XCOORD", "YCOORD", "ZCOORD");

    for $i (keys %{ $MOLS }) {
	$MOLDATA->{$i} = getAtoms($ATOMS, $MOLS->{$i});
	$CENTER = CoM($MOLDATA->{$i});
	$DIPOLE->{$i} = &calcDipoleMoment($MOLDATA->{$i}, $CENTER);
	if (! $prevTime) {
	    $DATA->{DIFFUSION_CONSTANT}{START_COM}{$i} = $CENTER;
	    %{ $DATA->{MOL_CORRELATION}{START}{$i} } = %{ $DIPOLE->{$i} };
	    %{ $DATA->{MOL_CORRELATION}{START_COORD}{$i} } = %{ $MOLDATA->{$i} };
	}
	$DISPLACEMENT->{$i} = &calcDisplacement($CENTER, $DATA->{DIFFUSION_CONSTANT}{START_COM}{$i});
    }
    $prevTime = $frameNum * $tStep;

    $DATA->{DENSITY}{DATA}{$prevTime}{T} = $totMass/($volume * 0.6023);
    $DATA->{DIFFUSION_CONSTANT}{DATA}{$prevTime} = &calcDiffusionConstant($DISPLACEMENT);
    $DATA->{MOL_CORRELATION}{DATA}{$prevTime} = &calcMolOrientCorrelation($DATA->{MOL_CORRELATION}, $DIPOLE, $ORIENT, $MOLDATA);

    $DIPOLE = ();
    $DIPOLE = &calcDipoleMoment($ATOMS, $CENTER);
    $DATA->{DIELECTRIC_CONSTANT}{DATA}{$prevTime} = &calcDielectricConstant(\%{ $DATA->{DIELECTRIC_CONSTANT} }, $DIPOLE, $volume);
    $DIPOLE = $CENTER = $DISPLACEMENT = ();
}

sub calcMolOrientCorrelation {
    my ($data, $dipole, $orientOpts, $molData) = @_;
    my ($i, $count, $corr, $j, $curr, $start, $type, $a1, $a2, $isSame);
    my (@tmp, $offset, $molMap, %MOLORIENT, $v1, $v2, $startCoord);

    $molMap->{0} = 0;
    $molMap->{"-1"} = -1;
    $startCoord = $data->{START_COORD};

    for $type (keys %{ $orientOpts }) {
	$count = 0;
	if ($type eq "T") {
	    for $i (keys %{ $molData }) {
		for $j ("X", "Y", "Z") {
		    $start->{"${j}COORD"} = $data->{START}{$i}{$j};
		    $curr->{"${j}COORD"} = $dipole->{$i}{$j};
		}
		$corr += DotProduct($curr, $start)/DotProduct($start, $start);
		$count++;
	    }
	} elsif ($orientOpts->{$type} =~ /(\d+)_(\d+)_x_(\d+)_(\d+)/) { #cross product
	    for $i (keys %{ $molData }) {
		@tmp = sort numerically keys %{ $molData->{$i} };
		for $j (0 .. $#tmp) {
		    $offset = $j + 1;
		    $molMap->{$offset} = $tmp[$j];
		}
		for $j ("X", "Y", "Z") {
		    $startCoord->{$i}{0}{"${j}COORD"}  = $data->{START}{$i}{$j};
		    $molData->{$i}{0}{"${j}COORD"}  = $dipole->{$i}{$j};
		    $startCoord->{$i}{-1}{"${j}COORD"} = 0;
		    $molData->{$i}{-1}{"${j}COORD"} = 0;
		}

		($v1, $v2) = getVecs($molMap, $startCoord->{$i}, $molData->{$i}, $1, $2);
		$start = CrossProduct($v1, $v2);

		($v1, $v2) = getVecs($molMap, $startCoord->{$i}, $molData->{$i}, $3, $4);
		$curr = CrossProduct($v1, $v2);

		$isSame = 1;
		for $j (keys %{ $curr }) {
		    if ($start->{$j} != $curr->{$j} or $start->{$j} > 0) {
			$isSame = 0;
			last;
		    }
		}
		if (! $isSame) {
		    $corr += DotProduct($curr, $start)/DotProduct($start, $start);
		} else {
		    $corr++;
		}
		$count++;				
	    }
	} elsif ($orientOpts->{$type} =~ /(\d+)_(\d+)/) { #regular vector
	    for $i (keys %{ $molData }) {
		@tmp = sort numerically keys %{ $molData->{$i} };
		for $j (0 .. $#tmp) {
		    $offset = $j + 1;
		    $molMap->{$offset} = $tmp[$j];
		}
		for $j ("X", "Y", "Z") {
		    $startCoord->{$i}{0}{"${j}COORD"}  = $data->{START}{$i}{$j};
		    $molData->{$i}{0}{"${j}COORD"}  = $dipole->{$i}{$j};
		    $startCoord->{$i}{-1}{"${j}COORD"} = 0;
		    $molData->{$i}{-1}{"${j}COORD"} = 0;
		}
		($start, $curr) = getVecs($molMap, $startCoord->{$i}, $molData->{$i}, $1, $2);
		$corr += DotProduct($curr, $start)/DotProduct($start, $start);
		$count++;
	    }		
	}
	
	$corr /= $count;
	$MOLORIENT{$type} = $corr;
    }

    return \%MOLORIENT;
}

sub calcDiffusionConstant {
    my ($displacement) = @_;
    my ($i, $count, $j, %diffuse);

    $count = 0;
    for $i (keys %{ $displacement }) {
	for $j (keys %{ $displacement->{$i} }) {
	    $diffuse{$j} += $displacement->{$i}{$j};
	}
	$count++;
    }
    
    for $i (keys %diffuse) {
	if ($count > 0) {
	    $diffuse{$i} /= $count;
	} else {
	    $diffuse{$i} = 0;
	}
    }
    return \%diffuse;
}

sub calcDielectricConstant {
    my ($data, $dipole, $V) = @_;
    my ($i, $delM, %dielectric, $kbT, $j, @types, $dat, $factor);

#   dielectric = (4 * Pi * (<M^2> - <avg M>^2))/(3 * e0 * kb * t * V)
#   M = C*m = 1.602176E-19 electron charge * 1E-10 Angstroms
#   e0 = 8.8541878176E-12 C^2 N^-1 m^-2
#   kb = 1.3806504E-23 J K^-1
#   V = m^3 = 1E-30 A^3
#   so have: 
#   4*PI*(1.602176E-29)^2/(3*300*1.306504 * 8.8541878176E-65) 
#   = (32.25747042608465465067 * 10E7)/(10411.22862040110336)
#	= 29319.41799957042610000000
    @types = keys %{ $dipole };

    $factor = 29319.41799957/$V;
    #(4 * PI)/(8.8541878176E-12 * 1.3806504E-23 * 300 * $V * 1E-30);

    if (! exists($data->{FrameCount})) {
	$data->{FrameCount} = 0;
	for $i (@types) {
	    $data->{M_squared}{$i} = 0;
	    $data->{M_bar}{$i} = 0;
	}
    }

    $data->{FrameCount}++;

    for $j (@types) {
	$dat = $dipole->{$j};
	$data->{M_squared}{$j} += $dipole->{$j}**2;
	if ($j ne "R") {
	    $data->{M_bar}{$j} += $dat;
	} else {
	    $data->{M_bar}{R} += $dat/(scalar keys %{ $MOLS });
	}
    }
    
    $data->{M_squared}{T} = ($data->{M_squared}{X} + $data->{M_squared}{Y} + $data->{M_squared}{Z});
    $data->{M_bar}{T} = ($data->{M_bar}{X} + $data->{M_bar}{Y} + $data->{M_bar}{Z});
    $dielectric{M_s} = $data->{M_squared}{T}/$data->{FrameCount};
    $dielectric{M_b} = $data->{M_bar}{T}/$data->{FrameCount};
    $dielectric{G} = $dielectric{M_s}/($totAtms * $dielectric{M_b}**2);

    push @types, ("T");
    for $i (@types) {
	$delM = ($data->{M_squared}{$i}/$data->{FrameCount}) -($data->{M_bar}{$i}/$data->{FrameCount})**2;
	$dielectric{$i} = 1 + $factor * $delM;
	#print "\n$i delM: $delM vol: $V dielectric: $dielectric{$i}\n";
    }

    $dielectric{R} = $data->{M_bar}{R} * 4.80320418270556774711/$data->{FrameCount};
    return \%dielectric;

}

sub calcDipoleMoment {
    my ($atoms, $com) = @_;
    my (@dim, $i, $j, %dipole, $dat);
    
    @dim = ("X", "Y", "Z");
    
    for $j (@dim) {
	$dipole{$j} = 0;
	for $i (keys %{ $atoms }) {
	    #dipole moment in atomic units
	    #$dipole{$j} += $atoms->{$i}{CHARGE} * ($atoms->{$i}{"${j}COORD"} - $com->{"${j}COORD"});
	    $dipole{$j} += $atoms->{$i}{CHARGE} * $atoms->{$i}{"${j}COORD"};#  * 1.602176E-29; #dipole moment in Coulomb-meter
	}
    }

    for $j (@dim) {
	$dipole{R} += $dipole{$j}**2;
    }
    $dipole{R} = sqrt($dipole{R});
    return \%dipole;
}

sub calcDisplacement {
    my ($com, $startCom) = @_;
    my ($i, %displacement, $totDist, $dist);

    $totDist = $dist = 0;
    for $i ("X", "Y", "Z") {
	$dist = ($com->{"${i}COORD"} - $startCom->{"${i}COORD"})**2;
	$totDist += $dist;
	$displacement{$i} = $dist
    }

    $displacement{R} = $totDist;
    $displacement{T} = sqrt($totDist);
    return \%displacement;
}

sub calcVolume {
    my ($box) = $_[0];
    my ($vol);

    $vol = 1;
    for ("XCOORD", "YCOORD", "ZCOORD") {
	$vol *= ($box->{$_}{hi} - $box->{$_}{lo});
    }

    return $vol;
}

sub getAtoms {
    my ($allAtoms, $atomList) = @_;
    my (%ATOMS, $i);

    for $i (keys %{ $atomList }) {
        $ATOMS{$i} = $allAtoms->{$i};
    }

    return \%ATOMS;
}

sub getTotMass {
    my ($atoms, $mols) = @_;
    my ($i, $j, $mass);

    for $i (keys %{ $mols }) {
	for $j (keys %{ $mols->{$i} }) {
	    $mass += $atoms->{$j}{MASS};
	}
    }

    return $mass;
}

sub determineNumFrames {
    my ($trjSelect) = $_[0];
    my ($originalFrames, @tmp, $power, $i, $tot);

    @tmp = sort numerically keys %{ $trjSelect };
    $originalFrames = scalar @tmp;
    $power = $tot = 2;
    die "ERROR: Need at least 2 frames to perform analysis!" if ($originalFrames == 1);
    while ($originalFrames >= ($power * 2)) {
	$power *= 2;
	$tot = $power;
    }

    if ($originalFrames > $tot) {
	print "NOTE: Using $tot frames for FFT: frames $tmp[0] to $tmp[$tot-1]\n";
	for $i ($tot .. $#tmp) {
	    delete $trjSelect->{$tmp[$i]};
	}
    }

    return $tot;
}
    
sub init {
    my (%OPTS, @tmp, $i, $j, $solvTmp, $count, $orientOpt);
    my ($atomSel, $bgfFile, $tSel, $FF, $ATOMSELECT);
    
    getopt('btalswfdo',\%OPTS);
    
    for ("b", "t", "a", "f", "d") {
	die &showUsage . "" if (! defined($OPTS{$_}));
    }
    
    print "Initializing...";
    ($atomSel,$trjFile,$bgfFile,$trjType,$tSel,$saveFile,$FF, $tStep, $orientOpt) =
        ($OPTS{a},$OPTS{t},$OPTS{b},$OPTS{l},$OPTS{s},$OPTS{w},$OPTS{f}, $OPTS{d}, $OPTS{o});
    
    FileTester($trjFile);
    FileTester($bgfFile);

    die "ERROR: Expected number for timestep, got \"$tStep\"\n"
	if ($tStep !~ /^\d+\.?\d*$/);

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
    $SELECT = TrjSelections($tSel);
    
    die "ERROR: No valid frames selected with selection $tSel!\n"
	if (! keys %{ $SELECT } and $tSel ne "*");
    
    if (-e $atomSel) {
	if (open(INFILE, $atomSel)) {
	    while (<INFILE>) {
		chomp;
		while ($_ =~ /(\S+)/g) {
		    push @tmp, $1;
		}
	    }
	    close INFILE;
	}
    } elsif ($atomSel =~ /\s+/) {
	@tmp = split /\s+/, $atomSel;
    } else {
	$tmp[0] = $atomSel;
    }
    $ATOMSELECT = GetSelections(\@tmp, 0);
    
    die "ERROR: No valid atoms selected for SOLUTE with selection $atomSel!\n"
	if (! keys %{ $ATOMSELECT });

    if (! defined($saveFile)) {
        $saveFile = basename($bgfFile);
        $saveFile =~ s/\.\w+$//;
	$saveFile .= ".rdf";
    }
    
    $saveFile =~ s/\.\w+$//;

    print "Done\nParsing BGF file $bgfFile...";
    ($BGF, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
    $ATOMSELECT = GetAtmList($ATOMSELECT, $BGF);
    $MOLS = SplitAtomsByMol($BGF, $ATOMSELECT);
    $totMols = scalar keys %{ $MOLS };
    $totAtms = scalar keys %{ $ATOMSELECT };
    die "ERROR: No molecule matches selection!\n" if (! $totMols);

    if (defined($orientOpt)) {
	$ORIENT = getOrientOpt($BGF, $MOLS, $orientOpt);
    } else {
	$ORIENT->{T} = "dM";
    }

    print "Done\n";

    if ($FF =~ /\s/) {
	@tmp = split /\s+/, $FF;
    } else {
	$tmp[0] = $FF;
    }
    $list = ();
    for $i (@tmp) {
	push @{ $list }, $i if (-e $i);
    }
}

sub showUsage {
    my ($usage) = "usage: $0 -b bgf file -t trajectory file -a atoms (see selection info)" . 
	"\n\t\t\t\t\t\t-f cerius2 forcefield -d trajectory dump timestep (ps)" .
        "\nOptional parameters:" .
	"\n\t-l traj type (amber(default)|lammps -w save name  -s trajectory selection (see selection info)\n" . 
	"\n\t-o orientation options (none default): For the molecule, expected 1_2:uM for correlation of 1->2 vector\n" .
	"\t\t labeled uM or {name1}_x_uM:p for the vector cross product of predefined name1 (e.g. 1->2) and uM called P\n" . 
	"\t\t Use dm for dipole moment. Enclose multiple in quotes\n" .
	&ShowSelectionInfo;
    return $usage;
}

sub getOrientOpt {
    my ($atoms, $mols, $orientOpts) = @_;
    my (%OPTS, $i, $currOpt, $molMap, $offset, $opt1, $opt2, %OPTLIST, $optName, @tmp);

    @tmp = sort numerically keys %{ $mols->{1} };
    for $i (@tmp) {
	$offset = $i - $tmp[0] + 1;
	$molMap->{$offset} = $i;
    }

    $OPTS{T} = "dM";
    $OPTLIST{"dm"} = "0_-1";
    $orientOpts = lc($orientOpts);

    while ($orientOpts =~ /(\S+)/g) {
	$currOpt = $1;
	if ($currOpt =~ /(\w+)_x_(\w+)\:(\w+)/) {
	    ($opt1, $opt2, $optName) = ($1, $2, $3);
	    for ($opt1, $opt2) {
		if ($_ =~ /(\d+)_(\d+)/) {
		    if (exists($molMap->{$1}) and exists($molMap->{$2}) and $1 != $2) {
			$OPTLIST{"${1}_${2}"} = 1;
		    }
		}
	    }
	    if (exists($OPTLIST{$opt1}) and exists($OPTLIST{$opt2})) {
		$opt1 = $OPTLIST{$opt1} if ($OPTLIST{$opt1} ne "1");
		$opt2 = $OPTLIST{$opt2} if ($OPTLIST{$opt2} ne "1");
		$OPTS{$optName} = "${opt1}_x_${opt2}";
	    }
	} elsif ($currOpt =~ /(\d+)_(\d+)\:(\w+)/) {
	    $optName = $3;
	    undef($opt1);
	    if (exists($molMap->{$1}) and exists($molMap->{$2}) and $1 != $2) {
		$OPTLIST{"${1}_${2}"} = 1;
		$OPTLIST{$optName} = "${1}_${2}";
		$opt1 = "${1}_${2}";
	    }
	    if (defined($opt1) and exists($OPTLIST{$opt1})) {
		$OPTS{$optName} = "${1}_${2}";
	    }
	}
    }

    return \%OPTS;
}

sub getVecs {
    my ($molMap, $startCom, $molData, $k, $l) = @_;
    my ($a1, $a2, $j, $v1, $v2);

    $a1 = $molMap->{$k};
    $a2 = $molMap->{$l};
    for $j ("XCOORD", "YCOORD", "ZCOORD") {
	$v1->{$j} = $startCom->{$a1}{$j} - $startCom->{$a2}{$j};
	$v2->{$j} = $molData->{$a1}{$j} - $molData->{$a2}{$j};
    }
    
    return ($v1, $v2);
}

#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use warnings;
no warnings "recursion";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo AddMass);
use General qw(FileTester CoM GetSelections TrjSelections 
			 ShowSelectionInfo DotProduct STDev LoadFFs);
use AMBER qw(ParseAmberTrj GetAmberByteOffset);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType);
use BOX qw(ConvertBox MakeBox);
use ManipAtoms qw(UnwrapAtoms ImageAtoms GetAtmList GetSolvent SplitAtomsByMol);
use Math::FFT;
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
sub AutoCorrelation;

my ($MOLS, $trjFile, $BGF, $BONDS, $printStr, $PARMS, $list, $i);
my ($getSnapshot, $getByteOffset, $trjType, $LAMMPSOPTS, $tStep, $bulkP);
my ($SELECT, $saveFile, $field, $DATA, $HEADERS, $totMass, $totFrames);

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
$getSnapshot->($BGF, $trjFile, $SELECT, $field, \&saveProperties, $printStr, $MOLS);

print "Writing data...\r";
$bulkP = calcProperties($DATA, $saveFile);

for $i ("DIELECTRIC_CONSTANT", "DIFFUSION_CONSTANT", "DENSITY") {
    printf "%-20s: %12.5G $bulkP->{$i}{UNITS} +/- %8.3G\n", $i, $bulkP->{$i}{AVG}, $bulkP->{$i}{STDEV};
}

sub AutoCorrelation {
    my ($data, $normFactor, $constant) = @_;
    my (@AutoCorrelationFunc, $i);

    $normFactor = 1 if (! defined($normFactor));
    $constant = 0 if (! defined($constant));
    for $i (0 .. $#{ $data }) {
	$AutoCorrelationFunc[$i] = DotProduct($data->[0], $data->[$i]) - $constant;
	$AutoCorrelationFunc[$i] /= $normFactor;
    }

    return \@AutoCorrelationFunc;
}

sub calcProperties {
    my ($data, $prefix) = @_;
    my ($i, $t, $dim, $M, $correlation, $spectrum, $fft, $kbT, $totMols, $tmp, $rFFT, $cFFT); 
    my (%options, $D, $avg, $stdev, $STATS, $M_squared, $M_bar, $delM, $V, $vals, $time, $norm);

    $totFrames--;
    $M_squared = $M_bar = 0;
    $kbT = 3.1668152e-06 * 300; # boltzmann constant in atomic units
    $STATS->{DENSITY}{UNITS} = "g/cm3";
    $STATS->{DIELECTRIC_CONSTANT}{UNITS} = "";
    $STATS->{DIFFUSION_CONSTANT}{UNITS} = "cm2/s";
    $totMols = scalar keys %{ $data->{DIPOLE_MOMENT} };

    print "Writing data...density\r";
    open DENSITY, "> ${prefix}_density.dat" or die "ERROR: Cannot create ${prefix}_density.dat: $!\n";
    for $i (0 .. $totFrames) {
        printf DENSITY "%8.5f %12.5f\n", ($i + 1), ($totMass/($data->{VOLUME}[$i] * 0.6023));
        $vals .= ($totMass/($data->{VOLUME}[$i] * 0.6023)) . " ";
        $V += $data->{VOLUME}[$i];
    }
    chop $vals;
    chop $V;
    ($STATS->{DENSITY}{AVG}, $STATS->{DENSITY}{STDEV}, undef) = STDev($vals);
    $vals = "";
    printf DENSITY "#%7s %12.5f %8.3f\n", "STATS", $STATS->{DENSITY}{AVG}, $STATS->{DENSITY}{STDEV};
    close DENSITY;
    $V /= $totFrames;

    print "Writing data...collecting data\r";

    for $i (keys %{ $data->{DIPOLE_MOMENT} }) {
        for $t (0 .. $totFrames) {
	    $D->[$t] .= "$data->{DISPLACEMENT}{$i}[$t] " if ($t < $totFrames);
	    $M_squared += $data->{DIPOLE_MOMENT}{$i}[$t]{TOT}**2;
	    $M_bar += $data->{DIPOLE_MOMENT}{$i}[$t]{TOT};
	    
	}
	$norm = DotProduct($data->{DIPOLE_MOMENT}{$i}[0],$data->{DIPOLE_MOMENT}{$i}[0]);
        $tmp = AutoCorrelation($data->{DIPOLE_MOMENT}{$i}, $norm);
	for $i (0 .. $#{ $tmp }) {
	    $correlation->[$i] += $tmp->[$i];
	}
    }

    print "Writing data...diffusion             \r";
    $vals = "";
    open DIFFUSION, "> ${prefix}_diffusion.dat" or die "ERROR: Cannot create ${prefix}_diffusion.dat: $!\n";
    for $i (0 .. ($totFrames - 1)) {
        chop $D->[$i];
        ($avg, $stdev, undef) = STDev($D->[$i]);
        $time = $tStep * ($data->{TIMESTEP}[$i + 1] - $data->{TIMESTEP}[0]);
        printf DIFFUSION "%8.5f %12.5f %8.3f\n", $time, $avg, $stdev;
        $avg /= $time * 2;
        $vals .= "$avg ";

    }
    chop $vals;
    ($STATS->{DIFFUSION_CONSTANT}{AVG}, $STATS->{DIFFUSION_CONSTANT}{STDEV}, undef) = STDev($vals);
    printf DIFFUSION "#%7s %12.5f %8.3f\n","STATS",$STATS->{DIFFUSION_CONSTANT}{AVG}, $STATS->{DIFFUSION_CONSTANT}{STDEV};
    $STATS->{DIFFUSION_CONSTANT}{AVG} *= 1e-04;
    $STATS->{DIFFUSION_CONSTANT}{STDEV} *= 1e-04;
    close DIFFUSION;
    $D = ();

    print "Writing data...dipole autocorrelation\r";

    $M_squared /= ($totFrames + 1);
    $M_bar /= $totMols * ($totFrames + 1);
    $delM = $M_squared - $M_bar**2;
    #$correlation = &AutoCorrelation($M, $M_bar**2, $delM);
    $V *= 1.889726**3; #volume in atomic units
    $STATS->{DIELECTRIC_CONSTANT}{AVG} = 4 * PI * $delM/(3 * $V * $kbT);
    $STATS->{DIELECTRIC_CONSTANT}{STDEV} = $M_bar;
    $M = ();

    open CORRFILE, "> ${prefix}_dipole_corr.dat" or die "ERROR: Cannot create ${prefix}_dipole_corr.dat: $!\n";
    for $i (0 .. $totFrames) {
        $time = $tStep * ($data->{TIMESTEP}[$i] - $data->{TIMESTEP}[0]);
	$correlation->[$i] /= $totMols;
	print CORRFILE "$time $correlation->[$i]\n";
        #$correlation->[$i] /= ($totMols * $delM);
        #$correlation->[$i] -= $M_bar**2;
    }
    close CORRFILE;

    print "Writing data...dipole fourier transform        \r";
    $fft = new Math::FFT($correlation);
    $spectrum = $fft->rdft();
    open rFFT, "> ${prefix}_rfft.dat" or die "ERROR: Cannot create ${prefix}_rfft.dat: $!\n"; # real part of fft
    open cFFT, "> ${prefix}_cfft.dat" or die "ERROR: Cannot create ${prefix}_cfft.dat: $!\n"; # complex part of fft
    for $i (0 .. ($#{ $spectrum } - 1)/2) {
	$rFFT->[$i] = $spectrum->[2*$i];
	$cFFT->[0] = 0;
	print rFFT "$i $rFFT->[$i]\n";
	next if ($i == 0);
	$cFFT->[$i] = $spectrum->[2*$i + 1];
	print cFFT "$i $cFFT->[$i]\n";
    }
    close rFFT;
    close cFFT;
    $spectrum = ();
    $correlation = ();

    print "Writing data...frequency dependent dipole \r";
    open DIPOLEFREQ, "> ${prefix}_dipole_freq.dat" or die "ERROR: Cannot create ${prefix}_dipole_freq.dat: $!\n";
    for $i (0 .. $#{ $cFFT }) {
	$vals = $STATS->{DIELECTRIC_CONSTANT}{AVG} - 
		($STATS->{DIELECTRIC_CONSTANT}{AVG} - 1) * 
		$i * ($rFFT->[$i] + $cFFT->[$i]);
        print DIPOLEFREQ "$i $vals\n";
    }
    close DIPOLEFREQ;	
    $cFFT = ();
    $rFFT = ();


    print "Writing data...Done                        \n";

    return $STATS;
}

sub saveProperties {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($CENTER, @tmp, $i, $j, $MOL, $dipole, $volume, $index);

    if ($trjType == 2) { #LAMMPS
        $frameNum = $ATOMS->{TIMESTEP}[0]/1000; #timestep in ps
        $BOX = $ATOMS->{"BOX BOUNDS"};
	$volume = ($BOX->[0]{hi} - $BOX->[0]{lo}) * ($BOX->[1]{hi} - $BOX->[1]{lo}) *
		  ($BOX->[2]{hi} - $BOX->[2]{lo});
        $ATOMS = $ATOMS->{ATOMS};
	$BOX = MakeBox(ConvertBox($BOX));
        for $i (keys %{ $ATOMS }) {
            for $j ("MASS", "CHARGE") {
                $ATOMS->{$i}{$j} = $_[1]->{$i}{$j};
            }
	    next if (! $LAMMPSOPTS->{scaled} and ! $LAMMPSOPTS->{imaged});
  	    for $j ("XCOORD", "YCOORD", "ZCOORD") {
	        $index = $j;
	        $index =~ s/COORD/INDEX/;
	        $index = $ATOMS->{$i}{$index};
		if ($LAMMPSOPTS->{scaled}) {
		    $ATOMS->{$i}{$j} *= $BOX->{$j}{len};
		    $ATOMS->{$i}{$j} += $BOX->{$j}{lo};
		}
		$ATOMS->{$i}{$j} += ($index * $BOX->{$j}{len}) if ($LAMMPSOPTS->{imaged});
	    }
        }
    } else {
        $volume = $BOX->{2}{DATA} * $BOX->{3}{DATA} * $BOX->{4}{DATA};
    }

    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    for $i (keys %{ $MOLS }) {
	$MOL = getAtoms($ATOMS, $MOLS->{$i});
	$CENTER = CoM($MOL);
	$dipole = &calcDipoleMoment($MOL, $CENTER);
	push @{ $DATA->{DIPOLE_MOMENT}{$i} }, $dipole;

	if (exists($DATA->{START_COM}{$i})) {
	    push @{ $DATA->{DISPLACEMENT}{$i} }, &calcDisplacement($CENTER, $DATA->{START_COM}{$i});
	} else {
	    $DATA->{START_COM}{$i} = $CENTER;
	}
    }
    
    push @{ $DATA->{VOLUME} }, $volume;
    push @{ $DATA->{TIMESTEP} }, $frameNum;
}

sub calcDipoleMoment {
    my ($atoms, $com) = @_;
    my (@dim, $i, $j, %dipole);
    
    @dim = ("XCOORD", "YCOORD", "ZCOORD");
    
    for $i (keys %{ $atoms }) {
	for $j (@dim) {
	    #$dipole{$j} += ($atoms->{$i}{$j} - $com->{$j}) * $atoms->{$i}{CHARGE};
	    $dipole{$j} += ($atoms->{$i}{$j}) * $atoms->{$i}{CHARGE};
	}
    }

    for $j (@dim) {
	$dipole{TOT} += $dipole{$j}**2;
    }
    $dipole{TOT} = sqrt($dipole{TOT});

    return \%dipole;
}

sub calcDisplacement {
    my ($com, $prevCom) = @_;
    my ($i, @displacement, $totDist);

    for $i ("XCOORD", "YCOORD", "ZCOORD") {
	$totDist += ($com->{$i} - $prevCom->{$i})**2;
    }
    return $totDist;
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
    my (%OPTS, @tmp, $i, $j, $solvTmp, $count);
    my ($atomSel, $bgfFile, $tSel, $FF, $ATOMSELECT);
    
    getopt('btalswfd',\%OPTS);
    
    for ("b", "t", "a", "f", "d") {
	die &showUsage . "" if (! defined($OPTS{$_}));
    }
    
    print "Initializing...";
    ($atomSel,$trjFile,$bgfFile,$trjType,$tSel,$saveFile,$FF, $tStep) =
        ($OPTS{a},$OPTS{t},$OPTS{b},$OPTS{l},$OPTS{s},$OPTS{w},$OPTS{f}, $OPTS{d});
    
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
    
    
    if ($atomSel =~ /\s+/) {
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
	&ShowSelectionInfo;
    return $usage;
}

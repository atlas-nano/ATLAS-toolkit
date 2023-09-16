#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use warnings;
no warnings "recursion";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo AddMass);
use General qw(FileTester CoM GetSelections TrjSelections ShowSelectionInfo 
			 DotProduct CrossProduct LoadFFs GetBondLength CenterText
			 GetStats GetEquilPoint);
use AMBER qw(ParseAmberTrj GetAmberByteOffset);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset 
			GetLammpsTrjType ParseLAMMPSLogFile);
use BOX qw(ConvertBox MakeBox);
use ManipAtoms qw(UnwrapAtoms GetAtmList GetSolvent SplitAtomsByMol);
use Bulk qw(CalcMolOrientCorrelation CalcDipole CalcDielectricConstant CalcDiffusionConstant CalcTExpand
		      GetOrientOpt CalcHVap CalcHCap CalcTCompress CalcTExpand GetFluct WriteStats WriteData SetOpts);
use constant PI => atan2(1,1) * 4;

sub saveTrjData;
sub saveLogData;
sub getAtoms;
sub getTotMass;
sub determineNumFrames;
sub init;
sub showUsage;
sub LoadBGF;
sub numerically { ($a<=>$b); }

my ($MOLS, $getByteOffset, $tStep, $saveFile, $logFile, $OPTS, $isVariable, $sysTemp);
my ($list, $trjType, $LAMMPSOPTS, $printStr, $getSnapshot, $BULK, $tStart);
my ($PARMS, $BGF, $totMass, $field, $SELECT, $trjFile, $totFrames, $logSelect);
my ($totMols, $totAtms, $ORIENT, $LOGDATA, $atomSelect, $bgfFile, $orientOpt);

$|++;
&SetOpts(\%{ $BULK });
$atomSelect = &init;
$PARMS = LoadFFs($list);
print "Parsing BGF file $bgfFile...";
&LoadBGF($bgfFile, $atomSelect, $orientOpt);
&AddMass($BGF, $PARMS);
$totMass = getTotMass($BGF, $MOLS);
$field = scalar keys %{ $BGF };
print "Done\n";
$getByteOffset->($SELECT, $trjFile, $field);
#$totFrames = &determineNumFrames($SELECT);
$totFrames = scalar keys %{ $SELECT };
if ($trjType == 2) {
    &GetLammpsTrjType($SELECT, $trjFile, "atom", \%{ $LAMMPSOPTS });
    $field = "atom";
}
$printStr = "Calculating bulk properties from $trjFile...";
$getSnapshot->($BGF, $trjFile, $SELECT, $field, \&saveTrjData, $printStr, $OPTS);
if (defined($logFile)) {
    print "Calculating bulk properties from log file $logFile...";
    undef($tStart);
    %{ $SELECT } = %{ $logSelect } if (defined($logSelect));
    &ParseLAMMPSLogFile($logFile, $SELECT, \&saveLogData, $OPTS);
    print "Done\n";
}
print "Writing data...";
&WriteData($BULK, $saveFile, 1);
print "Done\n";
&WriteStats($BULK, "${saveFile}.dat");

sub saveLogData {
    my ($data, $frameNum, $bulkOpts) = @_;
    my ($delT);

    $frameNum /= 1000;
    $tStart = $frameNum if (! defined($tStart));
    $delT = $tStep * ($frameNum - $tStart);

    $BULK->{THERMO_FLUCT} = () if (! exists($BULK->{THERMO_FLUCT}));
    &GetFluct($data, \%{ $BULK->{THERMO_FLUCT} });
	
    if (exists($bulkOpts->{HEAT_VAPORIZATION})) {
	$BULK->{HEAT_VAPORIZATION}{tStep}{$frameNum} = CalcHVap($BULK->{THERMO_FLUCT}, $totMols);
    }
	
    if (exists($bulkOpts->{HEAT_CAPACITY})) {
	$BULK->{HEAT_CAPACITY}{tStep}{$frameNum} = CalcHCap($BULK->{THERMO_FLUCT});
    }
	
    if (exists($bulkOpts->{THERMAL_COMPRESSIBILITY}) and exists($data->{volume})) {
	$BULK->{THERMAL_COMPRESSIBILITY}{tStep}{$frameNum} = CalcTCompress($BULK->{THERMO_FLUCT});
    }
	
    if (exists($bulkOpts->{THERMAL_EXPANSION}) and exists($data->{volume})) {
	$BULK->{THERMAL_EXPANSION}{tStep}{$frameNum} = CalcTExpand($BULK->{THERMO_FLUCT});
    }
}

sub saveTrjData {
    my ($ATOMS, $BOX, $frameNum, $bulkOpts) = @_;
    my ($CENTER, $i, $j, $volume, $MOLDATA);
    my ($DIPOLE, $DISPLACEMENT, $index, $delT);
    
    
    if ($trjType == 2) { #LAMMPS
	$bulkOpts = $frameNum;
        $frameNum = $ATOMS->{TIMESTEP}[0]; 
        $BOX = $ATOMS->{"BOX BOUNDS"};
	$volume = ($BOX->[0]{hi} - $BOX->[0]{lo}) * ($BOX->[1]{hi} - $BOX->[1]{lo}) *
		  ($BOX->[2]{hi} - $BOX->[2]{lo});
        $ATOMS = $ATOMS->{ATOMS};
	$BOX = MakeBox(ConvertBox($BOX));
        if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
            UnwrapAtoms($ATOMS, $BOX, $LAMMPSOPTS->{scaled});
        }

        for $i (keys %{ $ATOMS }) {
            for $j ("MOLECULE", "MOLECULEID", "MASS", "CHARGE") {
	 	next if(exists($ATOMS->{$i}{$j}));
                $ATOMS->{$i}{$j} = $_[1]->{$i}{$j};
            }
        }
    } else {
        $volume = $BOX->{2}{DATA} * $BOX->{3}{DATA} * $BOX->{4}{DATA};
    }

    if ($isVariable) {
	$MOLS = SplitAtomsByMol($ATOMS);
	$totMols = scalar keys %{ $MOLS };
	$totAtms = scalar keys %{ $ATOMS };
	$totMass = getTotMass($ATOMS, $MOLS);
    }

    $tStart = $frameNum if (! defined($tStart));
    $delT = $tStep * ($frameNum - $tStart); #time in ps
    $logSelect->{$frameNum} = 1;
    $BULK->{VOLUME}{tStep}{$delT} = $volume;
    $BULK->{TEMPERATURE}{tStep}{$delT} = $sysTemp;
    $BULK->{DENSITY}{tStep}{$delT}{T} = $totMass/(0.6023*$volume) if (exists($bulkOpts->{DENSITY}));
    if (exists($bulkOpts->{DIFFUSION_CONSTANT})) {
	$BULK->{DIFFUSION_CONSTANT}{DAT}{COM} = () if (! exists($BULK->{DIFFUSION_CONSTANT}{DAT}));
	$BULK->{DIFFUSION_CONSTANT}{tStep}{$delT} = 
	    &CalcDiffusionConstant($ATOMS, $MOLS, \%{ $BULK->{DIFFUSION_CONSTANT}{DAT} }, $delT, $isVariable);
    }
    if (exists($bulkOpts->{DIELECTRIC_CONSTANT})) {
	$BULK->{DIELECTRIC_CONSTANT}{DAT} = () if (! exists($BULK->{DIELECTRIC_CONSTANT}{DAT}));
	$BULK->{DIELECTRIC_CONSTANT}{tStep}{$frameNum} = 
	    &CalcDielectricConstant($ATOMS, $MOLS, \%{ $BULK->{DIELECTRIC_CONSTANT}{DAT} }, $BULK->{VOLUME}, $BULK->{TEMPERATURE});
	$BULK->{DIELECTRIC_CONSTANT}{tStep}{$frameNum}{u} =
	    &CalcDipole($ATOMS, $MOLS, \%{ $BULK->{DIELECTRIC_CONSTANT}{DAT} }, 1);
    }

    if (exists($bulkOpts->{MOL_CORRELATION})) {
	$BULK->{MOL_CORRELATION}{DAT}{DIPOLE} = () if (! exists($BULK->{MOL_CORRELATION}{DAT}));
	$BULK->{MOL_CORRELATION}{tStep}{$frameNum} = 
	    &CalcMolOrientCorrelation($ATOMS, $MOLS, \%{ $BULK->{MOL_CORRELATION}{DAT} }, $ORIENT, $isVariable);
    }
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
    
sub LoadBGF {
    my ($bgfFile, $ATOMSELECT, $orientOpt) = @_;
    my ($BONDS, $HEADERS);

    ($BGF, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
    $ATOMSELECT = GetAtmList($ATOMSELECT, $BGF);
    $MOLS = SplitAtomsByMol($BGF, $ATOMSELECT);
    $totMols = scalar keys %{ $MOLS };
    $totAtms = scalar keys %{ $ATOMSELECT };
    die "ERROR: No molecule matches selection!\n" if (! $totMols);

    #molecular orientation
    if (defined($orientOpt)) {
	$ORIENT = GetOrientOpt($BGF, $MOLS, $orientOpt);
	$OPTS->{MOL_CORRELATION} = 1 if ($ORIENT);
    }
}

sub init {
    my (%OPTS, @tmp, $i, $j, $solvTmp, $count, $tSel);
    my ($atomSel, $FF, $ATOMSELECT, $ensemble);

    getopt('btalrswfdoevk',\%OPTS);
    
    for ("b", "t", "f") {
	die &showUsage . "" if (! defined($OPTS{$_}));
    }
    print "Initializing...";
    ($atomSel,$trjFile,$bgfFile,$trjType,$tSel,$saveFile,$FF, $tStep, $orientOpt, $ensemble, $logFile, $isVariable,$sysTemp) =
        ($OPTS{a},$OPTS{t},$OPTS{b},$OPTS{r},$OPTS{s},$OPTS{w},$OPTS{f}, $OPTS{d}, $OPTS{o}, $OPTS{e}, $OPTS{l}, $OPTS{v}, $OPTS{k});
    
    #defaults
    FileTester($trjFile);
    FileTester($bgfFile);
    if ($FF =~ /\s/) {
	@tmp = split /\s+/, $FF;
    } else {
	$tmp[0] = $FF;
    }
    $list = ();
    for $i (@tmp) {
	push @{ $list }, $i if (-e $i);
    }
    die "ERROR: No valid forcefield found while searching \"$FF\"!\n" if (! @{ $list });

#   OPTIONALS
    
    #temperature
    $sysTemp = 300 if (! defined($sysTemp) or $sysTemp !~ /^\d+\.?\d*$/);
    $sysTemp = $1 if ($sysTemp =~ /^(\d+\.?\d*)$/ and $1 > 0);

    #variable atoms
    $isVariable = 0 if (! defined($isVariable));
    $isVariable = 1 if ($isVariable =~ /^(1|yes)/i);


    #timestep
    $tStep = 1 if (! defined($tStep));
    die "ERROR: Expected number for timestep, got \"$tStep\"\n" if ($tStep !~ /^\d+\.?\d*$/);

    # trajectory type
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
	$tStep  = 1;
    } else {
	$trjType = 2;
	$getSnapshot = \&ParseLAMMPSTrj;
	$getByteOffset = \&GetLammpsByteOffset;
    }

    #snapshot selection
    $tSel = "*" if (! defined($tSel));
    $SELECT = TrjSelections($tSel); 
    die "ERROR: No valid frames selected with selection $tSel!\n" if (! keys %{ $SELECT } and $tSel ne "*");
    
    #atom selection
    $atomSel = "*" if (! defined($atomSel));
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
    die "ERROR: No valid atoms selected for SOLUTE with selection $atomSel!\n" if (! keys %{ $ATOMSELECT });

    #save file
    if (! defined($saveFile)) {
        $saveFile = basename($bgfFile);
        $saveFile =~ s/\.\w+$//;
	$saveFile .= ".rdf";
    } 
    $saveFile =~ s/\.\w+$//;

    #log file
    undef($logFile) if (defined($logFile) and (! -e $logFile or ! -r $logFile or ! -T $logFile));

    #ensemble
    $ensemble = "npt" if (! defined($ensemble));
    $ensemble = lc($ensemble);
    if ($ensemble =~ /(npt|nvt|nve|all)/) {
	$ensemble = $1;
    } else {
	$ensemble = "npt";
    }
    if ($ensemble eq "nve") {
	$OPTS->{DIFFUSION_CONSTANT} = 1;
    } elsif ($ensemble eq "all") {
        $OPTS->{DIFFUSION_CONSTANT} = 1;
        $OPTS->{DENSITY} = 1;
        #$OPTS->{DIELECTRIC_CONSTANT} = 1;
        if (defined($logFile)) {
            $OPTS->{HEAT_VAPORIZATION} = 1;
            $OPTS->{HEAT_CAPACITY} = 1;
            $OPTS->{THERMAL_COMPRESSIBILITY}  = 1;
            $OPTS->{THERMAL_EXPANSION} = 1;
        }
    } else {
	$OPTS->{DENSITY} = 1;
	#$OPTS->{DIELECTRIC_CONSTANT} = 1;
	if (defined($logFile)) {
	    $OPTS->{HEAT_VAPORIZATION} = 1;
	    $OPTS->{HEAT_CAPACITY} = 1;
	    $OPTS->{THERMAL_COMPRESSIBILITY}  = 1;
	    $OPTS->{THERMAL_EXPANSION} = 1;
	}
    }

    print "Done\n";
    return $ATOMSELECT;
}

sub showUsage {
    my ($usage) = "usage: $0 -b bgf file -f force field -t trajectory file" . 
        "\nOptional parameters:" .
	"\n\t-a atoms (see selection info) -d md timestep (ps)" .	
	"\n\t-e ensemble (npt(default)|nvt|nve) -r traj type (amber(default)|lammps -w save name" .
	"\n\t-s trajectory selection (see selection info) -l lammps log file -k temperature (default 300)" . 
	"\n\t-o orientation options (none default): For the molecule, expected 1_2:uM for correlation of 1->2 vector" .
	"\n\t\t labeled uM or {name1}_x_uM:p for the vector cross product of predefined name1 (e.g. 1->2) and uM called P" . 
	"\n\t\t Use dm for dipole moment. Enclose multiple in quotes" .
	"\n\t-v is variable: Set to 1 if trajectory contains varying # atoms. Default 0\n" .
	&ShowSelectionInfo;
    return $usage;
}

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
use constant PI => atan2(1,1) * 4;

sub writeData;
sub saveTrjData;
sub saveLogData;
sub getAtoms;
sub getTotMass;
sub determineNumFrames;
sub init;
sub showUsage;
sub LoadBGF;
sub calcDielectricConstant;
sub calcDiffusionConstant;
sub calcMolOrientCorrelation;
sub calcDipole;
sub numerically { ($a<=>$b); }
sub calcHVap;
sub calcHCap;
sub calcTCompress;
sub calcTExpand;
sub getFluct;
sub writeStats;
sub setOpts;
sub getVecs; 
sub getOrientOpt;

my ($MOLS, $getByteOffset, $tStep, $saveFile, $logFile, $OPTS, $isVariable);
my ($list, $trjType, $LAMMPSOPTS, $printStr, $getSnapshot, $BULK, $tStart);
my ($PARMS, $BGF, $totMass, $field, $SELECT, $trjFile, $totFrames, $logSelect);
my ($totMols, $totAtms, $ORIENT, $LOGDATA, $atomSelect, $bgfFile, $orientOpt);

$|++;
&setOpts;
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
writeData($BULK, $saveFile);
print "Done\n";
&writeStats($BULK, "${saveFile}.dat");

sub writeStats {
    my ($stats, $saveName) = @_;
    my ($outStr, $i, $header);

    $outStr = sprintf("%-45s%21s%21s\n", "BULK PROPERTY",
		      CenterText("Calc",21),CenterText("Exp(300K 1atm)",21));

    for $i (1 .. 90) {
	$outStr .= "=";
    }
    $outStr .= "\n";

    for $i (keys %{ $stats }) {
	if ($stats->{$i}{UNITS}) {
	    $header = sprintf("%-45s", $stats->{$i}{STATS}{T}{CONVERGED} . "${i} (" . $stats->{$i}{UNITS} . ")");
	} else {
	    $header = sprintf("%-45s", $stats->{$i}{STATS}{T}{CONVERGED} . $i);
	}
	$outStr .= $header . 
	    sprintf("%8.3f +/- %-8.3f%8.3f +/- %-8.3f\n", $stats->{$i}{STATS}{T}{AVG}, 
		    $stats->{$i}{STATS}{T}{STDEV}, $stats->{$i}{EXP_A}, $stats->{$i}{EXP_D});
    }
    $outStr .= "note: * indicates converged result not obtained\n";
    print $outStr;
    open OUTDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    print OUTDATA $outStr;
    close OUTDATA;
}

sub writeData {
    my ($data, $savePrefix) = @_;
    my ($type, $fileName, @index, @headers, $startStats, $num_str);
    my ($i, $j, $k, @tmp, $avg_str, $stdev_str, $count, $sData);

    delete $data->{VOLUME};
    for $type (keys %{ $data }) {
	if (! keys %{ $data->{$type}{tStep} }) {
	    delete $data->{$type};
	    next;
	}
	$count = 0;
	@index = sort numerically keys %{ $data->{$type}{tStep} };
	@headers = sort {$a cmp $b } keys %{ $data->{$type}{tStep}{$index[0]} };
	$fileName = $savePrefix . "_" . lc($type) . ".dat";
        open DATA, "> $fileName" or die "ERROR: Cannot create $fileName: $!\n";
        printf DATA "#%8s ", "TIME(ps)";
	for $j (@headers) {
            printf DATA "%12s ", $j;
        }
        print DATA "\n";
        for $i (@index) {
	    $count++;
            printf DATA "%9.3f ", $i;
            for $j (@headers) {
                printf DATA "%12.5f ", $data->{$type}{tStep}{$i}{$j};
	    }
	    print DATA "\n";
	}

	$avg_str = sprintf("#%7s ", "AVG");
	$stdev_str = sprintf("#%7s ", "STDEV");
	$num_str = sprintf("#%7s ", "num_pts");
	for $j (@headers) {
	    ($k, $sData) = GetEquilPoint($data->{$type}{tStep}, 0.02, $j);
	    if (! $k) {
		$startStats = sprintf("%.0f", (0.66 * scalar(@index)));
	    } else {
		$startStats = 1;
	    }
	    $data->{$type}{STATS}{$j} = GetStats($sData, $startStats);
	    $data->{$type}{STATS}{$j}{CONVERGED} = "";
	    $data->{$type}{STATS}{$j}{CONVERGED} = "*" if (! $k);
	    $avg_str .= sprintf("%12.5G ",$data->{$type}{STATS}{$j}{AVG} );
	    $stdev_str .= sprintf("%12.8G ", $data->{$type}{STATS}{$j}{STDEV});
	    $num_str .= sprintf("%11d%1s ", $data->{$type}{STATS}{$j}{NUM}, $data->{$type}{STATS}{$j}{CONVERGED});
	}
	print DATA "${avg_str}\n${stdev_str}\n${num_str}\n";
	close DATA;
    }
}

sub saveLogData {
    my ($data, $frameNum, $bulkOpts) = @_;
    my ($delT);

    $frameNum /= 1000;
    $tStart = $frameNum if (! defined($tStart));
    $delT = $tStep * ($frameNum - $tStart);

    if (exists($data->{volume})) {
	
	$BULK->{THERMO_FLUCT} = () if (! exists($BULK->{THERMO_FLUCT}));
	&getFluct($data, \%{ $BULK->{THERMO_FLUCT} });
	
	if (exists($bulkOpts->{HEAT_VAPORIZATION})) {
	    $BULK->{HEAT_VAPORIZATION}{tStep}{$frameNum} = calcHVap($BULK->{THERMO_FLUCT});
	}
	
	if (exists($bulkOpts->{HEAT_CAPACITY})) {
	    $BULK->{HEAT_CAPACITY}{tStep}{$frameNum} = calcHCap($BULK->{THERMO_FLUCT});
	}
	
	if (exists($bulkOpts->{THERMAL_COMPRESSIBILITY})) {
	    $BULK->{THERMAL_COMPRESSIBILITY}{tStep}{$frameNum} = calcTCompress($BULK->{THERMO_FLUCT});
	}
	
	if (exists($bulkOpts->{THERMAL_EXPANSION})) {
	    $BULK->{THERMAL_EXPANSION}{tStep}{$frameNum} = calcTExpand($BULK->{THERMO_FLUCT});
	}
    }
}

sub saveTrjData {
    my ($ATOMS, $BOX, $frameNum, $bulkOpts) = @_;
    my ($CENTER, $i, $j, $volume, $MOLDATA);
    my ($DIPOLE, $DISPLACEMENT, $index, $delT);
    
    
    if ($trjType == 2) { #LAMMPS
	$bulkOpts = $frameNum;
        $frameNum = $ATOMS->{TIMESTEP}[0]/1000; #timestep in ps
	$tStep = 1;
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
    $delT = $tStep * ($frameNum - $tStart);
    $logSelect->{($frameNum *1000)} = 1;
    $BULK->{VOLUME}{tStep}{$delT} = $volume;
    $BULK->{DENSITY}{tStep}{$delT}{T} = $totMass/(0.6023*$volume) if (exists($bulkOpts->{DENSITY}));
    if (exists($bulkOpts->{DIFFUSION_CONSTANT})) {
	$BULK->{DIFFUSION_CONSTANT}{DAT}{COM} = () if (! exists($BULK->{DIFFUSION_CONSTANT}{DAT}));
	$BULK->{DIFFUSION_CONSTANT}{tStep}{$delT} = 
	    &calcDiffusionConstant($ATOMS, $MOLS, \%{ $BULK->{DIFFUSION_CONSTANT}{DAT} }, $delT);
    }
    if (exists($bulkOpts->{DIELECTRIC_CONSTANT})) {
	$BULK->{DIELECTRIC_CONSTANT}{DAT} = () if (! exists($BULK->{DIELECTRIC_CONSTANT}{DAT}));
	$BULK->{DIELECTRIC_CONSTANT}{tStep}{$frameNum} = 
	    &calcDielectricConstant($ATOMS, \%{ $BULK->{DIELECTRIC_CONSTANT}{DAT} }, $BULK->{VOLUME});
	$BULK->{DIELECTRIC_CONSTANT}{tStep}{$frameNum}{u} =
	    &calcDipole($ATOMS, $MOLS, \%{ $BULK->{DIELECTRIC_CONSTANT}{DAT} }, 1);
    }

    if (exists($bulkOpts->{MOL_CORRELATION})) {
	$BULK->{MOL_CORRELATION}{DAT}{DIPOLE} = () if (! exists($BULK->{MOL_CORRELATION}{DAT}));
	$BULK->{MOL_CORRELATION}{tStep}{$frameNum} = 
	    &calcMolOrientCorrelation($ATOMS, $MOLS, \%{ $BULK->{MOL_CORRELATION}{DAT} }, $ORIENT);
    }
}

sub calcTExpand {
    my ($STATS) = @_;
    my ($count, $factor, $thermalExpand, $results);

    $factor = 10E6/5.77668; # 1/kb; kb in kcal/mol
    $count = $STATS->{counter};

    $results->{vh} = $STATS->{VH}/$count;
    $results->{temp_2} = $STATS->{T_2}/$count;
    $results->{vol} = $STATS->{V}/$count;
    $results->{h} = $STATS->{H}/$count;

    $thermalExpand = ($results->{vh} - ($results->{vol}*$results->{h}))/($results->{temp_2} * $results->{vol});
    $results->{T} = $factor * $thermalExpand; # 10E-4 K-1
    return $results;
}
    
sub calcTCompress {
    my ($STATS) = @_;
    my ($count, $factor, $thermalCompress, $results);
    
    $factor = 10E3/1.380; # V/kb; kb in kj/mol
    $count = $STATS->{counter};

    $results->{vol_2} = $STATS->{V_2}/$count;
    $results->{vol} = $STATS->{V}/$count;
    $results->{temp} = $STATS->{T}/$count;

    $thermalCompress = ($results->{vol_2} - $results->{vol}**2)/($results->{temp} * $results->{vol});
    $results->{T} = $factor * $thermalCompress; # 10E-6 atm-1
    return $results;
}
    
sub calcHCap {
    my ($STATS) = @_;
    my ($count, $heatCapacity, $Gas_Constant, $results);

    $Gas_Constant = 1.987/1000; # kcal mol-1 K-1
    $count = $STATS->{counter};

    $results->{h_2} = $STATS->{H_2}/$count; #(kcal mol-1)**2
    $results->{h} = $STATS->{H}/$count;
    $results->{temp} = $STATS->{T}/$count;

    $heatCapacity = ($results->{h_2} - $results->{h}**2)/($Gas_Constant * $results->{temp}**2);
    $results->{T} = $heatCapacity; #kcal mol-1 K-1
    return $results;
}

sub calcHVap {
    my ( $STATS) = @_;
    my ($hVap, $counter, $results);
        
    $counter = $STATS->{counter};

    $results->{U} = $STATS->{E_int}/$counter;
    $results->{vol} = $STATS->{V}/$counter;
    $results->{press} = $STATS->{P}/$counter;
    $results->{temp} = $STATS->{T}/$counter;
    $results->{pV} = ($results->{press} * $results->{vol} * 0.6023)/(9.8692 * 1000); #kJ/mol;
    $results->{rT} = 8.314472 * $results->{temp}/1000; # kJ/mol                                            

    #$hVap = -1 * 1000 * $potEng/($pV);
    #$pV = $rt = 0 if (! $applyCorrection);
    $hVap = (-1*($results->{U}/$totMols) + $results->{rT} - $results->{pV});
    #$hVap = -1*$potEng/800;
    $results->{T} = $hVap;                                                                             
    #printf "HVap: %.3f (kJ/mol)...", $hVap;
    return $results;
}

sub getFluct {
    my ($data, $STATS) = @_;

    $STATS->{H} += $data->{poteng};
    $STATS->{H_2} += $data->{poteng}**2;
    $STATS->{E_int} += ($data->{poteng} - $data->{e_bond} - $data->{e_angle}); # kcal/mol
    #$rt += (8.314472 * $data->{$i}{Temp}/1000); # kJ/mol
    $STATS->{T} += $data->{temp};
    $STATS->{T_2} += $data->{temp}**2;
    $STATS->{P} += $data->{press};
    $STATS->{V} += $data->{volume};
    $STATS->{V_2} += $data->{volume}**2;
    $STATS->{VH} += $data->{volume}*$data->{poteng};
    #$STATS->{VH_2} += ($data->{volume}*$data->{poteng})**2;
    $STATS->{counter}++;
}

sub calcMolOrientCorrelation {
    my ($atoms, $mols, $data, $orientOpts) = @_;
    my ($i, $count, $corr, $j, $curr, $start, $type, $a1, $a2, $dipole);
    my (@tmp, $offset, $molMap, %MOLORIENT, $v1, $v2, $tmpMol, $molData);

    $molMap->{0} = 0;
    $molMap->{"-1"} = -1;
    
    for $type (keys %{ $orientOpts }) {
	$count = 0;
	if ($type eq "T") {
	    for $i (keys %{ $mols }) {
		$tmpMol->{$i} = $mols->{$i};
		if (! exists($data->{DIPOLE}{$i})) {
		    $data->{DIPOLE}{$i} = calcDipole($atoms, $tmpMol, undef, 0);
		    $corr++;
		} else {
		    $dipole->{$i} = calcDipole($atoms, $tmpMol, undef, 0);
		    for $j ("X", "Y", "Z") {
			$start->{"${j}COORD"} = $data->{DIPOLE}{$i}{$j};
			$curr->{"${j}COORD"} = $dipole->{$i}{$j};
		    }
		    $corr += DotProduct($curr, $start)/DotProduct($start, $start);
		}
		$count++;
		$tmpMol = ();
	    }
	} elsif ($orientOpts->{$type} =~ /(\d+)_(\d+)_x_(\d+)_(\d+)/) { #cross product
	    for $i (keys %{ $mols }) {
		@tmp = sort numerically keys %{ $mols->{$i} };
		for $j (0 .. $#tmp) {
		    $offset = $j + 1;
		    $molMap->{$offset} = $tmp[$j];
		}
		$a1 = $atoms->{ $molMap->{$1} };
		$a2 = $atoms->{ $molMap->{$2} };
		$v1 = getVecs($a1, $a2);

		$a1 = $atoms->{ $molMap->{$3} };
		$a2 = $atoms->{ $molMap->{$4} };
		$v2 = getVecs($a1, $a2);
		
		if (! exists($data->{$type}{$i})) { # if this is the first time
		    $data->{$type}{$i} = CrossProduct($v1, $v2);
		    $corr++;
		} else {
		    $start = $data->{$type}{$i};
		    $curr = CrossProduct($v1, $v2);
		    $curr = CrossProduct($v1, $v2);	    
		    $corr += DotProduct($curr, $start)/DotProduct($start, $start);
		}
		$count++;				
	    }
	} elsif ($orientOpts->{$type} =~ /(\d+)_(\d+)/) { #regular vector
	    for $i (keys %{ $mols }) {
		@tmp = sort numerically keys %{ $mols->{$i} };
		for $j (0 .. $#tmp) {
		    $offset = $j + 1;
		    $molMap->{$offset} = $tmp[$j];
		}
		$a1 = $atoms->{ $molMap->{$1} };
		$a2 = $atoms->{ $molMap->{$2} };

		if (! exists($data->{$type}{$i})) { # if this is the first time
		    $data->{$type}{$i} = getVecs($a1, $a2);
		    $corr++;
		} else {
		    $start = $data->{$type}{$i};
		    $curr = getVecs($a1, $a2);
		    $corr += DotProduct($curr, $start)/DotProduct($start, $start);
		}
		$count++;
	    }		
	}
	
	$corr /= $count;
	$MOLORIENT{$type} = $corr;
    }

    return \%MOLORIENT;
}

sub calcDipole {
    my ($atoms, $mols, $data, $doMoments) = @_;
    my ($i, %DIPOLE, $dMoment, %MOL_DIPOLE, $j, @tmp, $k);

    @tmp = ("X", "Y", "Z");

    for $j (@tmp) {
	$DIPOLE{$j} = 0;
	for $i (keys %{ $mols }) {
	    for $k (keys %{ $mols->{$i} }) {
		$DIPOLE{$j} += ($atoms->{$k}{CHARGE} * $atoms->{$k}{"${j}COORD"});
		$MOL_DIPOLE{$i}{$j} += ($atoms->{$k}{CHARGE} * $atoms->{$k}{"${j}COORD"});
	    }
	}
    }

    if ($doMoments) {
	for $i (keys %MOL_DIPOLE) {
	    $dMoment += sqrt($MOL_DIPOLE{$i}{X}**2 +  $MOL_DIPOLE{$i}{Y}**2 + $MOL_DIPOLE{$i}{Z}**2);
	}
	$dMoment /= (scalar keys %MOL_DIPOLE);
	$data->{DIPOLE_MOMENT} += $dMoment;
	$dMoment = ($data->{DIPOLE_MOMENT} * 4.80320418270556774711/$data->{Count});
	return $dMoment;
    } else {
	return \%DIPOLE;
    }
}

sub calcDielectricConstant {
    my ($atoms, $data, $volData) = @_;
    my ($i, $j, @tmp, %DIELECTRIC, $avgVol, $factor, $d, $d_2);

    @tmp = ("X", "Y", "Z");
    if (! defined($data)) {
	$data->{Count} = 1; 
	$data->{X_2} = $data->{X} = $data->{Y_2} = 
	    $data->{Y} = $data->{Z_2} = $data->{Z} = 0;
    } else {
	$data->{Count}++;
    }
    for $i (values %{ $volData->{tStep} }) {
	$avgVol += $i;
    }
    $avgVol /= (scalar keys %{ $volData->{tStep} });
    $factor = 2931.941799957/$avgVol;

    for $j (@tmp) {
	$d = $d_2 = 0;
	for $i (keys %{ $atoms }) {
	    $d += ($atoms->{$i}{CHARGE} * $atoms->{$i}{"${j}COORD"});
	    $d_2 += ($atoms->{$i}{CHARGE} * $atoms->{$i}{"${j}COORD"})**2;
	}
	$d_2 = $d**2;
	$data->{$j} += $d;
	$data->{"${j}_2"} += $d_2;
    }

    for $j (@tmp) {
	$DIELECTRIC{$j} = $factor * ($data->{"${j}_2"}/$data->{Count} - ($data->{$j}/$data->{Count})**2);
    }
    
    $DIELECTRIC{T} = $factor * ($data->{"X_2"} + $data->{"Y_2"} + $data->{"Z_2"})/$data->{Count} -
	(($data->{X}/$data->{Count})**2 + ($data->{Z}/$data->{Count})**2 + ($data->{Z}/$data->{Count})**2);
											  
    return \%DIELECTRIC;
}

sub calcDiffusionConstant {
    my ($atoms, $mols, $data, $delT) = @_;
    my ($i, $j, $d, $dist, $MOL, @tmp, $CENTER, $factor);

    $factor = 10/(6*$totMols);
    $d = ();
    @tmp = ("X", "Y", "Z");
    for $i (keys %{ $mols }) {
	$MOL = getAtoms($atoms, $mols->{$i});
	$CENTER = CoM($MOL);
	%{ $data->{COM}{$i} } = %{ $CENTER } if (! exists($data->{COM}{$i}));
	for $j (@tmp) {
	    $dist = $data->{COM}{$i}{"${j}COORD"} - $CENTER->{"${j}COORD"};
	    $d->{lc($j)} += $dist**2;
	}
    }
    
    $d->{r} = $d->{x} + $d->{y} + $d->{z};
    $d->{a} = sqrt($d->{r});
    for $i (keys %{ $d }) {
	if ($i !~ /r|a/) {	    
	    $d->{$i} *= 3;
	}
	$d->{$i} *= $factor;
    }
    $data->{r} = $d->{r} if (! exists($data->{r}));
    $d->{T}= 0;
    $d->{T} = ($d->{r} - $data->{r})/$delT if ($delT > 0);

    return $d;
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
	$ORIENT = getOrientOpt($BGF, $MOLS, $orientOpt);
	$OPTS->{MOL_CORRELATION} = 1 if ($ORIENT);
    }
}

sub getOrientOpt {
    my ($atoms, $mols, $orientOpts) = @_;
    my (%OPTS, $i, $currOpt, $molMap, $offset, $opt1, $opt2, %OPTLIST, $optName, @tmp, @tmp2);

    @tmp = keys %{ $mols };
    @tmp2 = sort numerically keys %{ $mols->{ $tmp[0] } };

    for $i (@tmp2) {
	$offset = $i - $tmp2[0] + 1;
	$molMap->{$offset} = $i;
    }

    $orientOpts = lc($orientOpts);

    while ($orientOpts =~ /(\S+)/g) {
	$currOpt = $1;
	if ($currOpt =~ /dm/) {
	    $OPTS{T} = "dM";
	    $OPTLIST{"dm"} = "0_-1";
	} elsif ($currOpt =~ /(\w+)_x_(\w+)\:(\w+)/) {
	    ($opt1, $opt2, $optName) = ($1, $2, $3);
	    next if (lc($opt1) eq lc($opt2));
	    for ($opt1, $opt2) {
		if ($_ =~ /(\d+)_(\d+)/) {
		    next if ($1 == $2);
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
	    next if ($1 == $2);
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
    my ($a1, $a2) = @_;
    my (%vec, $i);

    for $i ("XCOORD", "YCOORD", "ZCOORD") {
	$vec{$i} = $a1->{$i} - $a2->{$i};
    }

    return \%vec;
}

sub init {
    my (%OPTS, @tmp, $i, $j, $solvTmp, $count, $tSel);
    my ($atomSel, $FF, $ATOMSELECT, $ensemble);

    getopt('btalrswfdoev',\%OPTS);
    
    for ("b", "t", "f") {
	die &showUsage . "" if (! defined($OPTS{$_}));
    }
    print "Initializing...";
    ($atomSel,$trjFile,$bgfFile,$trjType,$tSel,$saveFile,$FF, $tStep, $orientOpt, $ensemble, $logFile, $isVariable) =
        ($OPTS{a},$OPTS{t},$OPTS{b},$OPTS{r},$OPTS{s},$OPTS{w},$OPTS{f}, $OPTS{d}, $OPTS{o}, $OPTS{e}, $OPTS{l}, $OPTS{v});
    
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
    if ($ensemble =~ /(npt|nvt|nve)/) {
	$ensemble = $1;
    } else {
	$ensemble = "npt";
    }
    if ($ensemble eq "nve") {
	$OPTS->{DIFFUSION_CONSTANT} = 1;
    } else {
	$OPTS->{DENSITY} = 1;
	$OPTS->{DIELECTRIC_CONSTANT} = 1;
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

sub setOpts {

    $BULK->{DIFFUSION_CONSTANT}{UNITS} = "x 10-9 m2 s-1";
    $BULK->{DIFFUSION_CONSTANT}{EXP_A} = 2.23;
    $BULK->{DIFFUSION_CONSTANT}{EXP_D} = 0.1;

    $BULK->{DENSITY}{UNITS} = "g/cm3";
    $BULK->{DENSITY}{EXP_A} = 0.99953;
    $BULK->{DENSITY}{EXP_D} = 0.0001;

    $BULK->{DIELECTRIC_CONSTANT}{UNITS} = "";
    $BULK->{DIELECTRIC_CONSTANT}{EXP_A} = 78.46;
    $BULK->{DIELECTRIC_CONSTANT}{EXP_D} = 0.04;

    $BULK->{HEAT_VAPORIZATION}{UNITS} = "kcal/mol";
    $BULK->{HEAT_VAPORIZATION}{EXP_A} = 10.5176;
    $BULK->{HEAT_VAPORIZATION}{EXP_D} = 0.03;

    $BULK->{HEAT_CAPACITY}{UNITS} = "cal mol-1 K-1";
    $BULK->{HEAT_CAPACITY}{EXP_A} = 18.004;
    $BULK->{HEAT_CAPACITY}{EXP_D} = 0.0006;

    $BULK->{THERMAL_COMPRESSIBILITY}{UNITS} = "x 10E-6 atm-1";
    $BULK->{THERMAL_COMPRESSIBILITY}{EXP_A} = 45.86;
    $BULK->{THERMAL_COMPRESSIBILITY}{EXP_D} = 0.2;

    $BULK->{THERMAL_EXPANSION}{UNITS} = "x 10E-4 K-1";
    $BULK->{THERMAL_EXPANSION}{EXP_A} = 2.558;
    $BULK->{THERMAL_EXPANSION}{EXP_D} = 0.2;
    
    $BULK->{MOL_CORRELATION}{UNITS} = "";
    $BULK->{MOL_CORRELATION}{EXP_A} = 0;
    $BULK->{MOL_CORRELATION}{EXP_D} = 0;
    
}
sub showUsage {
    my ($usage) = "usage: $0 -b bgf file -f force field -t trajectory file" . 
        "\nOptional parameters:" .
	"\n\t-a atoms (see selection info) -d trajectory dump timestep (ps)" .	
	"\n\t-e ensemble (npt(default)|nvt|nve) -r traj type (amber(default)|lammps -w save name" .
	"\n\t-s trajectory selection (see selection info) -l lammps log file" . 
	"\n\t-o orientation options (none default): For the molecule, expected 1_2:uM for correlation of 1->2 vector" .
	"\n\t\t labeled uM or {name1}_x_uM:p for the vector cross product of predefined name1 (e.g. 1->2) and uM called P" . 
	"\n\t\t Use dm for dipole moment. Enclose multiple in quotes\n" .
	"\n\t-v is variable: Set to 1 if trajectory contains varying # atoms. Default 0\n" .
	&ShowSelectionInfo;
    return $usage;
}

#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts/Packages";
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Chart::Graph::Gnuplot qw(gnuplot);
use Packages::General;
use Packages::LAMMPS;
use Packages::CERIUS2;
use Packages::FileFormats;
use File::Basename;
use Packages::Superimpose;
use Packages::Math::MatrixReal;

sub getEnergy;
sub plotEnergy;
sub setZero;
sub initialize;
sub numerically;
sub addFields;
sub calcCOM;
sub reImageAtoms;
sub comOpt;
sub engOpt;
sub bgfOpt;
sub getAvgCoords;
sub getMass;
sub unscaleAtoms;
sub alignMols;
sub rmsOpt;
sub getRMS;
sub densityOpt;

die "usage: $0 trajectory|directory searchParm atomSelection starting bgf [cerius2 ff]\n"
    if (! @ARGV or $#ARGV < 3);

my ($searchPath, $dataType, $selection, $bgfFile, $cerius2FF) = @ARGV;
my (@vFiles, $DATA, $saveName, $file, $SELECT, $avgBGF);
my ($startBGF, $bonds, $FFDATA, $scaled, $tmp, $atomsInCoM);

$|++;
initialize;

print "Parsing LAMMPS trajectories...";
for $file (@vFiles) {
    $scaled = 1;
    #print "File: " . $file->{"FILE"} . "\n";
    ($DATA->{$file->{"INDEX"}}, $tmp) = ParseLammpsTrajectoryFile($file->{"FILE"}, $dataType, $SELECT);
    $scaled = $tmp if ($scaled);
    if ($scaled) {
	unscaleAtoms($DATA->{$file->{"INDEX"}});
    }
    if ($startBGF) {
	addFields($DATA->{$file->{"INDEX"}}, $startBGF);
    }
    $DATA->{$file->{"INDEX"}}{"INDEX"} = $file->{"INDEX"};
}
print "Done\n";

if ($dataType =~ /com/) {
    print "Calculating COM...";
    comOpt($DATA, $startBGF, $FFDATA, $atomsInCoM);
    print "Done\n";
}

if ($dataType =~ /eng/) {
    print "Computing Energy profile...";
    engOpt($DATA, $selection);
    print "Done\n";
}

if ($dataType =~ /bgf/) {
    print "Creating BGF files for each timestep...";
    bgfOpt($DATA, $startBGF, $bonds, $saveName);
    print "Done\n";
}

if ($dataType =~ /avg/) {
    $avgBGF = avgOpt($DATA, $startBGF, $bonds, $saveName, $FFDATA);
}

if ($dataType =~ /rms/) {
    rmsOpt($DATA, $startBGF, $bonds, $saveName, $FFDATA, $avgBGF);
}

if ($dataType =~ /density/) {
    densityOpt($DATA, $startBGF, $FFDATA);
}

sub comOpt {
    my ($fData, $bgfData, $PARMS, $comAtoms) = @_;
    my ($i, $tStep, $HEADER);

    @{ $HEADER } = @{ $bgfData->{"HEADER"} };
    delete $bgfData->{"HEADER"};
    for $i (keys %{ $fData }) {
	for $tStep (keys %{ $fData->{$i} }) {
	    next
		if ($tStep !~ /\d+/);
	    AddMass($fData->{$i}{$tStep}{"ATOMS"}, $PARMS);
	    if (! exists($fData->{$i}{$tStep}{"ATOMS"})) {
		print "no atoms for tStep $tStep in file $i\n";
	    } else {
		calcCOM($fData->{$i}{$tStep}, $comAtoms);
	    }
	}
    }
    addHeader($bgfData, $HEADER);
}
   
sub engOpt {
    my ($fData, $tot) = @_;

    my ($ENG, $COMEng) = getEnergy($fData, $tot);
    print "Done\nPlotting and saving energy files...";
    plotEnergy($ENG, "");
    if ($dataType =~ /com/) {
	plotEnergy($COMEng, "com");
    }
}

sub bgfOpt {
    my ($fData, $bgfData, $bonds, $save) = @_;
    my ($i, $tStep);

    for $i (keys %{ $fData }) {
	for $tStep (keys %{ $fData->{$i} }) {
	    next if ($tStep !~ /\d+/);
            @{ $fData->{$i}{$tStep}{"ATOMS"}{"HEADER"} } = @{ $bgfData->{"HEADER"} };
	    createBGF(\%{ $fData->{$i}{$tStep}{"ATOMS"} }, $bonds, $save . $tStep . ".bgf");
	}
    }
}
    
sub avgOpt {
    my ($fData, $bgfData, $bonds, $save, $PARMS) = @_;
    my ($COORDS, $HEADER, $saveBGF, $index, $MASSES); 

    $MASSES = GetMass($PARMS);
    $save .= "_avg.bgf";
    $COORDS = alignMols(1, $fData, $bgfData, $MASSES);
    print "Creating average BGF file $save...";
    $saveBGF = getAvgCoords($COORDS);
    $COORDS = ();
    $HEADER = createHeaders(my $BOX, $save);
    addHeader($saveBGF, $HEADER);
    createBGF($saveBGF, $bonds, $save);
    delete $saveBGF->{"HEADER"};
    print "Done\n";
    return $saveBGF;
}

sub rmsOpt {
    my ($fData, $inputBGF, $bonds, $save, $PARMS, $avgBGF) = @_;
    my ($i, $COORDS, $fileOut, @tmp, $timeStart); 
    my ($tStep, $startBGF, $rmsd, $MASSES);

    delete $inputBGF->{"HEADER"};
    $MASSES = GetMass($PARMS);
    print "Calculating RMS using ${save}.bgf as reference...";
    $COORDS = alignMols(0, $fData, $inputBGF, $MASSES, $inputBGF);
    $rmsd->{"starting"} = getRMS($COORDS, $inputBGF, $MASSES);
    printf "avg rms: %5.3f +/- %5.3f ...Done\n", $rmsd->{starting}{stats}{avg}, $rmsd->{starting}{stats}{stdev};

    #@tmp = keys %{ $fData };
    #$i = 9999999999999;
    #for $tStep (keys %{ $fData->{ $tmp[0] } }) {
	#next if ($tStep !~ /\d+/);
	#$i = $tStep if ($tStep < $i);
    #}
    #$timeStart = "timestep $i";
    #$startBGF = $fData->{ $tmp[0] }{$i}{"ATOMS"};
    #print "Calculating RMS using ${timeStart} as reference...";
    #$COORDS = alignMols(0, $fData, $inputBGF, $MASSES, $startBGF);
    #$rmsd->{$timeStart} = getRMS($COORDS, $startBGF, $MASSES);
    #printf "avg rms: %5.3f +/- %5.3f ...Done\n", $rmsd->{$timeStart}{stats}{avg}, $rmsd->{$timeStart}{stats}{stdev};

    #if (defined($avgBGF)) {
	#print "Calculating RMS using ${save}_avg.bgf as reference...";
	#$COORDS = alignMols(0, $fData, $inputBGF, $MASSES, $avgBGF);
        #$rmsd->{"average"} = getRMS($COORDS, $avgBGF, $MASSES);
	#printf "avg rms: %5.3f +/- %5.3f ...Done\n", $rmsd->{average}{stats}{avg}, $rmsd->{average}{stats}{stdev};
    #}	

    $COORDS = ();
    print "Creating RMS data file ${save}_rms.dat...";
    open RMSDATA, "> ${save}_rms.dat" or die "ERROR: Cannot create RMS file ${save}_rms.dat: $!\n";
    printf RMSDATA "%-10s%20s", "#TSTEP", "STARTING";
    @tmp = ("starting");
    for $i (keys %{ $rmsd }) {
	next if ($i eq "starting");
	push @tmp, $i;
	printf RMSDATA "%20s", $i;
    }
    print RMSDATA "\n";

    delete $rmsd->{"starting"}{"stats"};
    for $tStep (sort numerically keys %{ $rmsd->{"starting"} }) {
	next if ($tStep !~ /\d+/);
	printf RMSDATA "%-10d", $tStep;
	for $i (@tmp) {
	    printf RMSDATA "%20.3f", $rmsd->{$i}{$tStep};
	}
	print RMSDATA "\n";
    }

    close RMSDATA or die "ERROR: Cannot close RMS data file ${save}_rms.dat: $!\n";
    print "Done\n";
}

sub densityOpt {
    my ($fDATA, $BGF, $PARMS) = @_;
    my ($i, $atomC, $totMass, $tStep, $vol, $dim, $density);
    my (%DENS, $HEADER);

    @{ $HEADER } = @{ $BGF->{"HEADER"} };
    delete $BGF->{"HEADER"};
    AddMass($BGF, $PARMS);
    print "Creating density file density.dat...";
    $totMass = 0;
    for $atomC (keys %{ $BGF }) {
	$totMass += $BGF->{$atomC}{"MASS"};
    }
    
    for $i (keys %{ $fDATA }) {
	for $tStep (keys %{ $fDATA->{$i} }) {
	    next if ($tStep !~ /\d+/);
	    $vol = 1;
	    for $dim (@{ $fDATA->{$i}{$tStep}{"BOX"}}) {
		$vol *= $dim->{"hi"} - $dim->{"lo"};
	    }
	    $density = $totMass/($vol * 0.6023);
	    $DENS{$tStep} = $density;
	}
    }

    open DENSFILE, "> density.dat" or die "ERROR: Cannot create density file density.dat: $!\n";
    for $i (sort numerically keys %DENS) {
	printf DENSFILE "%-15d%11.5f\n", $i, $DENS{$i};
    }
    close DENSFILE or die "ERROR: Cannot close file density.dat:$!\n";
    @{ $BGF->{"HEADER"} } = @{ $HEADER };
    print "Done\n";
}

sub getRMS {
   my ($COORDS, $refMol, $MASSES) = @_;
   my (%RMSDATA, $tStep, $rms, $atomC, $rStr); 
   my ($atom, $totMass, $counter, $mass, $dist, $dim);
   
   for $tStep (@{ $COORDS }) {
	$rms = 0;
	$totMass = 0;
	$counter = 0;
	for $atomC (keys %{ $tStep->{"DATA"} }) {
	    $atom = $tStep->{"DATA"}{$atomC};
	    $mass = $MASSES->{ $atom->{"FFTYPE"} };
	    $mass = 1;
	    $dist = GetBondLength($refMol->{$atomC}, $atom);
	    $dist = $dist * $dist;
	    $rms += $mass * $dist;
	    $totMass += $mass;
	    $counter++;
	}
	$RMSDATA{ $tStep->{"TSTEP"} } = sqrt($rms/$totMass);
	$rStr .= "$RMSDATA{ $tStep->{TSTEP} } ";
   }
   chop $rStr;
   ($RMSDATA{stats}{avg}, $RMSDATA{stats}{stdev}, $RMSDATA{stats}{total}) = STDev($rStr);
   return \%RMSDATA;
}

sub alignMols {
    my ($doPrint, $fData, $bgfData, $MASSES, $refBGF) = @_;
    my ($tStep, @COORDS, $i, $index, $rec, $matchBGF, $dim);

    $matchBGF = ();
    @COORDS = ();
    $rec = ();

    for $i (keys %{ $fData }) {
        for $tStep (keys %{ $fData->{$i} }) {
	    next if ($tStep !~ /\d+/);
	    #for $index (keys %{ $fData->{$i}{$tStep}{"ATOMS"} }) {
		#%{ $currATOMS->{$index} } = %{ $bgfData->{$index} };
		#for $dim ("XCOORD", "YCOORD", "ZCOORD") {
		    #$currATOMS->{$index}{$dim} = $fData->{$i}{$tStep}{"ATOMS"}{$index}{$dim};
		#}
	    #}
	    
	    %{ $matchBGF } = %{ $fData->{$i}{$tStep}{"ATOMS"} };
            if (! defined($refBGF)) {
		print "Using timestep $tStep as reference\n" if ($doPrint);
                $refBGF = $matchBGF;
            } else {
                print "Superimposing timestep $tStep onto ref..\r" if ($doPrint);
                SuperimposeAtoms($refBGF, $matchBGF, $MASSES);
            }
	    $rec = (
		    {
			"DATA"  => $matchBGF,
			"TSTEP" => $tStep,
		    }
		    );
	    push @COORDS, $rec;
	    $rec = ();
	    $matchBGF = ();
        }
    }
    return \@COORDS;
}

sub getAvgCoords {
    my ($FILES) = $_[0];
    my ($BGF, $atomC, $total, $i, $dim);

    %{ $BGF } = %{ $FILES->[0]{"DATA"} };

    $total = $#{ $FILES };
    print "total: " . ($total + 1) . " snapshots...";
    for $atomC (keys %{ $FILES->[0]{"DATA"} }) {
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    for $i (1 .. $total) {
		$BGF->{$atomC}{$dim} += $FILES->[$i]{"DATA"}{$atomC}{$dim};
	    }
	    $BGF->{$atomC}{$dim} /= ($total + 1);
	}
    }

    return $BGF;
}
	    

sub initialize {
    my ($accepted) = "bgf, avg, com, eng, rms, density";
    my (@tmp, $index, $HEADERS);

    $index = 0;

    if (! -d $searchPath) {
	FileTester($searchPath);
	push @vFiles, (
		       {
			   "FILE"  => $searchPath,
			   "INDEX" => 1,
		       }
		       );
	print "Using file $searchPath as input\n";
    } else {
	print "Searching for LAMMPS dump files...";
	opendir TRAJFILES, $searchPath or die "ERROR: Cannot access directory $searchPath: $!\n";
	@tmp = grep { /\.dump$/ && -f} map {"$searchPath/$_"} readdir TRAJFILES;
	closedir TRAJFILES or die "ERROR: Cannot close directory $searchPath: $!\n";
	for (@tmp) {
	    if ($_ =~ /(\d+)\.\w+$/) {
		$index = $1;
	    } else {
		$index++;
	    }
	    push @vFiles, (
			   {
			       "FILE"  => $_,
			       "INDEX" => $index,
			   }
			   );
	}
	die "ERROR: No valid files found in directory $searchPath\n"
	    if (! @vFiles);
	print "found " . ($#vFiles + 1) . " valid files\n";
				 
    }

    if (! $saveName) {
	$saveName = basename($vFiles[0]);
    }

    $selection = Trim($selection);
    if ($selection eq "*") {
	$SELECT = ();
	print "USING EVERY ATOM\n";
	undef($selection);
    } else {
	if ($selection =~ /\s/) {
	    @tmp = split /\s+/, $selection;
	} else {
	    @tmp = ($selection);
	}
	$selection = GetSelections(\@tmp);
    }

    FileTester($bgfFile);
    print "Obtaining atom information from $bgfFile...";
    ($startBGF, $bonds, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
    if (defined($selection)) {
        ($startBGF, $bonds, $SELECT) = GetBGFAtoms($selection, $startBGF, $bonds);
        die "ERROR: No atoms in BGF file corresponds to atom selection\n" if (! keys %{ $startBGF });
    }
    addHeader($startBGF, $HEADERS);
    $saveName = $bgfFile;
    print "Done\n";
    
    if ($dataType =~ /(com|avg|rms|bgf|density)/) {
	die "ERROR: A cerius2 force field is required when calculating $1\n"
	    if (! defined($cerius2FF));
	FileTester($cerius2FF);
	
        print "Parsing Cerius2 forcefield $cerius2FF...";
	$FFDATA = parseCerius2FF($cerius2FF);
	print "Done\n";
    }

    $saveName =~ s/\.\w+$//;
    $saveName = basename($saveName);
    $dataType = Trim($dataType);
    die "Invalid option for searchParm. Got $dataType\nExpected options are:\n$accepted\n"
	if ($dataType !~ /bgf|avg|eng|com|rms|density/i);
    if ($dataType =~ /com\:(\d+)/) {
	$atomsInCoM = $1;
	print "First $1 atoms selected as group 1... remainder is group 2\n";
    } elsif ($dataType =~ /com/) {
	print "Calculating CoM for entire selection\n";
	$atomsInCoM = 0;
    }

}

sub getEnergy {
    my ($data, $atomTot) = @_;
    my ($i, $tStep, $engTot, %ENG, $atomEng, $atom); 
    my ($stats, $fileTot, %COMEng, $comDist, $comList);
    
    for $i (keys %{ $data }) {
	$fileTot = "";
	$comList = "";
	for $tStep (keys %{ $data->{$i} }) {
	    next if ($tStep eq "INDEX");
	    $engTot = 0;
	    for $atom (keys %{ $data->{$i}{$tStep}{"ATOMS"} }) {
		#next if ($atom > $atomTot);
		$atomEng = $data->{$i}{$tStep}{"ATOMS"}{$atom}{"ENERGY"};
		$ENG{$data->{$i}{"INDEX"}}{"TIME"}{$tStep}{$atom} = $atomEng;
		$engTot += $atomEng;
	    }
	    $ENG{$data->{$i}{"INDEX"}}{"TIME"}{$tStep}{"TOTAL"} = $engTot;
	    $fileTot .= "$engTot ";
	    if (exists($data->{$i}{$tStep}{"COM"})) {
		$comDist = GetBondLength($data->{$i}{$tStep}{"COM"}{"1"},$data->{$i}{$tStep}{"COM"}{"2"});
		$COMEng{$data->{$i}{"INDEX"}}{"TIME"}{$comDist}{"TOTAL"} = $engTot;
		$ENG{$data->{$i}{"INDEX"}}{"COM"}{$tStep} = $comDist;
		$comList .= "$comDist ";
	    }
	}
	$stats = \%{ $ENG{$data->{$i}{"INDEX"}}{"STATISTICS"} };
	chop $fileTot;
	chop $comList;
	($stats->{"AVG"}, $stats->{"STDEV"}, $stats->{"TOTAL"}) = STDev($fileTot);
	if (exists($COMEng{$data->{$i}{"INDEX"}})) {
	    $stats = \%{ $COMEng{$data->{$i}{"INDEX"}}{"STATISTICS"} };
	    ($stats->{"AVG"}, $stats->{"STDEV"}, $stats->{"TOTAL"}) = STDev($fileTot);
	    $stats = \%{ $COMEng{$data->{$i}{"INDEX"}}{"COMSTATS"} };
	    ($stats->{"AVG"}, $stats->{"STDEV"}, $stats->{"TOTAL"}) = STDev($comList);	    
	}
    }
    
    return (\%ENG,\%COMEng);
}

sub plotEnergy {
    my ($eData, $save) = @_;
    my (@XDATA, @YDATA, @YHI, @YLO, $stdev, $eng);
    my ($i, $j, $fIndex, @fXDATA, @fYDATA, $tEngStr, $tStep); 
    my($avgEngStr, %TStepEng, %COMperTime, $comStr, $avgCom);
    
    for $fIndex (sort numerically keys %{ $eData }) {
	for $tStep (sort numerically keys %{ $eData->{$fIndex}{"TIME"} }) {
	    $eng  = $eData->{$fIndex}{"TIME"}{$tStep}{"TOTAL"};
	    push @XDATA, $tStep;
	    push @YDATA, $eng;
	    $TStepEng{$tStep} = $eng;
	    if (exists($eData->{$fIndex}{"COM"}{$tStep})) {
		$COMperTime{$tStep} = $eData->{$fIndex}{"COM"}{$tStep};
	    }
	}
	$eng = \%{ $eData->{$fIndex}{"STATISTICS"} };
	push @fXDATA, $fIndex;
	push @fYDATA, $eng->{"AVG"};
	push @YLO, ($eng->{"AVG"} - $eng->{"STDEV"});
	push @YHI, ($eng->{"AVG"} + $eng->{"STDEV"});
	if ($save ne "com") {
	    $avgEngStr .= sprintf("%11d%11.3f%8.3f\n",$fIndex, $eng->{"AVG"}, $eng->{"STDEV"});
	} else {
	    $avgCom = $eData->{$fIndex}{"COMSTATS"}{"AVG"};
	    $avgEngStr .= sprintf("%11.3f%11.3f%8.3f\n",$avgCom, $eng->{"AVG"}, $eng->{"STDEV"});
	}
    }
    
    for $tStep (sort numerically keys %TStepEng) {
	$tEngStr .= sprintf("%11.3f%11.4f\n", $tStep,$TStepEng{$tStep});
	if (exists($COMperTime{$tStep})) {
	    $comStr .= sprintf("%11.3f%11.4f\n", $tStep,$COMperTime{$tStep});
	}
    }

    if ($save) {
	$save .= "VS";
    }

    open OUTDATA, "> $save"  . "avgenergy.dat" or die 
	"ERROR:Cannot create file $save" . "avgenergy.dat: $!\n";
    print OUTDATA $avgEngStr;
    close OUTDATA;

    if (! $_[1]) {
	$save = "tStep";
    } else {
	$save = $_[1];
    }
    open OUTDATA, "> $save" . "VStotenery.dat" or die 
	"ERROR: Cannot create file $save" . "_tStepVStotenergy.dat: $!\n";
    print OUTDATA $tEngStr;
    close OUTDATA;
    if ($comStr) {
	open OUTDATA, "> $save" . "VScom.dat" or die 
	    "ERROR: Cannot create  file $save" . "_tstepVScom.dat: $!\n";
	print OUTDATA $comStr;
	close OUTDATA;
    }

    gnuplot(
	    {
		"title"       => "Plot of Avg Total Energy vs File Index",
		"xrange"      => "[1:10]",
		"yrange"      => "[50:125]",
		"output file" => $save . "_avgeng.png",
		"output type" => "png",
	    },
	    [
	     {
		 "title" => "Average Energy (kcal/mol)",
		 "style" => "yerrorbars",
		 "using" => "1:2:3:4",
		 "type" => "columns",
	     },
	     [ @fXDATA ], # x
	     [ @fYDATA ], # y
	     [ @YLO ], # ylow
	     [ @YHI ], # yhigh
	     ], 
	    );
}

sub addFields {
    my ($fDATA, $BGF) = @_;
    my ($atomC, $field, $tStep, $atom);

    for $tStep (keys %{ $fDATA }) {
	for $atomC (keys %{ $fDATA->{$tStep}{"ATOMS"} }) {
	    $atom = $fDATA->{$tStep}{"ATOMS"}{$atomC};
	    die "ERROR: No information for atom # $atom found in LAMMPS trajectory file. Cannot continue\n"
		if (! exists($BGF->{$atomC}));
	    for $field (keys %{ $BGF->{$atomC} }) {
		next if ($field =~ /COORD/);
		$atom->{$field}= $BGF->{$atomC}{$field};
	    }
	}
    }
}

sub reImageAtoms {
    my ($ATOMS, $BOX) = @_;
    my ($atom, $bLen, $dim);

    $bLen->{"X"} = $BOX->[0]{"lo"} - $BOX->[0]{"hi"};
    $bLen->{"Y"} = $BOX->[1]{"lo"} - $BOX->[1]{"hi"};
    $bLen->{"Z"} = $BOX->[2]{"lo"} - $BOX->[2]{"hi"};

    for $atom (keys %{ $ATOMS }) {
	next if (exists($ATOMS->{$atom}{"fixed"}));
	for $dim ("X", "Y", "Z") {
	    if ($ATOMS->{$atom}{$dim . "INDEX"}) {
		$ATOMS->{$atom}{$dim . "COORD"} += $ATOMS->{$atom}{$dim . "INDEX"} * $bLen->{$dim};
	    }
	}
	$ATOMS->{$atom}{"fixed"} = 1;
    }
}

sub unscaleAtoms {
   my ($fDATA) = @_;
   my ($atomC, $bLen, $dim, $tStep, $atom);

   for $tStep (keys %{ $fDATA }) {
       $bLen = ();
       $bLen->{"X"} = $fDATA->{$tStep}{"BOX"}[0]{"hi"} - $fDATA->{$tStep}{"BOX"}[0]{"lo"};
       $bLen->{"Y"} = $fDATA->{$tStep}{"BOX"}[1]{"hi"} - $fDATA->{$tStep}{"BOX"}[1]{"lo"};
       $bLen->{"Z"} = $fDATA->{$tStep}{"BOX"}[2]{"hi"} - $fDATA->{$tStep}{"BOX"}[2]{"lo"};

       for $atomC (keys %{ $fDATA->{$tStep}{"ATOMS"} }) {
	   $atom = $fDATA->{$tStep}{"ATOMS"}{$atomC};
	   for $dim ("X", "Y", "Z") {
	       $atom->{$dim . "COORD"} *= $bLen->{$dim};
	   }
       }
   }
}

sub calcCOM {
    my ($SNAPSHOT, $comAtoms) = @_;
    my ($atom, @tmp, $totAtoms, $myBGF, $comCounter);
    
    @tmp = sort numerically keys %{ $SNAPSHOT->{"ATOMS"} };
    $totAtoms = 0;
    $SNAPSHOT->{"COM"}{2} = 0;
    $SNAPSHOT->{"COM"}{1} = 0;
    $comCounter = 1;

    if ($comAtoms) {
	while ($totAtoms < $comAtoms and $#tmp > -1) {
	    $atom = shift @tmp;
	    %{ $myBGF->{$atom} } = %{ $SNAPSHOT->{"ATOMS"}{$atom} };
	    $totAtoms++;
	}
	$SNAPSHOT->{"COM"}{1} = CoM($myBGF);
	$myBGF = ();
	$comCounter = 2;
    }

    for $atom (@tmp) {
	%{ $myBGF->{$atom} } = %{ $SNAPSHOT->{"ATOMS"}{$atom} };
    }

    if ($myBGF) {
	$SNAPSHOT->{"COM"}{$comCounter} = CoM($myBGF);
    }
}


sub numerically {
    ($a<=>$b);
}


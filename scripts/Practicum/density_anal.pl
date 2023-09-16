#!/usr/bin/perl -w

use strict;
use constant PI => atan2(1,1) * 4;


#     This program generates cylindrical axial velocity profiles for
#     multiple components. Cylinder axis assumed to be x-axis.
#
#     This variant generates a block average of
#     the profiles every nblock frames
#
#     N.B. If you don't want to do block averages, this
#     kind of sucks. You have to set nblock = nsample
#
#     Usage density.exe [list of file_suffix separated by blanks]
#
#     If no file_suffix list is given, default output and input
#     filenames are assumed. If only one file_suffix is given, it is added
#     to default input and output filenames. If more than one file_suffix is
#     given, these are added to the default input filenames, in order, and
#     default output filenames are used.

sub Initialize();
sub ReadParmFile();
sub AssignParms(@);
sub GetFiles();
sub CalcCylVals();
sub ReadFile(@);
sub IsValidData(@);
sub CollectStats(@);
sub SymmetrizeData();
sub CalcBins(@);
sub GetTimeSteps(@);
sub InOut(@);
sub DoBlockAverages(@);
sub WriteBins(@);
sub ValidTimeStep(@);

die "usage: $0 parmfile maxatoms [file suffixes]\nUse -1 for maxatoms too use all atoms\n"
    if (! @ARGV or $#ARGV < 1);
	
my ($parmfile, $maxatoms, @SUFF) = @ARGV;
my (%PARMS, $ANG_cnv, %sigma, $fcounter, @DATA, $ntype, @delRsq, @RsqWall);
my ($bin, $VBin, $FBin, $rsqmax, @rf, @rsqAds);

Initialize();
ReadParmFile();
my ($FILES) = GetFiles();

CalcCylVals();

$fcounter = 0;
while ($fcounter <= $#{ $FILES }) {
    @DATA = ();
    $bin = $VBin = $FBin = {};

    print "Reading data in " . $FILES->[$fcounter]{"VELOCITIES"} . "...";
    $DATA[0]{"VELS"} = ReadFile(1, $FILES->[$fcounter]->{"VELOCITIES"}, 0);
    print "Sucess\n";
    
    print "Reading data in " . $FILES->[$fcounter]{"ATOMS"} . "...";
    $DATA[0]{"ATOMS"} = ReadFile(0, $FILES->[$fcounter]{"ATOMS"}, 0);
    print "Sucess\n";
    
    print "Reading data in " . $FILES->[$fcounter]{"FORCES"} . "...";
    $DATA[0]{"FORCES"} = ReadFile(2, $FILES->[$fcounter]->{"FORCES"}, 0);
    print "Sucess\n";

    next
        if (! SymmetrizeData());
  
    print "Opening File for output...";
    InOut(1, $fcounter);
    print "Sucess\n";
    print "Collecting Statistics...";
    CollectStats(0);
    print "Sucess\n\n";
    print "Closing file...";
    InOut(0, $fcounter);
    print "Sucess\n";

    $fcounter++;
}

sub GetFiles() {
    my (@ValidFiles, $rec, $suffix, $curr_header, $filenm, $counter);
    
    $rec = (
	    {
		"ATOMS"      => "ldump_atoms",
		"VELOCITIES" => "ldump_vels",
		"FORCES"     => "ldump_forces",
	    }
	    );

    $counter = 0;
    if (@SUFF) {
	for $suffix (@SUFF) {
	    for $curr_header (keys %{ $rec }) {
		$filenm = $rec->{$curr_header} . $suffix;
		if (! -e $filenm or ! -T $filenm) {
		    print "ERROR: File not found: $filenm\n";
		    delete $ValidFiles[$counter];
		    $counter--;
		} else {
		    $ValidFiles[$counter]{$curr_header} = $filenm;
		    $ValidFiles[$counter]{"suffix"} = $suffix;
		}
	    }
	    $counter++;
	}
    } else {
	for $curr_header (keys %{ $rec }) {
	    $filenm = $rec->{$curr_header};
	    if ( -e $filenm and  -T $filenm) {
		$ValidFiles[$counter]{$curr_header} = $filenm;
		$ValidFiles[$counter]{"suffix"} = "";
	    } else {
		die "ERROR: File not found when accessing $filenm. Aborting\n";
	    }
	}
    }

    if ($#ValidFiles == -1) {
	die "ERROR: No valid files found in directrory\n";
    }

    for $counter (@ValidFiles) {
	$counter->{"OUTFILE"}{"denfile"} = "density" . $counter->{"suffix"} . ".dat";
	$counter->{"OUTFILE"}{"adsfile"} = "adsfrac" . $counter->{"suffix"} . ".dat";
	$counter->{"OUTFILE"}{"velfile"} = "vel" . $counter->{"suffix"} . ".dat";
	$counter->{"OUTFILE"}{"vsqfile"} = "vsq" . $counter->{"suffix"} . ".dat";
	$counter->{"OUTFILE"}{"vavfile"} = "vav" . $counter->{"suffix"} . ".dat";
	$counter->{"OUTFILE"}{"frcfile"} = "frc" . $counter->{"suffix"} . ".dat"
	    if ($PARMS{"ANALYSIS"}{"LForce"});
    }


    return \@ValidFiles;
}

sub Initialize() {
    $ANG_cnv = 1/(1.6022E-19 * 1.0E20);    # convertion to Angstoms

    die "Cannot locate $parmfile: $!\n"
	if (! -e $parmfile);
    die "Error: $parmfile is not readable\n"
	if (! -r $parmfile or ! -T $parmfile);

    %PARMS = (
	      "SOLVENT"      => {
		                 "itypef"       => -9999.99,
				 "qf"           => -9999.99,
				 "rf"           => -9999.99,
				 "mf"           => -9999.99,
				 },
	      "IONS"          => {
				  "itypefion"    => -9999.99,
				  "qfion"        => -9999.99,
				  "rfion"        => -9999.99,
				  "mfion"        => -9999.99,
				  },
	      "SOLUTE"       => {
				 "incS"        => 0,
				 "itypeSol"    => -9999.99,
				 "qSol"        => -9999.99,
				 "rSol"        => -9999.99,
				 "buffer"      => -9999.99,
				 "mSol"        => -9999.99,
				 },
	      "CYLINDER"     => {
		                 "rhof_targ"   => -9999.99,
				 "sigma"       => -9999.99,
				 "rc"          => -9999.99,
				 "xc"          => -9999.99,
				 "qw"          => -9999.99,
			     },
	      "ANALYSIS"    => {
		                 "LForce"      => 0,
                                 "maxstep"     => -9999.99,
                                 "maxbin"      => -9999.99,
                                 "nblock"      => -9999.99,
				 "nskip"       => -9999.99,
	                     },
	      );
}

sub ReadParmFile() {
    my ($header, $pkeys, $invalid);

    open INFILE, $parmfile or die "Cannot read $parmfile: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^\s*(\w+)\s*=\s*(\-?\d+\.?\d*e*\d*)/) {
	    if (! AssignParms($1,$2)) {
		print "Invalid key: $1\n";
	    }
	}
    }
    close INFILE;
    
    $invalid = 0;
    for $header (keys %PARMS) {
	if ($header ne "SOLUTE" or ($PARMS{"SOLUTE"}{"incS"} == 1)) {
	    for $pkeys (keys %{ $PARMS{$header} }) {
		if ($PARMS{$header}{$pkeys} <= -9999) {
		    die "Error in parameter file. $pkeys is either missing or invalid\n";
		    $invalid = 1;
		    last;
		}
	    }
	}
	if ($invalid) {
	    last;
	}
    }
}

sub AssignParms(@) {
    my ($keynm, $keyval) = @_;
    my ($header, $pkeys, $isvalid);

    $isvalid = 0;
    for $header (keys %PARMS) {
	for $pkeys (keys %{ $PARMS{$header} }) {
	    if ($pkeys eq $keynm) {
		$PARMS{$header}{$pkeys} = $keyval;
		$isvalid = 1;
		last;
	    }
	}
	if ($isvalid) {
	    last;
	}
    }

    return $isvalid;
}

sub CalcCylVals() {
    my ($counter, $ibin);

    $ntype = 1;
#   SOLVENT
    $rf[0] = $PARMS{"CYLINDER"}{"rc"} - $PARMS{"SOLVENT"}{"rf"};
    $rsqAds[0] = ($rf[0] - $PARMS{"SOLVENT"}{"rf"}) **2;

#   IONS
    $rf[1] = $PARMS{"CYLINDER"}{"rc"} - $PARMS{"IONS"}{"rfion"};
    $rsqAds[1] = ($rf[1] - $PARMS{"IONS"}{"rfion"}) **2;

#   SOLUTE
    if ($PARMS{"SOLUTE"}{"incS"}) {
	$rf[2] = $PARMS{"CYLINDER"}{"rc"} - $PARMS{"SOLUTE"}{"rSol"};
	$rsqAds[2] = ($rf[2] - $PARMS{"SOLUTE"}{"rSol"}) **2;
	$ntype = 2;
    }	
    
    for $counter (0 .. $ntype) {
    	$rsqmax = $rf[$counter] ** 2;
	$delRsq[$counter] = $rsqmax / ($PARMS{"ANALYSIS"}{"maxbin"} - 1);
	$RsqWall[$counter] = $rf[$counter]**2 - 0.05*$delRsq[$counter];
    }

}

sub ReadFile(@) {
    my ($filetype, $filename, $fcounter) = @_;
    my ($isvalid, $start_record, $has_data, $rec, $counter, $get_box);
    my ($ATOMDATA, @BBOX, $searchstring, $OUTDATA, $timestep);

    if ($filetype == 0) {
	$searchstring = " ITEM: ATOMS";
    } elsif ($filetype == 1) {
	$searchstring = " ITEM: VELOCITIES";
    } else {
	$searchstring = " ITEM: FORCES";
    }	
    $isvalid = $start_record = $counter = $has_data = $get_box = 0;
    $timestep = -1;
    open INFILE, $filename or die "Cannot access file $filename: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^ ITEM: TIMESTEP/) {
	    if ($has_data > 0 and $timestep > -1) {
		push @{ $OUTDATA }, $ATOMDATA;
		$OUTDATA->[$counter]{"TIMESTEP"} = $timestep;
		$counter++;
	    }
	    $ATOMDATA = {};
	    @BBOX = ();
	    $has_data = $start_record = 0;
	    $timestep = -1;
	} elsif ($timestep == -1 and $_ =~ /^\s+(\d+)\s*$/) {
	    if (ValidTimeStep($filetype, $1, $fcounter)) {
		$timestep = $1;
	    }
	} elsif ($_ =~ /^ ITEM: BOX BOUNDS/ and $timestep > -1) {
	    $get_box = 1;
	} elsif ($_ =~ /$searchstring/ and $timestep > -1) {
	    $get_box = 0;
	    $start_record = 1;
	} elsif ($get_box and $timestep > -1) {
	    if ($_ =~ /^\s*(\-?\d+\.\d+E?\+?\d*)\s+(\-?\d+\.\d+E?\+?\d*)/) {
		push @BBOX, $1;
		push @BBOX, $2;
	    }
	} elsif ($start_record and $timestep > -1) {
	    if ($_ =~ /^\s+(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
		if (IsValidData($1, $2)) {
		    $rec = (
			    {
				"TYPEID" => $2,
				"XCOORD" => $3,
				"YCOORD" => $4,
				"ZCOORD" => $5,
				"NUM"    => $1,
			    }
			    );
		    $ATOMDATA->{$1} = $rec;
		    $has_data = $isvalid = 1;
		}else {
		    print "Error: Either Type $2 or Atom $1 is not valid. Ignoring.\n";
		}
	    }
	}
    }
    close INFILE;

    if ($has_data and $timestep > -1) {
	push @{ $OUTDATA }, $ATOMDATA;
	$OUTDATA->[$counter]{"TIMESTEP"} = $timestep;
    }

    die "Error: $filename is not a valid Lammps data file\n"
	if (! $isvalid);
	
    return $OUTDATA;
}

sub IsValidData(@) {
    my ($atomNum, $atomType) = @_;
    my ($isvalid);

    $isvalid = 0;

    if ($atomNum <= $maxatoms or $maxatoms == -1) {
	if ($PARMS{"SOLUTE"}{"incS"} and $atomType == $PARMS{"SOLUTE"}{"itypeSol"}) {
	    $isvalid = 1;
	} elsif ( $atomType == $PARMS{"SOLVENT"}{"itypef"}) {
	    $isvalid = 1;
	} elsif ( $atomType == $PARMS{"IONS"}{"itypefion"}) {
	    $isvalid = 1;
	}
    }
    
    return $isvalid;

}

sub CollectStats(@) {
    my ($fcounter) = $_[0];
    my ($itype, @nf, @nads, @vaSum, $timeCounter, $atom, $atmCounter);
    my ($rsq, $ibin, $nstep);
    my ($radnorminv, $vr, $vt, $nsample);

    $nsample = $nstep = 0;
    for $timeCounter (0 .. $#{ $DATA[$fcounter]{"ATOMS"} }) {
	
	$nstep = $DATA[$fcounter]{"ATOMS"}->[$timeCounter]{"TIMESTEP"};
	if ($nstep > $PARMS{"ANALYSIS"}{"nskip"}) {
	    $nsample++;
	}
	for $itype (0 ..$ntype) {
	    $nf[$itype] = $nads[$itype] = $vaSum[$itype] = 0;
	}
	for $atmCounter (keys %{ $DATA[$fcounter]{"ATOMS"}[$timeCounter] }) {
	    next
		if ($atmCounter eq "TIMESTEP");
	    $atom = \%{ $DATA[$fcounter]{"ATOMS"}[$timeCounter]{$atmCounter} };
	    $rsq = $atom->{"YCOORD"}**2 + $atom->{"ZCOORD"}**2;
	    $itype = $atom->{"TYPEID"} - 1;
	    $ibin = $rsq/$delRsq[$itype] + 1;
	    if ($ibin > ($PARMS{"ANALYSIS"}{"maxbin"} -1)) {
		printf "WARNING: Exceeded upper bin range: %12.6f%12.6f%12.6f%12.6f%12.6f",
		sqrt($rsq), $delRsq[$itype],  $atom->{"XCOORD"},  $atom->{"YCOORD"},
		$atom->{"ZCOORD"};
		printf "%6d%6d%6d\n", $atom->{"NUM"}, ($itype + 1), $nstep;
		$ibin = ($PARMS{"ANALYSIS"}{"maxbin"} -1);
	    } elsif ($ibin < 0) {
		$ibin = 1;
	    }
	    
	    ($bin, $VBin, $FBin) = CalcBins($itype, $ibin, $fcounter, $rsq, $timeCounter, 
					    $atmCounter, $atom, $bin, $VBin, $FBin);
	    
	    ($bin, $VBin, $FBin) = CalcBins($itype, $ibin, $fcounter, $rsq, $timeCounter, 
					    $atmCounter, $atom, $bin, $VBin, $FBin)
		if ($rsq > $RsqWall[$itype]);
	    
#     Check if absorbed
	    $nf[$itype] += 1;
	    $nads[$itype] += 1
		if ($rsq > $rsqAds[$itype]);
	    
#     Compute average axial velocity
	    $vaSum[$itype] += $DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"XCOORD"};
	    
	}
	printf ADSFILE "";
	printf VAVFILE "";

	printf ADSFILE "%6d%12.5f%6d%6d\n", $nstep, ($nads[$itype]/$nf[$itype]),1,$ntype;
	printf VAVFILE "%6d%12.5f%6d%6d\n", $nstep, ($vaSum[$itype]/$nf[$itype]),1,$ntype;
	
	DoBlockAverages($nsample)
	    if (($nsample % $PARMS{"ANALYSIS"}{"nblock"}) == 0);
    }
}
	
sub CalcBins(@) {
    my ($itype, $ibin, $fcounter, $rsq, $timeCounter, $atmCounter, $atom, $bin, $VBin, $FBin) = @_;
    my ($radnorminv, $vr, $vt);

    $bin->{$itype}{$ibin} += 1;
    $VBin->{"X"}{"BIN"}{$itype}{$ibin} += $DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"XCOORD"};
    $VBin->{"Y"}{"BIN"}{$itype}{$ibin} += $DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"YCOORD"};
    $VBin->{"Z"}{"BIN"}{$itype}{$ibin} += $DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"ZCOORD"};
    $VBin->{"X"}{"SQBIN"}{$itype}{$ibin} += $DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"XCOORD"}**2;
    $VBin->{"Y"}{"SQBIN"}{$itype}{$ibin} += $DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"YCOORD"}**2;
    $VBin->{"Z"}{"SQBIN"}{$itype}{$ibin} += $DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"ZCOORD"}**2;
    
    $radnorminv = 1/sqrt($rsq);
    $vr = ($DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"YCOORD"} * $atom->{"YCOORD"} +
	   $DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"ZCOORD"} * $atom->{"ZCOORD"})
	* $radnorminv;
    $vt = (-$DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"YCOORD"} * $atom->{"ZCOORD"} +
	   $DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"ZCOORD"} * $atom->{"YCOORD"})
	* $radnorminv;
    $VBin->{"X"}{"VRBIN"}{$itype}{$ibin} += $DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"XCOORD"} * $vr;
    $VBin->{"X"}{"VTBIN"}{$itype}{$ibin} += $DATA[$fcounter]{"VELS"}[$timeCounter]{$atmCounter}{"XCOORD"} * $vt;
    $VBin->{"R"}{"VTBIN"}{$itype}{$ibin} += $vt * $vr;
    
    
    if ($PARMS{"ANALYSIS"}{"Lforce"}) {
	$FBin->{"X"}{"BIN"}{$itype}{$ibin} += $DATA[$fcounter]{"FORCES"}[$timeCounter]{$atmCounter}{"XCOORD"};
	$FBin->{"Y"}{"BIN"}{$itype}{$ibin} += $DATA[$fcounter]{"FORCES"}[$timeCounter]{$atmCounter}{"YCOORD"};
	$FBin->{"Z"}{"BIN"}{$itype}{$ibin} += $DATA[$fcounter]{"FORCES"}[$timeCounter]{$atmCounter}{"ZCOORD"};
    }

    return ($bin, $VBin, $FBin);
}

sub SymmetrizeData() {
    my ($timeCounter, $fcounter, $Astep, $Vstep, $Fstep, $counter);
    my ($ACounter, $VCounter, $FCounter, $AHolder, $FHolder, $VHolder, $isValid);

    $counter = 0;
    for $fcounter (0 .. $#{ $FILES }) {
   
	$isValid = 0;
	$Astep = GetTimeSteps($fcounter, 1);
	$Vstep = GetTimeSteps($fcounter, 2);
	$Fstep = GetTimeSteps($fcounter, 3);

	for $VCounter (0 .. $#{ $Vstep }) {
	    for $FCounter (0 .. $#{ $Fstep }) {
		if ($Vstep->[$VCounter] == $Fstep->[$FCounter] and 
		    $Vstep->[$VCounter] <= $PARMS{"ANALYSIS"}{"maxstep"}) {
		    for $ACounter (0 .. $#{ $Astep }) {
			if ($Astep->[$ACounter] == $Fstep->[$FCounter]) {
			    push @{ $AHolder}, \%{ $DATA[$fcounter]{"ATOMS"}[$ACounter] };
			    push @{ $VHolder}, \%{ $DATA[$fcounter]{"VELS"}[$VCounter] };
			    push @{ $FHolder}, \%{ $DATA[$fcounter]{"FORCES"}[$FCounter] };
			    $counter++;
			    $isValid = 1;
			    last;
			}
		    }
		    last;
		}
	    }
	}
	
	if ($isValid) {
	    print "Found $counter compatible time steps\n";
	    $DATA[$fcounter]{"ATOMS"} = $AHolder;
	    $DATA[$fcounter]{"VELS"} = $VHolder;
	    $DATA[$fcounter]{"FORCES"} = $FHolder;
	    return 1;
	} else {
	    return 0;
	}
    }
}

sub GetTimeSteps(@) {
    my ($fcounter, $which_file) = @_;
    
    my ($chosen_file, $time_counter, $returnval);

    if ($which_file > 1) {
        if ($which_file == 3) {
	    $chosen_file = \@{ $DATA[$fcounter]{"FORCES"} };
	} else {
	    $chosen_file = \@{ $DATA[$fcounter]{"VELS"} };
	}
    } else {
	$chosen_file = \@{ $DATA[$fcounter]{"ATOMS"} };
    }
	
    for $time_counter (@{ $chosen_file }) {
	push @{ $returnval }, $time_counter->{"TIMESTEP"};
    }

    return $returnval;
}
   
sub InOut(@) {
    my ($isOpen, $fcounter) = @_;
    my ($file_key, $file_val, $eval_str);

    if ($isOpen) {
	for $file_key (keys %{ $FILES->[$fcounter]{"OUTFILE"} }) {
	    $file_val = $FILES->[$fcounter]{"OUTFILE"}{$file_key};
	    $file_key = uc($file_key);
	    $eval_str = qq(open $file_key, "> $file_val");
	    eval($eval_str);
	    if ($@) {
		die "ERROR: Cannot open $file_val: $@\n";
	    }
	}
    } else {
	for $file_key (keys %{ $FILES->[$fcounter]{"OUTFILE"} }) {
	    $file_key = uc($file_key);
	    $eval_str = "close " . $file_key;
	    eval($eval_str);
	}
    }
}

sub DoBlockAverages(@) {
    my ($nsample) = @_;
    my ($iblock, $sumavrhof, $itype, $vbin, $rhofac, $sumrhof, $nbin, $ibin);
    my ($r, $bini, $rhof, $avrhof);

    $nbin = $PARMS{"ANALYSIS"}{"maxbin"} -1;
    $iblock = $nsample/$PARMS{"ANALYSIS"}{"nblock"};
    $sumavrhof = 0;
    print "Completed Block # $iblock\n";
	
	for $itype (0 .. $ntype) {
	    print DENFILE "# Density Profile for Type $itype Block # $iblock\n";
		print VELFILE "# X, Y, Z Velocity Profile for Type $itype Block # $iblock\n";
		print VSQFILE "# Vx*Vx,Vy*Vy,Vz*Vz,Vx*Vr,Vx*Vt,Vr*Vt for Type $itype Block # $iblock\n";
		print FRCFILE "# X, Y, Z Force Profile for Type $itype Block # $iblock\n"
		if ($PARMS{"ANALYSIS"}{"Lforce"});
	    $vbin = PI * $delRsq[$itype] * $PARMS{"CYLINDER"}{"xc"};
	    $rhofac = 1.0/($vbin*$PARMS{"ANALYSIS"}{"nblock"});
	    $sumrhof = 0.0;
	    for $ibin (0 .. $nbin) {
		($ibin == 1) ?
		    $r = 0 :
		    $r = sqrt(($ibin - 0.5)* $delRsq[$itype]);
		$bini = $bin->{$itype}{$ibin};
		$rhof = $bini * $rhofac;
		print DENFILE "$r,$rhof\n";
		if ($bini > 0) {
		    WriteBins($r, $itype, $ibin, $bini);
		}
		
		$ibin++;
		$r = sqrt(($rf[$itype]**2 + $RsqWall[$itype])/2);
		$bini = $bin->{$itype}{$ibin};
		$vbin = PI * ($rf[$itype]**2 - $RsqWall[$itype]) * $PARMS{"CYLINDER"}{"xc"};
		$rhofac = 1.0/($vbin*$PARMS{"ANALYSIS"}{"nblock"});
		$rhof = $bini*$rhofac;
		print DENFILE "$r,$rhof\n";
		if ($bini > 0) {
		    WriteBins($r, $itype, $ibin, $bini);
		}
	    }
	    print DENFILE "End of set\n";
	    print VELFILE "End of set\n";
	    print VSQFILE "End of set\n";
	    print FRCFILE "End of set\n"
		if ($PARMS{"ANALYSIS"}{"Lforce"});
	    
	    $avrhof = $sumrhof/$nbin;
	    $sumavrhof = $sumavrhof+$avrhof;
	    print "Average density for type $itype $avrhof\n";
	}

    print "Average total density = $sumavrhof\n";
}

sub WriteBins(@) {
    my ($r, $itype, $ibin, $bini) = @_;
    
    printf VELFILE "%12.7f%12.7f%12.7f%12.7f\n", $r, 
    $VBin->{"X"}{"BIN"}{$itype}{$ibin}/$bini,$VBin->{"Y"}{"BIN"}{$itype}{$ibin}/$bini,
    $VBin->{"Z"}{"BIN"}{$itype}{$ibin}/$bini;
    
    printf VSQFILE "%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n", $r,
    $VBin->{"X"}{"SQBIN"}{$itype}{$ibin}/$bini, $VBin->{"Y"}{"SQBIN"}{$itype}{$ibin}/$bini,
    $VBin->{"Z"}{"SQBIN"}{$itype}{$ibin}/$bini, $VBin->{"X"}{"VRBIN"}{$itype}{$ibin}/$bini,
    $VBin->{"X"}{"VTBIN"}{$itype}{$ibin}/$bini, $VBin->{"R"}{"VTBIN"}{$itype}{$ibin}/$bini;
    
    if ($PARMS{"ANALYSIS"}{"Lforce"}) {
	print FRCFILE "%12.7f%12.7f%12.7f%12.7f\n", $r, 
	$FBin->{"X"}{"BIN"}{$itype}{$ibin}/$bini, $FBin->{"Y"}{"BIN"}{$itype}{$ibin}/$bini,
	$FBin->{"Z"}{"BIN"}{$itype}{$ibin}/$bini;
    }
}

sub ValidTimeStep(@) {
    my ($ftype, $tstep, $fcounter) = @_;
    my ($returnval, $isvalid);
    my ($timeCounter, $Vstep, $VCounter);

    if ($ftype == 1) {
	$returnval = 1;
    } else {
	$Vstep = GetTimeSteps($fcounter, 2);
	$VCounter = $returnval = 0;
	while ($VCounter <= $#{ $Vstep }) {
	    if ($tstep < $Vstep->[$VCounter]) {
		$returnval = 0;
		last;
	    } elsif ($tstep == $Vstep->[$VCounter]) {
		$returnval = 1;
		last;
	    } else {
		$VCounter++;
	    }
	}
    }
    return $returnval;
}

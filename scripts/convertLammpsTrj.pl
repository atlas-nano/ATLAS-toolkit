#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use strict;
use General qw(FileTester TrjSelections LoadElements Trim CenterOnMol CoM 
			GetSoluteAtoms GetBondLength AddElementField GetFileTypeStr);
use FileFormats qw(ParseStructFile createHeaders addHeader createBGF createPDB GetBGFAtoms);
use LAMMPS qw(ParseLAMMPSTrj CreateLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use AMBER qw(CreateAmberTrj);
use XYZ qw(CreateXYZTrj);
use ManipAtoms qw(UnwrapAtoms ScaleAtoms ImageAtoms SelectAtoms BuildAtomSelectionString  
								SplitAtomsByMol GetAtmData GetMols ReimageAtoms ReimageMols);
use BOX qw(MoveAtomsToOrigin InitBox);
use REPLICATE qw(GetBoxVol);
use Superimpose;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use Storable qw(dclone);

sub numerically { ($a<=>$b); }

my ($lammpsTrj, $struct_file, $selection, $reImageCoord, $saveName, $outType, $comTrj, $slabDim, $shellOpt, $si, $fixLongBonds);
my ($printAtoms, $SELECT, $LAMMPSOPTS, $pStr, $SOLUTEATMS, $count, $reScaleCoord, $MOLS, $atomSelectX);
my ($BGF, $BONDS, $soluteSelect, $isTrj, $molList, $comBGF, $imgPoint, $savePt, $com_center);
my ($AMBERVELHANDLE, $writeAMBERvels, $x, $y, $z, $atmCount,$frameNum, $ATOMSEL, $iHEADER);

$|++;
&init;
print "Parsing BGF file $struct_file...";
($BGF, $BONDS, $iHEADER) = ParseStructFile($struct_file, 1);
&checkAtomTypes($BGF);
$atmCount = scalar(keys %{ $BGF });
$MOLS = GetMols($BGF, $BONDS);
($SOLUTEATMS, $molList, $ATOMSEL) = getSoluteSolvent($BGF, $BONDS, $soluteSelect, $MOLS, $atomSelectX);
#die "ERROR: No valid solute atoms found!\n" if (! keys %{ $SOLUTEATMS });
#die "ERROR: No valid solvent atoms found!\n" if (! keys %{ $molList });
$comBGF = createCOMBGF($BGF) if ($comTrj);
print "Done\n";
&GetLammpsByteOffset($SELECT, $lammpsTrj, scalar keys %{ $BGF });
&GetLammpsTrjType($SELECT, $lammpsTrj, "atom", \%{ $LAMMPSOPTS });
$AMBERVELHANDLE = &initAMBERvels($saveName) if ($isTrj and $outType eq "AMBER" and $LAMMPSOPTS->{hasvel});
$pStr = "Parsing LAMMPS trajectory $lammpsTrj...";
$count = 0;
&ParseLAMMPSTrj($BGF, $lammpsTrj, $SELECT, "atom", \&saveCoords, $pStr, \*OUTDATA);
close OUTDATA or die "ERROR: Cannot close $saveName: $!\n" if ($isTrj);
close $AMBERVELHANDLE or die "ERROR: Cannot close AMBER velocity file: $!\n" if ($isTrj and $writeAMBERvels);
&writeCOMBGF($comBGF, $saveName) if ($comTrj);

sub initAMBERvels {
	my ($coordname) = @_;
	my ($velname) = $coordname;
	$velname =~ s/\.(\w)+$//;
	$velname .= ".mdvel";
	$velname = ".vel.mdcrd" if($1 eq "mdvel");

	open AMBERVEL, "> $velname" or die "ERROR: Cannot create AMBER velocity file $velname: $!\n";
	print AMBERVEL "TITLE: $outType Velocity file create by lammpsTrj2amber on " . scalar(localtime(time())) . "\n";

	$writeAMBERvels = 1;
	return \*AMBERVEL;
}

sub saveCoords {
	my ($DATA, $bgfInfo, $fileHandle) = @_;
	my ($BOX, $mID, $HEADERS, @tmp, $CENTER, $j, $tstep); 
	my ($MOLECULE, $i, $fName, $index, $totAtms, $currBonds);

	$tstep = $DATA->{"TIMESTEP"}[0];
	$frameNum++;
	$BOX = ConvertLammpsBox($DATA->{"BOX BOUNDS"});
	$BOX->{X}{angle}=$BOX->{ALPHA}; $BOX->{Y}{angle}=$BOX->{GAMMA}; $BOX->{Z}{angle}=$BOX->{GAMMA};
	&InitBox($BOX, $DATA->{ATOMS});
	&GetBoxVol($BOX);
	$totAtms = $DATA->{"NUMBER OF ATOMS"}[0];
	&addMissingAtomCoords($DATA->{ATOMS}, $BOX, $BGF, $frameNum, scalar(keys %{ $SELECT })) 
		if($totAtms < $atmCount);
	&UnwrapAtoms($DATA->{ATOMS}, $BOX, $LAMMPSOPTS->{scaled}) if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged});
	&fLB($DATA->{ATOMS}, $BONDS, $BOX) if ($fixLongBonds);

	#here we align the solute atoms if necessary
	if ($si->{Opt}) {
		for $i (keys %{ $BGF }) {
			for $j ("MASS", "FFTYPE", "ATMNAME") {
				$DATA->{ATOMS}{$i}{$j} = $BGF->{$i}{$j};
			}
		}
		if($frameNum==1) { #save the reference
			#$si->{REF} = GetAtmData($DATA->{ATOMS}, $SOLUTEATMS);
			$si->{REF} = GetAtmData($BGF, $SOLUTEATMS);
		}
		#} else {
			$si->{CURR} = GetAtmData($DATA->{ATOMS}, $SOLUTEATMS);
			&Kabash::SuperimposeAtoms($si->{REF}, $si->{CURR}, $DATA->{ATOMS});
		#}
	}
	if ($reImageCoord) {
		@tmp = ("XCOORD", "YCOORD", "ZCOORD");
		for $j (@tmp) {
			$CENTER->{$j} = $BOX->{$j}{CENTER} = $BOX->{$j}{len}/2;
		}
		if (defined($SOLUTEATMS)) {
			if(scalar(keys%{ $SOLUTEATMS }) > 1) {
				$MOLECULE = GetAtmData($DATA->{ATOMS}, $SOLUTEATMS);
				$CENTER = CoM($MOLECULE, $BOX);
			} else {
				for $j (@tmp) {
					$CENTER->{$j} = $DATA->{ATOMS}{ $SOLUTEATMS->{1} }{$j}
				}
			}
		} elsif (defined($imgPoint)) {
			$CENTER = $imgPoint;
		} elsif (defined($com_center)) {
			$CENTER = CoM($DATA->{ATOMS}, $BOX);
		}
		for $j (@tmp) {
			$BOX->{$j}{hi} = $BOX->{$j}{len};
			$BOX->{$j}{lo} = 0;
			$CENTER->{$j} = $BOX->{$j}{CENTER} if(defined($slabDim) and $j eq $slabDim);
			for $i (keys %{ $DATA->{ATOMS} }) {
				$DATA->{ATOMS}{$i}{$j} += ($BOX->{$j}{CENTER} - $CENTER->{$j});
			}
		}
		#for $index (keys %{ $MOLS }) {
		#	$MOLECULE = GetAtmData($DATA->{ATOMS}, $MOLS->{$index}{MEMBERS});
		#	$CENTER = CoM($MOLECULE, $BOX);
		#	ImageAtoms($MOLECULE, $CENTER, $BOX);
		#}
		&ReimageMols($DATA->{ATOMS}, $MOLS, $BOX);
		if ($savePt) {
			%{ $DATA->{ATOMS}{$totAtms+1} } = %{ $imgPoint };
		}
	}
	%{ $currBonds } = %{ $BONDS};
	&updateBGF($bgfInfo, $DATA->{ATOMS}, $BOX, $DATA->{TIMESTEP}[0]);
	&addShellAtoms($bgfInfo, $currBonds) if ($shellOpt);
	$bgfInfo = createCOMTrj($bgfInfo, $comBGF) if ($comTrj);
	($bgfInfo, $currBonds, undef) = GetBGFAtoms($ATOMSEL, $bgfInfo, $currBonds) if (defined($ATOMSEL));
	if ($outType eq "AMBER") {
		delete $bgfInfo->{HEADER};
		&CreateAmberTrj($bgfInfo,$BOX, $fileHandle);
		if($writeAMBERvels) {
		   for $i (keys %{ $bgfInfo }) {
			 for $j ("X", "Y", "Z") {
				$bgfInfo->{$i}{"${j}COORD"} = $bgfInfo->{$i}{"${j}VEL"} * 1000/20.455; #copy vel info into coord field
			 }
		   }
		   &CreateAmberTrj($bgfInfo,undef, $AMBERVELHANDLE);
		}
	} elsif ($outType eq "LAMMPS") {
		delete $bgfInfo->{HEADER};
		ScaleAtoms($bgfInfo, $BOX) if ($reScaleCoord);
		#$HEADERS = createLAMMPSHeader($DATA, $BOX);
		if ($reScaleCoord) {
			push @{ $HEADERS->{cType} }, "xs ys zs";
		} else {
			push @{ $HEADERS->{cType} }, "xu yu zu";
		}
		CreateLAMMPSTrj($bgfInfo, $DATA, $fileHandle);
	} elsif ($outType eq "XYZ") {
		delete $bgfInfo->{HEADER};
		&CreateXYZTrj($bgfInfo,$tstep,$fileHandle);
	} else { #bgf
		$count++;
		$fName = $saveName;
		@tmp = keys %{ $SELECT };
		if (defined($SELECT) && $#tmp > 0) {
			$fName =~ s/\.\w+$//;
			$fName .= ".${tstep}." . lc($outType);
		}
		for $i (keys %{ $bgfInfo }) {
			next if (ref($bgfInfo->{$i}) ne "HASH" or !exists($bgfInfo->{$i}{CHARGE}));
		}
		if($outType eq "BGF") {
			for $i (@{ $iHEADER }) {
				next if ($i !~ /^REMARK/);
				splice @{ $bgfInfo->{HEADER} }, 3, 0, $i;
			}	
			&createBGF($bgfInfo, $currBonds, $fName);
		} elsif ($outType eq "PDB") {
			delete $bgfInfo->{HEADER};
			&createPDB($bgfInfo, $currBonds, $fName);
		}
	}
}

sub fLB {
	my ($atoms, $bonds, $box) = @_;
	my ($i, $j, $k, $d, $ad, $bi, $bl);

    for $i (sort numerically keys %{ $atoms }) {
		next if (${ $BGF->{$i}{MOLSIZE} } != 3); #only the water molecules
        for $j (@{ $bonds->{$i} }) {
			next if ($j>$i);
            for $k ("XCOORD", "YCOORD", "ZCOORD") {
				$bl = $box->{$k}{len};
                $d = $atoms->{$i}{$k}-$atoms->{$j}{$k};
				$ad = abs($d);
				if($ad>5 and $ad>$bl/2) {
					$bi = sprintf("%.0f",($d/$bl));
					$atoms->{$j}{$k} += $bi*$bl;
				}
            }
        }
    }
}

sub addMissingAtomCoords {
	my ($atoms, $box, $oldCoords,$currFrame,$totFrame) = @_;
	my ($newCoords, $i, $j, $del, $valid);
	my ($MOLECULE, $CENTER,$inc);

	$del = 1.0;
	$newCoords = ();
	if (! defined($x)) {
	  $x = $box->{XCOORD}{lo};
	  $y = $box->{YCOORD}{lo};
	  $z = $box->{ZCOORD}{lo};
	}

	for $j (sort numerically keys %{ $MOLS }) {
	  $valid = 0;
	  if (!exists($MOLS->{$j}{FINAL})) {
		for $i (keys %{ $MOLS->{$j}{MEMBERS} }) {
		  if(keys %{ $atoms->{$i} }) {
			$valid = 1;
			last;
		  }
		}
		next if($valid);
		$x += $del;
		if($x > $box->{XCOORD}{hi}) {
		  $x = 0;
		  $y += $del;
		}
		if ($y > $box->{YCOORD}{hi}) {
			$y = 0;
		  $z += $del;
		}
		$inc = $totFrame - $currFrame + 2;
		$MOLS->{$j}{FINAL}{XCOORD} = $x;
		$MOLS->{$j}{FINAL}{YCOORD} = $z;
		$MOLS->{$j}{FINAL}{ZCOORD} = $y;
		$MOLECULE = GetAtmData($oldCoords, $MOLS->{$j}{MEMBERS});
		$CENTER = CoM($MOLECULE);
		$MOLS->{$j}{INITIAL}{XCOORD} = $CENTER->{XCOORD};
		$MOLS->{$j}{INITIAL}{YCOORD} = $CENTER->{YCOORD};
		$MOLS->{$j}{INITIAL}{ZCOORD} = 17;
		$MOLS->{$j}{TOTAL} = $inc;
		for ("XCOORD", "YCOORD", "ZCOORD") {
			$MOLS->{$j}{INC}{$_} = ($MOLS->{$j}{FINAL}{$_}-$MOLS->{$j}{INITIAL}{$_})/$inc;
		}
		for $i (keys %{ $MOLS->{$j}{MEMBERS} }) {
			for ("XCOORD", "YCOORD", "ZCOORD") {
			   $MOLS->{$j}{ATOMOFFSET}{$i}{$_} = $oldCoords->{$i}{$_} - $CENTER->{$_};
			}
		}
		$MOLS->{$j}{COUNT} = 1;
	  }
	  for $i (keys %{ $MOLS->{$j}{MEMBERS} }) {
		$atoms->{$i} = \%{ $oldCoords->{$i} };
		for ("XCOORD","YCOORD","ZCOORD") {
			$atoms->{$i}{$_} = $MOLS->{$j}{ATOMOFFSET}{$_} + $MOLS->{$j}{INITIAL}{$_} + $MOLS->{$j}{INC}{$_}*$MOLS->{$j}{COUNT};
		}
	  }
	  $MOLS->{$j}{COUNT}++;
	}
}

sub createCOMAtom {
	my ($index) = $_[0];
	my (%ATOM);

	$ATOM{"INDEX"} = $index;
	$ATOM{"ATMNAME"} = "Na";
	$ATOM{"RESNAME"} = "Na";
	$ATOM{"RESNUM"} = $index;
	$ATOM{"FFTYPE"} = "Na";
	$ATOM{"NUMBONDS"} = 0;
	$ATOM{"LONEPAIRS"} = 0;
	$ATOM{"CHARGE"} = 0;
	$ATOM{"RADII"} = 0;
	$ATOM{"LABEL"} = "HETATM";
	return \%ATOM;
}

sub createCOMBGF {
	my ($atoms) = $_[0];
	my ($mols, $i, $comAtoms, $comBGF, $j, $atomCoord);

	$mols = SplitAtomsByMol($BGF, $BGF);

	for $i (keys %{ $mols }) {
		$comBGF->{$i} = createCOMAtom($i);
		$comAtoms = GetAtmData($atoms, $mols->{$i});
		for $j (keys %{ $mols->{$i} }) {
			$comBGF->{$i}{MOLECULE}{MEMBERS}{$j} = $j;
		}
		$atomCoord = CoM($comAtoms);
		for $j (keys %{ $atomCoord }) {
			$comBGF->{$i}{$j} = $atomCoord->{$j};
		}
	}

	return $comBGF;
}

sub createCOMTrj {
	my ($atomCoords, $comAtoms) = @_;
	my ($i, $atomCoord, $j, $atoms);

	for $i (keys %{ $comAtoms }) {
		$atoms = GetAtmData($atomCoords, $comAtoms->{$i}{MOLECULE}{MEMBERS});
		$atomCoord = CoM($atoms);
		$comAtoms->{$i}{TYPE} = 1;
		for $j (keys %{ $atomCoord }) {
			$comAtoms->{$i}{$j} = $atomCoord->{$j};
		}
	}

	return $comAtoms;
}

sub writeCOMBGF {
	my ($comAtoms, $trjSaveName) = @_;
	my ($bgfSaveName, $i, $comBonds, $comHeaders);

	$bgfSaveName = $trjSaveName;
	$bgfSaveName =~ s/\.\w+$//;
	$bgfSaveName .= "_com.bgf";
	print "Creating COM BGF file $bgfSaveName...";
	for $i (keys %{ $comAtoms }) {
		$comBonds->{$i} = ();
	}
	$comHeaders = createHeaders(undef, $bgfSaveName);
	&addHeader($comAtoms, $comHeaders);
	&createBGF($comAtoms, $comBonds, $bgfSaveName);
	print "Done\n";
}

sub doResFix {
	my ($ATOMS, $BOX, $currAtm) = @_;
	my ($i, @dims, $j, $bondPos, $atomPos);

	$ATOMS->{$currAtm}{VISITED} = 1;
	@dims = ("X", "Y", "Z");
	for $j (@dims) {
		$atomPos->{$j} = $ATOMS->{$currAtm}{$j . "COORD"};
	}
	for $i ( @{ $BONDS->{$currAtm} }) {
		next if (exists($ATOMS->{$i}{VISITED}));
		for $j (@dims) {
			$bondPos = $ATOMS->{$i}{$j . "COORD"};
			if(abs($bondPos - $atomPos->{$j}) > 5) { #image problem
				if ($atomPos->{$j} > $bondPos) { # move in + direction
					while (($atomPos->{$j} - $bondPos) > 5) {
						$bondPos += ($BOX->{$j}{hi} - $BOX->{$j}{lo});
					}
				} else { # move in - direction
					while (($bondPos - $atomPos->{$j}) > 5) {
						$bondPos -= ($BOX->{$j}{hi} - $BOX->{$j}{lo});
					}
				}
				$ATOMS->{$i}{$j . "COORD"} = $bondPos;
				$ATOMS->{$i}{$j . "INDEX"} = $ATOMS->{$currAtm}{$j . "INDEX"};
			}
		}
		doResFix($ATOMS, $BOX, $i);
	}
}			 
			
sub updateBGF {
	my ($REF, $MOD, $BOX, $TSTEP) = @_;
	my (@tmp, $i, $HEADERS, $j);

	$BOX->{1}{DATA} = 90;
	@tmp = ("XCOORD", "YCOORD", "ZCOORD");
	for $i (2 .. 4) {
		$BOX->{$i}{DATA} = ($BOX->{$tmp[($i - 2)]}{hi} - $BOX->{$tmp[($i - 2)]}{lo});
	}

	$TSTEP = "ts_${TSTEP}";
	$HEADERS = createHeaders($BOX, $TSTEP);

	for $i (keys %{ $REF }) {
		next if (!exists($MOD->{$i}));
		for $j (keys %{ $MOD->{$i} }) {
			$REF->{$i}{$j} = $MOD->{$i}{$j};
		}
	}
	&addHeader($REF, $HEADERS);
}

sub addShellAtoms {
	my ($atoms, $bonds) = @_;
	my ($i, $flist, $tot, $shellAtm);

	$tot = scalar(keys %{ $atoms }) - 1;
	$shellAtm = (
					{
						"HEADER"	 => "HETATM",
						"INDEX"	     => -1,
						"ATMNAME"	 => "SATM",
						"RESNAME"	 => "RES",
						"RESNUM"	 => 999,
						"XCOORD"	 => -1,
						"YCOORD"	 => -1,
						"ZCOORD"	  => -1,
						"FFTYPE"	  => "He",
						"NUMBONDS"    => 1,
						"LONEPAIRS"   => 0,
						"CHARGE"	  => -1,
						"MOLECULEID"  => -1,
					}
				);
	for $i (keys %{ $atoms }) {
		next if ($i !~ /^\d+$/);
		next if (!exists($atoms->{$i}{RSX}) or ! exists($atoms->{$i}{RSY}) or ! exists($atoms->{$i}{RSZ}));
		$tot++;
		%{ $atoms->{$tot} } = %{ $shellAtm };
		$atoms->{$tot}{INDEX} = $tot;
		$atoms->{$tot}{RESNAME} = $atoms->{$i}{RESNAME};
		$atoms->{$tot}{RESNUM} = $atoms->{$i}{RESNUM};
		$atoms->{$tot}{XCOORD} = $atoms->{$i}{XCOORD}+$atoms->{$i}{RSX};
		$atoms->{$tot}{YCOORD} = $atoms->{$i}{YCOORD}+$atoms->{$i}{RSY};
		$atoms->{$tot}{ZCOORD} = $atoms->{$i}{ZCOORD}+$atoms->{$i}{RSZ};
		$atoms->{$tot}{MOLECULEID} = $atoms->{$i}{MOLECULEID};
		$atoms->{$tot}{CHARGE} = $atoms->{$i}{CHARGE} - 1;
		$atoms->{$i}{CHARGE} = 1;
		push @{ $bonds->{$i} }, $tot;
		push @{ $bonds->{$tot} }, $i;
	}
}

sub getSoluteSolvent {
	my ($atoms, $bonds, $selection, $mols, $extractAtmSel) = @_;
	my ($solute, $solventMols, $i, $solventList, $extractAtms);
	
	$solute = SelectAtoms($selection, $atoms) if (defined($selection));
	$extractAtms = SelectAtoms($extractAtmSel, $atoms) if (defined($extractAtmSel));
	for $i (keys %{ $atoms }) {
		next if (defined($solute) && exists($solute->{$i}));
		$solventMols->{${ $atoms->{$i}{MOLECULEID} }} = $atoms->{$i}{MOLECULE};
	}
	return ($solute, $solventMols, $extractAtms);
}


sub init {
	my ($i, %OPTS, $usage, $atomSelection, @tmp, $coordOpt, $pointStr, $extractAtmStr);

	getopt('lstboamcpqdnxhrf',\%OPTS);
	$usage = &showUsage;

	($lammpsTrj,$struct_file,$selection,$outType,$saveName,$coordOpt,$atomSelection,$comTrj,$pointStr,$savePt, $slabDim,$com_center, $shellOpt, $extractAtmStr, $si->{Opt}, $fixLongBonds) = 
	($OPTS{l},  $OPTS{b},$OPTS{t},  $OPTS{o},$OPTS{s}, $OPTS{m}, $OPTS{a},      $OPTS{c},$OPTS{p},$OPTS{q},$OPTS{d},$OPTS{n},    $OPTS{h},  $OPTS{x},       $OPTS{r},   $OPTS{f});

	for ($lammpsTrj, $struct_file, $selection) {
		die "$usage\n" if (! defined($_));
	}
	print "Initializing...";

	FileTester($lammpsTrj);
	
	$outType = "lammps" if (! defined($outType));
	if ($outType !~ /(lammps|amber|bgf|xyz|pdb)/i) {
		print "WARNING: Expected lammps, xyz, amber, pdb or bgf for out_type. Got $outType...\n" . 
			"Setting to default value: lammps\n";
		$outType = "lammps";
	} else {
		$outType = $1;
	}
	print "$outType traj ";

	$coordOpt = 0 if (! defined($coordOpt) or $coordOpt !~ /^(1|2|3)$/);
	if ($coordOpt == 0) {
		$reScaleCoord = 0;
		$reImageCoord = 0;
	} elsif($coordOpt==1) {
		$reScaleCoord = 0;
		$reImageCoord = 1;
	} elsif($coordOpt==2) {
		$reScaleCoord = 1;
		$reImageCoord = 1;
	} else {
		$reScaleCoord = 1;
		$reImageCoord = 0;
	}
	print "with coordinate imaging.." if ($reImageCoord);
	print "without coordinate imaging..." if (! $reImageCoord);

	$si->{Opt} = 0 if (! defined($si->{Opt}) or ! $coordOpt or $si->{Opt} !~ /1|yes/i);
	$si->{Opt} = 1 if ($si->{Opt} =~ /1|yes/i);
	print "superimposing solute..." if ($si->{Opt});

	if (! defined($saveName)) {
		$saveName = basename($lammpsTrj);
		$saveName =~ s/\.\w+$/_mod\.lammpstrj/;
	}

	$SELECT = TrjSelections($selection);
	$outType = uc($outType);
	
	$isTrj = 0;
	$isTrj = 1 if ($outType =~ /LAMMPS|AMBER|XYZ/i);
	if ($isTrj) {
		open OUTDATA, "> $saveName" or die "Cannot write to trajectory file $saveName: $!\n";
		$pStr = "Creating $outType trajectory $saveName...";
		print OUTDATA "TITLE: $outType Trajectory file create by lammpsTrj2amber on " . scalar(localtime(time())) . "\n"
				if ($outType !~ /XYZ/i);
	} else {
		$pStr = "Creating $outType files for each frame...";
	}

	if (defined($atomSelection)) {
		$soluteSelect = BuildAtomSelectionString($atomSelection);
	}
	if (defined($extractAtmStr)) {
		$atomSelectX = BuildAtomSelectionString($extractAtmStr);
	}
	if (defined($pointStr) and $pointStr =~ /(\-?\d+\.?\d*)\s+(\-?\d+\.?\d*)\s+(\-?\d+\.?\d*)/) {
		($imgPoint->{XCOORD}, $imgPoint->{YCOORD}, $imgPoint->{ZCOORD}) = ($1, $2, $3);
		if (defined($savePt) && $savePt =~ /^(1|yes)/) {
			$savePt = 1;
		} else {
			$savePt = 0;
		}
	}
	undef($soluteSelect) if defined($imgPoint);
	$comTrj = 0 if (! defined($comTrj) or $comTrj !~ /^\s*(1|yes)\s*$/i);
	$comTrj = 1 if ($comTrj =~ /^\s*(1|yes)\s*$/i);
	if(defined($slabDim) && $slabDim =~ /^(x|y|z)$/) {
		$slabDim = uc $1 . "COORD";
	}

	$com_center = 0 if (! defined($com_center) or $com_center !~ /1|yes/i);
	$com_center = 1 if ($com_center =~ /1|yes/i);

	$writeAMBERvels = 0;
	$frameNum = 0;

	$fixLongBonds = 0;
	$fixLongBonds = 1 if (exists($OPTS{f}) and $OPTS{f} =~ /1|yes/i);

	print "Done\n";
}

sub checkAtomTypes {
		my ($atoms) = @_;
		my ($i, $ELEMENTS, $ffType, $atmName, $eleList);
		
		$ELEMENTS = &LoadElements;
		&AddElementField($atoms, $ELEMENTS, undef);
		for $i (keys %{ $atoms }) {
			$atoms->{$i}{MASS} = $atoms->{$i}{ELEMENT}{MASS};
		}
#		for $i (keys %{ $ELEMENTS }) {
#						$eleList->{uc($ELEMENTS->{$i}{SYMBOL})} = $i;
#		}
#
#
#		for $i (keys %{ $atoms }) {
#				($atoms->{$i}{ELEMENT}, $atoms->{$i}{ELENUM}) = getElement($eleList, $ffType, $atoms->{$i}{ATMNAME});
#				if (! defined($atoms->{$i}{ELEMENT})) {
#						$ffType = substr($ffType,0,1);
#						($atoms->{$i}{ELEMENT}, $atoms->{$i}{ELENUM}) = getElement($eleList, $ffType, $atoms->{$i}{ATMNAME});
#						if (! defined($atoms->{$i}{ELEMENT})) {
#								die "ERROR: Cannot deerming the element type of fftype $ffType\n";
#						}
#				}
#		}
}

sub getElement {
		my ($eleList, $ffType, $atmName) = @_;
		my ($i, $element);

		for $i ($ffType, $atmName) {
				$element = $i;
				next if ($element !~ /^([A-Za-z])([a-z]?)/);
				$element = $1;
				$element .= $2 if ($2);
				return ($element, $eleList->{uc($element)}) if (exists($eleList->{uc($element)}));
				$element = $1;
				return ($element, $eleList->{uc($element)}) if (exists($eleList->{uc($element)}));
		}

}
sub showUsage {
	my ($fTypeStr) = &GetFileTypeStr;
	my ($usage) = <<DATA;
usage: $0 -b struct_file -l lammps trj -t trajectory selection -n [center com 2 box center = no] -m [reimage/rescale coordinates] -p [spacial point] -a [solute atoms] -o [output type] -s [save name] -c [create com trj] -q [save spacial point] -h [include shell atoms] -x [extract atoms] -r [rotate solute]
options:
		-b struct_file (Required): the location of the structure file
$fTypeStr		
		-l lammps trj (Required): the location of the lammps trajectory file
		-t trajectory selection: The number of frames to use. Can be a single integer, or several integers in quotes
				To specify a range specify it as :Ita-b:c, which will take frames a - b, every c. Specify multiple 
				ranges or a combination of ranges and single frames by enclosing them in quotes. \"*\" for all frames
		-m [reimage/rescale coordinates]: (Optional) Whether to reimage and rescale the coordinates back into the 
				unit cell. 0:No reimage or rescale (Default), 1:reimage no rescale, 2:rescale and rescale 
				3:rescale no reimage
		-p [spacial point]: (Optional) Image the coordinates wrt to point in space - or -
		-a [solute atoms]: (Optional) Image the coordinates wrt to the center of mass of certain atoms
				Any valid field (fftype, resname, xcoord etc) expression (e.g. "resname eq 'WAT'")
		-r [rotate solute]: (Optional) Rotate the system to superimpose the solute to the initial frame. Default no		
		-q [save spacial point]: (Optional) Save the coordinates of the spacial point. Defaul no
		-d [slab dimension]: (Optional) The non periodic dimension in a slab simulation
		-o [output type]: (Optional). Either lammps (default) , xyz or  amber or bgf. If bgf, will print file for every frame
		-h [shell flag]: (Optional). If the trajectory has relative shell positions (RSX,RSY,RSZ), the create dummy shell atoms, with core charge of +1
		-c [create com trj]: convert each molecule to an atom based on their center of mass. Default no
		-s [save name]: (Optional). The name the save the output trajectory
DATA

	return $usage;
}

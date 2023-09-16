#!/usr/bin/perl
#
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo GetBGFAtoms);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use General qw(FileTester TrjSelections CoM HasCell Rotate);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use BOX qw(GetBox);
use ManipAtoms qw(UnwrapAtoms ScaleAtoms BuildAtomSelectionString SelectAtoms GetMols GetAtmData ImageAtoms);
use Math::MatrixReal;

sub init;
sub getContactAngle;
sub calcContactAngle;
sub writeDistribution;
sub numerically { $a<=>$b }

my ($bgfFile, $selection, $trjFile, $SELECT, $saveFile, $savePrefix);
my ($field, $pStr, $LAMMPSOPTS, $getByteOffset, $getSnapshots, $trjType, $OUTFILE);
my ($atmList, $ATOMS, $BONDS, $HEADERS, $MOLS, $BOX, $tmp, $cumXZ_PLANE);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
$BOX = &GetBox($ATOMS, undef, $HEADERS) if (HasCell($HEADERS));
$MOLS = GetMols($ATOMS, $BONDS);
($atmList,undef) = SelectAtoms($selection, $ATOMS);
print "Done\n";
if (! defined($trjFile)) {
	($ATOMS, $BONDS, undef) = GetBGFAtoms($atmList, $ATOMS, $BONDS);
	print "==========================================================\n";
	print "Contact Angle Results\n";
	print "==========================================================\n";
	&calcContactAngle($ATOMS,$BOX,\%{ $tmp }, 1, 0);
	printf "%10s %10s %10s %10s %10s %10s %10s\n","Theta","+/-","R^2","Radius","CenterX","CenterY";
	printf "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$tmp->{angle},$tmp->{sigma},$tmp->{R2},$tmp->{R},$tmp->{xc},$tmp->{yc};
} else {
	$field = scalar keys %{ $ATOMS };
	$getByteOffset->($SELECT, $trjFile, $field);
	if ($trjType == 2) {
		&GetLammpsTrjType($SELECT, $trjFile, "coord", \%{ $LAMMPSOPTS });
		$field = "coord";
	}
	$pStr = "Calculating contact angle data from $trjFile...";
	open $OUTFILE, "> $saveFile" or die "ERROR: Cannot write to $saveFile: $!\n";
	printf $OUTFILE "%-10s %10s %10s %10s %10s %10s %10s %10s\n","Tstep","Theta","+/-","R^2","Radius","CenterX","CenterY";
	$getSnapshots->($ATOMS, $trjFile, $SELECT, $field, \&getContactAngle, $pStr, $OUTFILE);
	close $OUTFILE;
	print "==========================================================\n";
	print "Contact Angle Results\n";
	print "==========================================================\n";
	&calcContactAngle($tmp,$BOX,\%{ $tmp }, 1, 1);
	printf "%10s %10s %10s %10s %10s %10s %10s\n","Theta","+/-","R^2","Radius","CenterX","CenterY";
	printf "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$tmp->{angle},$tmp->{sigma},$tmp->{R2},$tmp->{R},$tmp->{xc},$tmp->{yc};
}

#From "A Simple approach for the Estimation of Circular Arc Center
#and Its radius", Thomas and Chan, 
#Computer vision, graphics and image processing 45, 362-370 (1989)
sub fit_Landau {
	my ($x, $y, $data) = @_;
	my ($N, $p, $i); 
	my ($a1, $a2, $b1, $b2, $c1, $c2);

	$N = $#{ $x };
	for $i (1 .. 9) {
		$p->[$i] = 0;
	}

	for $i (0 .. $N) {
		$p->[1] += $x->[$i];
		$p->[2] += $x->[$i]**2;
		$p->[3] += $x->[$i]*$y->[$i];
		$p->[4] += $y->[$i];
		$p->[5] += $y->[$i]**2;
		$p->[6] += $x->[$i]**3;
		$p->[7] += $x->[$i]*$y->[$i]**2;
		$p->[8] += $y->[$i]**3;
		$p->[9] += $x->[$i]**2*$y->[$i];
	}
	$N++;
	$a1 = 2 * ($p->[1]*$p->[1] - $N*$p->[2]);
	$b1 = 2 * ($p->[1]*$p->[4] - $N*$p->[3]);
	$a2 = $b1;
	$b2 = 2 * ($p->[4]*$p->[4] - $N*$p->[5]);
	$c1 = $p->[2]*$p->[1] - $N*$p->[6] + $p->[1]*$p->[5] - $N*$p->[7];
	$c2 = $p->[2]*$p->[4] - $N*$p->[8] + $p->[4]*$p->[5] - $N*$p->[9];
	$data->{xc} = ($c1*$b2-$c2*$b1)/($a1*$b2-$a2*$b1);
	$data->{yc} = ($a1*$c2-$a2*$c1)/($a1*$b2-$a2*$b1);
	$data->{R} = sqrt(($p->[2] - 2*$p->[1]*$data->{xc} + $N*$data->{xc}*$data->{xc} + $p->[5] - 2*$p->[4]*$data->{yc} + $N*$data->{yc}*$data->{yc})/$N);
}

sub MD_filter {
	my ($atoms, $cutmax, $cutmin) = @_;
	my ($cutmaxsq, $cutminsq);
	my ($i, $j, $n, $x1m, $x2m, $Xc, $Cx);
	my ($Cinv, $MDlensq, $tmp, $aList);

	$n = 0;
	$cutmaxsq = $cutmax**2;
	$cutminsq = $cutmin**2;
	@{ $aList } = sort numerically keys %{ $atoms };
	for $i (@{ $aList }) {
		$x1m += $atoms->{$i}{XCOORD};
		$x2m += $atoms->{$i}{ZCOORD};
		$n++;
	}
	$x1m /= $n;
	$x2m /= $n;
	$Xc = new Math::MatrixReal($n,2);
	$j = 1;
	for $i (@{ $aList }) {
		$Xc->assign($j,1,$atoms->{$i}{XCOORD}-$x1m);
		$Xc->assign($j,2,$atoms->{$i}{ZCOORD}-$x2m);
		$j++;
	}
	$Cx = (~$Xc*$Xc)/($n-1);
	$Cinv = $Cx->inverse();
	$j = 1;
	for $i (@{ $aList }) {
		$tmp = $Xc->row($j);
		$MDlensq = $tmp * $Cinv * ~$tmp; #Mahalanobis length squared
		$tmp = $MDlensq->element(1,1);
		if($tmp > $cutmaxsq || $tmp < $cutminsq) { #remove vapors
			delete $atoms->{$i};
		}
	}
}

sub calcContactAngle {
	my ($atoms, $box, $data, $sFlag, $cFlag) = @_;
	my ($i, $j, $L, $H, $B, $dz, $dx, $dtheta, $nrot);
	my ($norm_const, $pi, $h_rc, $h_rc2);
	my ($A1, $Nx, $Ny, $Nz, $XZ_PLANE, $rot_count);
	my ($MD_cutoffMAX, $MD_cutoffMIN, $angles); 
	my ($xy_count, $rz, $rx, $xdata, $ydata, $NxNz);
	my ($xi, $zi, $k1, $k2, $xp, $zp, $rsqr, $CoM);
	my ($ref_density, $threshold_MAX, $threshold_MIN);
	my ($cutoff, $y0, $x0, $x01, $x02, $arg1, $fi);
	
	$pi = 4*atan2(1,1);
	$dx = $dz = 1; # Angstroms
	$nrot = 18;
	$dtheta = $pi/$nrot;
	$h_rc = 2;
	$h_rc2 = $h_rc**2;
	$A1 = 1/$h_rc2;

	$CoM = CoM($atoms) if (!$cFlag);
	$box->{XCOORD}{len}=$box->{X}{len} if (!exists($box->{XCOORD}));
	$box->{YCOORD}{len}=$box->{Y}{len} if (!exists($box->{YCOORD}));
	$box->{ZCOORD}{len}=$box->{Z}{len} if (!exists($box->{ZCOORD}));
	$L = $box->{XCOORD}{len};
	$H = $box->{ZCOORD}{len};
	$B = $box->{YCOORD}{len};
	$norm_const = 2/($pi * $h_rc2 * $B);
	$Nx = int($L/$dx) + 1;
	$Nz = int($H/$dz) + 1;
	for $i (0 .. $Nx) {
		for $j (0 .. $Nz) {
			$XZ_PLANE->[$i][$j] = 0;
		}
	}
	$rot_count = 0;

	#1. remove initial noise
	$MD_cutoffMAX = 12;
	$MD_cutoffMIN = 0;
	&MD_filter($atoms, $MD_cutoffMAX, $MD_cutoffMIN) if (! $cFlag); #filter by Mahalonbis distance
	
	#2. smear molecules on grid
	@{ $angles } = (0,0,$dtheta);
	for $i (1 .. $nrot) {
		for $j (keys %{ $atoms }) {
			$atoms->{$j}{XCOORD} -= $CoM->{XCOORD};
			$atoms->{$j}{YCOORD} -= $CoM->{YCOORD};
			$atoms->{$j}{ZCOORD} -= $CoM->{ZCOORD};
		}
		&Rotate($atoms, $angles, 2) if ($i>1);
		for $j (keys %{ $atoms }) {
			$atoms->{$j}{XCOORD} += $CoM->{XCOORD};
			$atoms->{$j}{YCOORD} += $CoM->{YCOORD};
			$atoms->{$j}{ZCOORD} += $CoM->{ZCOORD};
			$xi = $atoms->{$j}{XCOORD};
			$zi = $atoms->{$j}{ZCOORD};
			#Hardy weight function
			for $k1 (0 .. $Nx) {
				$xp = $k1*$dx;
				$rx = abs($xi - $xp);
				next if ($rx > $h_rc);
				for $k2 (0 .. $Nz) {
					$zp = $k2*$dz;
					$rz = abs($zi - $zp);
					next if ($rz > $h_rc);
					$rsqr = $rx**2+$rz**2;
					if($rsqr < $h_rc2) {
						$XZ_PLANE->[$k1][$k2] += (1-$rsqr*$A1)*$norm_const;
						$cumXZ_PLANE->[$k1][$k2] += (1-$rsqr*$A1)*$norm_const;
					}
				}
			}
		}
		$rot_count++;
	}
	#normalize
	$ref_density = 0;
	$XZ_PLANE = $cumXZ_PLANE if ($cFlag);
	for $i (0 .. $Nx) {
		for $j (0 .. $Nz) {
			$XZ_PLANE->[$i][$j] /= $rot_count;
			printf OUTFILE "%10d %10d %10.5f\n",$i,$j,$XZ_PLANE->[$i][$j];
			$ref_density = $XZ_PLANE->[$i][$j] if ($XZ_PLANE->[$i][$j] > $ref_density);
		}
	}
	&plot_grid_data($XZ_PLANE,$Nx,$Nz,"${savePrefix}.grid_dens.dat") if ($sFlag);

	#3. select only probable droplet contour pts
	$threshold_MAX = 0.2 * $ref_density;
	$threshold_MIN = 0.1 * $ref_density;
	$NxNz = $Nx*$Nz;
	$xy_count=0;
	for $i (0 .. $Nx) {
		for $j (0 .. $Nz) {
			if($XZ_PLANE->[$i][$j] > $threshold_MIN && $XZ_PLANE->[$i][$j] < $threshold_MAX) {
				$xdata->[$xy_count] = $i;
				$ydata->[$xy_count] = $j;
				$y0 = $j if (!defined($y0) or $j < $y0);
				$xy_count++;
			}
		}
	}
	&plot_surface($xdata,$ydata,"${savePrefix}.grid_surface.dat") if ($sFlag);
	$y0 += $dz*2;

	#4. remove monolayer
	&fit_Landau($xdata,$ydata,\%{ $data }); #initial fit to circle for cleaning
	$i = 0;
	while ($i <= $#{ $xdata } ) {
		if ($ydata->[$i] <= $y0) {
			splice @{ $ydata }, $i, 1;
			splice @{ $xdata }, $i, 1;
		} else {
			$i++;
		}
	}
	&plot_surface($xdata,$ydata,"${savePrefix}.grid_surface_nomonolayer.dat") if ($sFlag);
	$xy_count = $i;

	#5. fit with circle
	&fit_Landau($xdata,$ydata,\%{ $data }); #final
	&plot_fit($data,"${savePrefix}") if ($sFlag);

	#6. Calculate contact angle
	$arg1 = $data->{R}**2-($y0-$data->{yc})**2;
	if ($data->{yc} < $y0) { #hydrophilic
		$x0 = $data->{xc} - sqrt($arg1);
	} else {
		$x01 = $data->{xc} + sqrt($arg1);
		$x02 = $data->{xc} - sqrt($arg1);
		$x0 = $x02;
		$x0 = $x01 if ($x01 < $data->{xc});
	}
	$data->{gradient} = ($data->{xc} - $x0)/($y0 - $data->{yc});
	$data->{angle} = atan2($data->{gradient},1)*180/$pi;
	$data->{angle} += 180 if ($data->{angle} < 0);
	#statistics
	&getStats($xdata,$ydata,$data);
}

sub getContactAngle {
	my ($atoms, $box, $frameNum, $fileHandle) = @_;
	my ($tot, $geomVal, $i, $j, $data, $CENTER, $MOLECULE);

	if ($trjType == 2) { #LAMMPS
		$frameNum = $atoms->{TIMESTEP}[0];
		$box = ConvertLammpsBox($atoms->{"BOX BOUNDS"});
		$tot = $ATOMS->{"NUMBER OF ATOMS"}[0];
		$atoms = $atoms->{ATOMS};
		if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
			UnwrapAtoms($atoms,  $box, $LAMMPSOPTS->{scaled});
		}
	} else {
		$box = ConvertAmberBox(\%{ $box });
		$tot = scalar(keys %{ $atoms });
	}

	$MOLECULE = GetAtmData($atoms, $atmList);
	$CENTER = CoM($MOLECULE);
	for $j ("XCOORD","YCOORD","ZCOORD") {
		$box->{$j}{hi} = $box->{$j}{len};
		$box->{$j}{lo} = 0;
		$box->{$j}{CENTER} = $box->{$j}{hi}/2;
		for $i (keys %{ $atoms }) {
			$atoms->{$i}{$j} += ($box->{$j}{CENTER} - $CENTER->{$j});
		}
	}
	for $i (keys %{ $MOLS }) {
		$MOLECULE = GetAtmData($atoms, $MOLS->{$i}{MEMBERS});
		$CENTER = CoM($MOLECULE);
		&ImageAtoms($MOLECULE, $CENTER, $box);
	}
	for $i (keys %{ $atoms }) {
		delete $atoms->{$i} if (!defined($atmList->{$i}));
	}
	&calcContactAngle($atoms,$box,\%{ $data },0, 0);
	printf $OUTFILE "%10d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$frameNum,$data->{angle},$data->{sigma},$data->{R2},$data->{R},$data->{xc},$data->{yc};
}

sub init {
	my (%OPTS, $atm, $tSel, $list, $i);

	getopt('btasow',\%OPTS);
	for ("b", "a") {
		&usage if (! exists($OPTS{$_}));
	}
	print "Initializing...";
	($bgfFile, $atm, $trjFile, $tSel, $saveFile, $trjType) = 
		($OPTS{b}, $OPTS{a}, $OPTS{t}, $OPTS{s}, $OPTS{w}, $OPTS{o});
	FileTester($bgfFile);
	$selection = BuildAtomSelectionString($atm);

	$savePrefix = basename($bgfFile);
	if (defined($trjFile) and -e $trjFile and -r $trjFile and -T $trjFile) {
		if (! defined($trjType)) {
			if ($trjFile =~ /\.lammps/) {
				$trjType = "lammps";
			} else {
				$trjType = "amber";
			}
		}

		if (lc($trjType) ne "lammps") {
			$trjType = 1;
			$getSnapshots = \&ParseAmberTrj;
			$getByteOffset = \&GetAmberByteOffset;
		} else {
			$trjType = 2;
			$getSnapshots = \&ParseLAMMPSTrj;
			$getByteOffset = \&GetLammpsByteOffset;
		}
		if (! defined($tSel)) {
			$tSel = "*";
		}
		$list = TrjSelections($tSel);
		for $i (keys %{ $list }) {
			$SELECT->{$i} = $list->{$i};
		}
		die "ERROR: No valid frames selected with selection $tSel!\n"
			if (! keys %{ $SELECT } and $tSel ne "*");
		if (! defined($saveFile)) {
		   $saveFile = basename($bgfFile);
		   $saveFile =~ s/\.\w+$//;
		   $saveFile .= ".ca.dat";
		}
		$savePrefix = $saveFile;
	}
	$savePrefix =~ s/\.\w+$//;
	print "Done\n";
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -a atom_selection -w (save_name) -t (traj_file) -o (traj_file_type) -s (traj_sel)
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save (optional)
  atom_selection:
	any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
  traj_file: name of trajectory file (LAMMPS or AMBER)
  traj_type: type of trajectory file (LAMMPS or AMBER)
  traj_sel: trajectory file selection: The number of frames to use. Can be a single integer, or several integers in quotes
	To specify a range specify it as :Ita-b:c, which will take frames a - b, every c. Specify multiple 
	ranges or a combination of ranges and single frames by enclosing them in quotes. \"*\" for all frames
DATA
die "\n";

}
sub sum {
	my ($arrayref) = @_;
	my $result;
	foreach(@$arrayref) { $result+= $_; }
	return $result;
}

sub Mean {
	my ($arrayref) = @_;
	my $result;
	foreach (@$arrayref) { $result += $_ }
	return $result / @$arrayref;
}

sub variance {
	my ($arrayref)= @_;
	my ($mean) = Mean($arrayref);
	my @a = map ( ($_ - $mean)**2, @$arrayref);
	return sum(\@a) / scalar @$arrayref;
}

sub getStats {
	my ($x,$y,$data) = @_;
	my ($R,$i,$xval,$yval,$yc,$xc, $val2, $val, $n, $errorsq,$mdsq,$ymean,$tmp);

	$R = $data->{R};
	$xc = $data->{xc};
	$yc = $data->{yc};
	$n = $#{ $x };
	for $i (@{ $y }) { $ymean += $i};
	$ymean /= ($n+1);
	for $i (0 .. $n) {
		$xval = $x->[$i];
		$tmp = $R**2-($xval-$xc)**2;
		next if ($tmp < 0);
		$yval = sqrt($tmp) + $yc;
		$errorsq = ($y->[$i]-$yval)**2;
		$yval -= $yc*2;
		$tmp = ($y->[$i]-$yval)**2;
		$errorsq = $tmp if ($tmp < $errorsq);
		$val2 += $errorsq;
		$val += sqrt($errorsq);
		$mdsq += ($y->[$i]-$ymean)**2;
	}
	$n++;
	$data->{R2} = 1-$val2/$mdsq;
	$val /= $n;
	$val2 /= $n;
	$data->{sigma} = sqrt($val2 - $val**2);
}

sub plot_fit {
	my ($data, $fprefix) = @_;
	my ($fname, $gd_name, $gs_name, $fs_name);

	$fname = "${fprefix}.plot_fit.plt";
	$gd_name = "${fprefix}.grid_dens.dat";
	$gs_name = "${fprefix}.grid_surface_nomonolayer.dat";
	$fs_name = "${fprefix}.fit_circle.dat";
	open OUTFILE, "> $fname" or die "ERROR: Cannot write to ${fname}: $!\n";
	print OUTFILE <<DATA;
set trange [0:2*pi]
set parametric
xc = $data->{xc}
yc = $data->{yc}
R = $data->{R}
fx(t) = R*cos(t)+xc
fy(t) = R*sin(t)+yc
set samples 500
set table '${fs_name}'
pl fx(t),fy(t)
unset table
unset parametric
set view map
set key top out horiz center
spl '${gd_name}' u 1:2:3 w pm3d not, '${fs_name}' u 1:2:(0.0) w l lt -1 lc rgb "red" dt 1 lw 1.5 t sprintf("fit: R = %f xc = %f yc = %f",R,xc,yc), '${gs_name}' u 1:2:(0.0) w p lt 2 pt 7 t "data", "<echo '$data->{xc} $data->{yc}'" u 1:2:(0.0) w p pt 9 ps 2 t "center"
DATA
	close OUTFILE;

}

sub plot_surface {
	my ($xdata, $ydata, $fname) = @_;
	my ($i);

	open OUTFILE, "> $fname" or die "ERROR: Cannot write to $fname: $!\n";
	$i = 0;
	while ($i <= $#{ $xdata } ) {
		printf OUTFILE "%10d %10d\n", $xdata->[$i], $ydata->[$i];
		$i++;
	}
	close OUTFILE;
}

sub plot_grid_data {
	my ($XZ_PLANE,$Nx,$Nz,$fname) = @_;
	my ($i, $j);

	open OUTFILE, "> $fname" or die "ERROR: Cannot write to $fname: $!\n";
	for $i (0 .. $Nx) {
		for $j (0 .. $Nz) {
			printf OUTFILE "%10d %10d %10.5f\n",$i,$j,$XZ_PLANE->[$i][$j];
		}
		printf OUTFILE "\n";
	}
	close OUTFILE;
}

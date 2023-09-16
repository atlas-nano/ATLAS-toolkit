#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use FileFormats qw(addHeader createBGF createMSI createMOL2 insertHeaderRemark 
				addBoxToHeader ParseStructFile AddMass createHeaders);	
use General qw(ReadFFs FileTester LoadFFs GetBondLength CoM Rotate 
			LoadElements AddElementField GetFileTypeStr);
use GRID;
use File::Basename qw(basename);
use BOX qw(InitBox Cart2Frac);
use Storable qw(dclone);

my ($DATA, $saveFile, $ffList, $mol2File, $msiFile, $createAuxFiles, $optFlag, $overlapFlag);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $FF, $FFILES, $pStr, $tdens, $cellL, $ELEMENTS);

$|++;
&init;
if(!defined($cellL)) {
	($FFILES, undef) = ReadFFs($ffList);
	$FF = LoadFFs($FFILES);
}
&readStructs($DATA, $FF, $cellL);
$pStr = "Creating system with $DATA->{TOT} monomers...";
($ATOMS, $BONDS, $BOX, $tdens) = &createSys($DATA, $pStr, $cellL, $overlapFlag);
print "Done\nUpdating fields...";
$HEADERS = createHeaders(undef, $saveFile);
&insertHeaderRemark($HEADERS, "REMARK Created by $0 at " . localtime);
&insertHeaderRemark($HEADERS, "REMARK amorphous_builder of $DATA->{N} units to density $tdens");
&addBoxToHeader($HEADERS, $BOX);
$ATOMS = &optimize($ATOMS, $BONDS, $HEADERS, $ffList) if ($optFlag);
if($createAuxFiles) {
	print "Done\nCreating MOL2 file ${mol2File}...";
	&createMOL2($ATOMS, $BONDS, $mol2File,1);
	print "Done\nCreating MSI file ${msiFile}...";
	&AddElementField($ATOMS, $ELEMENTS, $FF);
	&createMSI($ATOMS, $BONDS, $HEADERS, $BOX, $msiFile);
}
print "Done\nCreating BGF file ${saveFile}...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub optimize {
	my ($atoms, $bonds, $headers, $ffs) = @_;
	my ($execCmd, $prefix, $tot, $outStr, $inStr, $i, $tmp, $s1, $s2);

	print "Optimizing structure...Creating LAMMPS files\r";
	&addHeader($atoms, $headers);
	&createBGF($atoms, $bonds, "_tmp.bgf");
	$prefix = "_tmp";
	$tmp = "";
	$execCmd = "$Bin/createLammpsInput.pl -b ${prefix}.bgf -f \"$ffs\" -s $prefix -t 'min'";
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
	open LMPFILE, ">> in.${prefix}_singlepoint" or die "ERROR: Cannot write to in.${prefix}";
	print LMPFILE <<DATA;
timestep             0.5	
thermo_style         multi
thermo_modify        flush yes

#expand system to ~1/2 density	
change_box           all x scale 1.25 y scale 1.25 z scale 1.25 remap

#minimization
thermo               10
dump                 1 all custom 25 \${sname}.min.lammps id type xu yu zu vx vy vz
min_style            sd
minimize             1.0e-4 1.0e-4 500 5000
min_style            cg
minimize             1.0e-4 1.0e-4 500 5000
#now                 minimize the entire system
minimize             1.0e-4 1.0e-4 500 5000
undump               1	

thermo               100
dump                 1 all custom 100 \${sname}.nvt.lammps id type xu yu zu vx vy vz
reset_timestep       0
#high tempeature annealing
variable             htemp equal \${rtemp}*10
velocity             all create \${htemp} 12345678 dist gaussian	
fix                  2 all nvt temp \${htemp} \${htemp} 100.0
run                  10000 # run for 5 ns
#compression to original density (high temp)
fix                  1 all deform 100 x scale 0.8 y scale 0.8 z scale 0.8
run                  20000
unfix                1
#high -> regular temp
fix                  2 all nvt temp \${htemp} \${rtemp} 100.0
run                  10000
#regular temp
fix                  2 all nvt temp \${rtemp} \${rtemp} 100.0
run                  10000
unfix                2
undump               1	
DATA

	print "Optimizing structure...                                        \r";
	print "Optimizing structure...Running LAMMPS\r";
	system("mv in.${prefix}_singlepoint in.${prefix}");
	$execCmd = "module purge; module load cpu/0.17.3b  gcc/10.2.0/npcyll4 cmake/3.21.4/teqow32 openmpi/4.1.3/oq3qvsv fftw/3.3.10/bmvqsbr python/3.8.12/7zdjza7 intel-mkl/2020.4.304/ghfk3mu gsl/2.7/wtlsmyy gnuplot/5.4.2/mfinpvw sdsc slurm &> /dev/null 2>&1; mpirun -n 1 -mca btl vader,self $Bin/../codes/bin/lmp_expanse -in in.${prefix} -screen none -log ${prefix}.lammps.log -var rtemp 298";
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));

	#update coordinates
	print "Optimizing structure...                                        \r";
	print "Optimizing structure...updating coordinates\r";
	$execCmd = "$Bin/convertLammpsTrj.pl -b ${prefix}.bgf -o bgf -s ${prefix}.bgf -t last -l ${prefix}.nvt.lammps -m 1";
	die "ERROR: Cannot execute \"$execCmd\"\n" if(system("${execCmd} > /dev/null"));
	($atoms, undef) = ParseStructFile("${prefix}.bgf",0);

	system("rm -fr in.${prefix} in.${prefix}_singlepoint data.${prefix} ${prefix}.bgf"); 
	system("rm -fr ${prefix}.lammps.slurm ${prefix}.lammps.log log.cite");
	system("rm -fr ${prefix}.min.lammps _byte_offset_${prefix}.heat.lammpstrj ${prefix}.heat.lammpstrj ${prefix}.nvt.lammps");
	print "Optimizing structure...                                        \r";
	print "Optimizing structure...";

	return $atoms; 
}

sub createSys {
	my ($data, $pStr, $cell, $ov) = @_;
	my ($t, $o, $c, $i, $j, $k, $p, $l, $x); 
	my ($blen, $box, $atoms, $bonds);
	my ($tmass, $CoM, $vol, $tmp); 
	my ($gs, $gm, $grid, $bmax);
	my ($cMol, $angles);

	print "${pStr}\r";
	#initialze random number generator
	srand(time() ^($$ + ($$ <<15)));

	$data->{N}--;

	$vol = 0.0; #cell volume
	for $i (0 .. $data->{N}) {
		$tmass = 0.0; #total mass
		#determine center of mass
		$CoM = CoM($data->{$i}{ATOMS});

		for $j (keys %{ $data->{$i}{ATOMS} }) {
			if(defined($data->{$i}{ATOMS}{$j}{MASS})) {
				$tmass += $data->{$i}{ATOMS}{$j}{MASS} * $data->{$i}{N};
			} else {
				$tmass++;
			}
			for $k ("XCOORD", "YCOORD", "ZCOORD") {
				$data->{$i}{ATOMS}{$j}{$k} -= $CoM->{$k}; #place molecule center at origin
			}
		}
		$vol += $tmass / 0.6023 / $data->{$i}{DENS} if (! defined($cell));
	}

	#determine box size based on number of monomers n and density d
	if (! defined($cell)) {
		$blen->{XCOORD} = $blen->{YCOORD} = $blen->{ZCOORD} = $vol ** (1/3);
		$tdens = $tmass /0.6023/ $vol;
		for $i ("X", "Y", "Z") {
			$box->{$i}{hi} = $box->{$i}{len} = $blen->{"${i}COORD"};
			$box->{$i}{angle} = 90; $box->{$i}{lo} = 0;
		}
	} else {
		@{ $tmp } = split /\s+/, $cell;
		$box->{X}{hi} = $box->{X}{len} = $tmp->[0];
		$box->{Y}{hi} = $box->{Y}{len} = $tmp->[1];
		$box->{Z}{hi} = $box->{Z}{len} = $tmp->[2];
		$box->{X}{lo} = $box->{Y}{lo} = $box->{Z}{lo} = 0;
		$box->{X}{angle} = $box->{Y}{angle} = $box->{Z}{angle} = 90;
		($blen->{XCOORD},$blen->{YCOORD},$blen->{ZCOORD}) = @{ $tmp };
		$vol = $tmp->[0]*$tmp->[1]*$tmp->[2];
		$tdens = 1;
	}
	&InitBox($box, $data->{1}{ATOMS});

	#initialize grid
	$gs = 2.5; #grid size
	for $i ("X", "Y", "Z") {
		$gm->{$i} = int($blen->{"${i}COORD"}/$gs)+1;
	}
	$grid = new GRID ($gm->{X}, $gm->{Y}, $gm->{Z}, $gs);

	$o = 0;
	for $i (0 .. $data->{N}) {
		&Cart2Frac($data->{$i}{ATOMS},$box);

		#counters
		$t = scalar(keys %{ $data->{$i}{ATOMS} });

		#now insert molecule into box 1 by 1 
		$c = 0;
		while($c < $data->{$i}{N}) {
			$cMol = dclone($data->{$i}{ATOMS});
			#radomely rotate atoms
			@{ $angles } = (int(rand() * 180), int(rand() * 180), int(rand() * 180));
			for $j (0 .. $#{ $angles }) {
				&Rotate($cMol, $angles, $j);
			}
			$p = 0;
			for $l (1 .. 100) {
				#randomely place molecule in cell
				for $j ("XCOORD", "YCOORD", "ZCOORD") {
					$CoM->{$j} = int(rand() * $blen->{$j});
				}
				next if($ov == 0 and $c > 0 and $grid->molOverlap($cMol, $CoM)); # check overlap
				for $j (keys %{ $cMol }) {
					%{ $atoms->{$j+$o} } = %{ $cMol->{$j} };
					$atoms->{$j+$o}{RESNUM}++;
					for $k ("X", "Y", "Z") {
						$atoms->{$j+$o}{"${k}COORD"} += $CoM->{"${k}COORD"}; #update atom coordinate
					}
					$grid->store($atoms->{$j+$o}); #place atom on grid
					for $k (0 .. $#{ $data->{$i}{BONDS}{$j} } ) {
						$bonds->{$j+$o}[$k] = $data->{$i}{BONDS}{$j}[$k] + $o;
					}
				}
				$o += $t; #offset counter
				$p = 1; #placed flag
				last;
			}
			if ($p) {
				$c++; #molecule counter
				$x = 0;
				print "${pStr} ${c}/$data->{N}\r";
			} else {
				$x++;
			}
			last if ($x==100); #give up after trying to place 100 times
		}

		warn "WARNING: Cannot place molecule after 100 tries... Placed $c ... Giving up!\n"
			if($x == 100);
	}

	$data->{N}++;
	print "${pStr} ${c}/$data->{N}...";	
	return ($atoms, $bonds, $box, $tdens);
}
		
sub readStructs {
	my ($data, $ff, $cell) = @_;
	my ($i, $fle);

	for $i (0 .. $data->{N}-1) {
		$fle = $data->{$i}{FILE};
		print "Getting atom information from ${fle}...";
		($data->{$i}{ATOMS}, $data->{$i}{BONDS}, $data->{$i}{HEADERS}) = ParseStructFile($fle,1);
		print "Done\n";
		&AddMass($data->{$i}{ATOMS}, $ff) if (! defined($cell));
	}
}

sub init {
	my (%OPTS, $itype, $mList, $dList, $nList, $tmp, $i, $rec);

	getopt('mndsfcxov',\%OPTS);
	#($molFile, $ffList, $num, $dens, $saveFile, $itype) = ($OPTS{m},$OPTS{f},$OPTS{n},$OPTS{d},$OPTS{s},$OPTS{t});
	($mList,  $ffList, $nList,  $dList,  $saveFile, $itype,  $cellL,  $createAuxFiles, $optFlag, $overlapFlag) = 
	($OPTS{m},$OPTS{f},$OPTS{n},$OPTS{d},$OPTS{s},  $OPTS{t},$OPTS{c},$OPTS{x},        $OPTS{o}, $OPTS{v});

	for ($mList, $nList) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	$mList =~ s/^\s*//;
	$mList =~ s/\s*$//;
	@{ $tmp } = split /\s+/, $mList;
	for $i (0 .. $#{ $tmp }) {
		FileTester($tmp->[$i]);
		$DATA->{$i}->{FILE} = $tmp->[$i];
		$DATA->{$i}{DENS} = 1;
	}
	$DATA->{N} = scalar(@{ $tmp });
	$nList =~ s/^\s*//;
	$nList =~ s/\s*$//;
	@{ $tmp } = split /\s+/, $nList;
	for $i (0 .. $#{ $tmp }) {
		die "ERROR: Expected number for num_monomers. Got '$tmp->[$i]'!\n"
			if($tmp->[$i] !~ /^\d+$/);
		$DATA->{$i}{N} = $tmp->[$i];
		$DATA->{TOT}+= $tmp->[$i];
	}
	$DATA->{Nmon} = scalar(@{ $tmp });
	die "ERROR: Number of monomers $DATA->{Nmon} < monomer files $DATA->{N}\n"
		if ($DATA->{Nmon} < $DATA->{N});

	if (defined($dList)) {
		@{ $tmp } = split /\s+/, $dList;
		for $i (0 .. $#{$tmp}) {
			die "ERROR: Expected number/decimal for density. Got '$tmp->[$i]'\n"
				if ($tmp->[$i] !~ /^\d+\.?\d*$/);
			$DATA->{$i}{DENS} = $tmp->[$i];
		}
	}	

	if(!defined($saveFile)) {
		$saveFile = basename($DATA->{0}{FILE});
		$saveFile =~ s/\.\w+$//;
		$saveFile = "${saveFile}.bgf";
	}
	$mol2File = $saveFile;
	$mol2File =~ s/\.\w+$//;
	$mol2File .= ".mol2";
	$msiFile = $saveFile;
	$msiFile =~ s/\.\w+$//;
	$msiFile .= ".msi";
	$optFlag = 0 if (! defined($optFlag) or $optFlag !~ /1|yes/i);
	$optFlag = 1 if ($optFlag =~ /1|yes/i);
	$overlapFlag = 0 if (! defined($overlapFlag) or $overlapFlag !~ /1|yes/i);
	$overlapFlag = 1 if ($overlapFlag =~ /1|yes/i);

	if(defined($cellL)) {
		if($cellL =~ /(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)/) {
			$cellL = "$1 $2 $3";
		} else {
			undef $cellL;
		}
	} elsif (! defined($ffList)) {
		die "ERROR: Need to specify the forcefield if not specifying the cell length\n";
	}
	$ELEMENTS = &LoadElements;
	$createAuxFiles=0 if (! defined($createAuxFiles) or $createAuxFiles =~ /0|no/i);
	$createAuxFiles=1 if ($createAuxFiles =~ /1|yes/i);
}

sub usage {
	my ($fTypeStr) = GetFileTypeStr;
	print STDOUT <<DATA;
usage: $0 -m "monomer(s)_file"  -n "num_monomer(s)" -f force_field(s) -d "(den(s)=1.0)" -c ("cell_lengths") 
		-o (optimize) -s (save_name) -t (mol_file_type) -x (create_aux_files) -v (allow_overlap=no)
Arguments:
  	molecule_file: molecule structure file
$fTypeStr	
	num_monomers: number of monomers in final structure
	force_file: list of forcefields
	density: final density of simulation cell
	cell_length: a b c cell lengths
	optmize: Perform a LAMMPS minimization and short NVT (298K) dynamics of the final structure.
		Expected 1|yes or 0|no (default)
	save_name: name of file to save
	mol_file_type: type of input molecule.
	create_aux_files: Create accompanying MSI and MOL2 files. Default = 0
	allow_overlap: Specify whether to place the monomers do that atoms can be very close to each other. Default = 0
DATA
die "\n";

}

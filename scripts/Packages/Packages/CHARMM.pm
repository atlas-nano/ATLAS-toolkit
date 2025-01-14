package CHARMM;

require Exporter;
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use General qw(Permutate PrintProgress GetElementSymbolByMass LoadElements LoadConverter);
use CERIUS2 qw(getLammpsOpts findDuplicate);
use constant q2e => 18.2223;
use Math::Trig qw(pi);

use strict;

our (@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(parsePSFFile parseCharmmCoordFile parseCharmmFF);;
$VERSION = "1.00";

sub numerically { ($a<=>$b); }


sub parseCharmmCoordFile {
	my ($atoms, $cFile) = @_;
	my ($pattern, $count, $totAtoms);

	$totAtoms = scalar(keys %{ $atoms });
	$pattern = '^\s+(\d+)\s+(\d+)\s+(\w+)\s+(\w+)\s+(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s+(\w+)\s+(\d+)\s+(\-?\d+\.\d+)';
	open COORDFILE, $cFile or die "ERROR: Cannot open CHARMM coordinate file $cFile: $!\n";
	while (<COORDFILE>) {
		chomp;
		if ($_ =~ $pattern) {
			$count++;
			die "ERROR: Unequal number of atoms in psfFile and coordFile" if ($1 > $totAtoms || $count > $totAtoms);
			$atoms->{$1}{RESNUM} = $2;
			$atoms->{$1}{RESNAME} = $3;
			$atoms->{$1}{ATMNAME} = $4;
			$atoms->{$1}{XCOORD}  = $5;
			$atoms->{$1}{YCOORD}  = $6;
			$atoms->{$1}{ZCOORD}  = $7;
		}
	}
	close COORDFILE;
	die "ERROR: No valid data found while parsing $cFile\n"
		if ($totAtoms != $count);
}

sub parsePSFFile {
	my ($psfFile) = @_;
	my ($atoms, $bonds, $headers, $atm_pattern, $bnd_pattern, $rec, $count); 
	my ($i,$ j, $totAtoms, $recAtoms, $recBonds, @tmp);

	$totAtoms = scalar(keys %{ $atoms });
	$recAtoms = $recBonds = 0;
	$atm_pattern = '^\s+(\d+)\s+(\w+)\s+(\d+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\-?\d+\.\d+E?\-?\d*)\s+(\d+\.\d+)\s+(\d+)\s*(\-?\d*\.?\d*)\s*(\-?\d*\.?\d*E?\-?\d*)';
	$bnd_pattern = '^\s+(\d+)\s+(\d+\s*.*)$';
	open PSFFILE, $psfFile or die "ERROR: Cannot open CHARMM PSF file $psfFile: $!\n";
	while (<PSFFILE>) {
		chomp;
		if ($_ =~ /\!NATOM/) {
			$recAtoms = 1;
		} elsif ($_ =~ /(\d+)\s+\!NBOND/) {
			$recAtoms = 0;
			$recBonds = 1;
		} elsif ($_ =~ /^\s*$/) {
			$recAtoms = $recBonds = 0;
		} elsif ($recAtoms && $_ =~ $atm_pattern) {
			$rec = (
						{
							RESNAME => $2,
							RESNUM  => $3,
							ATMNAME => $5,
							FFTYPE  => $6,
							CHARGE  => $7,
							MASS    => $8,
							INDEX   => $1,
						}
			);
			%{ $atoms->{$1} } = %{ $rec };
		} elsif ($recBonds && $_ =~ $bnd_pattern) {
			@tmp = split /\s+/,$2;
			$i = $1; $j =  $tmp[0];
			push @{ $bonds->{$i} }, $j;
			push @{ $bonds->{$j} }, $i;
			$count = 0;
			while ($count < $#tmp) {
				$i = $tmp[++$count];
				$j = $tmp[++$count];
				push @{ $bonds->{$i} }, $j;
				push @{ $bonds->{$j} }, $i;
			}
		}
	}
	close PSFFILE;
	die "ERROR: No valid data found while parsing $psfFile\n"
		if (! defined($atoms) or ! defined($bonds));
	return ($atoms, $bonds);
}

sub parseCharmmFF {
	my ($ff, $alter, $parms) = @_;
	my ($i, $j, $sstr, $header, $ty_counter, $ele, $CNV, $rec, $aM); 
	my ($curr, $tmp, $t1, $t2, $inv, $cmap, $valid, $ty_index, $resN);

	$CNV = &General::LoadConverter();
	$ele = &General::LoadElements();
	$header = "";

	if(! keys %{ $parms } or !exists($parms->{PARMS})) {
		$parms->{PARMS} = (
							{
								'EXO_CYCLIC'							=> 1.00000,
								'HAS_HBONDS'							=> 0,
								'QEq'									=> 0,
								'SW'									=> 0,
								'USE_HBOND'								=> 0,
								'coul_accuracy'							=> 0.00001,
								'cut_coul'								=> 10.00000,
								'cut_vdw'								=> 10.00000,
								'dielectric'							=> 1.00000000000000000,
								'hbond_angle'							=> 90,
								'hbond_distance'						=> 2.5,
								'mix_rule'								=> 'arithmetic',
								'same_scale'							=> 1,
								'scale_cou'								=> 0.00000,
								'scale_cou_12'							=> 0,
								'scale_cou_13'							=> 0,
								'scale_cou_14'							=> 0.50000,
								'scale_torsions'						=> 0,
								'scale_torsions_by_n_defined_torsions'	=> 0,
								'scale_vdw'								=> 0.00000,
								'scale_vdw_12'							=> 0,
								'scale_vdw_13'							=> 0,
								'scale_vdw_14'							=> 0.50000,
								'single_inversion'						=> 1,
							}
		);
	}

    open CHARMMPAR, $ff->{FF} || die "ERROR: Cannot open CHARMM force field $ff: $!\n";
    while (<CHARMMPAR>) {
		chomp;
		next if ($_ !~ /^(\w+)/ and $header ne "CMAP");
		$sstr = $1;
		next if (! defined($sstr));
		if ($sstr =~ /^ATOMS/) {
			$header = "ATOMTYPES";
		} elsif ($sstr =~ /^(BONDS|ANGLES|DIHEDRALS|IMPROPER|NONBONDED|CMAP)/ or $sstr =~ /^RESI/i) {
			undef $ty_counter; 
			undef $curr;
			$ty_index = 0;
			$header = $1;
			$header = "TORSIONS" if ($_ =~ /DIHEDRALS/);
			$header = "INVERSIONS" if ($_ =~ /IMPROPER/);
			$header = "VDW" if ($_ =~ /NONBONDED/);
			if ($_ =~ /RESI\s+(\S+)/i) {
				$header = "POLARIZATION";
				$resN = $1; #residue name
			}
		#} elsif ($header eq "ATOMTYPES" and $_ =~ /^MASS\s+\-?(\d+)\s+(\w+)\s+(\d+\.\d+)\s+/) { #atomsa
		} elsif ($_ =~ /^MASS\s+\-?(\d+)\s+(\w+)\s+(\d+\.\d+)\s+/) { #atoms
			$header = "ATOMTYPES";
			$aM = $3;
			$aM = 1e-16 if ($aM == 0);
			$parms->{$header}{$2} = (
										{
											"TYPEID"     => $1,
											"LABEL"      => $2,
											"MASS"       => $aM,
											"USE_CHARGE" => 0,
											"USED"       => 0,
											"CHARGE_SET" => 0,
											"CHARGE"     => 0,
											"ATOM"       => General::GetElementSymbolByMass($ele, $3),
										}
			);
		} elsif ($header eq "POLARIZATION" and $_ =~ /^ATOM\s+(\w+)\s+(\w+)\s+(\-?\d+\.\d+).*\s+ALPHA\s+\-?(\d+\.\d+)\s*(THOLE)?\s*(\-?\d*\.?\d*)/) {
			$parms->{$header}{$resN}{$2} = (
											{
												m     => 0.4,
												k     => 1000, #default to 1000 kcal/mol for CHARMM
												q     => -sqrt($4*3.0114), #q_D(eV) = −(αK_D)^(1/2)
												alpha => $4,
												totQ  => $3,
											}
										);
			$parms->{$header}{$resN}{$2}{thole} = 2.6;
			$parms->{$header}{$resN}{$2}{thole} = $6*2 if (defined($6) and $6); #note that the thole in the CHARMM drude files is 1/2 the value in LAMMPS
		} elsif ($header eq "BONDS" and $_ =~ /^(\w+)\s+(\w+)\s+(\d+\.\d+)\s+(\d+\.\d+)/) { #bonds
			if(exists($parms->{$header}{$2}) and exists($parms->{$header}{$2}{$1})) {
				$curr = \%{ $parms->{$header}{$2}{$1} };
			} else {
				$curr = \%{ $parms->{$header}{$1}{$2} };
			}
			$ty_counter = findDuplicate($curr,"HARMONIC") if(defined($curr));
			$ty_counter = 1 if (!defined($ty_counter) or ! $ty_counter);

			$curr->{$ty_counter} = (
									{
										"INDEX"    => ++$ty_index,
										"TYPE"     => "HARMONIC",
										"Lammps"   => getLammpsOpts("HARMONIC","bond", $CNV),
										"VALS"     => [$3,$4],
										"USED"     => 0,
										"KEY"      => "$1 $2 ",
										"IGNORE"   => 0,
									}
			);
		} elsif ($header eq "ANGLES" and $_ =~ /^(\w+)\s+(\w+)\s+(\w+)\s+(\d+\.\d+\s+\d+\.\d+\s*[\-?0-9\. ]*)/) { #angles
			if(exists($parms->{$header}{$3}) and exists($parms->{$header}{$3}{$2}) and exists($parms->{$header}{$3}{$2}{$1})) {
				$curr = \%{ $parms->{$header}{$3}{$2}{$2} };
			} else {
				$curr = \%{ $parms->{$header}{$1}{$2}{$3} };
			}
			$ty_counter = findDuplicate($curr,"CHARMM") if(defined($curr));
			$ty_counter = 1 if (!defined($ty_counter) or ! $ty_counter);
			@{ $tmp } = split /\s+/,$4;
			if ($#{ $tmp } == 1) {
				push @{ $tmp }, 0;
				push @{ $tmp }, 0;
			}

			$curr->{$ty_counter++} = (
										{
											"INDEX"    => ++$ty_index,
											"TYPE"     => "CHARMM",
											"Lammps"   => getLammpsOpts("CHARMM","angle", $CNV),
											"VALS"     => [@{ $tmp }],
											"USED"     => 0,
											"KEY"      => "$1 $2 $3 ",
											"IGNORE"   => 0,
											"CTYPE"    => "",
										}
			);
		} elsif ($header eq "TORSIONS" and $_ =~ /^(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\d+\.\d+)\s+(\d+)\s+(\d+\.\d+)/) { #torsions
			if(exists($parms->{$header}{$4}) and exists($parms->{$header}{$4}{$3}) and exists($parms->{$header}{$4}{$3}{$2}) and exists($parms->{$header}{$4}{$3}{$2}{$1})) {
				$curr = \%{ $parms->{$header}{$4}{$3}{$2}{$1} };
			} else {
				$curr = \%{ $parms->{$header}{$1}{$2}{$3}{$4} };
			}
			if (! defined($curr) or ! findDuplicate($curr,"CHARMM")) {
				$ty_counter = 1;
				$curr->{$ty_counter} = (
										{
											"INDEX"    => ++$ty_index,
											"TYPE"     => "CHARMM",
											"Lammps"   => getLammpsOpts("CHARMM","dihedral", $CNV),
											"VALS"     => [$5,$6,$7],
											"USED"     => 0,
											"KEY"      => "$1 $2 $3 $4 ",
											"IGNORE"   => 0,
											"CTYPE"    => "",
											"NUM"      => 1,
											"COUL"     => "000",
											"VDW"      => "000",
											"PER"	   => 3,
											"do_scale" => 0,
											"FFTYPEID" => $ff->{FFTYPEID},
										}
				);
			} else {
				$ty_counter = findDuplicate($curr,"CHARMM");
				$curr->{$ty_counter}{NUM}++;
				push @{ $curr->{$ty_counter}{VALS} }, ($5,$6,$7);
			}
		} elsif ($header eq "INVERSIONS" and $_ =~ /^(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\-?\d+\.\d+)\s+(\d+)\s+(\d+\.\d+)/) { #improper
			@{ $tmp } = ($1,$2,$3,$4);
			@{ $tmp } = permute(\@{ $tmp });
			$rec = (
					{
						"TYPE"     => "IT_IJKL",
						"Lammps"   => getLammpsOpts("IT_IJKL","inversion", $CNV),
						"VALS"     => [$5, $7],
						"USED"     => 0,
						"IGNORE"   => 0,
						"CTYPE"    => "",
						"NUM"      => 1,
					}
			);
			++$ty_index;
			for $i (0, 3) { #make permutations
				$j = 5;
				($tmp->[$j][0], $tmp->[$j][$i]) = ($tmp->[$j][$i],$tmp->[$j][0]) if($i>0);
				$curr = \%{ $parms->{$header}{$tmp->[$j][0]}{$tmp->[$j][1]}{$tmp->[$j][2]}{$tmp->[$j][3]} };
				$rec->{KEY} = "@{ $tmp->[$j] } ";
				$ty_counter = findDuplicate($curr,"IT_IJKL") if(defined($curr));
				$ty_counter = 1 if (!defined($ty_counter) or ! $ty_counter);	
				$curr->{$ty_counter} = \%{ $rec }; 
				$curr->{$ty_counter}{INDEX} = $ty_index;
			}
		} elsif ($header eq "CMAP") { #CMAP: essentially a 5 body term 			
			if(! defined($curr) and $_ =~ /(\w+\s+){8}(\d+)/) { 
				undef $tmp; undef $j;
				while ($_ =~ /(\w+)\s+/g) {
					push @{ $tmp }, $1;
				}
				$j = $1 if ($_ =~ /(\d+)$/g);
				next if ($#{ $tmp } < 7 or ! defined($j));
				next if ($tmp->[1] ne $tmp->[4] or $tmp->[2] ne $tmp->[5] or $tmp->[3] ne $tmp->[6]); #check for correct torsion
				splice @{ $tmp }, 4, 3;
				#now search for already existing cmap
				@{ $cmap } = @{ $tmp };
				for (1 .. 1) {
					$curr = \%{ $parms->{$header} };
					$valid = 1;
					while ($#{ $cmap } > -1) {
						$i = pop @{ $cmap };
						$valid = 0 if (! exists($curr->{$i}));
						$curr = \%{ $curr->{$i} };
					}
					last if ($valid);
					@{ $cmap } = reverse @{ $tmp };
				}
				$ty_counter = scalar(keys %{ $curr }) + 1;
				$curr->{$ty_counter} = (
										{
											"N"     => $j, # number of entries
											"I"     => $j*5, #total number of integer lines (5 per entry)
											"C"     => 0, #current entry counter
											"USED"  => 1,
											"KEY"   => "@{ $cmap } ",
											"ATOMS" => "@{ $cmap } ",
											"INDEX" => ++${ty_index},
										}
				);
				$curr = \%{ $curr->{$ty_counter} };	
			} else {
				$_ =~ s/^\!/#/;
				$parms->{PARMS}{CMAP}{FILE_DATA} .= "$_\n";
				if ($_ =~ /^((\s+\-?\d+\.\d+){4})/) {
					$curr->{C}++;
					undef $curr if ($curr->{C} >= $curr->{I});
				}
			}
		} elsif ($header eq "VDW" and $_ =~ /CTOFNB\s+(\d+\\.\d+)\s+CTONNB\s+(\d+\.\d+)/) {
			$parms->{PARMS}{cut_coul} = $parms->{PARMS}{cut_vdw} = $1;
		} elsif ($header eq "VDW" and $_ =~ /^(\w+)\s+(\S+\s+\-?\d+\.\d+\s+\d+\.\d+\s*[\-?0-9\. ]*)/) { #vdw
			@{ $tmp } = split /\s+/, $2;
			@{ $tmp } = (@{ $tmp }, @{ $tmp }) if ($#{ $tmp } == 2);
			push @{ $tmp }, $1;
			if($tmp->[0] =~ /[a-z]/i) {
				next;
				$curr = \%{ $parms->{$header}{$tmp->[6]}{$tmp->[0]} }; #NBFIX off diag
				$t1 = $tmp->[6]; $t2 = $tmp->[0];
			} else {
				$curr = \%{ $parms->{$header}{$tmp->[6]}{$tmp->[6]} };
				$t1 = $t2 = $tmp->[6];
			}
			$ty_counter = findDuplicate($curr,"LJ_6_12");
			$ty_counter = 1 if (! defined($ty_counter) or ! $ty_counter);
			$tmp->[1] = abs($tmp->[1]); $tmp->[4] = abs($tmp->[4]);
			$tmp->[2] = 2*$tmp->[2]/2**(1/6); $tmp->[5] = 2*$tmp->[5]/2**(1/6);
			if(($tmp->[1] == $tmp->[4]) && ($tmp->[2] == $tmp->[5])) { #check for duplicate params
				@{ $tmp } = ($tmp->[1],$tmp->[2]);
			} else {
				@{ $tmp } = ($tmp->[1], $tmp->[2], $tmp->[4], $tmp->[5]);
			}
			$curr->{$ty_counter} = (
									{
										"INDEX"   => ++$ty_index,
										"TYPE"    => "LJ_6_12",
										"Lammps"  => getLammpsOpts("LJ_6_12","vdw", $CNV,0, $parms->{PARMS}),
										"KEY"     => "$t1 $t2 ",
										"VALS"    => [@{ $tmp }],
										"ATOM"    => "$t1 $t2 ",
										"USED"    => 0,
										"IGNORE"  => 0,
										"IT"      => "vdw",
									}
			);
		}
	}

	return $parms;
}

sub permute {
	my ($inArray) = @_;
	my (@PERMS, $firstAtm, $i);

	$firstAtm = shift @{ $inArray };
	@PERMS = General::Permutate($inArray, []);
	for $i (0 .. $#PERMS) {
		unshift @{ $PERMS[$i] }, $firstAtm;
	}
	return @PERMS;
};

1;

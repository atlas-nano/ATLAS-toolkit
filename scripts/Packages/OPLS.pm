package OPLS;

require Exporter;
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use constant q2e => 18.2223;
use Math::Trig qw(pi);
use General qw(FileTester FindElementByMass LoadConverter LoadElements);
use CERIUS2 qw(getLammpsOpts findDuplicate);

use strict;

our (@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(parseOplsItpFF);
$VERSION = "1.00";

sub numerically { ($a<=>$b); }

sub parseOplsItpFF {
	my ($ff_file, $alter, $oldFF) = @_;
	my ($i, $j, $k, $l, $in_data, $elements, $eleNum, %PARMS, $CNV, $curr); 
	my ($atmList, $EQUIV, $which_var, $validFF, @tmp, $type_counter, $rec);

	$validFF = $type_counter = $which_var = 0;
	$CNV = &General::LoadConverter();
	$elements = &General::LoadElements();
	if(defined($oldFF)) {
		%PARMS = %{ $oldFF };
	} else {
		$PARMS{PARMS} = (
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
								'mix_rule'								=> 'geometric',
								'same_scale'							=> 1,
								'scale_cou'								=> 0.50000,
								'scale_cou_12'							=> 0,
								'scale_cou_13'							=> 0,
								'scale_cou_14'							=> 0.50000,
								'scale_torsions'						=> 0,
								'scale_torsions_by_n_defined_torsions'	=> 0,
								'scale_vdw'								=> 0.50000,
								'scale_vdw_12'							=> 0,
								'scale_vdw_13'							=> 0,
								'scale_vdw_14'							=> 0.50000,
								'single_inversion'						=> 1,
							}
						);
	}

	open FORCEFIELD, $ff_file->{FF} or die "Cannot open force field file $ff_file->{FF}: $!\n";
	while (<FORCEFIELD>) {
		chomp;
		$in_data = $_;
		if($in_data =~ /\[ atomtypes \]/) {
			$which_var = 1;
			$type_counter = 0;
		} elsif($in_data =~ /\[ atoms \]/) {
			$which_var = 2;
			$type_counter = 0;
		} elsif($in_data =~ /\[ bonds \]/) {
			$which_var = 3;
			$type_counter = 1;
		} elsif($in_data =~ /\[ angles \]/) {
			$which_var = 4;
			$type_counter = 1;
		} elsif($in_data =~ /\[ dihedrals \]/) {
			$which_var = 5;
			$type_counter = 1;
		} elsif($which_var == 5 && $in_data =~ /IMPROPER DIHEDRAL ANGLES/) {
			$which_var = 6;
		} elsif($which_var == 1 && $in_data =~ /^\s*(opls_\S+)\s+(\S+)\s+(.*)/) {
			@tmp = split /\s+/, $3;
			next if ($#tmp < 4);
			$validFF = 1;
			$type_counter++;
			$eleNum = General::FindElementByMass($tmp[0], $elements);
			#$PARMS{EQUIVALENCE}{$1}{$2} = 1;
			$EQUIV->{$1} = $2;
			$rec = (
						{
							"TYPEID"        => $type_counter,
							"ATOM"          => $elements->{$eleNum}{SYMBOL},
							"MASS"          => $tmp[0],
							"USE_CHARGE"    => 1,
							"USED"          => 0,
							"OTHER"         => $tmp[1],
							"LABEL"         => $2,
							"CHARGE"        => $tmp[1],
							"CHARGE_SET"    => 0,
							}
					);
			$PARMS{ATOMTYPES}{$2} = $rec;
			@tmp = ($tmp[4]/4.184, $tmp[3]*10); #Angstroms and kcal
			$PARMS{VDW}{$2}{$2}{1} = (
										{
											"TYPE"    => "LJ_6_12",
											"Lammps"  => getLammpsOpts("LJ_6_12","vdw", $CNV,0,$PARMS{PARMS}),
											"KEY"     => "$2 $2 ",
											"VALS"    => [@tmp],
											"ATOM"    => "$2 $2 ",
											"USED"    => 0,
											"IGNORE"  => 0,
											"IT"      => "vdw",
										}
									);
		} elsif($which_var == 2 && $in_data =~ /^\s+(\d+)\s+(opls_\S+)\s+(.*)/) {
			@tmp = split /\s+/,$3;
			die "ERROR: Cannot find atom type entry for $2... Aborting\n"
				if(!exists($EQUIV->{$2}));
			$i = $EQUIV->{$2};
			$atmList->{$1}=$i;
			if(!$PARMS{ATOMTYPES}{$i}{CHARGE_SET}) {
				$PARMS{ATOMTYPES}{$i}{CHARGE_SET} = 1;
				$PARMS{ATOMTYPES}{$i}{CHARGE} = $tmp[4];
				$PARMS{OVERWRITE_FFTYPES}{$tmp[2]} = $i;
			} elsif($PARMS{ATOMTYPES}{$i}{CHARGE} != $tmp[4]) {
				print "";
			}
		} elsif($which_var == 3 && $in_data =~ /^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
			$i = $atmList->{$1}; $j = $atmList->{$2};
			if(exists($PARMS{BONDS}{$j}) and exists($PARMS{BONDS}{$j}{$i})) {
				$curr=\%{$PARMS{BONDS}{$j}{$i}};
			} else {
				$curr=\%{$PARMS{BONDS}{$i}{$j}};
			}
			$validFF = 1;
			@tmp = ($5/4.184/200,$4*10);
			$type_counter = findDuplicate($curr,"HARMONIC");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			$curr->{$type_counter} = (
										{
											"INDEX"   => $type_counter,
											"TYPE"    => "HARMONIC",
											"Lammps"  => getLammpsOpts("HARMONIC","bond", $CNV),
											"VALS"    => [@tmp],
											"USED"    => 0,
											"KEY"     => "$i $j ",
											"IGNORE"  => 0,
										}
									);
		}elsif($which_var == 4 && $in_data =~ /^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
			$i = $atmList->{$1}; $j = $atmList->{$2}; $k = $atmList->{$3};
			if(exists($PARMS{ANGLES}{$k}) and exists($PARMS{ANGLES}{$k}{$j}) and exists($PARMS{ANGLES}{$k}{$j}{$i})) {
				$curr=\%{ $PARMS{ANGLES}{$k}{$j}{$i} };
			} else {
				$curr=\%{ $PARMS{ANGLES}{$i}{$j}{$k} };
			}
			$validFF = 1;
			@tmp = ($6/4.184/2,$5);
			$type_counter = findDuplicate($curr,"HARMONIC");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			$curr->{$type_counter} = (
										{
											"INDEX"   => $type_counter,
											"TYPE"    => "THETA_HARM",
											"Lammps"  => getLammpsOpts("THETA_HARM","angle", $CNV),
											"VALS"    => [@tmp],
											"USED"    => 0,
											"KEY"     => "$i $j $k ",
											"IGNORE"  => 0,
											"CTYPE"   => "",
										}
									);
		}elsif($which_var == 5 && $in_data =~ /^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\-?\d+\.\d+\s+.*)/) {
			$i = $atmList->{$1}; $j = $atmList->{$2}; $k = $atmList->{$3}; $l = $atmList->{$4};
			@tmp = map { $_ /= 4.184 } split /\s+/,$6;
			@tmp = (-3*$tmp[3]/2-2*$tmp[1], -$tmp[4]-$tmp[2], -0.5*$tmp[3], -0.24*$tmp[4]);
			if(exists($PARMS{TORSIONS}{$l}) and exists($PARMS{TORSIONS}{$l}{$k}) and exists($PARMS{TORSIONS}{$l}{$k}{$j}) and exists($PARMS{TORSIONS}{$l}{$k}{$j}{$i}) and exists($PARMS{TORSIONS}{$l}{$k}{$j}{$i})) {
				$curr=\%{ $PARMS{TORSIONS}{$l}{$k}{$j}{$i} };
			} else {
				$curr=\%{ $PARMS{TORSIONS}{$i}{$j}{$k}{$l} };
			}
			$type_counter = findDuplicate($curr,"OPLS");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			$curr->{$type_counter} = (
										{
											"INDEX"   => $type_counter,
											"TYPE"    => "OPLS",
											"Lammps"  => getLammpsOpts("OPLS","dihedral",$CNV),
											"VALS"    => [@tmp],
											"USED"    => 0,
											"KEY"     => "$i $j $k $l ",
											"IGNORE"  => 0,
											"CTYPE"   => "",
											"NUM"	  => 1,
											"COUL"    => "0010.5",
											"VDW"     => "0010.5",
											"FFTYPEID" => $ff_file->{FFTYPEID},
										}
									);
		}elsif($which_var == 6 && $in_data =~ /^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+\s+.*)/) {
			$i = $atmList->{$3}; $j = $atmList->{$2}; $k = $atmList->{$4}; $l = $atmList->{$1};
			@tmp = split /\s+/,$6;
			if(exists($PARMS{INVERSIONS}{$j})) {
				$curr = \%{ $PARMS{INVERSION}{$j} };
				if(exists($curr->{$i}) and exists($curr->{$i}{$k}) and exists($curr->{$i}{$k}{$l})) {
					$curr = \%{ $curr->{$i}{$k}{$l} };
				}elsif(exists($curr->{$i}) and exists($curr->{$i}{$l}) and exists($curr->{$i}{$l}{$k})) {
					$curr = \%{ $curr->{$i}{$l}{$k} };
				}elsif(exists($curr->{$k}) and exists($curr->{$k}{$i}) and exists($curr->{$k}{$i}{$l})) {
					$curr = \%{ $curr->{$k}{$i}{$l} };
				}elsif(exists($curr->{$k}) and exists($curr->{$k}{$l}) and exists($curr->{$k}{$l}{$i})) {
					$curr = \%{ $curr->{$k}{$l}{$i} };
				}elsif(exists($curr->{$l}) and exists($curr->{$l}{$i}) and exists($curr->{$l}{$i}{$k})) {
					$curr = \%{ $curr->{$l}{$i}{$k} }
				}elsif(exists($curr->{$l}) and exists($curr->{$l}{$k}) and exists($curr->{$l}{$k}{$i})) {
					$curr = \%{ $curr->{$l}{$i}{$k} }
				}else {
					$curr = \%{ $curr->{$i}{$k}{$l} };
				}
			} else {
				$curr=\%{ $PARMS{INVERSIONS}{$j}{$i}{$k}{$l} };
			}
			$validFF = 1;
			$type_counter = findDuplicate($curr,"IT_IJKL");
			$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
			if($tmp[0] == 0) {
				$tmp[0] = 1;
			} else {
				$tmp[0] = -1;
			}
			@tmp = ($tmp[1]/4.184,$tmp[0],$tmp[2]);
			$curr->{$type_counter} = (
										{
											"INDEX"   => $type_counter,
											"TYPE"    => "IT_IJKL",
											"Lammps"  => getLammpsOpts("IT_JIKL","inversion", $CNV),
											"VALS"    => [@tmp],
											"USED"    => 0,
											"KEY"     => "$i $j $k $l ",
											"IGNORE"  => 0,
											"CTYPE"   => "",
										}
									);
		}

	}
	close FORCEFIELD;
	die "ERROR: Invalid OPLS-ITP force field file $ff_file->{FF}\n"
		if (! $validFF);
	return \%PARMS;	
}

1;

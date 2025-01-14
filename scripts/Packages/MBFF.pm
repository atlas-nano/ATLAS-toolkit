package MBFF;

require Exporter;
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use Math::Trig qw(pi);
use Storable qw(dclone); 
use General qw(FileTester FindElement LoadConverter LoadElements);
use CERIUS2 qw(getLammpsOpts findDuplicate);

use strict;

our (@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(parseSWFF);
$VERSION = "1.00";

sub numerically { ($a<=>$b); }

sub parseSWFF {
	my ($ff_file, $alter, $oldFF) = @_;
	my ($CNV, $elements, %PARMS, $valid, $eleNum, $rec, $id_counter);
	my ($in_data, @tmp, $curr, $type_counter, $i, $j, $k);

	$CNV = General::LoadConverter();
	$elements = General::LoadElements();

	if(defined($oldFF)) {
		%PARMS = %{ $oldFF };
	} else {
		$PARMS{PARMS} = (
							{
								'EXO_CYCLIC'							=> 1.00000,
								'HAS_HBONDS'							=> 0,
								'QEq'									=> 0,
								'SW'									=> 1,
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
								'scale_cou_12'							=> 0.50000,
								'scale_cou_13'							=> 0.50000,
								'scale_cou_14'							=> 0.50000,
								'scale_torsions'						=> 0,
								'scale_torsions_by_n_defined_torsions'	=> 0,
								'scale_vdw'								=> 0.50000,
								'scale_vdw_12'							=> 0.50000,
								'scale_vdw_13'							=> 0.50000,
								'scale_vdw_14'							=> 0.50000,
								'single_inversion'						=> 1,
							}
						);
	}

	$id_counter = 0;
	open FORCEFIELD, $ff_file->{FF} or die "Cannot open force field file $ff_file->{FF}: $!\n";
	while (<FORCEFIELD>) {
		chomp;
		$in_data = $_;
		#put data into array
		@tmp = split /\s+/, $in_data;
		#check that there are at least 14 columns
		next if ($#tmp < 13);
		#columns 4 - 14 should be numeric
		$valid = 1;
		for $i (3 .. 13) {
			if($tmp[$i] !~ /^\-?\d+\.?\d*e?\+?\-?\d*/) {
				$valid = 0;
				last;
			}
		}
		next if (! $valid);
		#1st 3 columns must be elements
		$valid = 1;
		for $i (0 .. 2) {
			$eleNum = General::FindElement($elements, $tmp[$i]);
			if(! $eleNum) {
				$valid = 0;
				print "ERROR: Cannot determine element of type $tmp[$i]...";
				last;
			}
			if(!exists($PARMS{ATOMTYPES}{$tmp[$i]})) {
				#record element data
				$rec = (
							{
								"TYPEID"        => $id_counter++,
								"ELEMENT"       => $eleNum,
								"ATOM"          => $elements->{$eleNum}{SYMBOL},
								"MASS"          => $elements->{$eleNum}{MASS},
								"USE_CHARGE"    => 0,
								"USED"          => 0,
								"LABEL"         => $tmp[$i],
								"CHARGE"        => 0,
							}
				);
				$PARMS{ATOMTYPES}{$tmp[$i]} = $rec;
			}
		}
		die "\n" if (! $valid);
		#check to see if all 3 elements are the same
		$i = $tmp[0]; $j = $tmp[1]; $k = $tmp[2];
		if($i == $j and $j == $k) {
			#record VDW entry
			$curr = ();
			$type_counter = 0;
			if (exists($PARMS{VDW}{$i}) and exists($PARMS{VDW}{$i}{$j})) {
				$curr = \%{ $PARMS{VDW}{$i}{$j} };
				$type_counter = scalar(%{ $curr });
			} elsif(exists($PARMS{VDW}{$j}) and exists($PARMS{VDW}{$j}{$i})) {
				$curr = \%{ $PARMS{VDW}{$j}{$i} };
				$type_counter = scalar(%{ $curr });
			} else {
				$curr = \%{ $PARMS{VDW}{$i}{$j} };
			}
			$curr->{$type_counter} = (
										{
											"TYPE"    => "SW",
											"Lammps"  => getLammpsOpts("SW","vdw", $CNV,0,$PARMS{PARMS}),
											"KEY"     => "$i $i ",
											"VALS"    => [$tmp[4], $tmp[5]],
											"ATOM"    => "$i $i ",
											"USED"    => 0,
											"IGNORE"  => 0,
											"IT"      => "vdw",
											"PAIRMIX" => 0,
										}
									);
		}
		#finally record angle entry
		$curr=\%{ $PARMS{ANGLES}{$i}{$j}{$k} };
		$type_counter = findDuplicate($curr,"SW");
		$type_counter = 1 if (! defined($type_counter) or ! $type_counter);
		#remove the first 3 elements of the tmp array, saving only the values
		splice @tmp, 0, 3;
		$curr->{$type_counter} = (
										{
											"INDEX"   => $type_counter,
											"TYPE"    => "SW",
											"Lammps"  => getLammpsOpts("SW","angle", $CNV),
											"USED"    => 0,
											"KEY"     => "$i $j $k ",
											"IGNORE"  => 0,
											"CTYPE"   => "",
											"VALS"    => [@tmp],
										}
								);
	}
	close FORCEFIELD;
	return \%PARMS;	
}

1;

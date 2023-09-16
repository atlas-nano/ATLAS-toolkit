#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use General qw(FileTester TrjSelections GetBondLength GetMinDist dPcmplx);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);
use BOX qw(GetBox GetNlist);
use Math::SphericalHarmonics qw(ylm);
use Math::Trig qw(acos atan);

sub init;
sub getOrderParm;
sub calcTrjOrderParm;
sub writeOrderParmData;
sub numerically { ($a<=>$b) }

my ($bgfFile, $selection, $trjType, $trjFile, $SELECT, $saveFile, $nQ, $qlist);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $SELECTIONS);
my ($OrderParmData, $field, $getByteOffset, $getSnapshots, $LAMMPSOPTS, $pStr);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
$BOX = GetBox($ATOMS, undef, $HEADERS);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
print "Done\n";
if (! defined($trjFile)) {
	printf "Calculating order parameter from $bgfFile...";
	$OrderParmData->[0] = getOrderParm($ATOMS, $BONDS, $SELECTIONS, $BOX);
	$OrderParmData->[0]{TSTEP} = 0;
	print "Done\n";
} else {
	$field = scalar keys %{ $ATOMS };
	$getByteOffset->($SELECT, $trjFile, $field);
	if ($trjType == 2) {
		&GetLammpsTrjType($SELECT, $trjFile, "coord", \%{ $LAMMPSOPTS });
		$field = "coord";
	}
	$pStr = "Calculating order parameter from $trjFile...";
	$getSnapshots->($ATOMS, $trjFile, $SELECT, $field, \&calcTrjOrderParmData, $pStr, undef);
}
print "Writing results to $saveFile...";
&writeOrderParmData($ATOMS, $OrderParmData, $saveFile);
print "Done\n";

sub writeOrderParmData {
	my ($atoms, $data, $outfile) = @_;
	my ($i, $j, $k, $tstep, $alist, $atm_label);

	@{ $alist } = sort numerically keys %{ $data->[0]{0} };
	shift @{ $alist };
	open OUTFILE, ">$outfile" or die "ERROR: Cannot create $outfile\n";
	printf OUTFILE "#%-11s ","Atom";
	for $i (0 .. $nQ) {
		printf OUTFILE "%12s ","Q = $qlist->[$i]";
	}
	printf OUTFILE "\n";
	for $i (@{ $data }) {
		printf OUTFILE "#TIMESTEP %s\n",$i->{TSTEP};
		for $j (@{ $alist }) {
			$atm_label = $atoms->{$j}{ATMNAME};
			$atm_label =~ s/\s+$//;
			$atm_label = "${j} ($atm_label)";
			printf OUTFILE "%-12s ",$atm_label;
			for $k (0 .. $nQ) {
				printf OUTFILE "%12.8f ",$i->{$k}{$j};
			}
			printf OUTFILE "\n";
		}
		printf OUTFILE "%-12s ","#AVG";
		for $k (0 .. $nQ) {
			printf OUTFILE "%12.8f ",$i->{$k}{AVG};
		}
		printf OUTFILE "\n\n";
	}

	close OUTFILE;
}

sub getOrderParm {
	my ($atoms, $bonds, $select, $box) = @_;
	my ($i, $j, $l, $m, $nlist, $n, $qlm, $qilm, $lval);
	my ($dotP, $n1, $n2, $theta, $phi, $ylm_val, $count);

	&GetNlist($atoms, $select, $box, 3.2, \%{ $nlist });
	$qlm = $qilm = ();
	for $i (keys %{ $nlist }) {
		$n = scalar keys %{ $nlist->{$i}  };
		for $j (keys %{ $nlist->{$i} }) {
			$theta = $nlist->{$i}{$j}{theta};
			$phi = $nlist->{$i}{$j}{phi};
			for $l (0 .. $nQ) {
				$lval = $qlist->[$l];
				for $m (-$lval .. $lval) {
					$ylm_val = ylm($lval,$m,$theta, $phi);
					$qilm->{$i}{$l}{$m}{re} += $ylm_val->{re}/$n;
					$qilm->{$i}{$l}{$m}{im} += $ylm_val->{im}/$n;
				}
			}
		}
	}
	for $i (keys %{ $nlist }) {
		$n = scalar keys %{ $nlist->{$i}  };
		for $l (0 .. $nQ) {
			$lval = $qlist->[$l];
			$qlm->{$l}{$i} = 0;
			for $j (keys %{ $nlist->{$i} }) {
				$dotP = dPcmplx($qilm->{$i}{$l},$qilm->{$j}{$l},-1);
				$n1 = sqrt(dPcmplx($qilm->{$i}{$l},$qilm->{$i}{$l},1));
				$n2 = sqrt(dPcmplx($qilm->{$j}{$l},$qilm->{$j}{$l},1));
				$qlm->{$l}{$i} += $dotP/$n1/$n2/$n; 
			}
		}
	}
	for $l (0 .. $nQ) {
		$qlm->{$l}{AVG} /= scalar(keys %{ $nlist });
	}

	return $qlm;

}

sub calcTrjOrderParmData {
	my ($atoms, $box, $frameNum, $fileHandle) = @_;
	my ($tot, $rec);

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

	for ("X", "Y", "Z") {
		$box->{$_} = $box->{"${_}COORD"};
	}
	$rec = getOrderParm($atoms, $BONDS, $SELECTIONS, $box);
	$rec->{TSTEP} = $frameNum;
	push @{ $OrderParmData }, $rec;
}

sub init {
	my (%OPTS, $atm, $tSel, $list, $i, $qtmp);

	getopt('btawsq',\%OPTS);
	die "usage: $0 -b bgf file -q [q list] -a (atom selection) -t (traj file) -s (traj selection) -w (save file)\n"
		if (! exists($OPTS{b}) or ! exists($OPTS{q}));
	print "Initializing...";
	($bgfFile, $qtmp, $atm, $trjFile, $tSel, $saveFile) = 
		($OPTS{b}, $OPTS{q}, $OPTS{a}, $OPTS{t}, $OPTS{s}, $OPTS{w});
	FileTester($bgfFile);
	$atm = "index>0" if (!defined($atm));
	$selection = BuildAtomSelectionString($atm);
	while($qtmp =~ /(\d+)/g) {
		push @{ $qlist }, $1;
	}
	die "ERROR: Expected integer(s) for q list. Got \"$qlist\"\n" if (! defined($qlist));
	$nQ = $#{ $qlist };
	
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
	}
	if (! defined($saveFile)) {
		$saveFile = basename($bgfFile);
		$saveFile =~ s/\.\w+$//;
		$saveFile .= ".orderParm.dat";
	}
	print "Done\n";
}

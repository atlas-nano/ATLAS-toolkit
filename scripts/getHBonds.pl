#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use General qw(FileTester TrjSelections LoadElements GetBondLength GetAngle AddElementField);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);
use BOX qw(GetBox);

sub init;
sub getHBdata;
sub calcTrjHBdata;
sub writeHBdata;
sub writeDonorInfo;
sub numerically { ($a <=> $b); };

my ($bgfFile, $selection, $trjType, $trjFile, $SELECT, $saveFile, $dInfo);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $SELECTIONS, $ELEMENTS, $MOLS);
my ($HBdata, $field, $getByteOffset, $getSnapshots, $LAMMPSOPTS, $pStr);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
$BOX = GetBox($ATOMS, undef, $HEADERS);
&AddElementField($ATOMS, $ELEMENTS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
$MOLS = GetMols($ATOMS, $BONDS, $SELECTIONS);
print "Done\n";
if (! defined($trjFile)) {
	$HBdata->{1} = getHBdata($ATOMS, $BONDS, $SELECTIONS, $BOX);
} else {
	$field = scalar keys %{ $ATOMS };
	$getByteOffset->($SELECT, $trjFile, $field);
	if ($trjType == 2) {
		&GetLammpsTrjType($SELECT, $trjFile, "coord", \%{ $LAMMPSOPTS });
		$field = "coord";
	}
	$pStr = "Calculating hbond stats data from $trjFile...";
	$getSnapshots->($ATOMS, $trjFile, $SELECT, $field, \&calcTrjHBdata, $pStr, undef);
}
&writeHBdata($ATOMS, $HBdata, $saveFile);
&writeDonorInfo($ATOMS, $MOLS, $HBdata, $saveFile) if ($dInfo);

sub writeDonorInfo {
	my ($atoms, $mols, $hbdata, $outFile) = @_;
	my ($i, $j, $k, $ndonor, $nacceptor, $mtype, $atmList);

	$outFile =~ s/\.\w+$//;
	$outFile .= ".donors.dat";

	print "Writing HB acceptor/donor stats to $outFile...";
	open OUTFILE, "> $outFile" or die "ERROR: Cannot create $outFile: $!\n";
	printf OUTFILE "#%-7s %8s %8s %8s %8s\n","MOLID","NDONOR","NACCEPTOR","CLASS", "ATOMS";
	for $i (sort numerically keys %{ $hbdata }) {
		printf OUTFILE "#TIMESTEP $i $hbdata->{$i}{COUNT} total hbonds\n";
		for $j (sort numerically keys %{ $mols }) {
			$ndonor = $nacceptor = 0;
			if (exists($hbdata->{$i}{MOLS}{$j})) {
				$ndonor = $hbdata->{$i}{MOLS}{$j}{ndonor};
				$nacceptor = $hbdata->{$i}{MOLS}{$j}{nacceptor};
			}
			$mtype = "dd";
			$mtype = "sd" if($ndonor == 1);
			$mtype = "nd" if ($ndonor == 0);
			$mtype = "nohb" if ($ndonor == 0 and $nacceptor == 0);
			$atmList = "";
			for $k (sort numerically keys %{ $mols->{$j}{MEMBERS} }) {
				$atmList .= "# $k ($atoms->{$k}{ATMNAME}) ";
			}
			printf OUTFILE "%-8d %8d %8d %8s #$atmList\n",$j, $ndonor, $nacceptor, $mtype;
		}
	}
	close OUTFILE or die "ERROR: Cannot complete writing to $outFile: $!\n";
	print "Done\n";
}

sub writeHBdata {
	my ($atoms, $data, $outfile) = @_;
	my ($i, $j, $k, $l, $curr);
	my ($iname, $jname, $kname);

	open OUTFILE, ">$outfile" or die "ERROR: Cannot create $outfile\n";
	printf OUTFILE "#%-12s %-12s %-12s %13s %12s %12s\n","Donor","Acceptor","Hydrogen","D-A dist","D-H-A Angle","Energy";
	for $l (sort numerically keys %{ $data }) {
		printf OUTFILE "#TIMESTEP $l $data->{$l}{COUNT} total hbonds\n";
		for $i (sort numerically keys %{ $data->{$l}{LIST} }) {
			$iname = $atoms->{$i}{ATMNAME};
			for $j (sort numerically keys %{ $data->{$l}{LIST}{$i} }) {
				$jname = $atoms->{$j}{ATMNAME};
				for $k (sort numerically keys %{ $data->{$l}{LIST}{$i}{$j} }) {
					$kname = $atoms->{$k}{ATMNAME};
					$curr = $data->{$l}{LIST}{$i}{$j}{$k};
					printf OUTFILE "%-4d %-8s %-4d %-8s %-4s %-8s %12.3f %12.3f %12.3f\n",
						$i,$iname,$j,$jname,$k,$kname,$curr->{dist},$curr->{angle},$curr->{eng};
				}
			}
		}
	}
	close OUTFILE;
}

sub getHBdata {
	my ($atoms, $bonds, $select, $box) = @_;
	my ($i, $j, $k, $alist, $hlist, $angle, $data, $rec);
	my ($dr, $de, $sigma, $mol1, $mol2);

	$de = 4.0;
	$sigma = 2.82;
	$alist = getDonorAcceptorList($atoms, $select, $box, 3.5);
	for $i (keys %{ $alist }) {
		$hlist = ();
		for $j (@{ $bonds->{$i} }) {
			next if ($atoms->{$j}{ELEMENT}{SYMBOL} ne "H"); #hydrogen check
			$hlist->{$j} = 1;
		}
		next if (! keys %{ $hlist });
		for $j (keys %{ $alist->{$i} }) {
			for $k (keys %{ $hlist }) {
				$angle = GetAngle($atoms->{$i}, $atoms->{$k}, $atoms->{$j}, $box, 1);
				#next if ($angle < 2.79 or $angle > 3.49); #angle = 180 +/-20
				next if ($angle < 2.618 or $angle > 3.666); #angle = 180 +/-30
				#next if ($angle < 2.53 or $angle > 3.753); #angle = 180 +/-35
				#next if ($angle < 2.44 or $angle > 3.84); #angle = 180 +/-40
				#next if ($angle < 2.356 or $angle > 3.927); #angle = 180 +/-45
				#next if ($angle < 1.57 or $angle > 4.72); #angle = 180 +/-90
				$rec = ();
				$rec->{dist} = $alist->{$i}{$j};
				$rec->{angle} = $angle*180/3.141592;
				$dr = $sigma/$rec->{dist};
				$rec->{eng} = $de*(5*$dr**12-6*$dr**6)*cos($angle)**4; #dreiding hbond
				next if ($rec->{eng} > 0);
				$data->{COUNT}++;
				$mol1 = ${ $atoms->{$i}{MOLECULEID} };
				$mol2 = ${ $atoms->{$j}{MOLECULEID} };
				$data->{MOLS}{$mol1}{ndonor}++;
				$data->{MOLS}{$mol2}{nacceptor}++;
				$data->{LIST}{ $atoms->{$i}{INDEX} }{ $atoms->{$j}{INDEX} }{ $atoms->{$k}{INDEX} } = $rec;
			}
		}
	}

	return $data;

}

sub getDonorAcceptorList {
	my ($atoms, $select, $box, $cut) = @_;
	my ($cutsq, $i, $j, $dist, $distsq, $adlist);

	$cutsq = $cut*$cut;

	for $i (keys %{ $select }) {
		next if ($atoms->{$i}{ELEMENT}{SYMBOL} !~ /O|N|F|S/); #donor/acceptor check
		for $j (keys %{ $select }) {
			next if ($i == $j);
			next if ($atoms->{$j}{ELEMENT}{SYMBOL} !~ /O|N|F|S/); #donor/acceptor check
			$dist = GetBondLength($atoms->{$i}, $atoms->{$j}, $box);
			$distsq = $dist * $dist;
			if ($distsq < $cutsq) {
				$adlist->{$i}{$j} = $adlist->{$j}{$i} = $dist;
			}
		}
	}

	return $adlist;
}

sub calcTrjHBdata {
	my ($atoms, $box, $frameNum, $fileHandle) = @_;
	my ($tot);

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

	$HBdata->{$frameNum} = getHBdata($atoms, $BONDS, $SELECTIONS, $box);
}


sub init {
	my (%OPTS, $atm, $tSel, $list, $i);

	getopt('btawds',\%OPTS);
	die "usage: $0 -b bgf file -a (atom selection) -t (traj file) -s (traj selection) -w (save file) -d (write donor info=no)\n"
		if (! exists($OPTS{b}));
	print "Initializing...";
	($bgfFile, $atm, $trjFile, $tSel, $saveFile, $dInfo) = 
		($OPTS{b}, $OPTS{a}, $OPTS{t}, $OPTS{s}, $OPTS{w}, $OPTS{d});
	FileTester($bgfFile);
	$atm = "index>0" if (!defined($atm));
	$selection = BuildAtomSelectionString($atm);
	
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
		$saveFile .= ".HBdata.dat";
	}
	$ELEMENTS = LoadElements();
	$dInfo = 0 if (! defined($dInfo) or $dInfo !~ /1|yes/i);
	$dInfo = 1 if ($dInfo =~ /1|yes/i);

	print "Done\n";
}

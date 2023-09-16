#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Storable qw(dclone);

use Math::Trig qw(pi);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester LoadElements);
use OPLS qw(parseOplsItpFF); 
use CERIUS2 qw(saveCeriusFF);

sub parseGromacsFF;
sub updateGroParams;
sub addFields;

my ($itpFile, $saveName);
my ($PARMS, $FF, $ELEMENTS, $itpType);

$|++;
&init;
print "Parsing OPLS ITP force field file $itpFile...";
$FF = parseOplsItpFF($itpFile, 0, undef);
print "Done\nCreating CERIUS2 formatted forcefield $saveName...";
&saveCeriusFF($FF, $saveName, $ELEMENTS);
if(-e "$Bin/../ff/OPLSFF_header.ff") {
	system("mv $saveName tmp");
	system("cat $Bin/../ff/OPLSFF_header.ff tmp > $saveName");
}
print "Done\n";

sub READITPFF {
	my ($parms) = $_[0];
	my ($i,$j, $k, $ffid, $header, $valType, $types, $ids, $tmp, $ele);

	if ($itpType == 1) {
		for $i (keys %{ $parms->{ATOMS} }) {
			$ffid = $parms->{ATOMS}{$i}{TYPE};
			$parms->{ATOMTYPES}{$ffid} = \%{ $parms->{ATOMS}{$i} };
		}
	} elsif ($itpType == 2) {
		$tmp = dclone($parms->{ATOMTYPES});
		$parms->{ATOMTYPES} = ();
		for $i (keys %{ $tmp }) {
			$ele = $ELEMENTS->{$tmp->{$i}{ELEMENT}}{SYMBOL};
			$ffid = $tmp->{$i}{FFTYPE};
			$ffid =~ s/opls_/$ele/;use Storable qw(dclone);
			$parms->{atoms}{$ffid}{VALS}[0]{r} = $tmp->{$i}{sigma}*2**(1/6)*5; #Angstroms and 2^(1/6) factor
			$parms->{atoms}{$ffid}{VALS}[0]{e} = $tmp->{$i}{epsilon}/4.184; #kcal
			$parms->{atoms}{$ffid}{VALS}[0]{element} = $tmp->{$i}{ELEMENT};
			$parms->{atoms}{$ffid}{VALS}[0]{mass} = $tmp->{$i}{MASS};
			$parms->{atoms}{$ffid}{VALS}[0]{hybrid} = 0;
			$parms->{atoms}{$ffid}{VALS}[0]{equivalence} = $tmp->{$i}{EQUIVALENCE};
		}
	}
	for $i ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS") {
		$header =  substr($i, 0, -1);
		$header .= "TYPES";
		for $j (keys %{ $parms->{$i} }) {
			if($itpType == 1) {
				@{ $ids } = ($parms->{$i}{$j}{I},$parms->{$i}{$j}{J},$parms->{$i}{$j}{K},$parms->{$i}{$j}{L});
				pop @{ $ids } if ($i =~ /ANGLES|BONDS/);
				pop @{ $ids } if ($i =~ /BONDS/);
				$types = getTypesFromId($ids, $parms->{ATOMS});
			} else {
				@{ $types } = @{ $parms->{$i}{$j}{ATOMS} };
			}
			$valType = \%{ $parms->{lc $i}{shift @{ $types }} };
			for $k (0 .. $#{ $types }) {
				$valType = \%{ $valType->{ $types->[$k] } };
			}
			&getC2Vals($parms->{$i}{$j},$header,$valType) if ($itpType == 1);
			&getNewC2Vals($parms->{$i}{$j},$i,$valType) if ($itpType == 2);
		}
	}
	if ($itpType == 2) { 
		$parms->{torsionOrders}{1}= "0 1 2 3";
		$parms->{inversionOrders}{1}= "1 0 2 3";
	}
}

sub getNewC2Vals {
	my ($parm, $header, $valType) = @_;
	my ($rec, $i, @vals);

	if($header =~ /BOND/) {
		if ($parm->{FUNCT} == 1) {
			$rec->{kb} = $parm->{VALS}[0]/418.4/2;
			$rec->{r0} = $parm->{VALS}[1]*10;
		}
	} elsif ($header =~ /ANGLE/) {
		if ($parm->{FUNCT} == 1) {
			$rec->{t0} = $parm->{VALS}[1]*pi/180;
			$rec->{kt} = $parm->{VALS}[0]/4.184/2;
		}
	} elsif ($header =~ /TORSION/) {
		if($parm->{FUNCT} == 3) { #Ryckaert-Bellemans (multi/harmonic in lammps)
			$valType->{TYPE} = "MHARMONIC";
			$valType->{counter} = 1;
			for $i (0 .. $#{ $parm->{VALS} }) {
				$rec->[$i] = $parm->{VALS}[$i]/4.184;
			}
		}
	} elsif ($header =~ /INVERSION/) {
		if (! exists($parm->{FUNCT})) { #old stype itp for opls
			$valType->{counter} = 1;
			$rec->{p0} = $parm->{VALS}[0];
			$rec->{kp} = $parm->{VALS}[1]/4.184;
			$rec->{n} =  $parm->{VALS}[2];
		}
	} 
	$valType->{VALS}[0] = $rec;
}

sub getC2Vals {
	my ($parm, $header, $valType) = @_;
	my ($rec, $j, @vals);

	if ($header =~ /BOND/) {
		if ($parm->{FUNCT} == 2) {
			$rec->{kb} = $parm->{C0}/418.4/2;
			$rec->{r0} = $parm->{C1}*10;
		}
	} elsif ($header =~ /ANGLE/) {
		if ($parm->{FUNCT} == 2) {
			$rec->{t0} = $parm->{THETA}*pi/180;
			$rec->{kt} = $parm->{K0}/4.184/2;
		}
	} elsif ($header =~ /DIHEDRAL/) {
		if ($parm->{FUNCT} == 1) {
			@vals = ();
#			$vals[0] = -2*$parms->{$i}{$j}{VALS}[1]-3*$parms->{$i}{$j}{VALS}[3]/2;
#			$vals[1] = -$parms->{$i}{$j}{VALS}[2] + $parms->{$i}{$j}{VALS}[4];
#			$vals[2] = -$parms->{$i}{$j}{VALS}[3]/2;
#				for (0 .. $#vals) {
#					next if ($vals[$_] < 0.0001);
#					$rec = ();
#					$rec->{n} = $_ + 1;
#					$rec->{kp} = $vals[$_]/4.184/4;
#					$rec->{p0} = pi;
#				}
#				$rec->{torsionOrders} = "0 1 2 3";
#			} elsif($i eq "INVERSIONS") {
#				$curr->{VALS}[0]{p0} = $parms->{$i}{$j}{VALS}[0];
#				$curr->{VALS}[0]{kp} = $parms->{$i}{$j}{VALS}[1]/4.184/2;
#				$curr->{VALS}[0]{n} = $parms->{$i}{$j}{VALS}[2];
#				$curr->{counter} = $c;
#				$ff->{inversionOrders}{$c} = "2 0 1 3";
#				$c++;
		}
	}

	$valType->{VALS}[0] = $rec;

}

sub getTypesFromId {
	my ($ids, $parms) = @_;
	my ($i, @fftypes, $str1, $str2);

	for $i (@{ $ids }) {
		die "ERROR: Index $i does not have a corresponding atom entry!\n"
			if (! exists($parms->{$i}));
		push @fftypes, $parms->{$i}{TYPE};
	}
	$str1 = $str2 = "";
	for $i (0 .. $#fftypes) {
		$str1 .= $fftypes[$i];
		$str2 .= $fftypes[$#fftypes-$i];
	}
	@fftypes = reverse @fftypes if ($str2 gt $str1);
	return \@fftypes;
}
		
sub init {
	my (%OPTS);

	getopt('is',\%OPTS);
	($itpFile,$saveName) = ($OPTS{i}, $OPTS{s});
	die "usage: $0 -i itp_file -s (cerius2_ff)\n"
		if (! exists($OPTS{i}));
	print "Initializing...";
	FileTester($itpFile) if (defined($itpFile));
	if (! defined($saveName)) {
		$saveName = basename($itpFile);
		$saveName =~ s/\.w+$//;
		$saveName .= ".opls.ff";
	}
	$ELEMENTS = LoadElements();
	print "Done\n";
}

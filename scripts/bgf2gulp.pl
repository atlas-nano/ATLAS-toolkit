#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename dirname);
use FileFormats qw(GetBGFFileInfo AddMass);
use General qw(FileTester LoadFFs ReadFFs);
use CERIUS2 qw(GenerateUFFParms);
use BOX qw(GetBox);

sub loadBGFS;
sub checkAtomTypes;
sub writeGulpData;
sub loadenergies;
sub getGulpType;

my ($BGF, $saveName, $FF, $opt_type, $efactor);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $fStr, $fflist, $cell_opt, $efile, $elast);

$|++;
&init;
($fflist, undef) = ReadFFs($fStr);
$FF = LoadFFs($fflist);
$cell_opt = loadBGFS($BGF, $FF);
&GenerateUFFParms($BGF->[0]{ATOMS}, $BGF->[0]{BONDS}, $FF) if exists($FF->{PARMS}{UFF});
$elast = loadenergies($BGF, $efile) if (defined($efile));
print "Writing GULP file $saveName...";
&pruneVDW($FF);
&writeGulpData($BGF, $FF, $cell_opt, $saveName, $elast);
print "Done\n";

sub writeGulpData {
	my ($fData, $parms, $has_cell, $outName, $elast) = @_;
	my ($i, $j, $k, $num_atoms, $masses, $charges, $vdwtype); 
	my($idx, $count, $opt_val, $bo, $data);

	open OUTDATA, "> $outName" or die "ERROR: Cannot write to $outName: $!\n";
	print OUTDATA "molmech conv noautobond nobonds preserve_q\n";
	#print OUTDATA "conp\n" if($has_cell);
	$idx = 0;
	for $i (sort { $a->{INDEX} <=> $b->{INDEX} } @{ $fData }) {
		next if(defined($efile) and ! exists($i->{OBSERVABLES}{energy}));
		$idx++;
		if (!defined($i->{DESCRP})) {
			$i->{DESCRP}[0] = "struct $idx";
		}
		print OUTDATA "name: $i->{DESCRP}[0] "; 
		print OUTDATA "index $i->{INDEX} " if(exists($i->{INDEX}));
		print OUTDATA "FILE: " . basename($i->{FILE}); 
		print OUTDATA " ENERGY: $i->{OBSERVABLES}{energy}{VAL}" if (defined($efile));
		print OUTDATA "\n";
		print OUTDATA "cell @{ $i->{CRYSTX} }\n" if (exists($i->{CRYSTX}) and $i->{CRYSTX});
		$num_atoms = scalar(keys %{ $i->{ATOMS} });
		print OUTDATA "cartesian ang $num_atoms\n";
		for $j (1 .. $num_atoms) {
			printf OUTDATA "%-10s core %11.5f %11.5f %11.5f %11.5f\n",$i->{ATOMS}{$j}{FFTYPE}, 
				$i->{ATOMS}{$j}{XCOORD},$i->{ATOMS}{$j}{YCOORD},$i->{ATOMS}{$j}{ZCOORD},$i->{ATOMS}{$j}{CHARGE};
			$masses->{ $i->{ATOMS}{$j}{FFTYPE} } = $i->{ATOMS}{$j}{MASS};
			$charges->{ $i->{ATOMS}{$j}{FFTYPE} } = $i->{ATOMS}{$j}{CHARGE};
		}
		if(defined($efile) and exists($i->{OBSERVABLES})) {
			print OUTDATA "observables\n";
			for $j (keys %{ $i->{OBSERVABLES} }) {
				print OUTDATA "$j $i->{OBSERVABLES}{$j}{UNIT}\n$i->{OBSERVABLES}{$j}{VAL} $i->{OBSERVABLES}{$j}{WEIGHT}\n";
			}
			print OUTDATA "end\n";
		}
		for $j (1 .. $num_atoms) {
			next if (!exists($i->{BONDS}{$j}) or $#{ $i->{BONDS}{$j} } == -1);
			for $k (0 .. $#{ $i->{BONDS}{$j} }) {
				next if ($i->{BONDS}{$j}[$k] < $j);
				$bo = "";
				if(exists($i->{ATOMS}{$j}{ORDER})) {
					$bo = $i->{ATOMS}{$j}{ORDER}[$k] if (exists($i->{ATOMS}{$j}{ORDER}));
					$bo = "single" if ($bo eq "1");
					$bo = "double" if ($bo eq "2");
					$bo = "triple" if ($bo eq "3");
					$bo = "resonant" if ($bo eq "ar");
				}
				print OUTDATA "connect $j $i->{BONDS}{$j}[$k] $bo\n";
			}
		}
		print OUTDATA "shift\n$elast eV\n" if(defined($efile) and $idx==1);
	}
	print OUTDATA "\nspecies\n";
	for $i (keys %{ $charges }) {
		print OUTDATA "$i core $charges->{$i} $i\n";
	}
	print OUTDATA "\nelement\n";
	for $i (keys %{ $masses }) {
		print OUTDATA "mass $i $masses->{$i}\n";
		$vdwtype = getGulpType($parms->{VDW}{$i}{$i}{1}{TYPE});
	}
	print OUTDATA "end\n";
	print OUTDATA "vary\nshift\nend\n" if (defined($efile));
	print OUTDATA lc($vdwtype) . "\n";
	for $i (keys %{ $parms->{VDW} }) {
		for $j (keys %{ $parms->{VDW} }) {
			next if (! exists($parms->{VDW}{$i}{$j}));
			$count = scalar(@{ $parms->{VDW}{$i}{$j}{1}{VALS} });
			print OUTDATA "$i $j @{ $parms->{VDW}{$i}{$j}{1}{VALS} }";
			print OUTDATA " 0 10";
			$opt_val = 1;
			$opt_val = 0 if ($opt_type == 2 and $i ne $j);
			$opt_val = 0 if ($opt_type == 3 and $i eq $j);
			for (1 .. $count) { print OUTDATA " $opt_val"; }
			print OUTDATA "\n";
		}
	}
	my $sname = basename($outName);
	$sname =~ s/\.\w+$//;
	print OUTDATA "dump ${sname}.grs";
	close OUTDATA;
}

sub loadenergies {
	my ($data, $infile) = @_;
	my ($i, $fileN, $str, $eng, $weight, $emin, $has_weights, $nfound, $elast);

	$has_weights = $nfound = 0;
	open INFILE, $infile or die "ERROR: Cannot open $infile: $!\n";
	while (<INFILE>) {
		chomp;
		if ($_ =~ /^(\S+)\s+(\-?\.?\d+\.?\d*E?\-?\+?\d*)\s*(\d*\.?\d*E?\-?\+?\d*)/i) {
			($str, $eng) = ($1, $2);
			$eng *= $efactor if(defined($efactor));
			$emin = $eng if(!defined($emin) or $eng < $emin);
			$weight = 1;
			if($3) {
				#$has_weights = 1;
				#$weight = $3;
			}
			for $i (@{ $data }) {
				$fileN = basename($i->{FILE});
				next if(exists($i->{OBSERVABLES}{energy}));
				if ($fileN =~ /(\-|\.)$str\./) {
					$i->{OBSERVABLES}{energy}{VAL} = $eng+0;
					$i->{OBSERVABLES}{energy}{WEIGHT} = $weight+0;
					$i->{OBSERVABLES}{energy}{UNIT} = "eV";
					$i->{INDEX}=$str;
					$nfound++;
					$elast = $eng;
					last;
				}
			}
		}
	}
	close INFILE;

	die "ERROR: No matching BGF files to entries in energy file!\n"
		if(!$nfound);

	for $i  (@{ $data }) {
		last if ($has_weights);
		next if (!exists($i->{OBSERVABLES}));
		$i->{OBSERVABLES}{energy}{WEIGHT} = exp(-($i->{OBSERVABLES}{energy}{VAL}-$emin)/.592126/23.06);
	}
	return $elast;
}

sub pruneVDW {
	my ($parms) = @_;
	my ($i, $j);

	for $i (keys %{ $parms->{VDW} }) {
		if (! exists($parms->{ATOMTYPES}{$i}{USED}) or ! $parms->{ATOMTYPES}{$i}{USED}) {
			delete $parms->{VDW}{$i};
			next;
		}
		for $j (keys %{ $parms->{VDW}{$i} }) {
			if (! exists($parms->{ATOMTYPES}{$j}{USED}) or ! $parms->{ATOMTYPES}{$j}{USED}) {
				delete $parms->{VDW}{$i}{$j};
			}
		}
	}
	for $i (keys %{ $parms->{VDW} }) {
		for $j (keys %{ $parms->{VDW} }) {
			next if ($i lt $j);
			if (! exists($parms->{VDW}{$i}{$j})) {
				$parms->{VDW}{$i}{$j}{1}{VALS} = getPairMix($parms->{VDW}{$i}{$i}{1}{VALS},
															$parms->{VDW}{$j}{$j}{1}{VALS},
															$parms->{VDW}{$i}{$i}{1}{TYPE});
				$parms->{VDW}{$i}{$j}{1}{TYPE} = $parms->{$i}{$i}{1}{TYPE};
			}
		}
   }
}

sub getPairMix {
	my ($iVals, $jVals, $pairType) = @_;
	my (@VALS, $mixType);

	$mixType = $FF->{PARMS}{mix_rule};

	if ($mixType eq "geometric") {
		$VALS[0] = sqrt($iVals->[0]*$jVals->[0]);
	} else {
		$VALS[0] = 0.5*($iVals->[0]+$jVals->[0]);

	}
	$VALS[1] = sqrt($iVals->[1]*$jVals->[1]);

	if ($pairType ne "LJ_12_10" and $#{ $iVals } == 2) {
		$VALS[2] = 0.5*($iVals->[2]+$jVals->[2]);
	}

	return \@VALS;
}

sub loadBGFS {
	my ($blist, $parms) = @_;
	my ($i, $count, $pStr, $pLen, $has_cell);

	$count = $has_cell = 0;
	for $i (@{ $blist }) {
		$pStr = "Getting atom information from $i->{FILE}...";
		print "${pStr}\r";
		($i->{ATOMS}, $i->{BONDS}, $i->{HEADERS}) = GetBGFFileInfo($i->{FILE},1);
		&checkAtomTypes($i->{ATOMS}, $parms);
		&AddMass($i->{ATOMS}, $parms);
		for (@{ $i->{HEADERS} }) {
			if ($_ =~ /^(DESCRP|CRYSTX)\s+(.*)$/) {
				@{ $i->{$1} } = split /\s+/,$2;
			}
		}
		$has_cell = 1 if($i->{CRYSTX});
		print "Done\r";
		$count++;
	}
	$pLen = length($pStr) + 10;
	printf "%-${pLen}s\n", "Getting atom information...read atom info from $count file(s)";
	return $has_cell;
}

sub checkAtomTypes {
	my ($atoms, $parms) = @_;
	my ($i, $ffType);

	for $i (keys %{ $atoms }) {
		$ffType = $atoms->{$i}{FFTYPE};
		die "ERROR: Force field type $ffType not found in forcefield(s). Aborting\n"
			if (! exists($parms->{ATOMTYPES}{$ffType}));
		$parms->{ATOMTYPES}{$ffType}{USED} = 1;
	}
}


sub init {
	my (%OPTS, $flist, $findCmd, $rec);

	getopt('bfseox',\%OPTS);
	($flist, $saveName, $efile, $fStr, $opt_type, $efactor) = ($OPTS{b},$OPTS{s},$OPTS{e}, $OPTS{f}, $OPTS{o}, $OPTS{x});
	for ($flist, $fStr) {
		die "usage: $0 -b bgf file(s) -f forcefield -e (eng file) -x (energy scaling factor to make eV) -s (save name) -o (opt_type:1=all 2=diagonal 3=offdiagonal)\n"
			if (! defined($_));
	}
	print "Initializing...";
	$findCmd = "find " . dirname($OPTS{b}) . " -name '" . basename($OPTS{b}) . "' -print | sort";
	if (-e $OPTS{b}) {
		$rec->{FILE} = $OPTS{b};
		push @{ $BGF }, $rec;
	} elsif (open(FINDCMD, "$findCmd |")) {
		while (<FINDCMD>) {
			$rec = ();
			chomp;
			$rec->{FILE} = $_;
			push @{ $BGF }, $rec;
		}
		close FINDCMD;
	} else {
		die "ERROR: No valid BGF files found while searching $OPTS{b}\n";
	}
	FileTester($efile) if (defined($efile));
	if (! defined($saveName)) {
		$saveName = basename($BGF->[0]{FILE});
		$saveName =~ s/\.\w+$//;
		$saveName .= "_gulp.gin";
	}
	$opt_type = 1 if (!defined($opt_type) or $opt_type !~ /(1|2|3)/);
	if ($opt_type =~ /(\d+)/) {
		$opt_type = $1;
	}
	print "Done\n";
}

sub getGulpType {
	my ($c2vdwType) = $_[0];
	my ($gulpType);

	if($c2vdwType =~ /LJ_6_12/) {
		$gulpType = "lennard epsilon kcal inter"
	}elsif($c2vdwType =~ /EXPO_6/) {
		$gulpType = "buckingham kcal";
	}elsif($c2vdwType =~ /MORSE/) {
		$gulpType = "morse kcal";
	}

   return $gulpType;
}

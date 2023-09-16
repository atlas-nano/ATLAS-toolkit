#!/usr/bin/perl -w
#
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use LAMMPS qw (ReadDataFile);
use General qw (FileTester LoadElements FindElementByMass);
use File::Basename;
use Getopt::Std qw(getopt);

sub numerically { ($a<=>$b); }

my ($lmp_dat, $lmp_ctl, $savePrefix);
my ($DATA, $PAIRS, %SHASH, $ELEMENTS);

$|++;
&init;
print "Parsing LAMMPS data file $lmp_dat...";
$DATA = ReadDataFile($lmp_dat);
&getTypesFromMasses($DATA);
print "Done\nParsing LAMMPS control file $lmp_ctl...";
&readCtlFile($DATA, $lmp_ctl, \%{ $PAIRS });
&fixPairCoeffs($DATA);
print "Done\nVerifying data...";
delete $DATA->{"ATOM COEFFS"};
&updateCoeffs($DATA);
print "Done\nCreating CP2K forcefield file ${savePrefix}.cp2k.ff...";
&writeCP2Kff($DATA, $savePrefix);
print "Done\nCreating AMBER crd file ${savePrefix}.crd...";
&createAMBERCrd($DATA, $savePrefix);
#print "Done\nCreating xPDB file ${savePrefix}.xpdb...";
#&writexPDB($DATA, $savePrefix);
print "Done\nCreating PSF connectivity file ${savePrefix}.psf...";
&createPSF($DATA, $savePrefix);
print "Done\nCreating CP2K control file ${savePrefix}.cp2k.inp...";
&createCP2Kinp($DATA, $savePrefix);
print "Done\n";

sub getMassStr {
	my ($masses) = $_[0];
	my ($i, $fftype, $massStr, $mass, $element);

	for $i (keys %{ $masses }) {
		$mass = $masses->{$i}{0};
		$fftype = $masses->{$i}{2};
		$element = FindElementByMass($mass);
		$massStr .= "    &KIND $fftype\n";
		$massStr .= "      ELEMENT $element\n";
		$massStr .= "      MASS $mass\n";
		$massStr .= "    &END KIND\n";
	}

	return $massStr;
}

sub createCP2Kinp {
	my ($data, $savePrefix) = @_;
	my ($xpdb_file, $psf_file, $lx, $ly, $lz, $ff_file, $val);
	my ($md_template, $inp_file, $instr, $masses, $ff, $crd_file);

	$masses = getMassStr($data->{MASS});
	$xpdb_file = "${savePrefix}.xpdb";
	$psf_file = "${savePrefix}.psf";
	$crd_file = "${savePrefix}.crd";
	$ff_file = "${savePrefix}.cp2k.ff";
	$lx = $data->{BOX}{X}{hi}-$data->{BOX}{X}{lo};
	$ly = $data->{BOX}{Y}{hi}-$data->{BOX}{Y}{lo};
	$lz = $data->{BOX}{Z}{hi}-$data->{BOX}{Z}{lo};

	$md_template = "$Bin/cp2k.fist.template.in";
	$inp_file = "${savePrefix}.inp";
	
	open FF, $ff_file or die "ERROR: Cannot open $ff_file: $!\n";
	$/ = undef;
	$ff = <FF>;
	close FF;

	open CP2Kinp, "> $inp_file" or die "ERROR: Cannot open $inp_file: $!\n";
	open TEMPLATE, "$md_template" or die "ERROR: Cannot open $md_template: $!\n";
	while (<TEMPLATE>) {
		chomp;
		$instr = $_;
		while ($instr =~ /(\S+)_here/g) {
			$val = eval('$' . $1);
			$instr =~ s/${1}_here/$val/;
		}

		print CP2Kinp "$instr\n";
	}
	close TEMPLATE or die "ERROR: Cannot close $md_template: $!\n";
	close CP2Kinp or die "ERROR: Cannot close $inp_file: $!\n";
}

sub createPSF {
	my ($data, $savePrefix) = @_;
	my ($i, $j, $k, $d, $c, $saveName, $header, $num, $tot);
	my ($ffid, $fftype, $molid, $charge, $mass, $nline, $resid);

	$saveName = "${savePrefix}.psf";
	open PSFFILE, "> $saveName" or die "ERROR: Cannot care $saveName: $!\n";
	printf PSFFILE "PSF EXT\n\n";
	printf PSFFILE "%10d !NTITLE\n    Created by $0\n\n",1;
	#atoms
	printf PSFFILE "%10d !NATOM\n", $data->{TOTAL}{ATOMS};
	for $i (1 .. $data->{TOTAL}{ATOMS}) {
		$ffid = $data->{ATOMS}{$i}{1};
		$fftype = "Du";
		$fftype = $data->{MASS}{$ffid}{2} if (exists($data->{MASS}{$ffid}{2}));
		$mass = $data->{MASS}{$ffid}{0};
		$molid = $data->{ATOMS}{$i}{0};
		$resid = "MOL${molid}";
		$charge = $data->{ATOMS}{$i}{2};
		printf PSFFILE "%10d %-9s%-5d%8s     %-4s%9s%11.6f%14.4f%12d\n",$i,$resid,$molid,"RES",$fftype,
															   $fftype,$charge,$mass,0;
	}
	printf PSFFILE "\n\n";

	#bonds
	for $j ("BONDS","ANGLES","DIHEDRALS","INVERSIONS") {
		$nline = 2;
		if ($j eq "BONDS") {
			$header = "NBOND";
			$num = 2;
			$nline = 4;
		} elsif($j eq "ANGLES") {
			$header = "NTHETA";
			$num = $nline = 3;
		} elsif ($j eq "DIHEDRALS") {
			$header = "NPHI";
			$num = 4;
		} elsif ($j eq "INVERSIONS") {
			$header = "NIMPHI";
			$num = 4;
		}
		$tot = 0;
		$tot = $data->{TOTAL}{$j} if (exists($data->{TOTAL}{$j}));
		printf PSFFILE "%10d !%s: %s\n", $tot, $header, lc $j;
		$c = 0;
		next if (!exists($data->{TOTAL}{$j}) or $data->{TOTAL}{$j} == 0);
		for	$i (1 .. $data->{TOTAL}{$j}) {
			$c ++;
			for $k (1 .. $num) {
				printf PSFFILE "%10d",$data->{$j}{$i}{$k};
			}
			printf PSFFILE "\n" if($c % $nline == 0);
		}
		printf PSFFILE "\n\n";
	}

	close PSFFILE or die "ERROR: Cannot close $saveName: $!\n";

}

sub writexPDB {
	my ($data, $savePrefix) = @_;
	my ($i, $saveName, $fftype, $ffid, $molid, $x, $y, $z, $charge);

	$saveName = "${savePrefix}.xpdb";
	open xPDBFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
	for $i (1 .. $data->{TOTAL}{ATOMS}) {
		$ffid = $data->{ATOMS}{$i}{1};
		$molid = $data->{ATOMS}{$i}{0};
		$fftype = "Du";
		$fftype = $data->{MASS}{$ffid}{2} if (exists($data->{MASS}{$ffid}{2}));
		$fftype = substr($fftype,0,3) if (length($fftype)>3);
		$x = $data->{ATOMS}{$i}{3};
		$y = $data->{ATOMS}{$i}{4};
		$z = $data->{ATOMS}{$i}{5};
		$charge = $data->{ATOMS}{$i}{2};
		printf xPDBFILE "%-6s%5d%4s  %3s  %4d %11.3f%8.3f%8.3f%6.2f%6.2f%4s%2s%8.5f\n","ATOM",$i,
										$fftype,"RES",$molid,$x,$y,$z,1,0,"Se","0",$charge;
		
	}
	close xPDBFILE;
}

sub createAMBERCrd {
	my ($data, $savePrefix) = @_;
	my ($i, $d, $c, $saveName);

	$saveName = "${savePrefix}.crd";

	open CRDFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
	printf CRDFILE "$savePrefix\n%5d\n", $data->{TOTAL}{ATOMS};
	$c = 0;
	for $i (1 .. $data->{TOTAL}{ATOMS}) {
		for $d (3 .. 5) {
			$c++;
			printf CRDFILE "%12.7f", $data->{ATOMS}{$i}{$d};
			printf CRDFILE "\n" if ($c % 6 == 0);
		}
	}
	printf CRDFILE "\n" if ($c % 6 > 0);
#	for $i (qw/X Y Z/) {
#		printf CRDFILE "%12.3f",$data->{BOX}{$i}{hi}-$data->{BOX}{$i}{lo};
#	}
#	for $i (1 .. 3) {
#		printf CRDFILE "%12.3f",90;
#	}
#	printf CRDFILE "\n";
	close CRDFILE or die "ERROR: Cannot close $saveName: $!\n";
}

sub writeCP2Kff {
	my ($data, $prefix) = @_;
	my ($i, $j, $k, $l, $outfile, $ns, $units, $c, $tot);

	$outfile = "${prefix}.cp2k.ff";
	open FF, "> $outfile" or die "ERROR: Cannot write to $outfile: $!\n";
	for $i (keys %{ $data } ) {
		next if ($i !~ / COEFFS/);
		$c = 0;
		$tot = scalar(keys %{ $data->{$i} });
		for $j (keys %{ $data->{$i} }) {
			$ns = 5;
			$c++;
			if ($c == 1 or !exists($data->{$i}{$j}{CP2K}{section})) {
				printf FF "\n%${ns}s&%-s","",$data->{$i}{$j}{CP2K}{header};
			}
			if (exists($data->{$i}{$j}{CP2K}{section})) {
				$ns = 7;
				printf FF "\n%${ns}s&%-s","",$data->{$i}{$j}{CP2K}{section};
			}
			$ns += 2;
			printf FF "\n%${ns}s%-12s @{ $data->{$i}{$j}{ATOMS} }","","ATOMS";
			printf FF "\n%${ns}s%-12s %-s","","KIND", $data->{$i}{$j}{CP2K}{kind}
				if(exists($data->{$i}{$j}{CP2K}{kind}));
			for $k (0 .. $data->{$i}{$j}{CP2K}{nparm}-1) {
				next if (!exists($data->{$i}{$j}{CP2K}{labels}{$k}));
				$units = "";
				$units = "[$data->{$i}{$j}{CP2K}{units}{$k}]" if ($data->{$i}{$j}{CP2K}{units}{$k});
				printf FF "\n%${ns}s%-12s $units %s","",$data->{$i}{$j}{CP2K}{labels}{$k},
													  $data->{$i}{$j}{$k};
				if(exists($data->{$i}{$j}{CP2K}{sameline}{$k})) {
					for $l (sort numerically keys %{ $data->{$i}{$j}{CP2K}{sameline}{$k} }) {
						$units = "";
						$units = "[$data->{$i}{$j}{CP2K}{units}{$l}]" if ($data->{$i}{$j}{CP2K}{units}{$l});
						printf FF " $units %s",$data->{$i}{$j}{$l};
					}
				}
			}
			$ns -= 2;
			if (exists($data->{$i}{$j}{CP2K}{section})) {
				printf FF "\n%${ns}s&END %-s","",$data->{$i}{$j}{CP2K}{section};
				$ns -=2;
			}
			if (!exists($data->{$i}{$j}{CP2K}{section}) or $c == $tot) {
				printf FF "\n%${ns}s&END %-s\n","",$data->{$i}{$j}{CP2K}{header};
			}
		}
	}
	close FF or die "ERROR: Cannot close $outfile: $!\n";
}

sub updateCoeffs {
	my ($data) = $_[0];
	my ($i, $j, $k, $l, $ptype, $vtype, $valid);

	for $i (keys %{ $data }) {
		next if ($i !~ /(\S+) COEFFS/);
		$ptype = lc $1;
		die "ERROR: No valid CP2K entry found for $ptype" if (!exists($SHASH{$ptype}));
		for $j (keys %{ $data->{$i} }) {
			$vtype = $data->{$i}{$j}{TYPE};
			$vtype =~ s/\//_/g;
			$valid = 0;
			for $k (keys %{ $SHASH{$ptype} }) {
				if ($vtype =~ /^$k/) {
					$valid = 1;
					$data->{$i}{$j}{CP2K} = $SHASH{$ptype}{$k};
					last;
				}
			}
			die "ERROR: Cannot find CP2K entry for $vtype while searching $ptype!\n"
				if (!$valid);
			$k = getLastIndex($data->{$i}{$j});
			$l = 0;
			$valid = 0;
			while($l < $k) {
				if ($data->{$i}{$j}{$l} eq "#") {
					$valid = 1;
					$l++;
					last;
				}
				$l++;
			}
			die "ERROR: Cannot determine atoms for $vtype while searching $ptype!\n"
				if (! $valid);
			for ($l .. $k) {
				push @{ $data->{$i}{$j}{ATOMS} }, $data->{$i}{$j}{$_};
			}
		}
	}
}

sub getTypesFromMasses {
	my ($data) = $_[0];
	my ($i);

	for $i (keys %{ $data->{MASS} }) {
		die "ERROR: No comments found in data file while searching for fftype of $i!\n"
			if (!exists($data->{MASS}{$i}{2}));
		$data->{TYPES}{$i} = $data->{MASS}{$i}{2};
	}
}

sub fixPairCoeffs {
	my ($data) = $_[0];
	my ($i, $j, $pairs, $vtype, $plist, $c, $t);

	$pairs = $data->{"PAIR COEFFS"};
	for $i (keys %{ $pairs }) {
		if(defined($pairs->{$i}{TYPE})) {
			$vtype = $pairs->{$i}{TYPE};
			last;
		}
	}
	die "ERROR: Cannot determine the pair type!\n" if (! defined($vtype));

	$pairs = $data->{"PAIR COEFFS"};
	for $i (keys %{ $pairs }) {
		if (!exists($pairs->{$i}{I})) {
			$pairs->{$i}{I} = $pairs->{$i}{J} = $i;
		}
		$pairs->{$i}{TYPE} = $vtype if (!defined($pairs->{$i}{TYPE}));
		$plist->{ $pairs->{$i}{I} }{ $pairs->{$i}{J} } = \%{ $pairs->{$i} };
	}

	$c = getLastIndex($pairs)+1;
	$t = getLastIndex($data->{MASS});

	for $i (1 .. $t) {
		for $j ($i .. $t) {
			next if (exists($plist->{$i}{$j}) or exists($plist->{$j}{$i}));
			$pairs->{$c+1} = pairMix($plist->{$i}{$i}, $plist->{$j}{$j});
		}
	}
	for $i (keys %{ $pairs }) {
		$c = getLastIndex($pairs->{$i});
		($pairs->{$i}{$c+1}, $pairs->{$i}{$c+2}, $pairs->{$i}{$c+3}) = 
			("#",$data->{TYPES}{$pairs->{$i}{I}},$data->{TYPES}{$pairs->{$i}{J}});
	}
}

sub pairMix {
	my ($ii, $jj) = @_;
	my ($i, $ij);

	for $i (grep {/^\d+$/} keys %{ $ii }) {
		last if ($ii->{$i} =~ /#/);
		$ij->{$i} = sqrt($ii->{$i}*$jj->{$i});
	}

	$ij->{I} = $ii->{I};
	$ij->{J} = $jj->{J};
	$ij->{TYPE} = $ii->{TYPE};

	return $ij;
}

sub readCtlFile {
	my ($data, $infile, $pair) = @_;
	my ($tmp, $i, $j, $vtype, $c);

	open CTLFILE, $infile or die "ERROR: Cannot open $infile: $!\n";
	while (<CTLFILE>) {
		chomp;
		next if ($_ !~ /(coeff|style)/i);
		if ($_ =~ /^(\S+)_coeff\s*(.+)$/) {
			$vtype = uc $1;
			$c = $2;
			next if (!exists($DATA->{"${vtype}S"}) and $vtype !~ /(PAIR|ATOM)/i);
			@{ $tmp } = split /\s+/, $c;
			if ($vtype =~ /PAIR/) {
				($i, $j) = ($tmp->[0], $tmp->[1]);
				$tmp->[0] = getLastIndex($data->{"$vtype COEFFS"})+1;
				undef $c;
				if ($tmp->[2] !~ /\-?\d+\.?\d*E?\-?\+?\d*/i) {
					$c = $tmp->[2]; 
					splice @{ $tmp }, 2, 1;
				}
				$tmp->[1] = $c;
				$data->{"$vtype COEFFS"}{$tmp->[0]}{J} = $j;
				$data->{"$vtype COEFFS"}{$tmp->[0]}{I} = $i;
			}
			$c = 0;
			for $i (2 .. $#{ $tmp }) {
				$data->{"$vtype COEFFS"}{$tmp->[0]}{$c} = $tmp->[$i];
				$c++;
			}
			$data->{"$vtype COEFFS"}{$tmp->[0]}{TYPE} = $tmp->[1];
		} elsif ($_ =~ /^(\S+)_style\s+(.+)$/) {
			$vtype = uc $1;
			$c = $2;
			next if (!exists($DATA->{"${vtype}S"}) and $vtype !~ /(PAIR|ATOM)/i);
			@{ $tmp } = split /\s+/, $c;
			next if ($tmp->[0] =~ /hybrid/i);
			for $j (keys %{ $data->{"$vtype COEFFS"} }) {
				$data->{"$vtype COEFFS"}{$j}{TYPE} = $tmp->[0];
			}
		}
	}
	close CTLFILE;
}

sub getLastIndex {
	my ($h) = $_[0];
	my ($tmp);

	@{ $tmp } = grep {/^\d+$/} keys %{ $h };
	@{ $tmp } = sort numerically @{ $tmp };
	return $tmp->[$#{ $tmp }];
}

sub init {
	my (%OPTS);

	getopt('dcs',\%OPTS);

	die "usage: $0 -d lmp_data_file -c lmp_control_file -s (cp2k_save_prefix)\n"
		if (!exists($OPTS{d}) or ! exists($OPTS{c}));
	($lmp_dat, $lmp_ctl, $savePrefix) = ($OPTS{d}, $OPTS{c}, $OPTS{s});

	print "Initializing...";
    FileTester($lmp_dat);
    FileTester($lmp_ctl);
    if (! defined($savePrefix)) {
		$savePrefix = basename($lmp_dat);
		$savePrefix =~ s/data\.//;
    }
	my ($converter) = "$Bin/dat/lammpsParms2CP2K.perldata";
	FileTester($converter);
	eval `cat $converter` or die "ERROR: Cannot recreate data in file $converter: $! $@\n";
	$ELEMENTS = &LoadElements;
	print "Done\n";
}


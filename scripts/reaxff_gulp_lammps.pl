#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester LoadElements);
use ManipAtoms qw(FindElement);

my ($inFile, $outFile, $infileType); 
my ($DATA, $EngVals, $readFF, $writeFF, $ELEMENTS);

sub usage;
sub init;
sub readLammpsReaxFF;
sub readGulpReaxFF;
sub writeLammpsReaxFF;
sub writeGulpReaxFF;
sub getMass;
sub printlines;
sub scaleEnergies;
sub getValParent;

&init;
print "Done\nReading " . uc($infileType) . " forcefield $inFile...";
$readFF->($inFile,\%{ $DATA });
&scaleEnergies($DATA,$EngVals);
print "Done\nConverting FF to $outFile...";
$writeFF->($DATA,$outFile);
print "Done\n";

sub readLammpsReaxFF {
	my ($fle, $data) = @_;
}

sub writeLammpsReaxFF {
	my ($data, $ofile) = @_;
	my ($i, $j, $k, $l, $vars, $pen);
	my ($idata, $types, $species, $curr, $scale, $index);

	open OUTFILE, "> $ofile" or die "ERROR: Cannot write to $ofile: $!\n";
	print OUTFILE "LAMMPS ReaxFF forcefield generated from $inFile\n 39\n";

	$data->{header}{reaxFFvdwcutoff}{VALS}[0] = 10 if (!exists($data->{header}{reaxFFvdwcutoff}));
	$data->{header}{reaxff0_vdw}{VALS}[0] = 1.5591 if (!exists($data->{header}{reaxff0_vdw}));
	@{ $idata } = ($data->{header}{reaxff0_bond}{VALS}[0],$data->{header}{reaxff0_bond}{VALS}[1],
					$data->{header}{reaxff0_val_angle},3,6.5,50,
					$data->{header}{reaxff0_over}{VALS}[2],15,$data->{header}{reaxff0_over}{VALS}[3],
					$data->{header}{reaxff0_over}{VALS}[4],-24.671,0,$data->{header}{reaxFFvdwcutoff}{VALS}[0],
					2.8793,$data->{header}{reaxff0_valence}{VALS}[0],$data->{header}{reaxff0_lonepair}{VALS}[0],
					$data->{header}{reaxff0_valence}{VALS}[2],$data->{header}{reaxff0_valence}{VALS}[3],6.1431,
					$data->{header}{reaxff0_penalty}{VALS}[0],$data->{header}{reaxff0_penalty}{VALS}[1],
					$data->{header}{reaxff0_penalty}{VALS}[2],-2.4837,$data->{header}{reaxff0_torsion}{VALS}[0],
					$data->{header}{reaxff0_torsion}{VALS}[1],$data->{header}{reaxff0_torsion}{VALS}[2],-1.2327,
					$data->{header}{reaxff0_torsion}{VALS}[3],$data->{header}{reaxff0_vdw}{VALS}[0],
					$data->{header}{reaxFFtol}{VALS}[0]/.01,$data->{header}{reaxff0_val_ang_con_par},
					$data->{header}{reaxff0_over}{VALS}[1],$data->{header}{reaxff0_over}{VALS}[0],
					$data->{header}{reaxff0_valence}{VALS}[1],.5,20,5,0,$data->{header}{reaxff0_val_ang_con_par2});
	&printlines(\*OUTFILE,$idata,"10.4f",1,0);

	#species
	@{ $types } = sort { $a cmp $b } keys %{ $data->{reaxff1_radii} };
	$i = 1;
	for (0 .. $#{ $types }) {
		if ($types->[$_] eq "name") {
			delete $types->[$_];
		} else {
			$data->{reaxff1_radii}{$_}{index} = $i;
			$species->{$types->[$_]}{index} = $i;
			$species->{$types->[$_]}{mass} = getMass($types->[$_],$ELEMENTS);
			$i++;
		}
	}
	$species->{X}{index} = 0;
	$i--;
	printf OUTFILE "%3d    ! Nr of atoms; cov.r; valency;a.m;Rvdw;Evdw;gammaEEM;cov.r2;#\n",$i;
	print OUTFILE "            alfa;gammavdW;valency;Eunder;Eover;chiEEM;etaEEM;n.u.\n";
	print OUTFILE "            cov r3;Elp;Heat inc.;n.u.;n.u.;n.u.;n.u.\n";
	print OUTFILE "            ov/un;val1;n.u.;val3,vval4\n";

	for $i (@{ $types }) {
		printf OUTFILE " %-2s",$i;

		$data->{reaxff1_under}{$i}{scale} = 23.06 if ($data->{reaxff1_under}{$i}{opts} !~ /kcal/);
		$data->{reaxff1_under}{$i}{scale} = 1 if ($data->{reaxff1_under}{$i}{opts} =~ /kcal/);
		@{ $idata } = ($data->{reaxff1_radii}{$i}{VALS}[0],$data->{reaxff1_valence}{$i}{VALS}[0],
						$species->{$i}{mass},$data->{reaxff1_morse}{$i}{VALS}[2],$data->{reaxff1_morse}{$i}{VALS}[1],
						$data->{reaxff_gamma}{$i}{VALS}[0],$data->{reaxff1_radii}{$i}{VALS}[1],
						$data->{reaxff1_valence}{$i}{VALS}[2],$data->{reaxff1_morse}{$i}{VALS}[0],
						$data->{reaxff1_morse}{$i}{VALS}[3],$data->{reaxff1_valence}{$i}{VALS}[3],
						$data->{reaxff1_under}{$i}{VALS}[0],100,$data->{reaxff_chi}{$i}{VALS}[0],
						$data->{reaxff_mu}{$i}{VALS}[0],2,$data->{reaxff1_radii}{$i}{VALS}[2],
						$data->{reaxff1_lonepair}{$i}{VALS}[1],-2.37,$data->{reaxff1_over}{$i}{VALS}[1],
						$data->{reaxff1_over}{$i}{VALS}[0],$data->{reaxff1_over}{$i}{VALS}[2],0,0,
						$data->{reaxff1_over}{$i}{VALS}[3],$data->{reaxff1_angle}{$i}{VALS}[0],1.0183,
						$data->{reaxff1_valence}{$i}{VALS}[1],$data->{reaxff1_angle}{$i}{VALS}[1],0,0,0);
		&printlines(\*OUTFILE,$idata,"9.4f",8,3);
	}

	push @{ $types }, "X";
	#bonds
	$i = 0;
	&getCount($data->{reaxff2_bond},\$i); #get bond count... recursive.. elegant!
	printf OUTFILE " %-2d       ! Nr of bonds; Edis1;LPpen;n.u.;pbe1;pbo5;13corr;pbo6\n",$i;
	print OUTFILE "                         pbe2;pbo3;pbo4;n.u.;pbo1;pbo2;ovcorr\n";
	for $i (@{ $types }) {
		for $j (@{ $types }) {
			next if (!exists($data->{reaxff2_bond}{$i}{$j}));
			@{ $index } = ($species->{$i}{index},$species->{$j}{index});
			@{ $index } = reverse @{ $index } if ($index->[0]>$index->[1]);
			@{ $vars } = (0,0);
			$curr = "";
			$curr = $data->{reaxff2_bo}{$i}{$j}{opts} if (exists($data->{reaxff2_bo}{$i}{$j}{opts}));;
			if ($curr =~ /bo13/ and $curr =~ /over/) {
				@{ $vars } = (1,1);
			} elsif ($curr =~ /bo13/) {
				@{ $vars } = (1,0);
			}
			$vars->[2] = 0;
			$vars->[2] = 1 if (exists($data->{reaxff2_pen}{$i}{$j}{opts}) and $data->{reaxff2_pen}{$i}{$j}{opts} =~ /kcal/);
			printf OUTFILE "%3d%3d",@{ $index };
			@{ $idata } = ($data->{reaxff2_bond}{$i}{$j}{VALS}[0],$data->{reaxff2_bond}{$i}{$j}{VALS}[1],
							$data->{reaxff2_bond}{$i}{$j}{VALS}[2],$data->{reaxff2_bond}{$i}{$j}{VALS}[3],
							$data->{reaxff2_bo}{$i}{$j}{VALS}[4],$vars->[0],
							$data->{reaxff2_bo}{$i}{$j}{VALS}[5],$data->{reaxff2_over}{$i}{$j}{VALS}[0],
							$data->{reaxff2_bond}{$i}{$j}{VALS}[4],$data->{reaxff2_bo}{$i}{$j}{VALS}[2],
							$data->{reaxff2_bo}{$i}{$j}{VALS}[3],1,$data->{reaxff2_bo}{$i}{$j}{VALS}[0],
							$data->{reaxff2_bo}{$i}{$j}{VALS}[1],$vars->[1],$vars->[2]);
			&printlines(\*OUTFILE,$idata,"9.4f",8,6);
		}
	}

	#off diagonal morse
	$i = 0;
	&getCount($data->{reaxff2_morse}, \$i);
	printf OUTFILE "%3d    ! Nr of off-diagonal terms; Ediss;Ro;gamma;rsigma;rpi;rpi2\n",$i;
	for $i (@{ $types }) {
		for $j (@{ $types }) {
			next if (!exists($data->{reaxff2_morse}{$i}{$j}));
			@{ $index } = ($species->{$i}{index},$species->{$j}{index} );
			@{ $index } = reverse @{ $index } if ($index->[0]>$index->[1]);
			printf OUTFILE "%3d%3d",@{ $index };
			@{ $idata } = ($data->{reaxff2_morse}{$i}{$j}{VALS}[0],$data->{reaxff2_morse}{$i}{$j}{VALS}[2],
							$data->{reaxff2_morse}{$i}{$j}{VALS}[1],$data->{reaxff2_morse}{$i}{$j}{VALS}[3],
							$data->{reaxff2_morse}{$i}{$j}{VALS}[4],$data->{reaxff2_morse}{$i}{$j}{VALS}[5]);
			&printlines(\*OUTFILE,$idata,"9.4f",6,6);
		}
	}

	#angles
	$i = 0;
	&getCount($data->{reaxff3_angle},\$i);
	printf OUTFILE " %-2d   ! Nr of angles;at1;at2;at3;Thetao,o;ka;kb;pv1;pv2;val(bo)\n",$i;
	for $i (@{ $types }) {
		for $j (@{ $types }) {
			for $k (@{ $types }) {
				next if (!exists($data->{reaxff3_angle}{$i}{$j}{$k}));
				@{ $index } = ($species->{$j}{index},$species->{$i}{index},$species->{$k}{index});
				@{ $index } = reverse @{ $index } if ($index->[0]>$index->[2]);
				$pen = 0;
				$pen = $data->{reaxff3_penalty}{$i}{$j}{$k}{VALS}[0] if (exists($data->{reaxff3_penalty}{$i}{$j}{$k}{VALS}));
				printf OUTFILE "%3d%3d%3d",@{ $index };
				$data->{reaxff3_conjugation}{$i}{$j}{$k}{VALS}[0] = 0 if (!exists($data->{reaxff3_conjugation}{$i}{$j}{$k}));
				@{ $idata } = ($data->{reaxff3_angle}{$i}{$j}{$k}{VALS}[0],$data->{reaxff3_angle}{$i}{$j}{$k}{VALS}[1],
								$data->{reaxff3_angle}{$i}{$j}{$k}{VALS}[2],$data->{reaxff3_conjugation}{$i}{$j}{$k}{VALS}[0],
								$data->{reaxff3_angle}{$i}{$j}{$k}{VALS}[4],$pen,
								$data->{reaxff3_angle}{$i}{$j}{$k}{VALS}[3]);
				&printlines(\*OUTFILE,$idata,"9.4f",7,9);
			}
		}
	}

	#torsion
	$i = 0;
	&getCount($data->{reaxff4_torsion},\$i);
	printf OUTFILE " %-2s    ! Nr of torsions;at1;at2;at3;at4;;V1;V2;V3;V2(BO);vconj;n.u;n\n",$i;
	for $i (@{ $types }) {
		for $j (@{ $types }) {
			for $k (@{ $types }) {
				for $l (@{ $types }) {
					next if (!exists($data->{reaxff4_torsion}{$i}{$j}{$k}{$l}));
					@{ $index } = ($species->{$i}{index},$species->{$j}{index},$species->{$k}{index},$species->{$l}{index});
					@{ $index } = reverse @{ $index } if ($index->[0]>$index->[3]);
					printf OUTFILE "%3d%3d%3d%3d",@{ $index };
					$curr = $data->{reaxff4_torsion}{$i}{$j}{$k}{$l}{VALS};
					@{ $idata } = ($curr->[0],$curr->[1],$curr->[2],$curr->[3],$curr->[4],0,0);
					&printlines(\*OUTFILE,$idata,"9.4f",7,12);
				}
			}
		}
	}

	#hbonds
	$i = 0;
	&getCount($data->{reaxff3_hbond},\$i);
	printf OUTFILE " %2d    ! Nr of hydrogen bonds;at1;at2;at3;Rhb;Dehb;vhb1\n",$i;
	for $i (@{ $types }) {
		for $j (@{ $types }) {
			for $k (@{ $types }) {
				next if (!exists($data->{reaxff3_hbond}{$i}{$j}{$k}));
				printf OUTFILE "%3d%3d%3d",$species->{$j}{index},$species->{$i}{index},$species->{$k}{index};
				$curr = $data->{reaxff3_hbond}{$i}{$j}{$k}{VALS};
				@{ $idata } = ($curr->[0],$curr->[1],$curr->[2],$curr->[3]);
				&printlines(\*OUTFILE,$idata,"9.4f",4,9);
			}
		}
	}
	close OUTFILE;
	print "";
}

sub writeGulpReaxFF {
}

sub readGulpReaxFF {
	my ($fle, $data) = @_;
	my ($header_parms, $type_headers, $curr, $opts); 
	my ($isval, $continue, $parent, $key, $val, $count, $list, $i);

	$header_parms->[0] = "reaxFFvdwcutoff|reaxFFqcutoff|reaxFFtol|reaxff0_bond|reaxff0_over";
	$header_parms->[1] = "reaxff0_valence|reaxff0_penalty|reaxff0_torsion|reaxff0_vdw|reaxff0_lonepair";
	$type_headers->[0] = "reaxff1_radii|reaxff1_valence|reaxff1_over|reaxff1_under|reaxff1_lonepair"; 
	$type_headers->[1] = "reaxff1_angle|reaxff1_morse|reaxff_chi|reaxff_mu|reaxff_gamma|reaxff2_bond";
	$type_headers->[2] = "reaxff2_bo|reaxff2_over|reaxff2_morse|reaxff3_angle|reaxff3_penalty"; 
	$type_headers->[3] = "reaxff3_conjugation|reaxff3_hbond|reaxff4_torsion";

	$isval = $continue = 0;

	open INFILE, $fle or die "ERROR: Cannot open $fle: $!\n";
	while (<INFILE>) {
		chomp;
		next if ($_ =~ /^#/);
		if ($_ =~ /($header_parms->[0]|$header_parms->[1])\s+(.+)$/) {
			$isval = 0;
			$curr = \%{ $data->{header}{$1} };
			$val = $2;
			$continue = 0;
			$continue = 1 if ($val =~ /&/);
			$val =~ s/&//;
			@{ $curr->{VALS} } = split /\s+/,$val;
			$curr->{parent} = "header";
			$curr->{opts} = "";
		}elsif(!$isval and $continue and $_ =~ /(\-?\d+.*)/) {
			push @{ $curr->{VALS} }, split /\s+/,$1;
			$continue = 0;
			$continue = 1 if ($1 =~ /&/);
		}elsif ($_ =~ /($type_headers->[0]|$type_headers->[1]|$type_headers->[2]|$type_headers->[3])\s*(.*)$/) {
			$parent = \%{ $data->{$1} };
			$parent->{name} = $1;
			$opts = "";
			$opts = $2 if ($2);
			$isval = 1;
		}elsif ($isval and ! $continue and $_ =~ /core/) {
			$curr = \%{ $parent };
			$key = "";
			$continue = 0;
			$continue = 1 if ($_ =~ /&/);
			$_ =~ s/&//;
			while ($_ =~ /(\S+)\s+core\s+/g) {
				$curr = \%{ $curr->{$1} };
				$key .= "$1 ";
			}
			chop $key;
			$_ =~ /.*core\s+(.+)$/;
			@{ $curr->{VALS} } = split /\s+/,$1;
			$curr->{opts} = $opts;
			$curr->{parent} = $parent->{name};
			$curr->{key} = $key;
		}elsif ($isval and $continue) {
			$continue = 0;
			$continue = 1 if ($_ =~ /&/);
			$_ =~ s/&//;
			while ($_ =~ /(\S+)/g) {
				push @{ $curr->{VALS} }, $1;
			}
		}

	}
	close INFILE;

	$count = 0;
	$list = ();
	&getValParent($data->{reaxff3_conjugation},\%{ $list },\$count);
	$data->{header}{reaxff0_val_angle} = $data->{header}{reaxff0_val_ang_con_par}=$data->{header}{reaxff0_val_ang_con_par2}=0;
	for $i (keys %{ $list }) {
		$curr = $list->{$i};
		$data->{header}{reaxff0_val_angle} = $curr->{VALS}[1] if ($curr->{VALS}[1]>0);
		$data->{header}{reaxff0_val_ang_con_par} = $curr->{VALS}[3] if ($curr->{VALS}[3]>0);
		$data->{header}{reaxff0_val_ang_con_par2} = $curr->{VALS}[2] if ($curr->{VALS}[2]>0);
	}
	$count = 0;
	$list = ();
	&getValParent($data->{reaxff1_valence},\%{ $list },\$count);
	$data->{header}{reaxff0_triple_bond_stab} = 0;
	for $i (keys %{ $list }) {
		$curr = $list->{$i};
		$data->{header}{reaxff0_triple_bond_stab} = $curr->{VALS}[0] if ($curr->{VALS}[0]>0);
	}
}

sub getCount {
	my ($data, $index) = @_;
	my ($i);

	for $i (keys %{ $data }) {
		if ($i eq "VALS") {
			$$index++;
		} elsif(ref($data->{$i}) eq "HASH") {
			&getCount($data->{$i},$index);
		}
	}
}

sub getMass {
	my ($ele, $list) = @_;
	my ($i,$mass);

	$mass = 0.0;
	for $i (keys %{ $list }) {
		if ($list->{$i}{SYMBOL} =~ /^${ele}$/i) {
			$mass = $list->{$i}{MASS};
		}
	}

	return $mass;
}

sub printlines {
	my ($outfile, $data, $fmt, $num, $loffset) = @_;
	my ($i, $j);

	$j = 0;
	while ($j <= $#{ $data }) {
		printf $outfile "%${loffset}s"," " if ($j > 1 and $loffset > 0);
		for $i (1 .. $num) {
			printf $outfile "%$fmt",$data->[$j];
			$j++;
		}
		printf $outfile "\n";
	}

}

sub scaleEnergies {
	my ($data, $Eng) = @_;
	my ($i, $j, $index, $list, $count);

	$list = ();
	$count = 0;
	&getValParent($data,\%{ $list }, \$count);
	for $i (keys %{ $list }) {
		next if (exists($list->{$i}{opts}) and $list->{$i}{opts} =~ /kcal/);
		next if (!exists($Eng->{ $list->{$i}{parent} }));
		@{ $index } = @{ $Eng->{ $list->{$i}{parent} } };
		for $j (@{ $index }) {
			$list->{$i}{VALS}[$j] /= 4.3364432032e-2;
		}
	}

}

sub getValParent {
	my ($valList, $VList, $counter) = @_;
	my ($i);

	for $i (keys %{ $valList }) {
		next if (ref($valList->{$i}) ne "HASH");
		if (exists($valList->{$i}{VALS})) {
			$VList->{$$counter} = \%{ $valList->{$i} };
			$$counter++;
		} else {
			&getValParent($valList->{$i}, $VList, $counter);
		}
	}
}

sub init {
	my (%OPTS);
	getopt('iot',\%OPTS);
	($inFile, $outFile, $infileType) = ($OPTS{i},$OPTS{o},$OPTS{t});
	die &usage if (! defined($inFile));

	print "Initializing...";
	FileTester($inFile);
	if(!defined($infileType) || $infileType !~ /lammps|gulp/) {
		$infileType = "lammps";
		my $tstr = `grep -c '^reaxFFtol' $inFile`;
		$infileType = "gulp" if($tstr =~ /1/);
	}
	$infileType =~ /(lammps|gulp)/;
	$infileType = lc $1;
	$readFF = \&readLammpsReaxFF;
	$writeFF = \&writeGulpReaxFF;
	$readFF = \&readGulpReaxFF if ($infileType eq "gulp");
	$writeFF = \&writeLammpsReaxFF if ($infileType eq "gulp");
	if (! defined($outFile)) {
		$outFile = basename($inFile);
		$outFile =~ s/\.\w+$//;
		$outFile .= ".${infileType}FF.dat"
	}
	$ELEMENTS = &LoadElements;
	$EngVals = {
		"reaxff1_under"			=> [0],
		"reaxff1_lonepair"		=> [1],
		"reaxff1_morse"			=> [1],
		"reaxff2_bond"			=> [0,1,2],
		"reaxff2_pen"			=> [0,1,2],
		"reaxff2_morse"			=> [0],
		"reaxff3_angle"			=> [1],
		"reaxff3_penalty"		=> [0],
		"reaxff3_conjugation"	=> [0,1,2,3],
		"reaxff3_hbond"			=> [1,2,3],
		"reaxff4_torsion"		=> [0,1,2,4],
	};

}

sub usage {
	print STDOUT <<DATA;
usage: $0 -i in_file -o (out_file) -t (in_file_type=lammps|gulp)
Arguments:
	-i in_file: name of input force field (required)
	-o out_file: name of output force field (optional)
	-t in_file_type: type of input file. Expected lammps or gulp (optional: will guess if not present)
DATA
die "\n";

}

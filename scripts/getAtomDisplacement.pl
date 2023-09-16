#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use warnings;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester CoM TrjSelections GetStats);
use ManipAtoms qw(SplitAtomsByMol SelectAtoms BuildAtomSelectionString GetMols GetAtmData ImageAtoms UnwrapAtoms);
use BOX qw (GetBox InitBox);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox CreateLAMMPSTrj);
use REPLICATE qw(GetBoxVol);

sub init;
sub showUsage;
sub groupAtomsListByField;
sub tabulate;
sub writeData;
sub calcAtomDisplacement;

my ($bgfFile, $prefix, $trj, $bin, $field, $frameNum);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $atmGrps, $totAtoms, $printStr, $DATA, $MOLS);
my ($bmin, $bmax);

$|++;
&init;
print "Parsing bgf file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$BOX = GetBox($ATOMS, undef, $HEADERS);
$MOLS = GetMols($ATOMS,$BONDS);
$trj->{REF_ATOMS} = SelectAtoms($trj->{REF}, $ATOMS)
    if (defined($trj) and exists($trj->{IMAGE_FLAG}) and $trj->{IMAGE_FLAG});
$atmGrps = groupAtomsByField($ATOMS, $field);
print "Done\n";
if (! defined($trj)) {
    print "Calculating displacement for $bgfFile...";
    &calcAtomDisplacement($ATOMS, $BOX, 1, undef);
    print "Done\n";
} else {
    $totAtoms = scalar keys %{ $ATOMS };
    $trj->{getByteOffset}($trj->{FRAMES}, $trj->{FILE}, $totAtoms);
    if ($trj->{TYPE} == 2) {
		&GetLammpsTrjType($trj->{FRAMES}, $trj->{FILE}, "atom", \%{ $trj->{LAMMPSOPTS} });
		$totAtoms = "atom";
    }
    $printStr = "Calulating atom displacement from $trj->{FILE}...";
    $trj->{getSnapshot}($ATOMS, $trj->{FILE}, $trj->{FRAMES}, $totAtoms, \&calcAtomDisplacement, $printStr, undef);
}
print "Tabulating data...";
&tabulate($ATOMS, $atmGrps);
print "Done\nWriting Data...";
&writeData($atmGrps,$prefix);
print "Done\n";

sub writeData {
	my ($grps, $sprefix) = @_;
	my ($i, $j, $k, $l, $val, $blen);

	$blen = ($bmax-$bmin)/$bin;
	for $i (keys %{ $grps }) {
		open OUTFILE, "> ${sprefix}.${field}${i}.distrib.dat" or die "ERROR: Cannot write to ${sprefix}.${field}${i}.distrib.dat:$!\n";
		printf OUTFILE "%-10s %12s %12s %12s\n","BIN","X","Y","Z";
		for $l (1 .. $blen) {
			$j = $l*$bin+$bmin;
			printf OUTFILE "%10.5f ",$j;
			for $k ("X", "Y", "Z") {
				$val = 0;
				$val = $grps->{$i}{HIST}{$k}{$j} if(exists($grps->{$i}{HIST}{$k}{$j}));
				printf OUTFILE "%12.5f ", $val;
			}
			printf OUTFILE "\n";
		}
		close OUTFILE;
	}
}

sub tabulate {
	my ($atoms, $grps) = @_;
	my ($i, $j, $k, $m, $curr);

	for $m (keys %{ $grps }) {
		for $i (keys %{ $grps->{$m} }) {
			next if ($i eq "HIST");
			for $j ("X", "Y", "Z") {
				$atoms->{$i}{STATS}{$j} = GetStats($atoms->{$i}{$j});
				for $k (0 .. $#{ $atoms->{$i}{$j} }) {
					$curr = int(($atoms->{$i}{$j}[$k] - $atoms->{$i}{STATS}{$j}{AVG})/$bin)*$bin;
					$grps->{$m}{HIST}{$j}{$curr}++;
					$bmin = $curr if (! defined($bmin) or $curr<$bmin);
					$bmax = $curr if (! defined($bmax) or $curr>$bmax);
				}
			}
		}
	}
}

sub calcAtomDisplacement {
    my ($DATA, $bgfInfo, $fileHandle) = @_;
	my ($CENTER, $tstep, $BOX, $totATms, $i, $MOLECULE, $j, $l);

	$tstep = $DATA->{"TIMESTEP"}[0];
	$frameNum++;
	$BOX = ConvertLammpsBox($DATA->{"BOX BOUNDS"});
	&InitBox($BOX, $DATA->{ATOMS});
	&GetBoxVol($BOX);
	if (exists($trj->{LAMMPSOPTS}) and $trj->{LAMMPSOPTS}{scaled} or $trj->{LAMMPSOPTS}{imaged}) {
		&UnwrapAtoms($DATA->{ATOMS}, $BOX, $trj->{LAMMPSOPTS}{scaled});
	}
	$CENTER = CoM($DATA->{ATOMS}, $BOX);
	for $i (keys %{ $MOLS }) {
		$MOLECULE = GetAtmData($DATA->{ATOMS}, $MOLS->{$i}{MEMBERS});
		$CENTER = CoM($MOLECULE);
		&ImageAtoms($MOLECULE, $CENTER, $BOX);
	}
	for $i (keys %{ $MOLS }) {
		$MOLECULE = GetAtmData($DATA->{ATOMS}, $MOLS->{$i}{MEMBERS});
		$CENTER = CoM($MOLECULE);
		for $j (keys %{ $MOLS->{$i}{MEMBERS} }) {
			for $l ("X", "Y", "Z") {
				push @{ $ATOMS->{$j}{$l} }, $DATA->{ATOMS}{$j}{"${l}COORD"} - $CENTER->{"${l}COORD"};
			}
		}
	}
}

sub groupAtomsByField {
    my ($atoms, $field) = @_;
    my ($list, $i, $val);

    for $i (keys %{ $atoms }) {
		$val = $atoms->{$i}{uc $field};
		$val = ${ $atoms->{$i}{uc $field} } if($field =~ /mol/i);
		$list->{$val}{$i} = 1;
    }

    return $list;
}

sub init {
    my (%OPTS, $list, $usage);

    &getopt('bdftsymaw',\%OPTS);
    $usage = showUsage();
    die "$usage\n" if (! exists($OPTS{b}));

    print "Initializing...";
    ($bgfFile, $bin, $field, $trj->{FILE}, $trj->{SELECTION}, $trj->{TYPE}, $trj->{IMAGE_FLAG}, $trj->{REF}, $prefix) = 
        ($OPTS{b}, $OPTS{d}, $OPTS{f}, $OPTS{t}, $OPTS{s}, $OPTS{y}, $OPTS{m}, $OPTS{a}, $OPTS{w});
    FileTester($bgfFile);
    if (! defined($prefix)) {
	$prefix = basename($bgfFile);
	$prefix =~ s/\.\w+//;
    }
    $field = "CHAIN" if (! defined($field) or $field !~ /^CHAIN|RESNUM|RESID|RES|MOL|MOLECULE|MOLID|MOLECULEID|MOLSIZE|NUMBONDS|FFTYPE$/i);
    $field = uc($field);
    $field = "RESNUM" if ($field =~ /RESID/);
    $field = "MOLECULEID" if ($field =~ /MOL/ && $field !~ /SIZE/);
    $field = "MOLSIZE" if ($field =~ /SIZE/);
    $bin = 0.01 if (! defined($bin) or $bin !~ /^\d+\.?\d*$/);
    if (defined($trj) and defined($trj->{FILE})) {
        if(! -e $trj->{FILE} or ! -r $trj->{FILE} or ! -T $trj->{FILE}) {
	    print "WARNING: Cannot acess $trj->{FILE}: $!...skipping...";
	    undef $trj;
        } else {
			$trj->{IMAGE_FLAG} = 0 if (! defined($trj->{IMAGE_FLAG}) or $trj->{IMAGE_FLAG} !~ /1|yes/i);
		    $trj->{IMAGE_FLAG} = 1 if ($trj->{IMAGE_FLAG} =~ /1|yes/i);
		    if ($trj->{IMAGE_FLAG}) {
				$trj->{REF} = "index>0" if (! defined($trj->{REF}));
				$trj->{REF} = BuildAtomSelectionString($trj->{REF});
			}
		    $trj->{SELECTION} = "all" if (! defined($trj->{SELECTION}));
		    $list = TrjSelections($trj->{SELECTION});
			die "ERROR: No valid frames found while searching \"$trj->{SELECTION}\"\n"
				if (! keys %{ $list } and $trj->{SELECTION} ne "*");
		    for (keys %{ $list }) {
		        $trj->{FRAMES}{$_} = $list->{$_};
			}
		    if (! defined($trj->{TYPE})) {
		        if ($trj->{FILE} =~ /\.lammps/) {
					$trj->{TYPE} = "lammps";
		        } else {
				    $trj->{TYPE} = "amber";
		        }
		    }
			if (lc($trj->{TYPE}) ne "lammps") {
		        $trj->{TYPE} = 1;
		        $trj->{getSnapshot} = \&ParseAmberTrj;
			    $trj->{getByteOffset} = \&GetAmberByteOffset;
		    } else {
			    $trj->{TYPE} = 2;
				$trj->{getSnapshot} = \&ParseLAMMPSTrj;
		        $trj->{getByteOffset} = \&GetLammpsByteOffset;
		    }
        }
		$frameNum = 0;
    } else {
		undef $trj;
	}

    print "Done\n"; 
}

sub showUsage {
    my ($usage) = <<DATA;
usage: $0 -b bgf_file -d (bin_size) -f (field) -t (trajectory) -s (trj_selection) -y (trj_type) -m (reimage coords) -a (reference) -w (save_prefix)
Options:
	-b bgf_file (Required): location of bgf file
	-d bin_ze (Optional): coordinate deviation bin size. Determines resolution of plot. Default 0.1 angstroms.
	-f field (Optional): Group atoms by specified field. Choices are molsize|moleculeid|resname|fftype|chain|index|resid.
			     If grouping corresponds to more than 1 atoms, center of mass is used. Default is index (i.e. no grouping).
	-t trajectory (Optional): Accumalate statics over trajectory file.
	-s trajectory selection: (Optional) Frames to use in trajectory. Expected single number, "*" for all or 
				 range :Ita-b:c for frames a to b every c. Default is all or "*"
	-y trj_type (Optional): Either lammps(default) or amber.
	-m reimage flag (Optional): Reimage each frame of the trajectory back into fundamental unit cell. Default is no.
	-a reference atom (Optional): If reimaging, place the atom(s) specified at the center of the cell. Default is moleculeid==1
	-w save prefix (Optional): Name of output data file. Default is {bgf_file} name without .bgf
DATA

    return $usage;
}

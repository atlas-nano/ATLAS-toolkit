#!/usr/bin/perl

# Will do molecule/molecule search and select molecules in group 1 with at least x number
# of group2 neighbors within a cutoff y

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use IO::Handle;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use General qw(TrjSelections FileTester GetBondLength CoM);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString GetAtmData);
use BOX qw(GetBox);

sub init;
sub numerically { $a<=>$b }
sub findNeighs;
sub updateBGF;
sub writeNewBGF;
sub writeDistribution;

my ($bgfFile, $selection, $trjFile, $SELECT, $saveFile, $groupA, $groupB, $nMin, $cutoff, $uField);
my ($field, $pStr, $LAMMPSOPTS, $getByteOffset, $getSnapshots, $trjType, $OUTFILE);
my ($ATOMS, $BONDS, $HEADERS, $DATA, $tot);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
die "ERROR: Invalid field \"$uField\"!\n" if (defined($uField) and ! exists($ATOMS->{1}{$uField}));
print "Done\nSelecting relevant Atoms...";
$groupA = SelectAtoms($groupA, $ATOMS);
die "ERROR: No valid atoms found for group A!\n" if (! keys %{ $groupA });
$groupA = GetMols($ATOMS, $BONDS, $groupA);
$groupB = SelectAtoms($groupB, $ATOMS);
die "ERROR: No valid atoms found for group B!\n" if (! keys %{ $groupB });
$groupB = GetMols($ATOMS, $BONDS, $groupB);
print "Done\n";
open $OUTFILE, "> $saveFile" or die "ERROR: Cannot write to $saveFile: $!\n";
$OUTFILE->autoflush(1);
printf $OUTFILE "%-10s% 10s\n","#Step","numMols";
if (! defined($trjFile)) {
    print "Finding Neighbors from $bgfFile...";
    &findNeighs($ATOMS, GetBox($ATOMS, undef, $HEADERS), 1, $OUTFILE);
    print "Done\n";
    &writeNewBGF($ATOMS,$BONDS,$HEADERS,$uField,$saveFile) if(defined($uField));
} else {
    $field = scalar keys %{ $ATOMS };
    $getByteOffset->($SELECT, $trjFile, $field);
    if ($trjType == 2) {
        &GetLammpsTrjType($SELECT, $trjFile, "coord", \%{ $LAMMPSOPTS });
        $field = "coord";
    }
    $pStr = "Neighbour Analysis using $trjFile...";
    $getSnapshots->($ATOMS, $trjFile, $SELECT, $field, \&findNeighs, $pStr, $OUTFILE);
    close $OUTFILE;
}
close $OUTFILE;
&writeDistribution($DATA,$saveFile,$tot);

sub writeDistribution {
    my ($data, $saveName, $count) = @_;
    my ($i, @tmp);

    $saveName =~ s/\.\w+$//;
    $saveName .= ".distrib.dat";

    print "Writing distribution data to $saveName...";
    @tmp = sort numerically keys %{ $data };
    open OUTFILE, "> $saveName" or die "ERROR: Cannot write to $saveName: :$!\n";
    printf OUTFILE "%-10s %10s\n","#Neighs","Prob";
    for $i (@tmp) {
	printf OUTFILE "%-10d %10.9f\n", $i, $data->{$i}/$count;
    }
    close OUTFILE;
    print "Done\n";
}
sub writeNewBGF {
    my ($atoms, $bonds, $headers, $field, $saveName) = @_;

    $saveName =~ s/\.\w+$//;
    $saveName .= ".bgf";
    print "Creating BGF file ${saveName}...";
    &updateBGF($atoms,$field);
    &addHeader($atoms,$headers);
    &createBGF($atoms,$bonds,$saveName);
    print "Done\n";
}  

sub findNeighs {
    my ($atoms, $box, $frameNum, $fileHandle) = @_;
    my ($count, $i, $j, $com, $mol1, $mol2, $dist);

    $count = 0;
    if (defined($trjFile) and $trjType == 2) { #LAMMPS
        $frameNum = $atoms->{TIMESTEP}[0];
	$box = ConvertLammpsBox($atoms->{"BOX BOUNDS"});
	%{ $box->{X} } = %{ $box->{XCOORD} };
	%{ $box->{Y} } = %{ $box->{YCOORD} };
	%{ $box->{Z} } = %{ $box->{ZCOORD} };
        $atoms = $atoms->{ATOMS};
        if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
            UnwrapAtoms($atoms,  $box, $LAMMPSOPTS->{scaled});
        }
    } elsif (defined($trjFile)) {
	$box = ConvertAmberBox(\%{ $box });
    }

    for $i (keys %{ $groupA }) {
	$groupA->{$i}{COM} = CoM(GetAtmData($atoms, $groupA->{$i}{MEMBERS}));
	$groupA->{$i}{NUM} = 0;
    }
    for $i (keys %{ $groupB }) {
	$groupB->{$i}{COM} = CoM(GetAtmData($atoms, $groupB->{$i}{MEMBERS}));
    }

    for $i (keys %{ $groupA }) {
	for $j (keys %{ $groupB }) {
	    $dist = GetBondLength($groupA->{$i}{COM}, $groupB->{$j}{COM},$box);
	    next if ($dist == 0 or $dist > $cutoff);
	    $groupA->{$i}{NUM}++;
	}
    }

    for $i (keys %{ $groupA }) {
	$tot++;
	$DATA->{ $groupA->{$i}{NUM} }++;
	next if ($groupA->{$i}{NUM} < $nMin);
	$count++;
	for $j (keys %{ $groupA->{$i}{MEMBERS} }) {
	   $atoms->{$j}{FLAG} = 1;
	}
    }

    printf $OUTFILE "%-10d %10d\n",$frameNum,$count;
}

sub updateBGF {
    my ($atoms, $field) = @_;
    my ($i);

    for $i (keys %{ $atoms }) {
	$atoms->{$i}{$field} = 0 if (exists($atoms->{$i}{FLAG}));
    }
}

sub showUsage {
    my ($usage) = <<DATA;
usage: $0 -s bgf_file -a group_a_mols -b group_b_mols -l (trj_file) -t (trj_type) -r (trj_range) -c (cutoff) -n (num_neighbors) -f (update_field) -o (output_file)
options:
	-s bgf_file: name and location of input bgf file
	-a group_a_mols: selection string for group A molecules
	-b group_b_mols: selection string for group B molecules
	-l trj_file: (Optional) LAMMPS or AMBER trjectory files
	-t trj_type: (Optional) Type of trjectory file. Default: LAMMPS
	-r trj_range: (Optional) Range of snapshots in trjectory. Either "*", "1" or ":It1-10:2". Default is all
	-c cutoff: (Optional) Cutoff distance to search for neighbors. Default 3.2A
	-n num_neighbors: (Optional) minimum number of neighbors. Default: 6
	-f update_field: (Optional) field to update. Can be either chain(defaul),resid,moleculeid or resname.
	-o output_file: (Optional) name of output data file
DATA

    return $usage;
}

sub init {
    my (%OPTS, $usage, $groupAStr, $groupBStr, $tSel, $list, $i);

    getopt('sabltrcnfo',\%OPTS);
    $usage = showUsage();
    for ("a","b","s") {
	die "$usage\n" if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($bgfFile, $groupAStr, $groupBStr, $trjFile, $trjType, $tSel, $cutoff, $nMin, $field, $saveFile) = 
	($OPTS{s}, $OPTS{a}, $OPTS{b}, $OPTS{l}, $OPTS{t}, $OPTS{r}, $OPTS{c}, $OPTS{n}, $OPTS{f}, $OPTS{o});

    FileTester($bgfFile);
    $groupA = BuildAtomSelectionString($groupAStr);
    $groupB = BuildAtomSelectionString($groupBStr);

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
	if (! defined($saveFile)) {
	   $saveFile = basename($trjFile);
	   $saveFile =~ s/\.\w+$//;
	   $saveFile .= ".neighAnal.dat";
	}
    }
    $cutoff = 3.2 if (! defined($cutoff) or $cutoff !~ /^\d+\.?\d*$/);
    $nMin = 6 if (! defined($nMin) or $nMin !~ /^\d+$/);
    $uField = "chain" if (! defined($uField) or $uField !~ /^(chain|moleculeid|resid|resnum|moleculenum|molecule|resname)/i);
    $uField =~ s/num/id/i;
    $uField = "moleculeid" if ($field =~ /molecule/i);
    $uField = uc $uField;
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= ".neighAnal.dat";
    }
    print "Done\n";
}

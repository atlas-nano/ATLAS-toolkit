#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";

use Getopt::Std qw(getopt);
use File::Basename qw(basename);

use FileFormats qw(ParseStructFile);
use General qw(FileTester TrjSelections CoM GetFileTypeStr);
use ManipAtoms qw(SelectAtoms BuildAtomSelectionString GetMols GetAtmData ImageAtoms);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);

sub usage;
sub getBounds;
sub numerically { ($a<=>$b); }
sub getTrjMinMax;

my ($structFile, $saveName, $selection, $trj);
my ($getSnapshots, $getByteOffset, $field, $pStr, $LAMMPSOPTS, $OUTFILE, $saveFile);
my ($ATOMS, $BONDS, $SELECTIONS, $BOUNDS, $dim, $RMOLS);

$|++;
&init;
print "Getting atom information from $structFile...";
($ATOMS, $BONDS) = ParseStructFile($structFile,0);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
print "Done\n";
if (!defined($trj)) {
	print "ATOM Bounds\n";
	$BOUNDS = getBounds($SELECTIONS, $ATOMS);
	for (qw/X Y Z /) {
		$dim = $_ . "COORD";
	    print "$_ $BOUNDS->{$dim}{min} $BOUNDS->{$dim}{max}\n";
	}
} else {
	if($trj->{reimage}) {
		print "Selecting reimage atoms...";
		$RMOLS = GetMols($ATOMS, $BONDS, SelectAtoms($trj->{ref},$ATOMS));
		print "Done\n";
	}
	$field = scalar keys %{ $ATOMS };
	$getByteOffset->($trj->{frames}, $trj->{name}, $field);
	if ($trj->{typeID} == 2) {
		&GetLammpsTrjType($trj->{frames}, $trj->{name}, "coord", \%{ $LAMMPSOPTS });
		$field = "coord";
	}
	$pStr = "Calculating minmax data from $trj->{name}...";
	open $OUTFILE, "> $trj->{save}" or die "ERROR: Cannot write to $trj->{save}: $!\n";
	printf $OUTFILE "%-10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n",
	"#Tstep","xmin","xmax","xlen","ymin","ymax","ylen","zmin","zmax","zlen";
	$getSnapshots->($ATOMS, $trj->{name}, $trj->{frames}, $field, \&getTrjMinMax, $pStr, $OUTFILE);
	close $OUTFILE; 
}

sub getTrjMinMax {
	my ($atoms, $box, $frameNum, $fileHandle) = @_;
	my ($tot, $bounds, $mol, $com, $i);

    if ($trj->{typeID} == 1) { #LAMMPS
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
	if($trj->{reimage}) {
		for $i (keys %{ $RMOLS }) {
			$mol = GetAtmData($atoms, $RMOLS->{$i}{MEMBERS});
			$com = CoM($mol);
			&ImageAtoms($mol, $com, $box);
		}
	}
	$bounds = getBounds($SELECTIONS, $atoms);
	printf $OUTFILE "%-10d%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n",$frameNum,
	$bounds->{XCOORD}{min},$bounds->{XCOORD}{max},$bounds->{XCOORD}{max}-$bounds->{XCOORD}{min},
	$bounds->{YCOORD}{min},$bounds->{YCOORD}{max},$bounds->{YCOORD}{max}-$bounds->{YCOORD}{min},
	$bounds->{ZCOORD}{min},$bounds->{ZCOORD}{max},$bounds->{ZCOORD}{max}-$bounds->{ZCOORD}{min};
}

sub getBounds {
    my ($select, $atoms) = @_;
    my ($i, $j, $bounds, $coords, $minmax);

	for $j (keys %{ $select }) {
		for $i (qw/XCOORD YCOORD ZCOORD/) {
			push @{ $coords->{$i} }, $atoms->{$j}{$i};
		}
	}
	for $i (qw/XCOORD YCOORD ZCOORD/) {
		@{ $minmax } = sort numerically @{ $coords->{$i} };
	    $bounds->{$i}{min} = shift @{ $minmax };
	    $bounds->{$i}{max} = pop @{ $minmax };
	}
	undef $coords;
    return $bounds;
}

sub init {
    my (%OPTS, $atomSel, $tlist, $i, $trjName, $tmp);

    getopt('bolntrms',\%OPTS);
    ($structFile,$atomSel,$trjName,$tlist,$trj->{type},$trj->{reimage},$trj->{ref},$trj->{save}) = 
		($OPTS{b},$OPTS{o},$OPTS{l},$OPTS{n},$OPTS{t},$OPTS{r},$OPTS{m},$OPTS{s});
    &usage if (! defined($structFile));
    print "Initializing...";
	$atomSel = "index>0" if (!defined($atomSel));
    $selection = BuildAtomSelectionString($atomSel);
	#trajectory options
	if(defined($trjName)) {
		if (!-e $trjName or ! -r $trjName or ! -T $trjName) {
			undef $trj;
		} else {
			$trj->{name} = $trjName;
			if (!defined($trj->{save})) {
				$trj->{save} = basename($trj->{name});
				$trj->{save} =~ s/\.\w+$//;
				$trj->{save} .= ".minmax.dat";
			}
			$tlist = "*" if (!defined($tlist));
			$tmp = TrjSelections($tlist);
			for $i (keys %{ $tmp }) {
				$trj->{frames}{$i} = $tmp->{$i};
			}
			die "ERROR: No valid frames selected with selection '$OPTS{n}'!\n"
				if (! keys %{ $trj->{frames} } and $tlist ne "*");

			if (! defined($trj->{type})) {
				$trj->{type} = "amber";		
				$trj->{type} = "lammps" if ($trj->{name} =~ /\.lammps/);
			} else {
				$trj->{type} = "lammps" if ($trj->{type} !~ /lammps|amber/i);
			}
			$trj->{typeID} = 2;
			$getSnapshots = \&ParseAmberTrj;
			$getByteOffset = \&GetAmberByteOffset;
			if($trj->{type} =~ /lammps/i) {
				$trj->{typeID} = 1;
				$getSnapshots = \&ParseLAMMPSTrj;
				$getByteOffset = \&GetLammpsByteOffset;
			}
		}
		$trj->{reimage} = 0 if (!defined($trj->{reimage}) or $trj->{reimage} !~ /1|yes/i);
		$trj->{reimage} = 1 if ($trj->{reimage} =~ /1|yes/i); 
		if($trj->{reimage}) {
			$trj->{ref} = $atomSel if (! defined($trj->{ref}));
			$trj->{ref} = BuildAtomSelectionString($trj->{ref});
		}
	} else {
		undef $trj;
	}
}

sub usage {
	my ($fTypeStr) = &GetFileTypeStr;
    print STDOUT <<DATA;
usage: $0 -b struct_file -o (atom_selection) -l (trj_name) -n (trj_sel) -t (trj_type) -r (trj_reimage) -m (trj_ref) -s (trj_save)
Arguments:
  struct_file: (REQUIRED) name of structure file
$fTypeStr  
  trj_name: (OPTIONAL) name of trajectory file
  trj_sel: (OPTIONAL) frames from trajectory to select. Either single number or range :Ita-b:c from a to b every c. Default if all "*"
  trj_type: (OPTIONAL) type of trajectory. Either LAMMPS (default) or amber
  trj_reimage: (OPTIONAL) reimage the trajectory. Default is no
  trj_ref: (OPTIONAL) reference atom(s) to reimage trajectory. Default is atom selection (see below)
  trj_save: (OPTIONAL) name of file to write trajectory minmax data too
  atom selection (OPTIONAL):
    any valid bgf field expression. E.g. resname eq 'WAT' will select
    all the "WAT" residues while index > 10 will select all indices > 10.
    combine multiple expressions to make complicated selections: e.g.
    (xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
	Default is all, i.e. "index>0"
DATA
die "\n";

}

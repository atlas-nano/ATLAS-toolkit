#!/usr/bin/perl -w
#
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo createPDB);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use MESO qw(GetMesoParms);

# This script will take a coarse grain file and convert it to it's atomistic
# counter part, using the vectors in the vecs file

sub init;
sub getVectors;
sub numerically { ($a<=>$b); }
sub mapAtoms;
sub getClosestBead;
sub mapPoint;
sub applySugerFix;
sub createAtomModel;
sub getResName;

my ($mesoBGF, $vecsFile, $parmFile, $saveName);
my ($PARMS, $MESO_INFO, $ATOMS, $VECTORS, $Connections);

$|++;
&init;
print "Obtaining Bead parameters ....";
$PARMS = GetMesoParms($parmFile);
print "Done\nReading Vector file $vecsFile...";
$VECTORS = getVectors($PARMS, $vecsFile);
print "Done\nReading Meso file $mesoBGF...";
($MESO_INFO, $Connections) = GetBGFFileInfo($mesoBGF,0);
print "Done\nMapping Meso->Atoms...";
$ATOMS = createAtomModel($PARMS, $MESO_INFO, $Connections, $VECTORS);
print "Done\nCreating Atomistic PDB file $saveName...";
createPDB($ATOMS, undef, $saveName);
print "Done\n";

sub createAtomModel {
    my ($beadData, $mesoData, $mesoBonds, $vecList) = @_;
    my ($counter, $index, $tmp, $res, $old_val, $not_found, $i);
    my (%PDBINFO);

    $index = $res = $not_found = 0;
    $old_val = "";
    for $counter (sort numerically keys %{ $mesoData }) {
	if ($mesoData->{$counter}{"RESNUM"} ne $old_val) {
	    $old_val = $mesoData->{$counter}{"RESNUM"};
	    $res++;
	}
	$mesoData->{$counter}{"INDEX"} = $counter;
	($tmp, $index) = mapAtoms($beadData, $mesoData->{$counter}, $mesoBonds, $vecList, $res, $index);
	for $i (keys %{ $tmp }) {
	    %{ $PDBINFO{$i} } = %{ $tmp->{$i} };
	}
	undef($tmp);
    }

    return \%PDBINFO;
}

sub init {
    my (%OPTS);

    getopt('bpvs', \%OPTS);

    for ("b", "p", "v") {
	die "usage: $0 -b meso bgf file -p meso parm file -v atom vector file -s [save_name(optional)]\n" 
	    if (! defined($OPTS{$_}));
    }

    print "Initializing...";
    ($mesoBGF, $parmFile, $vecsFile, $saveName) = ($OPTS{b}, $OPTS{p}, $OPTS{v}, $OPTS{s});
    FileTester($parmFile);
    FileTester($mesoBGF);
    FileTester($vecsFile);

    if (! defined($saveName)) {
	$saveName = basename($mesoBGF);
	$saveName =~ s/.\w+//;
	$saveName .= "_atom.pdb";
    }

    print "Done\n";
}


sub getVectors {
    my ($beads, $in_file) = @_;
    my (%ATOMS_INFO, $name, $res, @header, $counter, $rec, $patern);
    
    $patern = '(\S+)\s+(\w+)\s+(\S+)\s+(\S+)\s+(\-?\d+\.\d+)\s+\d+\.\d+\s+(\-?\d+\.\d+)\s+\d+\.\d+' .
	'\s+(\-?\d+\.\d+)\s+\d+\.\d+';
    open INFILE, $in_file or die "Cannot open vectors file $in_file: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /$patern/) {
	    $name = $1;
	    $res = $2;
	    
	    $rec = (
		    {
			"P1"     => $5,
			"P2"     => $6,
			"P3"     => $7,
			"BEAD1"  => getResName($beads, $3),
			"BEAD2"  => getResName($beads, $4),
		    }
		    );
	    
	    if ($rec) {
		$ATOMS_INFO{$res}{$name} = $rec;
	    }
	}
    }
    close INFILE;
    
    return \%ATOMS_INFO;
}

sub getResName {
    my ($beads, $atom) = @_;
    my ($i);

    for $i  (keys %{ $beads->{BEADS} }) {
	if ($beads->{BEADS}{$i}{ELEMENT} eq $atom) {
	    return $beads->{BEADS}{$i}{NAME};
	}
    }
    $atom =~ s/_//;
    return $atom;
}

sub mapAtoms {
    my ($BEADS, $curr_bead, $CON, $vect_info, $res, $index) = @_;
    my ($b1, $b2, %RES, $resName, $curr_atm);
    my ($atom_counter, $xpos, $ypos, $zpos, $found_atm);
    
    $resName = $curr_bead->{FFTYPE};
    $found_atm = 0;
    
    
    for $atom_counter (keys %{ $vect_info->{$resName} }) {
	$curr_atm = \%{ $vect_info->{$resName}{$atom_counter} };
	
	if ($curr_bead->{"ATMNAME"} =~ /SUO|SUE/i) {
	    applySugarFix(\%{ $curr_atm });
	}
	$curr_atm->{BEAD1} .= ",PHO,PHE" if ($curr_atm->{BEAD1} =~ /PH?/);
	$curr_atm->{BEAD1} .= ",SUO,SUE" if ($curr_atm->{BEAD1} =~ /SU?/);
	$b1 = getClosestBead($curr_atm->{"BEAD1"}, $curr_bead, 0, $CON->{$curr_bead->{"INDEX"}});

	if (defined($b1)) {
	    $curr_atm->{BEAD2} .= ",SUO,SUE" if ($curr_atm->{BEAD2} =~ /SU?/);
	    $curr_atm->{BEAD2} .= ",PHO,PHE" if ($curr_atm->{BEAD2} =~ /PH?/);
	    $b2 = getClosestBead($curr_atm->{"BEAD2"}, $curr_bead, $b1->{"INDEX"}, $CON->{$curr_bead->{"INDEX"}});
	}
		
	$index++;
	$found_atm = 1;
	$xpos = mapPoint($curr_atm, $b1, $curr_bead, $b2, "XCOORD");
	$ypos = mapPoint($curr_atm, $b1, $curr_bead, $b2, "YCOORD");
	$zpos = mapPoint($curr_atm, $b1, $curr_bead, $b2, "ZCOORD");
	$RES{$index} = (
			{
			    "LABEL" => "ATOM",
			    "RESNAME" => $curr_bead->{RESNAME},
			    "RESNUM" => $res,
			    "ATMNAME" => $atom_counter,
			    "ATMNUM" => $index,
			    "XCOORD" => mapPoint($curr_atm, $b1, $curr_bead, $b2, "XCOORD"),
			    "YCOORD" => mapPoint($curr_atm, $b1, $curr_bead, $b2, "YCOORD"),
			    "ZCOORD" => mapPoint($curr_atm, $b1, $curr_bead, $b2, "ZCOORD"),
			}
			);
    }
        
    return (\%RES, $index);
}

sub getClosestBead {
# Finds the closest bead of a specified name to the current bead
    my ($search_name, $curr_bead, $skip, $Connections) = @_;
    my ($rec, $counter);

    for $counter (@{ $Connections }) {
	next if ($counter == $skip);
	
	if ($search_name =~ /$MESO_INFO->{$counter}{ATMNAME}/) {
	    $rec = \%{ $MESO_INFO->{$counter} };
	    $rec->{"INDEX"} = $counter;
	    last;
	}
    }
    
    return $rec;
}

sub mapPoint {
    my ($atom, $bead1, $bead2, $bead3, $dimension) = @_;
    my ($result);
    
    $result = $atom->{"P1"} * $bead1->{$dimension} +
	$atom->{"P2"} * $bead2->{$dimension} +
	$atom->{"P3"} * $bead3->{$dimension};
    
    return $result;
}

sub applySugarFix {
    my ($bead) = $_[0];
    
    if ($bead->{"BEAD1"} =~ /PHO|PHE/i) {
	$bead->{"BEAD2"} = "GUA,CYT,ADE,THY";
    } else {
	$bead->{"BEAD1"} = "GUA,CYT,ADE,THY";
    }
}


		    
		    

#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}

use strict;
use Packages::FileFormats;
use Packages::General;

# This script will take a coarse grain file and convert it to it's atomistic
# counter part, using the vectors in the vecs file

sub Initialize();
sub GetVectors($);
sub Numercially;
sub MapAtoms($$$$$);
sub WriteFile($$);
sub GetClosestBead($$$$);
sub MapPoint($$$$$);
sub ApplySugerFix($);

die "usage: $0 meso_file vecs_file [save_name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($meso, $vecs, $save_name) = @ARGV;
my ($PARMS, $MESO_INFO, %ATOM_INFO, $VECTORS, $Connections);
my ($counter, $out_string, $index, $tmp, $res, $old_val, $not_found);

print "Initializing...";
Initialize();

print "Done\nReading Vector file $vecs...";
$VECTORS = GetVectors($vecs);

print "Done\nReading Meso file $meso...";
($MESO_INFO, $Connections) = GetBGFFileInfo($meso,0);
print "Done\nMapping Meso->Atoms...";

$index = $res = $not_found = 0;
$old_val = "";
for $counter (sort Numerically keys %{ $MESO_INFO }) {
    if ($MESO_INFO->{$counter}{"RESNUM"} ne $old_val) {
		$old_val = $MESO_INFO->{$counter}{"RESNUM"};
		$res++;
    }
	$MESO_INFO->{$counter}{"INDEX"} = $counter;
    ($tmp, $index) = MapAtoms(\%{ $MESO_INFO->{$counter} }, $Connections, $VECTORS, $res, $index);
    $out_string .= $tmp;
    $tmp = "";
}
print "Done\nCreating Atomistic file $save_name...";
WriteFile($out_string, $save_name);
print "Done\n";

if ($not_found > 0) {
	print "WARNING: Didn't Apply mappings for $not_found atoms\n";
}

sub WriteFile($$) {
	my ($outtext, $outfile) = @_;
	
	open OUTFILE, "> $outfile" or die "Cannot write to file $outfile: $!\n";
	print OUTFILE $outtext;
	close OUTFILE;
	
}

sub Initialize() {

    die "ERROR: $! for file $meso\n"
	if (! -e $meso or ! -T $meso or ! -r $meso);
    
    die "ERROR: $! for file $vecs\n"
	if (! -e $vecs or ! -T $vecs or ! -r $vecs);
    
    
    if (! $save_name) {
		$save_name = "out.pdb";
    }

}


sub Numerically {
    ($a<=>$b);
}

sub GetVectors($) {
    my ($in_file) = $_[0];
    my (%ATOMS_INFO, $name, $res, @header, $counter, $rec, $patern);
    
	$patern = '(\S+)\s+(\w+)\s+(\-?\d+\.\d+)\s+\d+\.\d+\s+(\-?\d+\.\d+)\s+\d+\.\d+' .
	'\s+(\-?\d+\.\d+)\s+\d+\.\d+\s+(\w+)\s+(\w+)';
    open INFILE, $in_file or die "Cannot open vectors file $in_file: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /$patern/) {
	    $name = $1;
	    $res = $2;
		
		$rec = (
				{
				 "P1"     => $3,
				 "P2"     => $4,
				 "P3"     => $5,
				 "BEAD1"  => $6,
				 "BEAD2"  => $7,
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

sub MapAtoms($$$$$) {
    my ($curr_bead, $CON, $vect_info, $res, $index) = @_;
	my ($b1, $b2, $outstring, $res_counter, $curr_atm);
	my ($atom_counter, $xpos, $ypos, $zpos, $found_atm);

    $outstring = "";
	
	$res_counter = $curr_bead->{"ATMNAME"};
	$found_atm = 0;
	
	
	for $atom_counter (keys %{ $vect_info->{$res_counter} }) {
		$curr_atm = \%{ $vect_info->{$res_counter}{$atom_counter} };

		if ($curr_bead->{"ATMNAME"} eq "SUG") {
			ApplySugarFix(\%{ $curr_atm });
		}
		$b1 = GetClosestBead($curr_atm->{"BEAD1"}, $curr_bead, 0, \@{ $CON->{$curr_bead->{"INDEX"}} });
		if (defined($b1)) {
			$b2 = GetClosestBead($curr_atm->{"BEAD2"}, $curr_bead, $b1->{"INDEX"}, \@{ $CON->{$curr_bead->{"INDEX"}} });
		}
		
		if (! defined($b2)) {
			print "";
			$not_found++;
			next;
		}
		
		$index++;
		$found_atm = 1;
		$xpos = MapPoint($curr_atm, $b1, $curr_bead, $b2, "XCOORD");
		$ypos = MapPoint($curr_atm, $b1, $curr_bead, $b2, "YCOORD");
		$zpos = MapPoint($curr_atm, $b1, $curr_bead, $b2, "ZCOORD");
		$outstring .= sprintf("ATOM%7d%5s%4s%6d%12.3f%8.3f%8.3f\n",
				      $index, $atom_counter, $curr_bead->{"RESNAME"},
				      $res, $xpos, $ypos, $zpos);
	}
	
	if (! $found_atm) {
			$index++;
			$xpos = $curr_bead->{"XCOORD"};
			$ypos = $curr_bead->{"YCOORD"};
			$zpos = $curr_bead->{"ZCOORD"};
			$outstring .= sprintf("ATOM%7d%5s%4s%6d%12.3f%8.3f%8.3f\n",
								  $index, $curr_bead->{"ATMNAME"}, $curr_bead->{"RESNAME"},
								  $res, $xpos, $ypos, $zpos);
	}
			
    return ($outstring, $index);
}

sub GetClosestBead($$$$) {
# Finds the closest bead of a specified name to the current bead
    my ($search_name, $curr_bead, $skip, $Connections) = @_;
    my ($rec, $counter);

    for $counter (@{ $Connections }) {
		if ($counter == $skip) {
			next;
		}
	
		if ($search_name =~ /$MESO_INFO->{$counter}{"ATMNAME"}/) {
			$rec = \%{ $MESO_INFO->{$counter} };
			$rec->{"INDEX"} = $counter;
			last;
		}
	}
    
    return $rec;
}

sub MapPoint($$$$$) {
	my ($atom, $bead1, $bead2, $bead3, $dimension) = @_;
	my ($result);
	
	$result = $atom->{"P1"} * $bead1->{$dimension} +
		$atom->{"P2"} * $bead2->{$dimension} +
		$atom->{"P3"} * $bead3->{$dimension};
	
	return $result;
	
}

sub ApplySugarFix($) {
	my ($bead) = $_[0];
	
	if ($bead->{"BEAD1"} eq "PHO") {
		$bead->{"BEAD2"} = "GUA,CYT,ADE,THY";
	} else {
		$bead->{"BEAD1"} = "GUA,CYT,ADE,THY";
	}
}


		    
		    

#!/usr/bin/perl -w

use strict;

#  This script will read in a lammps trajectory file, and "param" file
#  And create PQR (similar to PDB but with charge and radii info) files
#  That can be visualized with vmd

sub Initialize();
sub ReadParm();
sub ParseTrajectory();
sub WriteFile(@);
sub Numerically;

die "usage: $0 lammps_trajectory lammps_param [interval]\n"
    if (! @ARGV or $#ARGV < 1);

my ($traj, $myparm, $interval) = @ARGV;
my (%ATOMTYPE);

Initialize();
ReadParm();
ParseTrajectory();

sub Initialize() {
    die "Error: Cannot find trajectory file $traj: $!\n"
	if (! -e $traj);
    die "Error: Cannot read trajectory file $traj: $!\n"
	if (! -r $traj or ! -T $traj);

    die "Error: Cannot find parameter file $myparm: $!\n"
	if (! -e $myparm);
    die "Error: Cannot read parameter file $myparm: $!\n"
	if (! -r $myparm or ! -T $myparm);
    $interval = -1
	if (! $interval);
}

sub ReadParm() {

    my ($isvalid);

    $isvalid = 0;
    open PARMFILE, $myparm or die "Error while accessing parameter file $myparm: $!\n";

    while (<PARMFILE>) {
	chomp;
	if ($_ =~ /^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\-?\d+\.\d+)\s+(\d+\.\d+)/) {
	    $ATOMTYPE{$1} = (
			     {
				 "ATOMNAME"    => $2,
				 "RESIDUENAME" => $3,
				 "ATOMCHARGE"  => $4,
				 "ATOMRADII"   => $5,
			     }
			     );
	    $isvalid = 1;

	}
    }
    close PARMFILE;

    die "Error: Invalid parameter file: $myparm\n"
	if (! $isvalid);

}

sub ParseTrajectory() {
    my ($isvalid, $start_record, $has_data, $rec, $counter, $get_box);
    my (%ATOMDATA, @BBOX, $interval_counter);

    $isvalid = $start_record = $counter = $has_data = $get_box = 0;
    $interval_counter = 0;
    open TRAJFILE, $traj or die "Cannot access trajectory file $traj: $!\n";
    while (<TRAJFILE>) {
	chomp;
	if ($_ =~ /^\s*ITEM: TIMESTEP/) {
	    if ($has_data and $interval_counter >= $interval) {
		WriteFile(\%ATOMDATA, \@BBOX, $counter);
		$counter++;
		$interval_counter = 0;
	    }
	    %ATOMDATA = ();
	    @BBOX = ();
	    $has_data = $start_record = 0;
	    $interval_counter += 1;
	}
	if ($interval_counter >= $interval) {
	    if ($_ =~ /^\s*ITEM: BOX BOUNDS/) {
		$get_box = 1;
	    } elsif ($_ =~ /^\s*ITEM: ATOMS/) {
		$get_box = 0;
		$start_record = 1;
	    } elsif ($get_box) {
		if ($_ =~ /^\s*(\-?\d+\.\d+)\S+\s+(\-?\d+\.\d+)/) {
		    push @BBOX, $1;
		    push @BBOX, $2;
		}
	    } elsif ($start_record) {
		if ($_ =~ /^\s*(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
		    if ($ATOMTYPE{$2}) {
			$rec = (
			        {
				    "TYPEID" => $2,
				    "XCOORD" => $3,
				    "YCOORD" => $4,
				    "ZCOORD" => $5,
				    "ATOMNAME"    => $ATOMTYPE{$2}{"ATOMNAME"},
				    "RESIDUENAME" => $ATOMTYPE{$2}{"RESIDUENAME"},
				    "CHARGE"  => $ATOMTYPE{$2}{"ATOMCHARGE"},
				    "RADII"   => $ATOMTYPE{$2}{"ATOMRADII"},
				}
				);
			$ATOMDATA{$1} = $rec;
			$has_data = $isvalid = 1;
		    }else {
			print "Error: $2 does not correspond to any type in the parm file. Ignoring.\n";
		    }
		}
	    }
        }
    }
    close TRAJFILE;

    die "Error: $traj is not a valid Lammps atom trajectory file\n"
	if (! $isvalid);
}

sub WriteFile(@) {
    my ($atm_data, $bbox, $file_counter) = @_;
    my ($out_file_name, $fmt, $key_counter, $out_string);
    my ($xcoord, $ycoord, $zcoord);

    $out_file_name = "out_atoms_" . $file_counter . ".pqr";

    print "Creating $out_file_name...";
    $fmt = "%-6s%5d  %4s %3s  %4d    %8.3f%8.3f%8.3f %5.1f %5.1f\n";

    for $key_counter (sort Numerically keys %{ $atm_data }) {
	$xcoord = $atm_data->{$key_counter}->{"XCOORD"} * ($bbox->[1] - $bbox->[0]) + $bbox->[0];
	$ycoord = $atm_data->{$key_counter}->{"YCOORD"} * ($bbox->[3] - $bbox->[2]) + $bbox->[2];
	$zcoord = $atm_data->{$key_counter}->{"ZCOORD"} * ($bbox->[5] - $bbox->[4]) + $bbox->[4];

	$out_string .= sprintf($fmt, "ATOM", $key_counter, $atm_data->{$key_counter}->{"ATOMNAME"},
			       $atm_data->{$key_counter}->{"RESIDUENAME"}, $atm_data->{$key_counter}->{"TYPEID"},
			       $xcoord, $ycoord, $zcoord, $atm_data->{$key_counter}->{"CHARGE"},
			       $atm_data->{$key_counter}->{"RADII"});
    }

    open OUTDATA, "> $out_file_name" or die "Cannot create $out_file_name: $!\n";
    print OUTDATA $out_string;
    close OUTDATA;
    print "Sucess\n";
}

sub Numerically {
    ($a <=> $b);
}

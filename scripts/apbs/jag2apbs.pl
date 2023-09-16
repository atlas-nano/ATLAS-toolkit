#!/usr/bin/perl -w
use strict;
use File::Basename;

# Reads a jaguar output file, and creates a pqr file for use in apbs

sub VerifyInput(@);
sub ReadJagOut(@);
sub CreateOutput();
sub WriteOutput(@);
sub Numerically;

my ($jag_out, $use_esp, $save_name) = @ARGV;
my ($Jag_Info, $ENG_Summary);

die "usage: $0 jaguar_output use_esp [save_name]\n"
    if (! @ARGV or $#ARGV < 1);

$save_name = $jag_out
    if (! $save_name);

VerifyInput($jag_out);

($Jag_Info, $ENG_Summary) = ReadJagOut($jag_out);

WriteOutput($save_name);

sub VerifyInput(@) {
    die "Cannot locate $_[0]: $!\n"
	if (! -e $_[0]);
    die "$_[0] is not a text file!\n"
	if (! -T $_[0]);
    die "$_[0] is not readable!\n"
	if (! -r $_[0]);
    
}

sub ReadJagOut(@) {
    my ($jag_out) = $_[0];
    my ($inline, $is_radii, $hash_key, $is_geometry, $is_charge);
    my (@atom_id, @atom_charge, $counter, $is_valid, $curr_id);
    my (%Jaguar_Data, %Energy_Summary);

    $is_radii = $is_geometry = $is_charge = $is_valid = 0;
    open JAGOUT, $jag_out or die "Cannot read from $jag_out: $!\n";
    while (<JAGOUT>) {
	chomp;
	$inline = $_;
	if ($inline =~ /van der Waals radii/) {
	    $is_radii = 1;
	} elsif ($inline =~ /Number of Lewis structures found/) {
	    $is_radii = 0;
	} elsif ($is_radii) {
	    if ($inline =~ /([a-zA-Z]+)(\d+)\s+(\d+\.\d+)/) {
		$hash_key = $2;
		$Jaguar_Data{$hash_key}{"RADII"} = $3;
		$Jaguar_Data{$hash_key}{"ATOMNO"} = $2;
		$Jaguar_Data{$hash_key}{"ATOMLABEL"} = $1 . $hash_key;
		$is_valid = 1
	    }
	}elsif ($inline =~ /final geometry/) {
	    $is_geometry = 1;
	} elsif ($inline =~ /bond lengths \(angstroms\)/) {
	    $is_geometry = 0;
	}elsif ($is_geometry && $is_valid) {
	    if ($inline =~ /([A-Z]+)(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
		$hash_key = $2;
		if ($Jaguar_Data{$hash_key}) {
		    $Jaguar_Data{$hash_key}{"XPOS"} = $3;
		    $Jaguar_Data{$hash_key}{"YPOS"} = $4;
		    $Jaguar_Data{$hash_key}{"ZPOS"} = $5;
		}
	    }
	} elsif ($inline =~ /Atomic charges from electrostatic potential/ and $use_esp) {
	    $is_charge = 1;
	} elsif ($inline =~ /Atomic charges from Mulliken population analysis/ and ! $use_esp) {
	    $is_charge = 1;
	} elsif ($inline =~ /sum of atomic charges/) {
	    $is_charge = 0;
	}elsif ($is_charge && $is_valid) {
	    print "Got Charge\n";
	    if ($inline =~ /Atom\s+(.+)$/) {
		@atom_id = split /\s+/, $1;
	    }elsif ($inline =~ /Charge\s+(.+)$/) {
		@atom_charge = split /\s+/, $1;
		for $counter (0 .. $#atom_charge) {
		    $atom_id[$counter] =~ /(\d+)/;
		    $curr_id = $1;
		    if ($Jaguar_Data{$curr_id}) {
#			print "Wrote Charge\n";
			$Jaguar_Data{$curr_id}{"CHARGE"} = $atom_charge[$counter];
			$is_valid = 1;
		    }
		}
	    }
	}elsif ($inline =~ /Molecular surface area:\s+(\d+\.\d+\s+\[\w+\])/) {
	    $Energy_Summary{"MOL_AREA"} = $1;
	}elsif ($inline =~ /Reaction field energy: \s+(\-\d+\.\d+)\s+\[\w+\]/) {
	    $Energy_Summary{"ELECTROSTATICS"} = sprintf("%3.8f", ($1 / 2.478)) . " kcal/mol";
	}elsif ($inline =~ /Solvent-accessible surface area:\s+(\d+\.\d+\s+\[\w+\])/) {
	    $Energy_Summary{"SASA"} = $1;
	}elsif ($inline =~ /Cavity energy:\s+(\-?\d+\.\d+)\s+\[\w+\]/) {
	    $Energy_Summary{"CAVITY"} = sprintf("%3.8f", ($1 / 2.478)) . "  kcal/mol";
	}
    }
    close JAGOUT;

    die "Error reading Jaguar Output $jag_out\n"
	if (! $is_valid);

    return (\%Jaguar_Data, \%Energy_Summary);
}

sub WriteOutput(@) {
    my ($pqr_file, $f_extension);
    my ($counter, $out_string, $total_charge);

    $f_extension = basename($_[0]);
    if ($f_extension =~ /(\.\w{3})$/) {
	$f_extension = $1;
	$pqr_file = $_[0];
	$pqr_file =~ s/$f_extension/\.pqr/;
   }else {
	$pqr_file = $_[0];
	$pqr_file .= ".pqr";
    }


    open OUTPQR, "> $pqr_file" or die "Cannot create $pqr_file: $!\n";

   for $counter (sort Numerically keys %{ $Jag_Info } ) {
	printf OUTPQR "ATOM%7d%5s%5s%5d%12.4f%9.4f%9.4f%9.5f%8.4f\n", $counter, 
	$Jag_Info->{$counter}{"ATOMLABEL"},"RES", "444", $Jag_Info->{$counter}{"XPOS"},
	$Jag_Info->{$counter}{"YPOS"}, $Jag_Info->{$counter}{"ZPOS"},
	$Jag_Info->{$counter}{"CHARGE"}, $Jag_Info->{$counter}{"RADII"}, "";
	$total_charge += $Jag_Info->{$counter}{"CHARGE"};
    }
    close OUTPQR;

    if ($ENG_Summary) {
	
	print "-== SUMMARY -==\n";
	for $counter (keys %{ $ENG_Summary }) {
	    printf "%18s %16s\n", $counter . ":", $ENG_Summary->{$counter};
	}
#	printf "%18s %16s\n", "CHARGE:", $total_charge . " ev";
    }

    print "Create PQR file: $pqr_file\n";

}

sub Numerically {

    ($a<=>$b);

}

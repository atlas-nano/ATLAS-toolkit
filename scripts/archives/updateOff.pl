#!/usr/bin/perl -w
use strict;
use File::Basename;

#  UpdateOff.pl - This script will update the position and charge of an
#  Amber Off file with the output from a Jaguar Gemoetry optimization on that
#  Structure.
#  usage: updateOff.pl Amber_Off_File Jaguar_Out_File [save_name]

sub Initialize();
sub ParseOffFile();
sub ParseJagFile();
sub WriteOutput();

die "usage: $0  Amber_Off_File Jaguar_Out_File [save_name]\n"
    if (! @ARGV or $#ARGV < 1);
my ($amber_file, $jag_file, $save_name) = @ARGV;
my ($Jag_data, $out_data);

Initialize();

print "Parsing Jaguar Output file $jag_file...";
$Jag_data = ParseJagFile();
print "Done\n";
print "Parsing Amber Off File $amber_file...";
$out_data = ParseOffFile();
print "Done\n";

print "Writing changes to $save_name...";
WriteOutput();
print "Done\n";

sub Initialize() {
    my ($out_cmd);

    die "Cannot access $amber_file: $!\n"
	if (! -e $amber_file or ! -r $amber_file or ! -T $amber_file);

    die "Cannot access $jag_file: $!\n"
	if (! -e $jag_file or ! -r $jag_file or ! -T $jag_file);

    if ($save_name) {
	if (-e $save_name) {
	    print "NOTE: File already exists: $save_name. Renaming to $save_name";
	    print ".old\n";
	    $out_cmd = "cp $save_name $save_name" . ".old";
	    system($out_cmd);
	}
    }else {
	$out_cmd = "cp $amber_file $amber_file" . ".old";
	system($out_cmd);
	$save_name = $amber_file;
    }
}

sub ParseJagFile() {

    my ($inline, $hash_key, $is_geometry, $is_charge);
    my (@atom_id, @atom_charge, $counter, $curr_id, $is_valid);
    my (%Jaguar_Data);

    $is_geometry = $is_charge = $is_valid = 0;

    open JAGOUT, $jag_file or die "Cannot open $jag_file: $!\n";
    while (<JAGOUT>) {
	chomp;
	$inline = $_;
	
        if ($inline =~ /final geometry/) {
            $is_geometry = 1;
        } elsif ($inline =~ /bond lengths \(angstroms\)/) {
            $is_geometry = 0;
        }elsif ($is_geometry) {
            if ($inline =~ /([A-Z]+)(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/i) {
                $hash_key = $1 . $2;
		$Jaguar_Data{$hash_key}{"XCOORD"} = $3;
		$Jaguar_Data{$hash_key}{"YCOORD"} = $4;
		$Jaguar_Data{$hash_key}{"ZCOORD"} = $5;
		$Jaguar_Data{$hash_key}{"LABEL"} = $1;
            }
        } elsif ($inline =~ /Atomic charges from Mulliken population analysis/) {
            $is_charge = 1;
        }elsif ($is_charge) {
#           print "Got Charge\n";
            if ($inline =~ /Atom\s+(.+)$/) {
                @atom_id = split /\s+/, $1;
            }elsif ($inline =~ /Charge\s+(.+)$/) {
                @atom_charge = split /\s+/, $1;
                for $counter (0 .. $#atom_charge) {
                    $curr_id = $atom_id[$counter];
                    if ($Jaguar_Data{$curr_id}) {
 #                      print "Wrote Charge\n";
                        $Jaguar_Data{$curr_id}{"CHARGE"} = $atom_charge[$counter];
                        $is_valid = 1;
                    }
                }
            }
	}
    }
    close JAGOUT;

    die "Error processing Jaguar output file: $jag_file\n"
	if (! $is_valid);

    return \%Jaguar_Data;

}

sub ParseOffFile() {
    my ($out_data, $in_line);
    my ($charge_start, $position_start, $counter, @Atoms, $rec);
    my (@holder, $h_count, $my_label, $atm_counter, $Jag_label);

    $charge_start = $position_start = $counter = $atm_counter = 0;

    open OffFile, $amber_file or die "Cannot open $amber_file: $!\n";
    while (<OffFile>) {
	chomp;
	$in_line = $_;

	if ($in_line =~ /unit.atoms /) {
	    $charge_start = 1;
	    $position_start = 0;
	    $out_data .= $in_line . "\n";
	} elsif ($in_line =~ /unit.atomspertinfo /) {
	    $charge_start = 0;
	    $position_start = 0;
	    $out_data .= $in_line . "\n";
	} elsif ($in_line =~ /unit.positions /) {
	    $charge_start = 0;
	    $position_start = 1;
	    $out_data .= $in_line . "\n";
	} elsif ($in_line =~ /unit.residueconnect /) {
	    $charge_start = 0;
	    $position_start = 0;
	    $out_data .= $in_line . "\n";
	} elsif ($charge_start and ! $position_start) {
	    @holder = split /\s+/, $in_line;
	    $atm_counter = $holder[6];
	    if ($#holder ==8 and $atm_counter =~ /^\d+/) {
		$holder[1] =~ /([A-Z][a-z]?)/;
		$my_label = $1 . ($#Atoms + 2);
		$Jag_label = "";

		if (defined($Jag_data->{$my_label})) {
	  	    $Jag_label = $my_label;
		}
		$rec = (
			{
				"LABEL"    => $my_label,
				"JAGLABEL" => $Jag_label,
			}
		      );
		push @Atoms, $rec;

		if ($Jag_label) {
		    for $h_count (0 .. 7) {
			$out_data .= $holder[$h_count] . " ";
		    }
		    $out_data .= $Jag_data->{$Jag_label}->{"CHARGE"} . "\n";
		} else {
		    print "FOUND: $my_label but didn't find any corresponding Atom Label in Jaguar\n";
		    $out_data .= $in_line . "\n";
		}
	    } else {
		print "Invalid line found when looking for charge: \n$in_line\n";
		$out_data .= $in_line . "\n";
	    }
	} elsif (! $charge_start and $position_start) {
	    for $counter (@Atoms) {
		if ($counter->{"JAGLABEL"} ne "") {
		    $out_data .= sprintf(" %8.6f %8.6f %8.6f\n", $Jag_data->{$counter->{"JAGLABEL"}}{"XCOORD"},
					$Jag_data->{$counter->{"JAGLABEL"}}{"YCOORD"}, 
					$Jag_data->{$counter->{"JAGLABEL"}}{"ZCOORD"});
		} else {
		    print "Got to an atom that is not in the Jaguar file: " . $counter->{"LABEL"} . "\n";
		    $out_data .= $in_line . "\n";
		}
	    }
	} else {
	    $counter = 0;
	    $out_data .= $in_line . "\n";
	}
    }
    close OffFile;

    return $out_data;
}

sub WriteOutput() {
    open OUTDATA, "> $save_name" or die "Cannot write to $save_name: $!\n";
    print OUTDATA $out_data;
    close OUTDATA;
}

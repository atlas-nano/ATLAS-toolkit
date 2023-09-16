#!/usr/bin/perl
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}

use strict;
use File::Basename;
use Packages::General;

sub IsAcceptableRange(@);
sub ValidateInput();
sub GetInfo();
sub ShowTotals(@);
sub GetRelevantInfo(@);
sub GetHeading(@);
sub GetUnitName(@);
sub GetParms(@);

if (!@ARGV or $#ARGV < 2) {
    die "usage: $0 filebase start_value end_value\n";
}


my ($fle_nm, $start_no, $end_no) = @ARGV;
my (@Backbone, @Torsions, @valid_files, @strand_info, @inter_strand_info);
my ($total_strands, $total_units);
my ($curr_fle, $i, $j, $k, $counter);
my (@Params) = (
		"Shift_(Dx)_+/-", 
		"Slide_(Dy)_+/-", 
		"Rise_(Dz)_+/-", 
		"Tilt_(tau)_+/-", 
		"Roll_(rho)_+/-", 
		"Twist_(omega)_+/-", 
		);

my (@Torsion_Headers) = (
			"Chi_C1'-N_+/-",
			"Gamma_C5'-C4'_+/-",
			"Delta_C4'-C3'_+/-",
			"Epsil_C3'-O3'_+/-",
			"Zeta_O3'-P_+/-",
			"Alpha_P-O5'_+/-",
			"Beta_O5'-C5'",
			);
my (@Backbone_Headers) = (
			  "C1'-C2'_ _+/-", 
			  "C2'-C3'_ _+/-",
			  "Phase_ _+/-",
			  "Ampli_ _+/-",
			  "Pucker_ _+/-",
			  "C1'_ _+/-",
			  "C2'_ _+/-",
			  "C3'_ _+/-",
			  );
print "\n\n--=== START GET_HELICAL_PARMS ===--\n\n";
ValidateInput();
GetInfo();

print "Processing $counter files...\n";
for $i (@valid_files) {
    GetRelevantInfo($i);
}
print "Done\n\n";

print "Step 1: Writing Global Inter-Base Parameters to files...";
ShowTotals(1, \@strand_info);
print "Completed\n\n";

print "Step 2: Writing Global Inter-Base pair Parameters to files...";
ShowTotals(0, \@inter_strand_info);
print "Completed\n\n";

print "Step 3: Writing Backbone Parameters to files...";
ShowTotals(2, \@Backbone);
print "Completed\n\n";

print "Step 4: Writing Backbone Torsions to files...";
ShowTotals(3, \@Torsions);
print "Completed\n\n";

sub GetRelevantInfo(@) {
    my $curr_file = $_[0];
    
    my ($i, $in_str, @base_base_info, @inter_base_info, $my_strand);
    my ($is_part_E, $is_part_F, $curr_unit, $curr_header, $linfo);
    my ($is_Backbone, $is_torsion, @backbone_parms);

    ($is_part_E, $is_part_F) = (0,0);
    $is_Backbone = 0;
    if(open INFILE, $curr_file) {
	while (<INFILE>) {
	    chomp;
	    if ($_ =~ /Global Inter-Base pair Parameters/) {
		print "CHANGING\n";
		$is_part_E = 0;
		$is_part_F = 1;
	    } elsif ($_ =~ /Global Inter-Base Parameters/) {
		print "STARTING\n";
		$is_part_E = 1;
		$is_part_F = 0;
	    } elsif ($_ =~ /Local Inter-Base Parameters/) {
		print "ENDING\n";
		$is_part_E = 0;
		$is_part_F = 0;
	    } elsif ($_ =~ /Backbone Parameters/) {
		$is_Backbone = 1;
	    } elsif ($_ =~ /Groove parameters/) {
		$is_Backbone = 0;
	    } elsif ($_ =~ /^[\s+\d+\)]|[Strand]/ && ! $is_Backbone) {
		if ($is_part_F) { push @inter_base_info, $_ }
		elsif ($is_part_E) { push @base_base_info, $_ };
	    } elsif ($is_Backbone && $_ =~ /\S+/) {
		push @backbone_parms, $_;
	    }
		
 	}
	close INFILE;

#	Now extract the information from the Units
	$my_strand = "";
	for $linfo (@base_base_info) {
	    if ($linfo =~ /^\s+(\d+)\w+\s+strand/) {
		$my_strand = $1 - 1;
		print "NEW STRAND: $my_strand\n";
	    } elsif ($linfo =~ /^\s+(\d+)\)\s+(\w)\s*(\d+).(\w)\s*(\d+)/ && $my_strand ne "") {
		$curr_unit = "$1_$2-$3/$4-$5";
		$i = 0;
		while ($linfo =~ m/(\-?\d+\.\d+)/g) {
		    $strand_info[$my_strand]->{$curr_unit}->{$Params[$i]} .= "$1 ";
		    $i++;
		}
	    }
	}

#	Now get the inter base info

	for $linfo (@inter_base_info) {
	    if ($linfo =~ /^\s+(\d+)\)\s+(\w)\s*(\d+).(\w)\s*(\d+)/) {
		$curr_unit = "$1_$2-$3/$4-$5";
		$i = 0;
		while ($linfo =~ m/(\-?\d+\.\d+)/g) {
		    $inter_strand_info[0]->{$curr_unit}->{$Params[$i]} .= "$1 ";
		    $i++;
		}
	    }
	}

#	Now get the backbone and the Torsion Info
	for $linfo (@backbone_parms) {
	    if ($linfo =~ /^\s+(\d+)\w+\s+strand/) {
		$my_strand = $1 - 1;
		print "NEW STRAND: $my_strand\n";
		$is_torsion = 0;
	    } elsif ($linfo =~ /^\s+Torsions/) {
		$is_torsion = 1;
		print "Getting Torsions for $my_strand\n";
	    } elsif ($my_strand ne "" && $_ =~ /^\s+(\d+)\)\s+(\w)\s+(\d+)\s+/) {
		$curr_unit = "$1_$2-$3";
		$i = 0;
		while ($linfo =~ m/(\-?\d+\.\d+)/g) {
		    if ($is_torsion) {
			@Torsions[$my_strand]->{$curr_unit}->{$Torsion_Headers[$i]} .= "$1 ";
		    }else {
			@Backbone[$my_strand]->{$curr_unit}->{$Backbone_Headers[$i]} .= "$1 ";
		    }
		    $i++;
		}
	    }
	}
    } else { 
	print "Unable to open $curr_file: $!\n";
    }
    
}


sub ShowTotals(@) {
    my ($is_E, $CurrArray) = @_;
    my ($curr_strand, $strand_counter, $unit_counter, $parms_counter, $parms_holder);
    my (%Parm_Avg, $outText, $avg, $stdev, $outFile);

    $outText = "\n";
    print "IS_E: $is_E\n";
    $strand_counter = 0;
    for $curr_strand (@{ $CurrArray }) {
	$strand_counter++;
	%Parm_Avg = ();
	$outText .= GetHeading($is_E, $strand_counter);
	for $unit_counter (sort { $a <=> $b } keys %{ $curr_strand }) {
	    $parms_holder = \%{ $curr_strand->{$unit_counter} };
	    $outText .= GetUnitName($unit_counter);
	    for $parms_counter (@Params) {
		chop $parms_holder->{$parms_counter};
		($avg, $stdev) = STDev($parms_holder->{$parms_counter});
		$outText .= sprintf("%10.2f%6.2f", $avg, $stdev);
		$Parm_Avg{$parms_counter} .= "$avg ";
	    }
	    $outText .= "\n";
	}
	$outText .= sprintf("%-16s\n", "#Average");
	$outText .= sprintf("%16s", " ");
	for $parms_counter (@Params) {
	    chop $Parm_Avg{$parms_counter};
	    ($avg, $stdev) = STDev($Parm_Avg{$parms_counter});
	    $outText .= sprintf("%10.2f %5.2f", $avg, $stdev);
	}
	$outText .= "\n#\n";
    }

    if ($is_E == 1) { 
	$outFile = basename($fle_nm) . "_StrandInfo.dat";
    }elsif ($is_E == 0) {
	$outFile = basename($fle_nm) . "_InterStrandInfo.dat";
    } elsif ($is_E == 2) {
	$outFile = basename($fle_nm) . "_Backbone.dat";
    } elsif ($is_E == 3) {
	$outFile = basename($fle_nm) . "_Torsions.dat";
    }	
    
    print "Creating $outFile...";
    open OUTFILE, "> $outFile" || die "Cannot write to $outFile: $!\n";
    print OUTFILE $outText;
    close OUTFILE;
    print "Sucess";
}


sub ValidateInput() {
    if (! $start_no =~ /^\d+$/ || ! $end_no =~ /^\d+$/) {
	die "Invalid start and end values. Expected integers\n";
    }
    
    if ($start_no > $end_no) {
	($start_no, $end_no) = ($end_no, $start_no);
    }
    
    for $i ($start_no .. $end_no) {
	$curr_fle = $fle_nm . "." . $i;
	if (-e $curr_fle) {
	    push @valid_files, $curr_fle;
	    $counter++;
	} else {
	    print "Warning Cannot find file $curr_fle\n";
	}
    }
    
    $counter ?
	print "Found $counter valid files. Beginning Extraction\n" :
	die "Cannot find any valid files\n" if (! $counter);
}

sub GetInfo() {
# Get the info from the first file
    my ($valid_file, @frags);

    $curr_fle = $valid_files[0];
    $valid_file = 0;
    open INFILE, $curr_fle or die "Cannot open $curr_fle: $!\n";

    while (<INFILE>) {
	chomp;
	if ($_ =~ /=/) {
	    @frags = split /=/, $_;
	    $valid_file = 1;
	    if ($frags[1] =~ /(\d+)/) {
		$total_strands = $1;
	    }
	    if ($frags[4] =~ /(\d+)/) {
		$total_units = $1;
	    }
	    last;
	}
    }
    
    close INFILE;

    $valid_file ? 
	print "Generating files for $total_strands strands ", 
	"and $total_units units (" . ($total_units/$total_strands),
	" units per strand)\n" :
	warn "Invalid format for file $curr_fle\n";
}

sub GetHeading(@) {

    my ($is_E, $curr_strand) = @_;
    my ($result, $counter, $tmp1, $tmp2);

    if ($is_E) {
	$result = sprintf("%-16s", "# Strand $curr_strand");
    } else {
	$result = "# Strand 1 with Strand 2\n";
	$result .= sprintf("%-16s", "     #Duplex");
	}
    
    for $counter (@Params) {
	$counter =~ /^(.+)_(.+)_(.+)$/;
	$tmp1 .= sprintf("%16s", CenterText($1, 15));
	$tmp2 .= sprintf("%10s%6s", CenterText($2, 9), CenterText($3, 6));
    }

    $result .= "$tmp1\n" . sprintf("%-16s", "# ") . "$tmp2\n";

    return $result;

}

sub GetUnitName(@) {
    my ($unit_counter, $unit_name) = split /_/, $_[0];
    my ($returnval);

    $returnval = sprintf("%4d%12s", $unit_counter, $unit_name);
    return $returnval;

}

sub GetParms(@) {
    my ($which_index) = $_[0];
    my ($counter, @row1, @row2, @row3);

    for $counter (@Params) {
	$counter =~ /^(.+)_(.+)_(.+)$/;
	push @row1, $1;
	push @row2, $2;
	push @row3, $3;
    }

    if ($which_index == 0) {
	return @row1;
    }elsif ($which_index == 1) {
	return @row2;
    } else {
	return @row3;
    }
}

sub IsAcceptableRange(@) {
    my ($tilt_angle) = $_[0];
    my ($roll_angle) = $_[1];

    my ($isValidRange);

    $isValidRange = 0;

    if ((-60.00) < $tilt_angle and $tilt_angle < 60.0) {
	if ((-60.00) < $roll_angle and $roll_angle < 60.0) {
	    $isValidRange = 1;
	}
    }

    return $isValidRange;
}


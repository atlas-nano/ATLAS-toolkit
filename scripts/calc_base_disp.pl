#!/usr/bin/perl -w
use strict;
use File::Basename;

sub DetermineLegitInt(@);
sub GetValidFiles();
sub ProcessFile(@);
sub CalculateStats();
sub PrintOutput();
sub AddLine(@);
sub STDev(@);
sub numerically { ($a<=>$b); }

#--- Main ---

if (! @ARGV or $#ARGV < 2) {
    die "usage: $0 filebase start end\n";
}

my ($filebase, $startnum, $endnum) = @ARGV;
my (@valid_files, $i, $processedfiles, %DataHash);

$processedfiles = 0;
DetermineLegitInt($startnum, "Starting Number");
DetermineLegitInt($endnum, "Ending Number");
GetValidFiles();

print "Processing files...";
for $i (@valid_files) {
    ProcessFile($i);
}

print "Done\n";
CalculateStats();
PrintOutput();

#--- End Main ---

#--- subroutine declarations

sub ProcessFile(@) {
    my ($in_file) = $_[0];
    my ($is_valid, $file_valid);

    $is_valid  = 0;
    $file_valid = 0;
    open INFILE, $in_file or die "Cannot open $in_file: $!\n";
    while (<INFILE>) {
	if ($_ =~ /Global Base pair-Axis Parameters/) {
	    $is_valid = 1;
	} elsif ($_ =~ /Global Base-Base Parameters/) {
	    $is_valid = 0;
	} elsif ($is_valid) {
	    if (AddLine($_)) {
		$file_valid = 1;
		$processedfiles++;
	    }
	}
    }
    close INFILE;
    if (! $file_valid) {
	print "WARNING: $in_file is invalid\n";
    }
}

sub AddLine(@) {
    my ($in_text) = $_[0];
    my ($hashkey, $rec, $returnval);
    if ($in_text =~ /\s+(\d+)\)\s+(\w\s*\d+\-\w\s*\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
	$hashkey = sprintf("%03d", $1);
	$DataHash{$hashkey}->{"BASE_ID"} = $2;
	$DataHash{$hashkey}->{"X_DISP"}->{"DATA"} .= "$3 ";
	$DataHash{$hashkey}->{"Y_DISP"}->{"DATA"} .= "$4 ";
	$DataHash{$hashkey}->{"INCLIN"}->{"DATA"} .= "$5 ";
	$DataHash{$hashkey}->{"TIP"}->{"DATA"} .= "$6 ";
	$DataHash{$hashkey}->{"IsValid"} = 1;
	$returnval = 1;
    } else {
	$returnval = 0;
    }

    return $returnval;
}

sub PrintOutput() {
    my ($hkeys, %outdata, $comp_key, $outfile);

    print "Creating data file...";
    for $hkeys (sort numerically keys %DataHash) {
	if ( $DataHash{$hkeys}->{"IsValid"} ) {
	    $DataHash{$hkeys}->{"BASE_ID"} =~ s/\s+/_/g;
	    for $comp_key (sort numerically keys %{ $DataHash{$hkeys} }) {
		if ($comp_key ne "BASE_ID" && $comp_key ne "IsValid") {
		    $outdata{$comp_key} .= sprintf("%5d", int($hkeys));
		    $outdata{$comp_key} .= sprintf("%11s", $DataHash{$hkeys}->{"BASE_ID"});
		    $outdata{$comp_key} .= sprintf("%9.2f",  $DataHash{$hkeys}->{$comp_key}->{"AVG"});
		    $outdata{$comp_key} .= sprintf("%6.2f",  $DataHash{$hkeys}->{$comp_key}->{"STDev"});
		    $outdata{$comp_key} .= "\n";
		}
	    }
	}
    }

    for $hkeys (sort numerically keys %outdata) {
#	$outfile = basename($filebase) . "_" . $hkeys . ".dat";
	$outfile = $hkeys . ".dat";
	if (open OUTFILE, "> $outfile") {
	    print OUTFILE $outdata{$hkeys};
	    close OUTFILE;
	} else {
	    print "WARNING: Cannot create $outfile: $!\n";
	}
    }
    print "Done\n\nAll Tasks completed... Exiting\n";
}

sub CalculateStats() {
    my ($hkeys, $avg, $total, $stdev, $hdata, $base_comp);
    print "Calculating Statistics...";

    for $hkeys (sort numerically keys %DataHash) {
	if ( $DataHash{$hkeys}->{"IsValid"} ) {
	    for $base_comp (keys %{ $DataHash{$hkeys} }) {
		if ( $base_comp ne "BASE_ID" && $base_comp ne "IsValid") {
		    $hdata = $DataHash{$hkeys}->{$base_comp}->{"DATA"};
		    chop $hdata;
		    ($avg, $stdev, $total) = STDev($hdata);
#		    print "$base_comp\n";
		    $DataHash{$hkeys}->{$base_comp}->{"AVG"} = $avg;
		    $DataHash{$hkeys}->{$base_comp}->{"STDev"} = $stdev;
		    $DataHash{$hkeys}->{$base_comp}->{"TOTAL"} = $total;
		}
	    }
	}
    }
    print "Done\n";
}

sub STDev(@) {
    my (@datavalues, $n_total, $avg, $result, $i);

    @datavalues = split / /, $_[0];
    $n_total = $#datavalues;
    $avg = 0.0;
    $result = 0.0;

    foreach $i (@datavalues) {
        $avg += $i;
    }

    $avg = $avg/($n_total + 1);

    foreach (@datavalues) {
        $result += ($_ - $avg) **2;
    }

    if ($n_total ==0) {
        $n_total = 1;
    }


    $result = sqrt($result/$n_total);
    return ($avg, $result, ($avg * $n_total));
}

sub DetermineLegitInt(@) {
    my ($inval, $e_message) = @_;

    if (! $inval =~ /^\d+$/) {
	die "Invalid value in $e_message. Expected Integer\n";
    }
}

sub GetValidFiles() {
    my ($i, $counter, $curr_fle, $line_in);

    $counter = 0;
    for $i ($startnum .. $endnum) {
	$curr_fle = $filebase . "_" . $i . ".lis";
	if (-e $curr_fle) {
	    push @valid_files, $curr_fle;
	    $counter +=1;
	} else {
	    $curr_fle = $filebase . ".$i";
	    if (-e $curr_fle) {
		push @valid_files, $curr_fle;
		$counter +=1;
	    } else {
		print "Warning: Cannot find file $curr_fle\n";
	    }
	}
    }

    if ($counter ==0) {
	die "Cannot find any valid files\n";
    }
    
    print "Found $counter valid file(s)...\n";


}


#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename;
use General;
use FileFormats qw(GetBGFFileInfo);

sub ValidateInput();
sub GetNOEList();
sub ProcessFile(@);
sub CalcNOEDist(@);
sub ResetData();
sub WriteOutput();
sub DoViolations();
sub GetAvgs();
sub CreateGaussianPlots();
sub GetGaussianData();
sub GetGaussianBar(@);
sub GetBinSize(@);
sub DetermineGoodModels();
sub WriteTracker();
sub numerically { ($a<=>$b); }
sub getSelection;
sub fixName;

my (@validfiles, $i, %NOEList,  $valid_data, %filedata, $has_avgs, %avgs, @headings);
my ($noelistfile, $templatefile, $selection, $dist_file, $avg_file, $vio) = @ARGV;
my (@NOE_pairs, $SELECT, $NOEATOMS);

$|++;

if (!@ARGV or $#ARGV < 2) {
    die "usage: $0 noelistfile templatefile selection ", 
    "[noe distance] [avg_file] [violation_limit]\n";
}

#system "clear";

$has_avgs = 0;

open OUTFILE, "> scan_results" or die "Cannot write file scan_results: $!\n";

print "\n-==Starting Program==-\n";
print "Step 1: Validating Input..";
ValidateInput();
print "Sucess\nStep 2: Getting NOE Distance Information...";
GetNOEList();
print "Sucess\nStep 3: Resetting Data..."; 
ResetData();
$NOEATOMS = getNOEAtmMap(\@validfiles, \%NOEList);
print "Sucess\nStep 4: Processing Files...";
$valid_data = 0;

for $i (@validfiles) {
    if (ProcessFile($i, $NOEATOMS)) {
	$valid_data++;
    }
}
close OUTFILE;

if ($valid_data > 0) {
    print "Sucess\nStep 5: Creating Tables...";
    WriteOutput();
    if ($has_avgs) {
	DoViolations();
    }
    CreateGaussianPlots();

    if ($dist_file && -e $dist_file) {
	DetermineGoodModels();
    }
    print "\n-==End Program==-\n";
} else {
    die "ERROR: No PDB files were found with the NOE information\n\nAborting...\n";
}

sub CreateGaussianPlots() {

    my ($h_keys, $gauss_keys, $outText, $counter, $labels); 
    my ($curr_rec, $curr_val, $total_datpoints);
    print "Creating Plot of Gaussian Distributions...";

    &GetGaussianData;

    $counter = 1;
    $labels = sprintf("%4s%25s%10s%10s%10s", "#", "NOE Couple", "MIN", "MAX", "BIN");
    for $h_keys (sort numerically keys %{ $NOEList{"GAUSSDATA"} }) {
	$outText .= sprintf("%4d%25s", $counter, $h_keys);
	$labels .= sprintf("%4d%25s", $counter, $h_keys);

	$curr_rec = $NOEList{"GAUSSDATA"}->{$h_keys};

	$outText .= sprintf("%10.3f%10.3f%10.3f", $curr_rec->{"MIN"}, 
			   $curr_rec->{"MAX"}, $curr_rec->{"BINSIZE"});
	$curr_val = $curr_rec->{"MIN"};

#	print "KEY: $h_keys, min: " . $curr_rec->{"MIN"} . " max: " . $curr_rec->{"MAX"};
#	print " bin: " . $curr_rec->{"BINSIZE"} . "\n";
	$total_datpoints = 0;
	while ($curr_val <= ($curr_rec->{"MAX"} + $curr_rec->{"BINSIZE"})) {
#	for $curr_val (sort { $a <=> $b } keys %{ $curr_rec->{"RESULTS"} } ) {
	    if ($curr_rec->{"RESULTS"}->{$curr_val}) {
		$outText .= sprintf("%10d", $curr_rec->{"RESULTS"}->{$curr_val});
		$total_datpoints += $curr_rec->{"RESULTS"}->{$curr_val};
	    } else {
		$outText .= sprintf("%10d", 0);
	    }
#	    $labels .= sprintf("%5.3f-%5.3f", $curr_val, 
#			       ($curr_val + (2 * $curr_rec->{"BINSIZE"})) );
	    $curr_val += $curr_rec->{"BINSIZE"};
	}
	$labels .= "\n";
	$outText .= sprintf("%10d\n", $total_datpoints);
	$counter++;
    }

    open OUTDATA, "> NOE_Gaussian.dat" or die "Cannot create NOE_Gaussian.dat: $!\n";
    print OUTDATA "$outText";
    close OUTDATA;
    print "Sucess\nCreated NOE_Gaussian.dat\n";
}

sub GetGaussianData() {
    my ($i, $rec1, $rec2);
    my ($bin_size, $data_val, $min_val, $max_val);

    my (@file_listing) = sort numerically keys %{ $NOEList{"DATA"} };

#    Calc Spread

    $bin_size = 0.05;

    for $i (@NOE_pairs) {
	$NOEList{"GAUSSDATA"}->{$i}->{"MAX"} = 0;
	$NOEList{"GAUSSDATA"}->{$i}->{"MIN"} = 9999;
	$data_val = "";
	for (@file_listing) {
	    $min_val = $max_val = $NOEList{"DATA"}->{$_}->{$i};
	    $data_val .= $NOEList{"DATA"}->{$_}->{$i} . " ";
	    if ($NOEList{"GAUSSDATA"}->{$i}->{"MAX"} < $max_val) {
		$NOEList{"GAUSSDATA"}->{$i}->{"MAX"} = $max_val;
	    }
	    if ($NOEList{"GAUSSDATA"}->{$i}->{"MIN"} > $min_val) {
		$NOEList{"GAUSSDATA"}->{$i}->{"MIN"} = $min_val;
	    }
	}
	chop $data_val;
#		    print "GETBINSIZE: unit: $i";
	$NOEList{"GAUSSDATA"}->{$i}->{"BINSIZE"} = GetBinSize($data_val);
    }
		
	
    for $i (@NOE_pairs) {
	$min_val = $NOEList{"GAUSSDATA"}->{$i}->{"MIN"};
	$max_val = $NOEList{"GAUSSDATA"}->{$i}->{"MAX"};
	$bin_size = $NOEList{"GAUSSDATA"}->{$i}->{"BINSIZE"};
#		    print "$i: $min_val, $max_val, $bin_size\n\t"; 
	for (@file_listing) {
	    $data_val = GetGaussianBar($NOEList{"DATA"}->{$_}->{$i}, 
				       $min_val, $bin_size);
#			print "$_: " . $NOEList{"DATA"}->{$_}->{$i}. ", $data_val, ";
	    $NOEList{"GAUSSDATA"}->{$i}->{"RESULTS"}->{$data_val}++;
	}
    }
}

sub GetGaussianBar(@) {

    my ($in_data, $min_val, $bin_size) = @_;
    my ($curr_pos,$returnval);

    $curr_pos = int( ($in_data - $min_val)/$bin_size );
    return $min_val + ($bin_size * $curr_pos);
	
}

sub GetBinSize(@) {
# This will implement the nonparametric density estimation method too calculate
# the bin size.
# This result was also obtained by Freedman and Diaconis sumarize in: 
# Izenman, A. J. 1991. "Recent developments in nonparametric density estimation."
#                       Journal of the American Statistical Association, 86(413):205-224. 
#
# W = 2 (IQR) N**-1/3 
#    where  IQR is the interquartile range 
#    (the 75th percentile minus the 25th percentile). 

    my (@datavals) = split /\s+/, $_[0];
    @datavals = sort numerically @datavals;

#    print " data vals: $#datavals";
    my ($IQR, $first_percentile, $third_percentile);

    $first_percentile = sprintf("%.0f", .25 * ($#datavals + 1));
    $first_percentile = ($datavals[$first_percentile + 1] + $datavals[$first_percentile])/2;

    $third_percentile = sprintf("%.0f", .75 * ($#datavals + 1));
    $third_percentile = ($datavals[$third_percentile + 1] + $datavals[$third_percentile])/2;

    $IQR = $third_percentile - $first_percentile;
 #   print " 1: $first_percentile, 3: $third_percentile, IQR: $IQR ";
    my ($result) = 2 * ($IQR) / (($#datavals +1) ** (1/3));

 #   print "result: $result\n";

#    $result = 0.25;
    return $result;


}
sub DoViolations() {
    my ($hkey, $outdata, $curr_avg, $i, @violation_list, $framekey, @largest_v);
    my ($vio_data, $avg, $stdev, $max_vio, $junk, $max_index, $max_vio_name);

    print "Obtaining Violations List using infomation in $avg_file....";
    my ($result, $head_len, $vio_len, $num_data) = GetAvgs();

    $head_len += 5;
    $vio_len += 5;

    $max_vio = 0;
    $max_index = 0;
    if ($result) {
	$outdata = sprintf("%-" . $head_len . "s", "NOE Couple");
	for $i (0 .. $num_data) {
	    $outdata .= sprintf("%" . $vio_len . "s", $headings[$i]);
	}
	$outdata .= "\n";
	for $hkey (sort numerically keys %avgs) {
	    if ($#{ $avgs{$hkey} } > -1) {
		if ($#{ $avgs{$hkey} } > $max_index) {
		    $max_index = $#{ $avgs{$hkey} };
		}
		for $i (0 .. $#{ $avgs{$hkey} }) {
		    $violation_list[$i] = 0;
		    $curr_avg = $avgs{$hkey}->[$i];
		    for $framekey (sort numerically keys %{ $NOEList{"DATA"} }) {
			$result = $NOEList{"DATA"}->{$framekey}->{$hkey} - $curr_avg;
			$NOEList{"TRACKER"}->{"REL"}->[$i]->{$framekey}->{$hkey} = $result;
			$result = abs($result);
			if ($result > $vio) {
			    $violation_list[$i]++;
			    $NOEList{"VIOLATIONS"}->{$framekey}->[$i]++;
			}
			if ($#largest_v >= $i) {
			    if ($result > $largest_v[$i]) {
				$largest_v[$i] = $result;
			    }
			}else {
			    $largest_v[$i] = $result;
			}
		    }
		}
	    }
	    $outdata .= sprintf("%-" . $head_len . "s", $hkey);
	    for $i (@violation_list) {
		$outdata .= sprintf("%" . $vio_len . "d", $i);
	    }
	    $outdata .= "\n";
	}


	$outdata .= sprintf("%-" . $head_len . "s", "Violation Max"); 
	for $i (@largest_v) {
	    $outdata .= sprintf("%" . $vio_len . ".3f", $i);
	}
	$outdata .= "\n";

	$outdata .= sprintf("%-" . $head_len . "s", "Violation/Model");
	for $i (0 .. $max_index) {
	    $vio_data = "";
	    for $framekey (sort numerically keys %{ $NOEList{"VIOLATIONS"} }) {
		if ($NOEList{"VIOLATIONS"}->{$framekey}->[$i]) {
		    $vio_data .= $NOEList{"VIOLATIONS"}->{$framekey}->[$i] . " ";
		    if ($NOEList{"VIOLATIONS"}->{$framekey}->[$i] > $max_vio) {
			$max_vio = $NOEList{"VIOLATIONS"}->{$framekey}->[$i];
			$max_vio_name = $framekey;
		    }
		} else {
		    $vio_data .= "0 ";
		}
	    }
	    if ($vio_data) {
		chop($vio_data);
#		print "DATA: $vio_data\n";
		($avg, $stdev, $junk) = STDev($vio_data);
	    } else {
		$avg = $stdev = 0;
	    }
	    $outdata .= sprintf("%" . $vio_len . ".3f%" . $vio_len . ".3f", 
				$avg, $stdev);
	}
	$outdata .= "\n";
	$outdata .= sprintf("%-" . $head_len . "s", "Worst Model");
	$outdata .= sprintf("%" . $vio_len . "s%" . $vio_len . "d\n", 
			    $max_vio_name, $max_vio);
	open OUTFILE, "> Violation_List.txt" || 
	    die "Cannot create Violation_List.txt: $!\n";
	print OUTFILE "$outdata";
	close OUTFILE;
	print "Sucess\nCreate file Violation_List.txt\n";
	WriteTracker();
    } else {
	die "Invalid data in $avg_file\n";
    }
}


sub GetAvgs() {

    my (@f_avg, $i, $valid_data, $valid_heading); 
    my ($largest_head_len, $max_head_len, @tmp, $hkey);

    $valid_data = 0;
    $valid_heading = 0;
    $largest_head_len = 0;
    $max_head_len = 0;

    open AVGFILE, $avg_file or die "Cannot open $avg_file: $!\n";
    while (<AVGFILE>) {
	chomp;
	if ($_ =~ /^(\w+\s+\w+\*?\s+\-\s+\w+\s+\w+\*?)\s+(.+)$/) {
	    $hkey = $1;
	    if ($NOEList{$hkey} ne "") {
#		print "Data: $_\n";
		@f_avg = (split /\s+/, $2);
		if ($#f_avg > -1) {
		    @tmp = ();
		    for $i (@f_avg) {
			if ($i =~ /^\d+\.\d+$/) {
			    push @{ $avgs{$hkey} }, $i;
			    $valid_data = 1;
			}
		    }
		}
	    }

	    if (length($1) > $largest_head_len) {
		$largest_head_len  = length($1);
	    }

	}elsif ($_ =~ /^NOE Pair\s+(.+)$/) {
	    @headings = (split /\s+/, $1);
#	    print "HEAD: $_ \nlen:$#headings\n";
#	    push  @{ $avgs{"HEADINGS"} }, @headings;
	    for $i (@headings) {
		if (length($i) > $max_head_len) {
		    $max_head_len = length($i);
		}
	    }
	    $valid_heading = 1;
	}
    }
    close AVGFILE;
#    print "DATA: $valid_data, $largest_head_len, $max_head_len, $#headings\n";
    if (! $valid_heading) {
	$valid_data = 0;
    }

    return ($valid_data, $largest_head_len, $max_head_len, $#headings);
}

sub WriteOutput() {
    my ($counter, $duplicate) = (0, 0);
    my ($i, $avg, $stdev, $outdata);
    my ($rec1, $rec2, $total);

    $outdata = "";

    for $i (sort numerically keys %NOEList) {
	($rec1, $rec2)  = split /\s\-\s/, $i;
	if ($rec1 && $rec2) {
	    if ($rec1 gt $rec2) {
		if ($NOEList{$i}) {
		    push @NOE_pairs, $i;
		    chop $NOEList{$i};
		    if ($NOEList{$i} =~ /\s/) {
			($avg, $stdev, $total) = STDev($NOEList{$i});
		    } else {
			($avg, $stdev, $total) = ($NOEList{$i}, 0.00, 1);
		    }
		    $outdata .= sprintf(
					"%5d%8.3f%8.3f%25s%8d",
					($counter +1), 
					$avg, 
					$stdev, 
					$i, 
					$total/$avg
					);
		    $outdata .= "\n";
		    $counter++;
		}
	    } else {
		$duplicate++;
	    }
	}
    }
    open OUTFILE, "> NOE_data_Table.txt" || 
	die "Cannot create  NOE_data_Table.txt: $!\n";
    print OUTFILE $outdata;
    close OUTFILE;
    print "Sucess\nWrote $counter NOE distance(s) to NOE_data_Table.txt.\n";
    print "Ignored $duplicate duplicate distance(s)\n";
    
	
    
}
sub CalcNOEDist(@) {

    my (@fkeys) = (keys %NOEList);
    
    my ($res1, $res2, $i, $distance, $hkey, $found, $notfound, @scratch);
    my ($fdata, $f_info) = @_;

    $found = $notfound = 0;
    for $i (@fkeys) {
	($res1, $res2) = split /\s\-\s/, $i;
	if ($res1 && $res2) {
#	print "IN: $res1, $res2\n";
#	if ($res2 gt $res1) {
#	    ($res1, $res2) = Swap($res1, $res2);
#	}
#	print "OUT: $res1, $res2\n";
	    
	    if ($fdata->{$res1} && $fdata->{$res2}) {
		$distance = ($fdata->{$res1}->{"XPOS"} - $fdata->{$res2}->{"XPOS"}) ** 2;
		$distance += ($fdata->{$res1}->{"YPOS"} - $fdata->{$res2}->{"YPOS"}) ** 2;
		$distance += ($fdata->{$res1}->{"ZPOS"} - $fdata->{$res2}->{"ZPOS"}) ** 2;
		$distance = sqrt($distance);
		$hkey = $res1 . " - " . $res2;
		$NOEList{$hkey} .= "$distance ";
		$NOEList{"DATA"}->{$f_info}->{$hkey} = $distance;
#	    print "distance for $res1 - $res2: $distance\n";
		$found++;
	    } else {
		print OUTFILE "Could not find key for $res1 or $res2\n";
		$notfound++;
	    }
	}
    }
#    print "Found $found NOE distances... $notfound not found\n";
}

sub getNOEAtmMap {
    my ($v_files, $NList) = @_;
    my ($curr_file, $i, @junk, $ATOM, $BONDS, %MAP, %NOE_info, $hkey);

    my (@fkeys) = (keys %NOEList);
    for $i (@fkeys) {
	push @junk, (split /\s\-\s/, $i);
	for (@junk) {
	    $NOE_info{$_} = 1;
	}
    }
    $curr_file = $v_files->[0];

    ($ATOM, $BONDS) = GetBGFFileInfo($curr_file, 0, 0);
    
    for $i (keys %{ $ATOM }) {
	$ATOM->{$i}{RESNAME} =~ s/D//g;
	$ATOM->{$i}{ATMNAME} = fixName($ATOM->{$i}{ATMNAME}) if ($ATOM->{$i}{ATMNAME} =~ /\'/);
	$hkey = $ATOM->{$i}{RESNAME} . $ATOM->{$i}{RESNUM} . " " . $ATOM->{$i}{ATMNAME};
	if (exists($NOE_info{$hkey})) {
	    $MAP{$hkey} = $i;
	}
    }

    die "ERROR: BGF file does not contain any valid NOE atoms!\n" if (! keys %MAP);
    return \%MAP;
}

sub ProcessFile(@) {

    my (@fkeys) = (keys %NOEList);
    my ($i, $hkey, $mcounter, $is_noe_point, $rec, %tmp);
    my ($curr_file, $MAP) = @_;
    my (%NOE_info, $indata, $fkey, $ATOM, $BONDS, @junk);

    for $i (@fkeys) {
	push @junk, (split /\s\-\s/, $i);
	for (@junk) {
	    $NOE_info{$_} = 1;
	}
    }
    $mcounter = 0;
    ($ATOM, $BONDS) = GetBGFFileInfo($curr_file, 0, 0);

    for $hkey (keys %{ $MAP }) {
	$i = $MAP->{$hkey};
	next if (! exists($ATOM->{$i}));
	$ATOM->{$i}{RESNAME} =~ s/D//g;
	$ATOM->{$i}{ATMNAME} = fixName($ATOM->{$i}{ATMNAME}) if ($ATOM->{$i}{ATMNAME} =~ /\'/);
	$rec = (
		{
		    "XPOS" => $ATOM->{$i}{XCOORD},
		    "YPOS" => $ATOM->{$i}{YCOORD},
		    "ZPOS" => $ATOM->{$i}{ZCOORD},
		}
		);
	if (! $tmp{$hkey}) {
	    push @{ $filedata{$hkey} }, $rec;
	    $tmp{$hkey} = $rec;
	    $mcounter++;
	}
    }
    
    if ($mcounter == 0) {
	print OUTFILE "WARNING: $_[0] is not a valid Xleap PDB file\n";
	return 0;
    } else {
	$fkey = basename($curr_file);
	$fkey =~ s/\.pdb//;
	CalcNOEDist(\%tmp, $fkey);
	return 1;
    }
       
}

sub GetNOEList() {

    my ($curr_res, $hash_key, $counter, $myflag);

    $curr_res = "";
    $counter = 0;
    $myflag = 0;

    open NOEFILE, $noelistfile or die "Cannot open $noelistfile: $!\n";
    while (<NOEFILE>) {
	chomp;
	if ($_ =~ /^RES: (\w+\d+)/) {
	    $curr_res = $1;
	}elsif ($_ =~ /^(\w+\*?)\s+\-\s+(\w+\*?)\s+(\w+\d+)/) {
	    if ($curr_res) {
		$hash_key = "$3 $2 - $curr_res $1";

		$NOEList{$hash_key} = 1;
		$counter++;

	    } else {
		print OUTFILE "ERROR: Found valid data but no residue: $_ $curr_res...";
	    }
	} else {
#	    print "ERROR: Invalid data: $_";
	}
    }
    close NOEFILE;

    print "$counter valid NOE distances found...";
    if ($counter == 0) {
	die "Error reading $noelistfile: No valid data found\n";
    }
}
sub ValidateInput() {
    my ($i, $currfile);

    -e $noelistfile or die "Cannot locate $noelistfile: $!\n";
    
    if ($avg_file && -e $avg_file) {
	if ($vio) {
	    if ($vio  =~ /^\d+\.?\d?$/) {
		$has_avgs = 1;
	    }
	}
    }

    $SELECT = getSelections($selection);
    for $i (sort numerically keys %{ $SELECT }) {
	$currfile = $templatefile . $i . ".bgf";
	if (! -e $currfile) {
	    print OUTFILE "WARNING: Cannot locate $currfile: $!\n";
	} else {
	    push @validfiles, $currfile;
	}
    }

    if ($#validfiles == -1) {
	die "ERROR: No valid PDB file found. Aborting Execution\n";
    }
}

sub getSelections {
    my ($select) = $_[0];
    my (@tmp, $i, %SELECTION, $j, $counter);

    @tmp = split /\s+/, $select;

    $counter = 0;
    for $i (@tmp) {
        if ($i =~ /^\:It(\d+)\-(\d+)\:(\d+)$/) {
            $j = $1;
            while ($j < $2) {
                $SELECTION{$j} = 1;
                $counter++;
                $j += $3;
            }
        } elsif ($i =~ /^(\d+)$/) {
            $SELECTION{$1} = 1;
            $counter++;
        }
    }

    die "ERROR: No valid selection found. See help file. Got $select\n" if ($select ne "*" && (! %SELECTION || ! keys %SELECTION));

    print "trajectory: ";
    if ($select eq "*") {
        print "using all frames...";
    } else {
        print "using $counter frames...";
    }
    return \%SELECTION;
}

sub ResetData() {

    my ($res1, $res2, $hkey1, $keylist, $counter);

    for $i (keys %NOEList) {
	$keylist .= " $i ";
    }

    $counter  = 0;
    print OUTFILE "\nNOE List:\n";
    for $i (keys %NOEList) {
	($res1, $res2) = split /\s\-\s/, $i;
	$hkey1 = $res2 . " - " . $res1;
	if (! $keylist =~ /$hkey1/) {
	    print OUTFILE "WARNING: Duplicate NOE spec not found for $hkey1\n";
	    $counter++;
	    delete $NOEList{$i};
	} else {
	    $NOEList{$i} = "";
	    print OUTFILE "$i\n";
	}
    }

    if ($counter > 0) {
	print "" . ($counter + 1) . " NOE specs were not duplicated..";
    }
}

sub DetermineGoodModels() {
    my ($inData, $curr_NOE, %Dist, $has_NOE, $curr_Model, $counter);
    my ($Dist_data, $outHash, $i, $should_proceed);
    $has_NOE = 0;

    print "Step 6: Determining Best Models...";
    open INFILE, $dist_file or die "Cannot open $dist_file: $!\n";
    while (<INFILE>) {
	chomp;
   	$inData = $_;
	if ($inData =~ /^(\w+\s+\w+\*?\s+\-\s+\w+\s+\w+\*?)\s+(\-?\d+\.\d+)/) {
	    $curr_NOE = $1;
	    if ($NOEList{$curr_NOE} ne "") {
		$Dist{$curr_NOE}{"DATA"} = $2;
		$has_NOE = 1;
	    }
	}
    }

    close INFILE;

    if (! $has_NOE) {
	die "No valid NOE Information obtained\n";
    }

    for $curr_Model (keys %{ $NOEList{"DATA"} }) {
	for $curr_NOE(keys %Dist) {
	    $counter = 0;
	    $Dist_data  = $Dist{$curr_NOE}{"DATA"};
	    if (( $NOEList{"DATA"}{$curr_Model}{$curr_NOE} <= ($Dist_data - $vio) ) or
		( $NOEList{"DATA"}{$curr_Model}{$curr_NOE} >= ($Dist_data + $vio) ) ) {
		$should_proceed = 0;
		last;
	    }
	    $should_proceed = 1;

	}
	if ($should_proceed) {
	    $outHash .= "$curr_Model ";
	}
    }

    if ($outHash) {
	open OUTFILE, "> Good_Models.txt" or die "Cannot create Good_Models: $!\n";
	print OUTFILE $outHash . "\n";
	close OUTFILE;
    }

    print "Done\n";
}

sub WriteTracker() {
    my ($counter, $index, $outstring, $my_counter, $i, $NOE_names, @tmp);

    @tmp = sort numerically keys %{ $SELECT };
    for $i (0 .. $#{ $NOEList{"TRACKER"}{"REL"} }) {
	$NOE_names = "";
	$my_counter = 1;
	$outstring = sprintf("%-15s", "TIMESTEP");
	for $counter (keys %{ $NOEList{"TRACKER"}{"REL"}->[$i] }) {
	    for $index (keys %{ $NOEList{"TRACKER"}{"REL"}->[$i]->{$counter} }) {
		$outstring .= sprintf("%10s", $my_counter);
		$NOE_names .= sprintf("%10d-%30s\n", $my_counter, $index);
		$my_counter++;
	    }
	    last;
	}	 
	
	$outstring .= "\n";
	
	$my_counter = pop @tmp;
	for $counter (keys %{ $NOEList{"TRACKER"}{"REL"}->[$i] }) {
	    $outstring .= sprintf("%-15s", $my_counter);
	    for $index (keys %{ $NOEList{"TRACKER"}{"REL"}->[$i]->{$counter} }) {
		$outstring .= sprintf("%10.4G", $NOEList{"TRACKER"}{"REL"}->[$i]->{$counter}{$index});
	    }
	    $outstring .= "\n";
	    $my_counter++;
	}
	
	open OUTFILE, "> NOE_tracker_" . $i . ".dat" or die "Cannot create NOE_tracker.dat: $!\n";
	print OUTFILE $outstring;
	close OUTFILE;

#	print NOE NAMES
	open OUTFILE, "> NOE_names.txt" or die "Cannot create NOE_name.txt:$!\n";
	print OUTFILE $NOE_names;
	close OUTFILE;
    }
}

sub fixName {
    my ($atmname) = $_[0];
    my ($retName);
    
    $atmname =~ s/\'/\*/g;
    if ($atmname =~ /^([a-zA-Z]+)(\d+)\*(\d+)/) {
	$retName = $3 . $1 . $2 . "*";
    } else {
	$retName = $atmname;
    }

    return $retName;
}

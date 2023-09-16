#!/usr/bin/perl -w

sub GetLinkInfo();
sub GetExtension(@);
sub DetermineRelevantFiles();
sub FixVals(@);

my ($fle_nm, $strand_length, $helix_no);
my (@valid_files, @link_array);

if (!@ARGV or $#ARGV < 4) {
    die "usage: do_curve_anal.pl filebase startnum endnumber whichhelix units_per_strand\n";
}

$fle_nm = $ARGV[0];

if (! $ARGV[1] =~ /^\d+$/ or ! $ARGV[2] =~ /^\d+$/) {
    die "Invalid start and end values. Expected integers\n";
}

if ($ARGV[1] > $ARGV[2]) {
    die "Invalid start and end values: start value exceedes end value\n";
}


if (! $ARGV[4] =~ /^\d+$/) {
    die "Invalid strand length. Expected integer\n";
} else { $strand_length = $ARGV[4]; }

if (! $ARGV[3] =~ /^([1|2])$/) {
    die "Invalid helix specification. Expected 1 or 2\n";
}

$fle_nm .= "_Helix" . $ARGV[3];

DetermineRelevantFiles();
GetLinkInfo();

sub DetermineRelevantFiles() {

    my $i = 0;
    my $j = 0;
    my $k = 0;
    my $counter = 0;
    my $curr_fle = "";
    my $line_in = "";

    for ($i = $ARGV[1]; $i<=$ARGV[2];$i++) {
	$curr_fle = $fle_nm . "_" . GetExtension($i) . ".pdb";
	if (-e $curr_fle) {
	    push @valid_files, $curr_fle;
	    $counter +=1;
	} else {
	    warn "Warning Cannot find file $curr_fle\n";
	}
    }

    if ($counter ==0) {
	die "Cannot find any valid files\n";
    }
    
    print "Found $counter valid file(s). Running Curve....\n";

}

sub GetExtension(@) {
    my $inval = $_[0];
    return $inval;
}

sub GetLinkInfo() {

    my (@test_cmd);
    my ($i);

    $test_cmd[0] = "/ul/maiti/Curve_linux/Cur5_s <<\!";
    $test_cmd[1] = $valid_files[0];
    $test_cmd[1] = "\&inp file=" . $test_cmd[1] . ",comb=.t.,lis=" . $test_cmd[1] . "_temp,";
    $test_cmd[2] = " pdb=" . $valid_files[0] . "_1,grv=.t., \&";
 
    if ($ARGV[3] == 1) {
	$test_cmd[3] = "2 $strand_length -" . $strand_length . " 0 0";
    } else {
	$test_cmd[3] = "2 -" . $strand_length . " $strand_length 0 0";
    }

    for $i (1 .. $strand_length) {
	$test_cmd[4] .= sprintf("%3d", $i);
    }
    
#    my ($startpoint) = $strand_length +1;
#    my ($endpoint) = $strand_length * 2;
    my ($counter) = $strand_length *2;
 
    while ($counter >= ($strand_length +1)) {
        $test_cmd[5] .= sprintf("%3d", $counter);
        $counter--;
#        print "Got here!!\n";
    }
	   
#    for ($i = $endpoint; $i>= $startpoint; $i--) {
#	$test_cmd[5] .= sprintf("%3d", $i);
#    } 

   $test_cmd[6] = "0. 0. 0. 0.";
    $test_cmd[7] = "\!";

    open TMPFILE, "> junk" or die "Cannot write to file junk: $!\n";

    for $i (0 .. 7) {
	print TMPFILE $test_cmd[$i] . "\n";
    }

    close TMPFILE;

    system "chmod +x junk";
    system "junk";

    open TMPFILE, $valid_files[0] . "_temp.lis" or die "Cannot locate file junkresults: $!\n";
    while (<TMPFILE>) {
	chomp;
	if ($_ =~ /linkage from atom O3. \(\s*(\d+)\)/) {
	    push @{$link_array[0]}, $1;
#	    print "$1 " . ($1 +1) . "\n";
	}
    }

    close TMPFILE;

    my ($tmp_file) = $valid_files[0] . "_temp.lis";
    my ($tmp_pdb_file) = $valid_files[0] . "_1.pdb";

    if ($#link_array <= -1) {
	die "Could not parse $tmp_file to obtain linking information. The file is not valid\n";
    } else {
	CreateCurveFile();
    }

    system "rm -f $tmp_file $tmp_pdb_file";


}

sub CreateCurveFile() {
    my ($i, $out_file);
    my ($helix_no) = $_[0];

    my (@tmp);

    $out_file = "helix" . $ARGV[3] . "_ana";

    if (! -d "helix" . $ARGV[3] . "_results") {
	system "mkdir helix" . $ARGV[3] . "_results";  
    }


    system "cp /ul/tpascal/curve_anal/helix_ana $out_file";

    open INFILE, "$out_file" or die "Cannot open $out_file: $!\n";

    while (<INFILE>) {
	chomp;
	push @tmp, $_;
    }
    close INFILE;

    for $i (0 .. $#tmp) {
	my ($inval) =  $tmp[$i];
	my ($outval) = FixVals($inval,1);
	$tmp[$i] = $outval;
    }


    open OUTFILE,"> $out_file" or die "Cannot open $out_file: $!\n";

    for $i (0 .. $#tmp) {
	print OUTFILE "$tmp[$i]\n";
    }
    close OUTFILE;

    system "chmod +x $out_file";
    system "$out_file";
#    system "rm -f $out_file";
    system "mv $fle_nm" . "_lis* helix" . $ARGV[3] . "_results";
    system "mv $fle_nm" . "_curve* helix" . $ARGV[3] . "_results";
}

sub FixVals(@) {
    my ($instr) = $_[0];
    my ($isHelix1) = $_[1];
    my ($i, $tmpstr);

    for $i (0 .. 5) {
	my ($filevar) = INFO . $i;
	my ($replacevar) = $link_array[0]->[$i]  . " " . ($link_array[0]->[$i] +1);
	$instr =~ s/$filevar/$replacevar/g;
    }
    my ($fle_name) = "../" . $fle_nm;

    $instr =~ s/molname/$fle_nm/g;
    $instr =~ s/startval/$ARGV[1]/g;
    $instr =~ s/endval/$ARGV[2]/g;

    if ($instr eq "0. 0. 0. 0.") {
	$tmpstr = $instr;
	if ($ARGV[3] == 1) {
	    $instr = "2 $strand_length -" . $strand_length . " 0 0\n"; 
	} else {
	    $instr = "2 -" . $strand_length . " $strand_length 0 0\n"; 
	}
	for $i (1 .. $strand_length) {
	    $instr .= sprintf("%3d", $i);
	}
	$instr .= "\n";
        
        my ($counter) = $strand_length *2;	
        while ($counter >= ($strand_length +1)) {
             $instr .= sprintf("%3d", $counter);
             $counter--;
#        print "Got here!!\n";
        }

#	my ($startpoint) = $strand_length +1;
#	my ($endpoint) = $strand_length * 2;
	
#	for $i ($startpoint .. $endpoint) {
#	    $instr .= sprintf("%3d", $i);
#	}
	$instr .= "\n" . $tmpstr;
    }
    

#    print "$instr\n";
    return $instr;
}

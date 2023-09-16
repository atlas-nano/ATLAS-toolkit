#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename;
use General;

sub GetLinkInfo();
sub DetermineRelevantFiles();
sub FixVals(@);
sub GenerateCurveInput(@);
sub CreateCurveFile(@);

if (!@ARGV or $#ARGV < 4) {
    die "usage: do_curve_anal.pl filebase startnum endnumber ", 
    "whichhelix num_basepairs\n";
}

my ($fle_nm, $start_no, $end_no, $helix_no, $strand_length) = @ARGV;
my (@valid_files, @link_info, $i);

# --== START ==--

# test input
die "Invalid start and end values. Expected integers\n" 
    if (! IsInteger($start_no) || ! IsInteger($end_no)) ;

if ($start_no > $end_no) {
    ($start_no, $end_no) = Swap($start_no, $end_no);
}

die "Invalid strand length. Expected integer\n" 
    if (! IsInteger($strand_length) );
die "Invalid helix specification. Expected 1 or 2\n" 
    if (! $helix_no =~ /^[1|2]$/);

$fle_nm .= "_Helix" . $helix_no;

#Call Subroutines
DetermineRelevantFiles();
GetLinkInfo();

for $i (@valid_files) {
    CreateCurveFile($i);
}
# --=== END ==---

sub DetermineRelevantFiles() {
    my ($i, $j, $k, $curr_fle, $line_in);

    for $i ($start_no .. $end_no) {
	$curr_fle = $fle_nm . "_" . $i . ".pdb";
	-e $curr_fle ? 
	    push @valid_files, $curr_fle  : 
	    warn "Warning Cannot find file $curr_fle\n";
    }

    $#valid_files > -1 ? 
	print "Found " . ($#valid_files + 1) ." valid file(s). Running Curve....\n" : 
	die "Cannot find any valid files\n";
}

sub GenerateCurveInput(@) {

    my ($curr_file, $is_temp) = @_;
    my ($returnval, $i, $filebase, $counter);

    if (basename($curr_file) =~ /(.+)_(\d+)\.pdb/) {
        ($filebase, $counter) = ($1, $2);
    } else {
	($filebase, $counter) = (basename($curr_file =~ s/\.pdb//), 1);
    }

    $curr_file = $_[0];

    $returnval = "/ul/maiti/Curve_linux/Cur5_s <<\!\n&inp ";
    if (! $is_temp && $#link_info > -1) {
	$returnval .= "ibond=" . ($#link_info + 1) . ", ";
    }
    $returnval .= "file=$curr_file, comb=.t., lis=" . $filebase;
    $returnval .= "_" . $counter . ", \n pdb=" . $filebase . "_curve_";
    $returnval .= $counter . ", grv=.t., &\n";

    if (! $is_temp) {
	for $i (@link_info) {
	    $returnval .= $i . " " . ($i + 1) . "\n";
	}
    }

    $returnval .= "2 $strand_length -" . $strand_length . " 0 0";
    $returnval .= "\n";

    for (1 .. $strand_length) {
	$returnval .= sprintf("%3d", $_);
    }
    $returnval .= "\n";

    my ($start, $end) = ( ($strand_length + 1), ($strand_length * 2) );

    for (reverse $start .. $end) {
	$returnval .= sprintf("%3d", $_);
    } 
    $returnval .= "\n";

    $returnval .= "0. 0. 0. 0.\n\!\n";

    return $returnval;

}

sub GetLinkInfo() {
    my ($test_cmd, $curr_lisfile, $tmp_pdb_file, $i, $counter);

    my ($file_base) = basename($valid_files[0]);

    $file_base =~ s/\.pdb//;

    $curr_lisfile = $tmp_pdb_file = $file_base;
    $curr_lisfile .= ".lis";
    $tmp_pdb_file =~ s/(\d+)$/curve_$1/;
    $tmp_pdb_file .= ".pdb";

    $test_cmd = GenerateCurveInput($valid_files[0], 1);

    open TMPFILE, "> getlinkinfo" or die "Cannot write to file getlinkinfo: $!\n";
    print TMPFILE $test_cmd;
    close TMPFILE;

    system "chmod +x getlinkinfo";
    if (open CURVEFILE, "getlinkinfo |") {
	while (<CURVEFILE>) {
	    chomp;
	    if ($_ =~ /error/i) {
		die "ERROR while creating $curr_lisfile: $_\n";
		last;
	    }
	}
	close CURVEFILE;
    } else {
	print "ERROR creating $curr_lisfile: $!\n";
    }

    $counter = 0;
    open TMPFILE, $curr_lisfile or die "Cannot locate file $curr_lisfile: $!\n";
    while (<TMPFILE>) {
	chomp;
	if ($_ =~ /linkage from atom O3. \(\s*(\d+)\)/) {
	    push @link_info, $1;
	    $counter++;
	}
    }
    close TMPFILE;


    system "rm -f $curr_lisfile $tmp_pdb_file ";

#    if (! $counter) {
#	die "Could not parse $curr_lisfile to obtain linking information. ", 
#	"The file is not valid\n";
#    }
}

sub CreateCurveFile(@) {

    my ($curr_file) = $_[0];
    my ($curr_lisfile, $new_lisfile, $outText, $out_file);

    my ($file_base) = basename($curr_file);

    $file_base =~ s/\.pdb//;

    $curr_lisfile = $file_base;
    $curr_lisfile .= ".lis";

    $outText = GenerateCurveInput($curr_file, 0);

    $out_file = "helix" . $helix_no . "_ana";

    open OUTFILE,"> $out_file" or die "Cannot open $out_file: $!\n";
    print OUTFILE "$outText";
    close OUTFILE;

    system "chmod +x $out_file";
    if (open CURVEFILE, "$out_file |") {
	while (<CURVEFILE>) {
	    chomp;
	    if ($_ =~ /error/i) {
		print "ERROR while creating $curr_lisfile: $_\n";
		last;
	    }
	}
	close CURVEFILE;
    } else {
	print "ERROR creating $curr_lisfile: $!\n";
    }

# Fix name of listfile
    $new_lisfile = $curr_lisfile;
    $new_lisfile =~ s/(\d+)\.lis/lis\.$1/;
    system "mv $curr_lisfile $new_lisfile";

    system "rm -f $out_file";
}

sub FixVals(@) {
    my ($instr) = $_[0];

    my ($i, $tmpstr, $filevar, $replacevar);

    my ($num_link) = ($#link_info) + 1;

    $filevar = "LINKINFO";
    for $i (0 .. ($num_link -1)) {
	$replacevar .= $link_info[$i]  . " " . ($link_info[$i] +1) . "\n";
    }
    chomp $replacevar;
    $instr =~ s/$filevar/$replacevar/g;

    $instr =~ s/molname/$fle_nm/;
    $instr =~ s/startval/$start_no/;
    $instr =~ s/endval/$end_no/;
    $instr =~ s/NUMLINK/$num_link/;

    if ($instr eq "0. 0. 0. 0.") {
	$tmpstr = $instr;
	$instr = "2 $strand_length -" . $strand_length . " 0 0";

	$instr .= "\n";

	for $i (1 .. $strand_length) {
	    $instr .= sprintf("%3d", $i);
	}

	$instr .= "\n";
	
	for  ($i = ($strand_length * 2); $i >= ($strand_length + 1); $i-- ) {
	    $instr .= sprintf("%3d", $i);
	}

	$instr .= "\n" . $tmpstr;
    }
    

    return $instr;
}


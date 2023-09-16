#!/usr/bin/perl -w

use strict;
use File::Basename;

#   get_valence_coefficients: This script will use xmgrace to obtain the valence
#   coeffiecients (spring constants, equilibrium positions) using data in the input file
#   (assumed to be gaussian)
sub Initialize(@);
sub CreateBatchFile(@);
sub CreateFiles(@);
sub CreatePics(@);
sub SaveData(@);
sub PrintData(@);
sub GetMesh(@);
sub ExtractName(@);
sub AppendSubtitle(@);
sub sortString;

die "usage: $0 file|directory VALENCE_TERM template_file grace_cmd\n"
    if (! @ARGV or $#ARGV < 3);

my ($datafile, $valence, $template, $grace_cmd) = @ARGV;
my ($curr_file, @FILES, $grace_batch, $exec_string, $status, $PARMS, %DATA);

print "Initializing data...";
($valence, @FILES) = Initialize($datafile, $valence);

print "Done\n";

system "mkdir -p pics";

for $curr_file (@FILES) {
    print "$curr_file...";
    $grace_batch = CreateBatchFile($curr_file, $valence);
    if ($grace_batch eq "") {
	print "too little data...skipping\n";
	next;
    }
    $exec_string = "$grace_cmd -nosafe -hardcopy -legend load -batch $grace_batch";
    print ".";

    $PARMS = ();
    print "\...";
    ($PARMS) = CreateFiles($exec_string);
    print ".";
    system("rm -fr $grace_batch");
    if ($PARMS->{"WASSUCESSFUL"}) {
	print "...";
	CreatePics(\%{ $PARMS->{"FILES"} }, $curr_file);
	print ".";
	SaveData(\%DATA, $PARMS, $curr_file, $valence);
	print "Done\n";
    } else {
	print "error\n";
    }
}
PrintData(\%DATA);

sub Initialize(@) {
    my ($search_loc, $v_data) = @_;
    my (@FILES, $tmp_data);

    @FILES = ();
    if (-e $search_loc) { # if the file/direcotry exists
	if (! -f $search_loc and -d $search_loc) { # is not a file then must be a directory
	    opendir DATFILES, $search_loc or die "ERROR: Cannot access directory $search_loc: $!\n";
	    @FILES = grep { /\.dat$/ && -f} map { "$search_loc/$_"} readdir DATFILES;
	    closedir DATFILES;
	} else {
	    $FILES[0] = $search_loc;
	}
    }

    die "ERROR: No valid files found while searching $search_loc\n"
	if ($#FILES == -1);

    if ($v_data =~ /^\s*(\d+)\s*$/) {
	if ($1 >= 1 and $1 <= 3) {
	    $v_data = $1;
	} else {
	    $v_data = "";
	}
    } elsif ($v_data =~ /^\s*(\w+)\s*$/) {
	$tmp_data = $1;
	if ($tmp_data =~ /([A-Za-z])/) {
	    if ($1 =~ /BOND/) {
		$v_data = 1;
	    } elsif ($1 =~ /ANGLE/) {
		$v_data = 2;
	    } elsif ($1 =~ /TORSION/) {
		$v_data = 3;
	    } else {
		$v_data = 0;
	    }
	} else {
	    $v_data = 0;
	}
    } else {
	$v_data = 0;
    }


    die "INVALID option for valence_term: EXPECTED BOND/ANGLE/TORSION or 1-3, GOT $_[1]\n"
	if (! $v_data);

    die "Cannot find grace at $grace_cmd\n"
	if (! -e $grace_cmd);

    return ($v_data, @FILES);
}

sub CreateBatchFile(@) {
    my ($cfile, $valence_term) = @_;
    my ($instring, $out_string, $TITLE, $AXISLABEL, $FORMULAR);
    my ($GAS_CONST, $TEMP, $DATA_TITLE, $MESHDATA, $tmp, $min, $max, $constraint);
    my ($sub_title, $rt2, $counter);    
    $GAS_CONST = 0.0019872;
    $TEMP = 298.15;
    $rt2 = 2 * $GAS_CONST * $TEMP;
    
    if ($valence_term == 1) {
	$TITLE       = "BONDS";
	$AXISLABEL   = "BOND LENGTH (A)";
	$FORMULAR    = "y = a0*exp(-a1*(x-a2)^2/$rt2)";
    } elsif ($valence_term == 2) {
	$TITLE       = "ANGLES";
	$AXISLABEL   = "BENDING ANGLE (degrees)";
	$FORMULAR    = "y = a0*exp(-a1*(cos(x deg) - cos(a2 deg))^2/$rt2)";
    } else {
	$TITLE       = "TORSIONS";
	$AXISLABEL   = "SHIFTED DIHEDRAL (degrees)";
	$FORMULAR    = "y = a0*exp(-a1*(1 + cos(a2*x - a3 deg))/$rt2)";
	
    }

    $DATA_TITLE = basename($cfile);
    $DATA_TITLE =~ s/\.\w+$//g;
    $DATA_TITLE =~ s/\-/_/g;
 
    $out_string = "";
    ($MESHDATA, $min, $max, $counter) = GetMesh($cfile, $valence_term);
    return "" if ($counter < 20);
    if ($valence_term != 3) {
	$constraint = "a2 = " . ($min + ($max - $min)/2) . "\n";
	$constraint .= "a2min = " . sprintf("%.1f", $min) . "\n";
	$constraint .= "a2max = " . sprintf("%.1f", $max) . "\n";
    } else {
	$constraint = "a2 = 1\na2 constraints off";
    }
    open INFILE, "$template" or die "Cannot open grace template $template :$!\n";
    while (<INFILE>) {
	chomp;
	$instring = $_;
	$instring =~ s/INFILENAME/$cfile/g;
	$instring =~ s/DATATITLE/$DATA_TITLE/g;
	$instring =~ s/VALENCE/$TITLE/g;
	$instring =~ s/AXISLABEL/$AXISLABEL/g;
	$instring =~ s/FORMULAR/$FORMULAR/g;
	$instring =~ s/MESHDATA/$MESHDATA/g;
	$instring =~ s/a2options/$constraint/g;
	$out_string .= "$instring\n";
    }

    close INFILE;

    $tmp = basename($cfile);
    $tmp =~ s/\.w+$//g;
    $tmp = "_" . $tmp . "_batch.tmp";

    open OUTFILE, "> $tmp" or die "Cannot write to $tmp: $!\n";
    print OUTFILE $out_string;
    close OUTFILE;
    
    return $tmp;
}    

sub CreateFiles(@) {
    my ($exec_string) = $_[0];
    my (%RETURNDATA, $subtitle);

    if (system("$exec_string >& _junk ")) {
	die "Error while executing command $exec_string\n";
    }

    $exec_string = "_junk";

    $RETURNDATA{"WASSUCESSFUL"} = 0;
 
    open PIPEPROG, $exec_string or die "Cannot open $exec_string: $!\n";
    while (<PIPEPROG>) {
	if ($_ =~ /^Computed values:/) {
	    $RETURNDATA{"WASSUCESSFUL"} = 1;
	} elsif ($RETURNDATA{"WASSUCESSFUL"}) {
	    if ($_ =~ /^\s+(a\d+) = (\-?\d+\.\d+)/) {
		if ($1 ne "a0") {
		    $RETURNDATA{"PARMS"}{$1} = $2;
		    $subtitle .= "$1 - $2, ";
		}
	    } elsif ($_ =~ /Chi-square: (\d+\.\d+)/) {
		$RETURNDATA{"STATS"}{"Chi-square"} = $1;
	    } elsif ($_ =~ /Correlation coefficient: (\d+\.\d+)/) {
		$RETURNDATA{"STATS"}{"Correlation"} = $1;
	    } elsif ($_ =~ /Theil U coefficent: (\d+\.\d+)/) {
		$RETURNDATA{"STATS"}{"Theil_U"} = $1;
	    }
	}
    }
    close PIPEPROG;
    
    if (defined($subtitle)) {
	chop $subtitle;
	chop $subtitle;
	AppendSubTitle("test2.agr", $subtitle);
    }
    system "rm -fr _junk";
    $RETURNDATA{"FILES"}{"1"} = "test.agr";
    $RETURNDATA{"FILES"}{"2"} = "test2.agr";
    return \%RETURNDATA;
}

sub CreatePics(@) {
    my ($pic_list, $myfile) = @_;

    my ($counter, $pic_file, $index);

    $index = 1;
    for $counter (values %{ $pic_list}) {
	$pic_file = basename($myfile);
	$pic_file =~ s/\.\w+$//g;
	$pic_file .= "_" . $index . ".jpg";
	$pic_file = "pics/" . $pic_file;
	$exec_string = "$grace_cmd $counter -hardcopy -hdevice 'JPEG' -printfile '" . $pic_file . "'";
	system($exec_string);
	$index++;
	system ("rm -fr $counter");
    }
	
}

sub SaveData(@) {
    my ($data_ptr, $parms_ptr, $curr_file, $valence) = @_;
    my ($counter);
    
    $curr_file = basename($curr_file);
    $curr_file =~ s/\.\w+$//g;
    
    for $counter (keys %{ $parms_ptr->{"PARMS"} }) {
	$data_ptr->{$valence}{$curr_file}{$counter} = $parms_ptr->{"PARMS"}{$counter};
    }

    for $counter (keys %{ $parms_ptr->{"STATS"} }) {
	$data_ptr->{$valence}{$curr_file}{$counter} = $parms_ptr->{"STATS"}{$counter};
    }


}
sub ExtractName(@) {
    my ($instr) = $_[0];

    my (@names) = split /\-/, $instr;
    my ($out_string);

    $out_string = " ";
    for (@names) {
	$out_string .= "$_  ";
    }

    return $out_string;



}
sub PrintData(@) {
    my ($data) = $_[0];
    my ($valence, $file_nm, $parm, $out_text, $tmp, $tmp1);

    print "STATS\n";
    for $parm (0 .. 120) {
	print "-";
    }

    for $valence (keys %{ $data }) {
	if ($valence == 1) {
	    $tmp = "BOND_STRETCH";
	    $tmp1 = "HARMONIC";
	} elsif ($valence == 2) {
	    $tmp = "ANGLE_BEND";
	    $tmp1 = "THETA_HARM";
		
	} else {
	    $tmp = "TORSIONS";
	    $tmp1 = "DIHEDRAL";
	}

	print "\n$tmp\n";

	print "\n";
	printf "%-20s%10s%15s%15s%15s%15s%15s%15s\n", "Atoms", "Potential", "Chi-square", "Correlation", "Theil_U",
	"a1", "a2", "a3";
	for $file_nm (sort sortString keys %{ $data->{$valence} }) {
	    printf "%-20s%10s", ExtractName($file_nm), $tmp1;
	    for $parm("Chi-square", "Correlation", "Theil_U") {
		if (defined($data->{$valence}{$file_nm}{$parm})) {
		    printf "%15.4G", $data->{$valence}{$file_nm}{$parm};
		} else {
		    printf "%15s", "-";
		}
	    }
	    for $parm ("a1", "a2", "a3") {
		if (defined($data->{$valence}{$file_nm}{$parm})) {
		    printf "%15.4f", $data->{$valence}{$file_nm}{$parm};
		}
	    }
	    print "\n";
	}
    }

    for $parm (0 .. 120) {
	print "-";
    }
    print "\n";
}

sub GetMesh(@) {
    my ($curr_file, $type) = @_;
    my ($curr_val, $max_val, $min_val, $num_bins, $returnstr, $counter);
    my ($spacing);

    $max_val = -999;
    $min_val = 9999;
    $counter = 0;
    open INFILE, $curr_file or die "Cannot open file $curr_file: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^\s+\d+\s+(\d+\.\d+)/) {
	    $counter++;
	    $curr_val = $1;
	    $max_val = $curr_val
		if ($curr_val > $max_val);
	    $min_val = $curr_val
		if ($curr_val < $min_val);
	}
    }
    close INFILE;

    $num_bins = int($counter/4);
    $num_bins++;
    if ($counter > 0) {
	$spacing = ($max_val - $min_val)/$counter;
    }

    if ($num_bins > 20) {
	$num_bins = 20;
    }
    if ($type == 1) {
	$min_val = $min_val - ($spacing);
	$max_val = $max_val + ($spacing);
    } else {
	$min_val = int($min_val/10) * 10;
	$max_val = (int($max_val/10) + 1) * 10;
    }

    $returnstr = "MESH(" . $min_val . "," . $max_val . "," . $num_bins . ")";
    return ($returnstr, $min_val, $max_val, $counter);

}

sub AppendSubTitle(@) {
    my ($file, $stitle) = @_;
    my ($counter, $cmd, $inStr, $out_string);

    if (open INFILE, $file) {
	while (<INFILE>) {
	    chomp;
	    $inStr = $_;
	    $inStr =~ s/STITLE/$stitle/g;
	    $out_string .= "$inStr\n";
	}
	close INFILE;
	open OUTFILE, "> " . $file;
	print OUTFILE $out_string;
	close OUTFILE;
    }
}
		
sub sortString {
    { lc($a) cmp lc($b) }
}

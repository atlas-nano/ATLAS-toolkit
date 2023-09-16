#!/usr/bin/perl -w

use strict;
sub DisplayHelpFile();
sub DetermineLegitInt(@);
sub GetValidFiles();
sub GetParameters();
sub CreatePdbfiles();
sub MovePdbFiles();
sub MoveLisFiles();
sub RenamePDBFiles(@);


my ($i, $j, $cmd_str);
my (@valid_files);

if (!@ARGV) {
    die "usage: $0 paramfile filebase startnum endnumber interval\n--h for the help file\n";
}

if ($ARGV[0] eq "--h") {
    DisplayHelpFile();
} else {
    if ($#ARGV < 4) {
	die "usage: $0 paramfile filebase startnum endnumber interval\n--h for the help file\n";
    }
}

my ($paramfile, $filebase, $startnum, $endnum, $interval) = @ARGV;
my ($groove_desc, $is3prime, $numEndBases, $strand_length); 
my(@crossovers, $topfile, $trajfile, $firstfile);

GetParameters();

DetermineLegitInt($startnum, "Starting Number");
DetermineLegitInt($endnum, "Ending Number");
DetermineLegitInt($numEndBases, "Bases at the end");
DetermineLegitInt($strand_length, "Length of strand");
DetermineLegitInt($interval, "Time step");

for $i (0 .. $#crossovers) {
    DetermineLegitInt($crossovers[$i], "Crossover specification");
}

if (! $is3prime =~ /^[1|2]$/) {
    die "Invalid value in is3prime. Expected Integer\n";
}

CreatePdbfiles();
print "\n";


#GetValidFiles();

#die "$#valid_files total files. Crossovers: $#crossovers $groove_desc $is3prime $numEndBases\n";
for $i (0 .. $#valid_files) {
    my ($cmd_str) = "/ul/tpascal/scripts/helix_from_strands.pl ";
    $cmd_str .= $valid_files[$i] . " $groove_desc $is3prime $numEndBases ";
    for $j (0 .. $#crossovers) {
	$cmd_str .= $crossovers[$j] . " ";
    }

    system "$cmd_str";
}

# now execute curve for helix1
print "\n";
$cmd_str = "/ul/tpascal/scripts/do_curve_anal_newtest.pl";
$cmd_str .= " $filebase $startnum $endnum 1 $strand_length";
system "$cmd_str";
print "\n";

$cmd_str = "/ul/tpascal/scripts/get_helical_parms.pl";
$cmd_str .= " helix1_results/" . $filebase . "_Helix1_lis $startnum $endnum";
system "$cmd_str";
print "\n";

chdir "helix1_results";
MovePdbFiles();
MoveLisFiles();
chdir "../";

$cmd_str = "/ul/tpascal/scripts/getmoreparms.pl helix1_results/lisfiles $filebase";
$cmd_str .= "_Helix1_lis $startnum $endnum 1 " . ($strand_length * 2);
system "$cmd_str";
print "\n";

chdir "helix1_results";
MovePdbFiles();
MoveLisFiles();
chdir "../";

$cmd_str = "mv *Helix1*.pdb helix1_results/pdbfiles";
system "$cmd_str";

chdir "helix1_results";
print "\nCalculating strand shortening...\n";
print "--------------------------------\n";
RenamePDBFiles(0, "_Helix1");
$cmd_str = "/ul/maiti/bgfstr/linux/elastic " . ($strand_length * 2) . " 1 $firstfile 1 " . ($endnum - $startnum + 1);
chdir "pdbfiles";
if (system "$cmd_str > junk") {
    print "Failure!!\n";
} else {
    print "Sucess!!\n";
    system "mv junk ../$filebase" . "_Helix1_strandshortening.dat";
}
print "\n";

chdir "../";

RenamePDBFiles(1, "_Helix1");

chdir "../";



#helix2
$cmd_str = "/ul/tpascal/scripts/do_curve_anal_newtest.pl";
$cmd_str .= " $filebase $startnum $endnum 2 $strand_length";
system "$cmd_str";
print "\n";

$cmd_str = "/ul/tpascal/scripts/get_helical_parms.pl";
$cmd_str .= " helix2_results/" . $filebase . "_Helix2_lis $startnum $endnum";
system "$cmd_str";
print "\n";

chdir "helix2_results";
MovePdbFiles();
MoveLisFiles();
chdir "../";

$cmd_str = "/ul/tpascal/scripts/getmoreparms.pl helix2_results/lisfiles $filebase";
$cmd_str .= "_Helix2_lis $startnum $endnum 1 " . ($strand_length * 2);
system "$cmd_str";
print "\n";

chdir "helix2_results";
MovePdbFiles();
MoveLisFiles();

chdir "../";

$cmd_str = "mv *Helix2*.pdb helix2_results/pdbfiles";
system "$cmd_str";

chdir "helix2_results";
print "\nCalculating strand shortening...\n";
print "--------------------------------\n";
RenamePDBFiles(0, "_Helix2");
$cmd_str = "/ul/maiti/bgfstr/linux/elastic " . ($strand_length * 2) . " 1 $firstfile 1 " . ($endnum - $startnum + 1);
chdir "pdbfiles";
if (system "$cmd_str > junk") {
    print "Failure!!\n";
} else {
    print "Sucess!!\n";
    system "mv junk ../$filebase" . "_Helix2_strandshortening.dat";
}
print "\n";

chdir "../";

RenamePDBFiles(1, "_Helix2");

chdir "../";

MovePdbFiles();

print "\nCalculating RMSd values\n";
print "-----------------------\n";
print "Using $startnum to $ARGV[3] to calculate average structure\n";
$cmd_str = "/ul/tpascal/scripts/ptraj_rms.pl $filebase $trajfile $topfile $startnum " . $ARGV[3] . " 1 $strand_length";
#system "$cmd_str";


sub DisplayHelpFile() {
    my ($error_str) = "usage: px_anal.pl paramfile filebase startnum endnumber\n";

    $error_str .= "\nArguments:\n";
    $error_str .= sprintf("%-30s%-50s\n", "filebase", "The base name of the file");
    $error_str .= sprintf("%-30s%-70s\n", "startnum", "The lowest index of the files in the series");
    $error_str .= sprintf("%-30s%-70s\n", "endnum", "The highest index of the files");
    $error_str .= sprintf("%-30s%-70s\n", "majorgroove:minorgroove", "The ratio of the bases in the grooves, e.g. 6:5");
    $error_str .= sprintf("%-30s%-70s\n", "is3primeIn", "Is the helices arranged as 5' 3' 3' 5' from top to bottom");
    $error_str .= sprintf("%-30s%-70s\n", "numBasesAtEnd" , "The number of bases at the ends of the crossover moleculs");
    $error_str .= sprintf("%-30s%-70s\n", "strandLength", "The number of bases in each strand of the crossover molecule");
    $error_str .= sprintf("%-30s%-70s\n", "crossovers", "the base number of the crossovers, delimited by spaces");

    die "$error_str";
}

sub GetParameters() {

    open PARAMFILE, $paramfile or die "Cannot open $paramfile: $!\n";
    my ($in_data, $c_list);

    $is3prime = -1;

    while (<PARAMFILE>) {
	chomp;
        $in_data = $_;
     
	if ($in_data =~ /Mol:\s(\d+).(\d+)/) {
	    $groove_desc = $1 . ":" . $2;
	} elsif ($in_data =~ /Total bases:\s(\d+)/) {
	    $strand_length = int($1/4);
	} elsif ($in_data =~ /Bases at end:\s(\d+)/) {
	    $numEndBases = $1;
	} elsif ($in_data =~ /Crossovers:\s(.+)$/) {
	    $c_list = $1;
	    @crossovers = split(/ /, $c_list);
	} elsif ($in_data =~ /Topology file:\s(.+)$/) {
	    $topfile = $1;
	} elsif ($in_data =~ /Trajectory file:\s(.+)$/) {
	    $trajfile = $1;
	} elsif ($in_data =~ /3 prime in:\s([1|0])$/) {
	    $is3prime = $1;
	}
    }

    close PARAMFILE;
    if (! $numEndBases or ! $strand_length or ! $groove_desc or ($is3prime == -1) ) {
	die "Invalid Parameter file\n";
    }
}

sub CreatePdbfiles() {
    my ($out_cmd, $i, $j, $curr_pdb_file, $save_nm);
    my ($starttm, $num_precision);

    $num_precision = $endnum . "";
    $num_precision = length($num_precision);

    $starttm = sprintf("%0" . $num_precision . "d", $startnum);


    $out_cmd = "> tmp_gen_rst_file";

    open OUTCMD, $out_cmd or die "Cannot create file tmp_gen_rst_file: $!\n";
    print OUTCMD "trajin $trajfile $startnum $endnum $interval\n";

    $i = 1;
    $j = 1;
    my ($total_bases) = $strand_length * 4;

    while ($i < $total_bases) {

	print OUTCMD "center :1-" . ($j * $strand_length) . " mass origin\n";
	print OUTCMD "image origin center\n";
	$i +=$strand_length;
	$j++;
    }
    print OUTCMD "trajout " . $filebase . "_pdb pdb\n";

    close OUTCMD;

    $out_cmd = "ptraj $topfile < tmp_gen_rst_file > junk";
    
    system $out_cmd;
    system "rm -f junk tmp_gen_rst_file";

    $out_cmd = "/ul/tpascal/scripts/fix_pdb_name.pl $filebase" . "_pdb $startnum $endnum";
    system "$out_cmd";

    $i = $startnum;
    $j = $startnum;

    system "mkdir -p pdbfiles/completefiles";
    while ( $i <= $endnum) {
	$curr_pdb_file = $filebase . "_$i" . ".pdb";
	$save_nm = $filebase . "_$j" . ".pdb";
	push @valid_files, $save_nm;
	system "cp $curr_pdb_file pdbfiles/completefiles/$save_nm";
	system "/ul/tpascal/scripts/stripNaH20.pl $curr_pdb_file $save_nm";
	system "/ul/maiti/src/util/scripts/fixcurvepdb.pl $save_nm > junk";
	$i += $interval;
	$j++;
    }

    $endnum = $startnum + int( ($endnum - $startnum) / $interval );
    system "rm -f junk";

}

sub DetermineLegitInt(@) {
    my ($inval) = $_[0];
    my ($e_message) = $_[1];

    if (! $inval =~ /^\d+$/) {
	die "Invalid value in $e_message. Expected Integer\n";
    }
}

sub GetValidFiles() {
    my $i = 0;
    my $j = 0;
    my $k = 0;
    my $counter = 0;
    my $curr_fle = "";
    my $line_in = "";

    for ($i = $startnum; $i<=$endnum;$i++) {
	$curr_fle = $filebase . "_" . GetExtension($i) . ".pdb";
	if (-e $curr_fle) {
	    push @valid_files, $curr_fle;
	    $counter +=1;
	} else {
	    warn "Warning: Cannot find file $curr_fle\n";
	}
    }

    if ($counter ==0) {
	die "Cannot find any valid files\n";
    }
    
    print "Found $counter valid file(s)...\n";


}

sub GetExtension(@) {
    my $inval = $_[0];
    return $inval;
}

sub MovePdbFiles() {
    if (! -d "pdbfiles") {
	mkdir "pdbfiles", 0777;
	system "mv *.pdb pdbfiles";
    }
}

sub MoveLisFiles() {
    if (! -d "lisfiles") {
	mkdir "lisfiles", 0777;
	system "mv *_lis.* lisfiles";
    }
}

sub RenamePDBFiles(@) {
    my ($is_num_b4_pdb, $addinfo) = @_;
    my ($i);
    my ($curr_file, $new_file);


    if ($is_num_b4_pdb) {
	$i = $startnum;

	while( $i <= $endnum) {
	    $curr_file = "pdbfiles/" . $filebase . $addinfo . "_pdb." . $i;
	    $new_file =  "pdbfiles/" . $filebase . $addinfo . "_" . $i . ".pdb";
	    if (-e $curr_file) {
		system "mv $curr_file $new_file";
	    } else {
		system "pwd";
		print "ERROR: $curr_file does not exist\n";
	    }
	    $i += 1;
	}
    } else {
	$i = $startnum;

	while( $i <= $endnum) {
	    $new_file = "pdbfiles/" . $filebase . $addinfo . "_pdb." . $i;
	    $curr_file =  "pdbfiles/" . $filebase . $addinfo . "_" . $i . ".pdb";
	    if (-e $curr_file) {
		system "mv $curr_file $new_file";
	    } else {
		print "ERROR: $curr_file does not exist\n";
	    }
	    $i += 1;
	}
	$firstfile = $filebase . $addinfo . "_pdb." . $startnum;

    }
}

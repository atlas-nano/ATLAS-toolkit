#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename;
use General qw(FileTester);
use GetParms;
use HelixLayout;

sub DisplayHelpFile();
sub DetermineLegitInt(@);
sub GetValidFiles();
sub GetParameters();
sub CreatePdbfiles;
sub MovePdbFiles();
sub MoveLisFiles();
sub RenamePDBFiles(@);
sub GetHelixData(@);
sub CreateAvgStructure(@);
sub CreateNoCrossFile();
sub CreateHelix1Files();
sub init;
sub execCmd;

my ($i, $j, $cmd_str, $run_carnal);
my (@valid_files, $SELECT);

die "usage: $0 parmfile filebase trajSelect [multi-helix=yes]\n--help for the help file\n"
    if (! @ARGV || ($#ARGV < 2 && $ARGV[0] ne "--h"));

if ($ARGV[0] eq "--h") {
    DisplayHelpFile();
}

my ($paramfile, $filebase, $selection, $multi_helix) = @ARGV;
my ($groove_desc, $is3prime, $numEndBases, $strand_length); 
my(@crossovers, $topfile, $trajfile, $firstfile, $create_pdb);
my ($PARMS, $VFILES, $cmd_str);

$|++;

print "\n--===Begin STRUCTUREANAL===--\n\n";
print "Initializing...";
&init;
print "Done\nCreating PDB files...";
$VFILES = CreatePdbfiles($PARMS, $SELECT, $filebase);
print "Done\n";

if ($multi_helix) {
    print "Splitting Crossover structure into seperate helices...";
    for $i (0 .. $#{ $VFILES }) {
	$cmd_str = "/ul/tpascal/scripts/helix_from_strands.pl ";
	$cmd_str .= $VFILES->[$i] . " $PARMS->{Molecule}{major_groove}:$PARMS->{Molecule}{minor_groove} "; 
	$cmd_str .= $PARMS->{Molecule}{is3PrimeIn} . " " . $PARMS->{Molecule}{bases_at_end};
	for $j (0 .. $#{ $PARMS->{Molecule}{crossovers} }) {
	    $cmd_str .= " $PARMS->{Molecule}{crossovers}[$j]";
	}
	$cmd_str .= " > structure_anal_dump.txt";
	execCmd($cmd_str);
    }
    print "Done\n";
} else {
    CreateHelix1Files();
}

# Helix1
print "Writing Helix1 Data...";
system "echo  >> structure_anal_dump.txt";
system "echo --==Helix1==-- >> structure_anal_dump.txt";
system "echo  >> structure_anal_dump.txt";
GetHelixData(1);

print "Done\n";

if ($multi_helix) {
# Helix2
    print "Writing Helix2 Data...";
    system "echo  >> structure_anal_dump.txt";
    system "echo --==Helix2==-- >> structure_anal_dump.txt";
    system "echo  >> structure_anal_dump.txt";
    GetHelixData(2);
    
    print "Done\n";
    print "Creating Average Structure with no crossovers...";
    CreateNoCrossFile();
}


MovePdbFiles();

print "Done\n\nAll Tasks Completed\n";
print "\n--===End STRUCTUREANAL===--\n";

sub init {
    FileTester($parmFile);

    if ( ! defined($multi_helix) or $multi_helix ne "0") {
	$multi_helix = 1;
    } else {
	if ($multi_helix =~ /(\d+)/) {
	    if ($1 == 0) {
		$multi_helix = 0;
	    } else {
		$multi_helix = 1;
	    }
	} else { 
	    $multi_helix = 1;
	}
    }

    print "Getting parms from $parmFile...";
    $PARMS = GetParms->new();
    if (! $PARMS->IsValidParams($parmFile)) {
	die "Error in Paramater file\n";
    }
}

sub GetHelixData(@) {

    my ($current_helix) = $_[0];
    my ($curr_fbase)  = $filebase;
    my ($avg_stru);

    system "mkdir -p helix" . $current_helix . "_results/pdbfiles";
    $cmd_str = "mv *Helix" . $current_helix . "*.pdb helix";
    $cmd_str .= $current_helix . "_results/pdbfiles";
    system "$cmd_str";

    chdir "helix" . $current_helix . "_results";

    $cmd_str = "/ul/tpascal/scripts/do_curve_anal.pl";
    $cmd_str .= " pdbfiles/$filebase $startnum $endnum ";
    $cmd_str .= $current_helix . " $strand_length >> ../structure_anal_dump.txt";
    system "$cmd_str";

    $cmd_str = "/ul/tpascal/scripts/get_helical_parms_new.pl ";
    $cmd_str .= $filebase . "_Helix" . $current_helix;
    $cmd_str .= "_lis $startnum $endnum >> ../structure_anal_dump.txt";
    system "$cmd_str";
    
    MovePdbFiles();
    MoveLisFiles();
    
#   Calculate the base displacement
    $cmd_str = "/ul/tpascal/scripts/calc_base_disp.pl ";
    $cmd_str .= "lisfiles/$filebase" . "_Helix" . $current_helix . "_lis";
    $cmd_str .= " $startnum $endnum >> ../structure_anal_dump.txt";
    system "$cmd_str";

    $cmd_str = "/ul/tpascal/scripts/getmoreparms.pl ";
    $cmd_str .= "lisfiles $filebase";
    $cmd_str .= "_Helix" . $current_helix ."_lis $startnum $endnum ";
    $cmd_str .= "1 " . ($strand_length * 2) . " >> ../structure_anal_dump.txt";
    system "$cmd_str";
    
#    RenamePDBFiles(0, "_Helix" . $current_helix);
#    $cmd_str = "/ul/maiti/bgfstr/linux/elastic " . ($strand_length * 2);
#    $cmd_str .= " 1 $firstfile 1 " . ($endnum - $startnum + 1);
#    chdir "pdbfiles";
#    if (system "$cmd_str >& junk1") {
#	print "Failure while calculating strand shortening ";
#    } else {
#	system "mv junk1 ../$filebase" . "_Helix1_strandshortening.dat";
#    }
#    chdir "../";
#    RenamePDBFiles(1, "_Helix" . $current_helix);

#   Calculate the CRMS per base
    $avg_stru = CreateAvgStructure($current_helix);
    $cmd_str = "/ul/tpascal/scripts/carnal_rms.pl ";
    $cmd_str .= " ../$paramfile $startnum $endnum";
    if ($run_carnal) {
	system "$cmd_str >> ../structure_anal_dump.txt";
    }

    chdir "../";
}

sub CreateNoCrossFile() {
    my ($avg_file) =  basename($filebase) . "_avg.pdb";
    my ($dna_string, $solute_string);

    die "Cannot locate file $avg_file" if (! -e $avg_file);

    system "mkdir -p avg_structure";
    system "cp $avg_file avg_structure/";

    open INFILE, $avg_file or die "Cannot open $avg_file: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /Na|WAT|Mg/i) {
	    $solute_string .= "$_\n";
	} else {
	    $dna_string .= "$_\n";
	}
    }
    
    close INFILE;

    open OUTFILE, "> tmp_dna_file" or die "Cannot create tmp_dna_file: $!\n";
    print OUTFILE $dna_string;
    close OUTFILE;

#    system "/ul/maiti/src/util/scripts/fixcurvepdb.pl tmp_dna_file>> structure_anal_dump.txt";
    my ($cmd_str) = "/ul/tpascal/scripts/helix_from_strands.pl ";
    $cmd_str .= "tmp_dna_file $groove_desc $is3prime $numEndBases ";
    for $j (0 .. $#crossovers) {
	$cmd_str .= $crossovers[$j] . " ";
    }
    system "$cmd_str >> structure_anal_dump.txt";
    system "rm -f tmp_dna_file*";
    $cmd_str = "mv pdbfiles/nocrossovers/tmp_dna_file avg_structure/";
    $cmd_str .= basename($filebase) . "_nocross.pdb";
    system $cmd_str;

    $avg_file = "avg_structure/" . basename($filebase) . "_solvent.pdb";
    open OUTFILE, ">> $avg_file" or die "Cannot write to $avg_file: $!\n";
    print OUTFILE $solute_string;
    close OUTFILE;
#    system "/ul/tpascal/scripts/fixxleappdb.pl $avg_file >> structure_anal_dump.txt";

}
sub CreateAvgStructure(@) {
    my ($curr_helix) = $_[0];
    my ($counter, $avg_structure, $out_cmd);

    $avg_structure = basename($filebase) . ".pdb";
    $out_cmd = "> tmp_gen_rst_file";

    open OUTCMD, $out_cmd or die "Cannot create file tmp_gen_rst_file: $!\n";
    print OUTCMD "trajin $trajfile $startnum $endnum $interval\n";

    for $counter (1 .. 4) {
        print OUTCMD "center :1-" . ($counter * $strand_length) . " mass origin\n";
        print OUTCMD "image origin center\n";
    }

    print OUTCMD "average $avg_structure pdb\n";

    close OUTCMD;

    $out_cmd = "ptraj $topfile < tmp_gen_rst_file >& junk";

    if ( system "$out_cmd") {
	system "rm -f tmp_gen_rst_file junk";
	die "Failure!! Error executing ptraj\n";
    } else {
	system "rm -f tmp_gen_rst_file junk";
	system "cp $avg_structure ../" . basename($filebase) . "_avg.pdb";
	system "/ul/tpascal/scripts/stripNaH20.pl $avg_structure";
	return $avg_structure;
    }

}
sub DisplayHelpFile() {
    my ($error_str) = "usage: $0 parmfile filebase trajSelection [multi-helix=yes]\n";

    $error_str .= "\nArguments:\n";
    $error_str .= sprintf("%-30s%-50s\n", "parmfile", "The parameter file");    
    $error_str .= sprintf("%-30s%-50s\n", "filebase", "The prefix for the saved files");
    $error_str .= sprintf("%-30s%-50s\n", "trajSelection", "trajectory selection. ex \":Ii1-10:2\" frames 1 - 10 every 2\n");
    $error_str .= sprintf("%-30s%-50s\n", " ", "can have multiple selections within the quotes. for all use \"*\"\n");
    $error_str .= sprintf("%-30s%-70s\n", "multi-helix", "Optional. Specify whether to split structure into different helices. Default=yes\n");

    die "$error_str";
}

sub GetParameters() {

    my ($has_reference, $has_chains);

    open PARAMFILE, $paramfile or die "Cannot open $paramfile: $!\n";
    my ($in_data, $c_list);

    $is3prime = -1;

    while (<PARAMFILE>) {
	chomp;
        $in_data = $_;
     
	if ($in_data =~ /Mol:\s(\d+).(\d+)/) {
	    $groove_desc = $1 . ":" . $2;
	} elsif ($in_data =~ /Total bases:\s(\d+)/) {
	    if ($multi_helix == 1) { 
		$strand_length = int($1/4); 
            } else { 
		$strand_length = int($1/2);
	    }
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
	} elsif ($in_data =~ /Create PDB: ([1|0])/) {
	    $create_pdb = $1;
	} elsif ($in_data =~ /Reference: (\w+)/) {
	    if (-e $1) {
		$has_reference = 1;
	    }
	} elsif ($in_data =~ /Chains: (\d+)/) {
	    $has_chains = 1;
	}
    }

    close PARAMFILE;

    if ($has_chains && $has_reference) {
	$run_carnal = 1;
    } else {
	$run_carnal = 0;
    }

    if (! $numEndBases or ! $strand_length or ! $groove_desc or ($is3prime == -1) ) {
	die "Invalid Parameter file\n $numEndBases, $strand_length, $groove_desc, $is3prime\n";
    }

    die "Cannot find files: $!\n"
	if (! $topfile or ! $trajfile);
    if (! -e $topfile or ! -e $trajfile or ($#crossovers == -1 and $multi_helix == 1)) {
	die "Invalid Parameter file\n$topfile $trajfile $#crossovers\n";
    }
}

sub CreatePdbfiles {    
    my ($parms, $select, $filebase) = @_;
    my (@tmp, $start, $end, $tmp_file, $out_cmd, @FILES);

    my ($out_cmd, $i, $j, $curr_pdb_file, $save_nm);
    my ($starttm, $num_precision);

    @tmp = sort numerically keys %{ $SELECT };
    $start = $tmp[0];
    $end = $tmp[$#tmp];

    if ($parms->{Amber_Options}{run_ptraj} == 0) {
	print "Assuming PDB files already exists\n";
    } else {
	$precision =  sprintf("%0" . length($end) . "d", $start);
	$tmp_file = "tmp_gen_rst_file";
    
	open OUTCMD, "> $tmp_file" || die "ERROR: Cannot create file $tmp_file: $!\n";
	print OUTCMD "trajin $trajfile $startnum $endnum $interval\n";
	print OUTCMD "center :1-" . $parms->{Molecule}{total_bases} . " mass origin\n";
	print OUTCMD "image origin center\n";
	print OUTCMD "trajout " . $filebase . "_pdb pdb\n"; 
	close OUTCMD;
	
	$out_cmd = "ptraj $topfile < tmp_gen_rst_file >& structure_anal_dump.txt";
	execCmd($out_cmd, $tmp_file);
	
	$out_cmd = "/ul/tpascal/scripts/fix_pdb_name.pl $filebase";
	$out_cmd .= "_pdb $start $end > structure_anal_dump.txt";
	execCmd($out_cmd);
	
    }	

    $j = $start;
    for $i (@tmp) {
	$currFile = $filebase . "_$i" . ".pdb";
	$save_nm = $filebase . "_$j" . ".pdb";
	next (if (! -e $currFile || ! -r $currFile || ! -T $currFile));
	push @FILES, $save_nm;
	if ($parms->{Amber_Options}{run_ptraj} != 0) {
	    $out_cmd =  "/ul/tpascal/scripts/stripNaH20.pl $curr_pdb_file $save_nm";
	    execCmd($out_cmd);

	    $outCmd = "/ul/maiti/src/util/scripts/fixcurvepdb.pl $save_nm > structure_anal_dump.txt";
	    execCmd($out_cmd);
	}
	$j++;
    }    
    die "ERROR: No valid PDB files found!\n" if (! @FILES);

    return \@FILES;
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
	die "Cannoting find any valid files\n";
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
    }
    system "mv *.pdb pdbfiles";   
}

sub MoveLisFiles() {
    if (! -d "lisfiles") {
	mkdir "lisfiles", 0777;
    }
    system "mv *_lis.* lisfiles";
}

sub RenamePDBFiles(@) {
    my ($is_num_b4_pdb, $addinfo) = @_;
    my ($i, $curr_file, $new_file);


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

sub CreateHelix1Files() {
    my ($old_name, $new_name);

    for $old_name (@valid_files) {
	if ($old_name =~ /^(.+)_(\d+)\.pdb/) {
	    $new_name = $1 . "_Helix1_" . $2 . ".pdb";	    
	    system "cp $old_name $new_name";
	}
    }

}

sub execCmd {
    my ($cmd, $fileToDelete) = @_;
    
    die "Error while executing $cmd\n" if (system($cmd));
    system "rm -f $fileToDelete" if (defined($fileToDelete);
}

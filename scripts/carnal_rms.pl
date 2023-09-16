#!/usr/bin/perl -w
use strict;
use File::Basename;

sub ValidateFields();
sub CreateCarnalFile();
sub CreatePtrajFile();
sub RunCmd();
sub CalculateStats();
sub GetParms();
sub CreateTable();

if (! @ARGV) {
    die "usage: $0 parmfile\n"; #trajfile topfile no_base_pairs no_chains reference\n";
}

my ($parmfile) = @ARGV;
my ($trajfile, $topfile, $no_chains, $reference, @restraints);
my ($chain_length, $rms_template, $no_atoms, %RMS_ARRAY, $no_bases, $other_trj_files);
my ($start, $end) = ("", "");

$other_trj_files = "";
print "Getting Paramaters...";
GetParms();
ValidateFields();
print "Sucess\nCreating Carnal File...";
#CreateCarnalFile();
CreatePtrajFile();
print "Sucess\nRunning Carnal...";
RunCmd();
print "Sucess\nFinalizing Data...";
CalculateStats();
print "Sucess\nCreating Table...";
CreateTable();
print "Sucess\n\nAll tasks completed - Ending Execution\n";

system "rm -f out tmp_ptraj_file";
sub CreateTable() {
    my ($i, $j, $currkey, $outdata, $cfile, $avgdata);

    for $i (1 .. $no_bases) {
	$currkey = sprintf("%0" . length($no_bases) . "d", $i);
	if ($RMS_ARRAY{$currkey}->{"HASDATA"}) {
	    $outdata .= sprintf("%-3s", "$i");
	    $avgdata .= sprintf("%-3s", "$i");
	    for $j (keys %{ $RMS_ARRAY{$currkey}->{"RMSDATA"} }) {
		$outdata .= sprintf("%7.3f", $RMS_ARRAY{$currkey}->{"RMSDATA"}->{$j});
	    }
	    $outdata .= sprintf("%7.3f", $RMS_ARRAY{$currkey}->{"AVERAGE"});
	    $outdata .= sprintf("%8.4f", $RMS_ARRAY{$currkey}->{"STDEV"});
	    $outdata .="\n";

	    $avgdata .= sprintf("%7.3f", $RMS_ARRAY{$currkey}->{"AVERAGE"});
	    $avgdata .= sprintf("%8.4f", $RMS_ARRAY{$currkey}->{"STDEV"});
	    $avgdata .="\n";
	}
    }

    $cfile = basename($topfile) . ".dat";
    open TABLEFILE,"> $cfile" or die "Cannot create datafile $cfile: $!\n";
    print TABLEFILE "$outdata";
    close TABLEFILE;

    $cfile = basename($topfile) . "_avg_stdev.dat";
    open TABLEFILE,"> $cfile" or die "Cannot create datafile $cfile: $!\n";
    print TABLEFILE "$avgdata";
    close TABLEFILE;
}

sub CreatePtrajFile() {
    my ($i, $fname, $counter, $isfound, $bases_restrained, $outInfo);

    system "mkdir -p datfiles";
    $fname = "datfiles/" . basename($topfile); 
    open PTRAJFILE, "> tmp_ptraj_file" or die "Cannot create temporary file: $!\n";
    print PTRAJFILE "$other_trj_files";
    print PTRAJFILE "trajin $trajfile $start $end\n";
    print PTRAJFILE "reference $reference\n";
    print PTRAJFILE "strip :Na+\nstrip :WAT\n";

    for $i (1 .. $no_bases) {
	print PTRAJFILE "rms reference mass out $fname" . "_";
	print PTRAJFILE $i . "bp.dat :" . $i . " \n";
    }

# Create listing for the restrained, unrestrained and all residues

    $outInfo = "";
    print PTRAJFILE "rms reference mass out datfiles/All_res.dat :1-" . $no_bases . "\n";
    if ($#restraints > -1) {
	print PTRAJFILE "rms reference mass out datfiles/Restraints.dat :";
	for $bases_restrained (@restraints) {
	    $outInfo .= "$bases_restrained,";
	}
	chop $outInfo;
	print PTRAJFILE "$outInfo\n";
	$outInfo = "";
	print PTRAJFILE "rms reference mass out datfiles/UnRestrained.dat :";
	for $counter (1 .. $no_bases) {
	    $isfound = 0;
	    for $bases_restrained (@restraints) {
		if ($bases_restrained == $counter) {
		    $isfound = 1;
		    last;
		}
	    }
	    if (! $isfound) {
		$outInfo .= "$counter,";
	    }
	}
	chop $outInfo;
	print PTRAJFILE "$outInfo\n";

    }

    close PTRAJFILE;
}


sub GetParms() {
    -e $parmfile or die "Cannot locate parameter file $parmfile: $!\n";
    open INFILE, $parmfile or die "Cannot open $parmfile: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /Trajectory file: (.+)$/) {
	    if ($trajfile) {
        	$other_trj_files .= "trajin $trajfile\n";
            }
	    $trajfile = $1;
	} elsif ($_ =~ /Topology file: (.+)$/) {
	    $topfile = $1;
	} elsif ($_ =~ /Reference file: (.+)$/) {
	    $reference = $1;
	} elsif ($_ =~ /Total bases: (.+)/) {
	    $no_bases = $1;
	} elsif ($_ =~ /Chains: (.+)/) {
	    $no_chains = $1;
	} elsif ($_ =~ /Start Step: (.+)/) {
	    $start = $1;
	} elsif ($_ =~ /End Step: (.+)/) {
	    $end = $1;
	} elsif ($_ =~/Restraints: /) {
	    while ($_ =~ /(\d+)/g) {
		push @restraints, $1;
	    }
	}	
    }
    close INFILE;

    if (! $trajfile || ! $topfile || ! $reference || ! $no_bases || ! $no_chains) {
	die "ERROR: Could not parse paramater file $ARGV[0]\n";
    }
}


sub CalculateStats() {
    my ($i, $bname, $fname, $currbase, $counter, $isvalid);
    my ($avg, $stdev, $validbase);

    $isvalid = 0;
    $bname = "datfiles/" . basename($topfile); 
    for $i (1 .. $no_bases) {
	$fname = $bname . "_" .  $i . "bp.dat";
	$currbase = sprintf("%0" . length($no_bases) . "d", $i);
	if ($start) {
	    $counter = $start;
	} else {
	    $counter =1;
	}
	$validbase = 0;
	if (open INFILE, $fname) {
	    while (<INFILE>) {
		if ($_ =~ /^\s+(\d+)\.\d+\s+(\d+\.\d+)/) {
#		    if ($1 == $counter) {
			$RMS_ARRAY{$currbase}->{"RMSDATA"}->{$1} = $2;
			$RMS_ARRAY{$currbase}->{"DATA"} .= "$2 ";
			$isvalid = 1;
			$validbase = 1;
#		    }
		    $counter++;
		}
	    }
	    close INFILE;
	    $RMS_ARRAY{$currbase}->{"HASDATA"} = $validbase;
	    if ($validbase) {
		($avg, $stdev) = STDev($RMS_ARRAY{$currbase}->{"DATA"});
		$RMS_ARRAY{$currbase}->{"AVERAGE"} = $avg;
		$RMS_ARRAY{$currbase}->{"STDEV"} = $stdev;
	    }
	} else {
	    print "WARNING: Cannot locate datafile for base pair $i in $fname: $!\n";
	}
    }
    if (! $isvalid) {
	die "ERROR: Could not find any valid data for the bases\n";
    }
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
    return ($avg, $result);
}

sub RunCmd() {
    -e "tmp_ptraj_file" or die "Cannot locate temporary file tmp_rms: $!\n";
    if (system "/ul/maiti/ptraj-6.3/linux/ptraj $topfile < tmp_ptraj_file >& out ") {
	die "Cannot execute ptraj\n";
    }
}


sub CreateCarnalFile() {
# Modification - well it turns out that Carnal is too slow too do this
# so use Ptraj which is faster
    my (@tempdata, $in_text, $tblname, $groupinfo, $tblinfo, $i, $counter);
    my ($trajtext);

    if (! $start or ! $end) {
	$trajtext = $trajfile;
    } else {
	$trajtext = "WIN " . $start . " ";
	$trajtext .= abs($end - $start) . " $trajfile";
    }
    
    $tblname = basename($topfile) . "_rms";
    $tblinfo = " TABLE tbl";
    $groupinfo = " GROUP grp1 (RES 1 - $no_bases);\n";
    $groupinfo .= " RMS fit1 FIT grp1";
    $groupinfo .= " s1 ref_set RES";
    $tblinfo .= " fit1";

    open RMSFILE, $rms_template or die "Cannot open $rms_template: $!\n";
    while (<RMSFILE>) {
	chomp;
	$in_text = $_;
	$in_text =~ s/parmfile/$topfile/;
	$in_text =~ s/trajfile/"ATOM $no_atoms $trajtext"/;
	$in_text =~ s/reffile/$reference/;
	$in_text =~ s/tablename/$tblname/;
	$in_text =~ s/groupinfo/$groupinfo/;
	$in_text =~ s/tableinfo/$tblinfo/;
	push @tempdata, $in_text;
    }
    close RMSFILE;

    open OUTFILE, "> tmp_rms" or die "Error creating temp file: $!\n";
    for $i (@tempdata) {
	print OUTFILE "$i\n";
    }
    close OUTFILE;
}

sub ValidateFields() {

    $no_atoms = 0;
    $rms_template = "/ul/tpascal/scripts/carnal_rms_template";
    -e $rms_template or die "Cannot locate template " . \
	"file $rms_template: $!\n";
    -e $trajfile or die "Cannot locate trajectory file $trajfile: $!\n";
    -e $topfile or die "Cannot locate topology file $topfile: $!\n";
    -e $reference or die "Cannot locate reference file $reference: $!\n";
    
    if (! $no_bases =~ /\d+/) {
	die "Expected integer for no_bases, got: $no_bases\n";
    }
    
    if ($no_bases < 2) {
	die "Invalid no_bases: Must be larger than 2. Got: $no_bases\n";
    }

    if (! $no_chains =~ /(2|4)/) {
	die "Invalid number in no_chains. Expected 2 or 4 got $no_chains\n";
    }

    if ($no_bases % $no_chains > 0) {
	die "Error in no_bases: Must be integer of $no_chains\n";
    }
}

#!/usr/bin/perl -w
use strict;

if (!@ARGV) {
    die "Usage: $0 inputfile [avg_structure] [outfolder]\n";
}

my ($infile, $avg_file, $outfolder) = @ARGV;
my ($isfirstmodel, $curr_model, %modelsinfo, $headerinfo);
my ($models_key, $arr_data, $counter, $precision, @validfiles);
my (@rms_data, $out_line, $start_i, $end_i);

sub FixLine(@);
sub CalcRMS();
sub numerically { ($a<=>$b); }

-e $infile or die "$infile: $!\n";

open INFILE, $infile or die "Cannot open $infile: $!\n";
print "Examining $infile...";

$isfirstmodel = 0;
$counter = 0;
$precision = 1;

while (<INFILE>) {
    chomp;
    if (! $isfirstmodel) {
	if ($_ =~ /^MODEL/) {
	    $isfirstmodel = 1;
	    if ($_ =~ /^MODEL\s+(\d+)/) {
		$curr_model = $1;
	    } else {
		$curr_model = 1;
	    }
	    $counter = 1;
	} else {
	    $headerinfo .= "$_\n";
	}
    } elsif ($_ =~ /^MASTER/) {
	last;
    } else {
	if ($_ =~ /^MODEL/) {
	    if ($_ =~ /^MODEL\s+(\d+)/) {
		$curr_model = $1;
	    } else {
		$curr_model += 1;
	    }
	    $counter++;
	    if (length($curr_model) > $precision) {
		$precision = length($curr_model);
	    }
	} else {
	    push @{ $modelsinfo{$curr_model} }, $_;
	}
    }
}

close INFILE;

print "Done\nFound $counter unique molecules.\n";

if ($outfolder) {
    if (! -d $outfolder) {
	mkdir $outfolder, 0777;
    }
    chdir $outfolder;
} else {
    system "mkdir -p Models";
    chdir "Models";
}

print "Creating Models...";
$counter = 0;

for $models_key (sort numerically keys %modelsinfo) {
    $curr_model = "Model_";
    $curr_model .= sprintf("%0" . $precision . "d", $models_key);
    $curr_model .= ".pdb";

    if (open OUTFILE, "> $curr_model") {
	print OUTFILE "$headerinfo";
	for $arr_data (@{ $modelsinfo{$models_key} }) {
	    if ($arr_data ne "") {
		print OUTFILE "$arr_data\n";
	    }
	}
	$counter++;
	push @validfiles, substr($curr_model,0, -4);
	close OUTFILE;
    } else {
	print "WARNING: Could not create $curr_model: $!\n";
    }
}
if ($counter > 0) {
    print "Done. Created $counter models\n";
} else {
    print "Error. No files written\n";
}

chdir "../";
CalcRMS();

print "Creating table...";
$precision += 2;
$out_line = sprintf("%".  $precision . "s", " ");
for $arr_data (@validfiles) {
    $out_line .= sprintf("%" . $precision . "s", $arr_data);
}
$out_line .= "\n";

for $models_key (0 .. $#validfiles - 1) {
    $out_line .= sprintf("%" . $precision . "s", $validfiles[$models_key]);
    $start_i = ($models_key * $#validfiles);
    $end_i = $start_i + $#validfiles;
    for $arr_data ($start_i .. $end_i) {
	$out_line .= sprintf("%" . $precision . "s", $rms_data[$arr_data]);
    }
    $out_line .= "\n";
}


open RMS_FILE, "> rms_data" or die "Cannot create rms_data: $!\n";
print RMS_FILE "$out_line";
print "Done\n";

sub FixLine(@) {
    my ($indata) = $_[0];
    my ($returnval) = "";

    if ($indata =~ /^ATOM/) {
	if ($indata =~ /^(ATOM\s+\d+\s+[C|N|O|P].{2}\s+)\s(\w\s\w\s+\d+\s+.+)/){
	    $returnval = $1 . "D" . $2;
	}
    } else {
	$returnval = $indata;
    }

    return $returnval;

}

sub CalcRMS() {
    my ($i, $j);

    $precision = 1;

    if ($avg_file) {
	system "cp $avg_file Model_avg.pdb";
	$avg_file = "Model_avg";
	push @validfiles, $avg_file;
    }

    print "Calculating RMSD...";

    open TMPFILE, "> tmpfile" or die "Cannot create tmpfile: $!\n";
    for $i (@validfiles) {
	if (length($i) > $precision) {
	    $precision = length($i);
	}

	for $j (@validfiles) {
	    print TMPFILE "_GUI/MODAL 1\n";
	    print TMPFILE "_SYSTEM/REINIT_ALL\n";
	    print TMPFILE "FILES/LOAD_FORMAT   PDB\n";
	    print TMPFILE "FILES/LOAD  \"Models/$i" . ".pdb\"\n";
	    print TMPFILE "FILES/LOAD  \"Models/$j" . ".pdb\"\n";
	    print TMPFILE "MOVE/MATCH_WEIGHT  MASS\n";
	    print TMPFILE "MOVE/REFERENCE_MODEL  \"$i\"\n";
	    print TMPFILE "MOVE/FIT_MODEL  \"$j\"\n";
	    print TMPFILE "MOVE/MATCH_MODELS\n";
	}
    }
    close TMPFILE;
#    run command
    if (open RMS_FILE, "cerius2 -n tmpfile |") {
	while (<RMS_FILE>) {
	    chomp;
#		    print "$_\n";
	    if ($_ =~ /^RMS fit of \d+ atoms: (.+)/) {
		push @rms_data, $1;
#		print "$1\n";
		if (length($1) > $precision) {
		    $precision = length($1);
		}
	    }
	}
	print "Done\n";
    } else {
	die "ERROR: Cannot run cerius2: $!\n";
    }
}

#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}
                                                                                  

use strict;
use File::Basename;
use Packages::General;

# This script will take a jaguar output, create a pqr file and run apbs
# It will then collect and tabulate the differences between the the two

sub GetAllFiles(@);
sub ValidateInput();
sub ExtractName(@);
sub RunApbs();
sub RunJaguar(@);
sub WriteOutput(@);
sub Numercially;
sub ComputeACC(@);
sub GetMolNames(@);

die "Usage: $0 directory|file is_jaguar_input use_esp calculate_epsilon [mol_name] [use fine grids]\n"
    if (!@ARGV or $#ARGV < 3);

my ($location, $is_input, $use_esp, $calc_epsilon, $mol_name, $use_fine) = @ARGV;
my (%MOL_INFO, %file_list);
my ($outfile);


ValidateInput();
RunApbs();
print "Wrote data $outfile\n";
#system "rm -f *.log *apbs* *.in jagstatus io.mc";

sub ValidateInput() {
    if (-d $location) {
	%file_list = GetAllFiles($location);
    } elsif (-e $location) {
	$file_list{"0"} = (
			   {
			       "FileLocation" => $location,
			       "FileName" => ExtractName($location),
			   }
			   );
    }else {
	die "ERROR: Please specify a valid directory or Jaguar output file!\n";
    }
    $calc_epsilon = 0
	if ($calc_epsilon !~ /^[1|0]$/);
    $use_esp = 1
	if ($use_esp !~ /^[1|0]$/);
    if ($#ARGV == 2) {
	$use_fine = 1;
    }
    $use_fine = 1
	if (! $use_fine or $use_fine !~ /^[1|0]$/);

}

sub GetAllFiles(@) {
    my ($curr_dir) = $_[0];
    my ($counter, %file_list);

    opendir THISDIR, $curr_dir or die "Cannot read $curr_dir: $!\n";

    my (@allfiles) = grep { /^\./ && -f } map { "$curr_dir/$_" } readdir THISDIR;
    closedir THISDIR;

    $counter = 0;
    for (@allfiles) {
	if ($_ =~ /\.[out|in]/) {
	    $file_list{$counter} = (
				    {
					"FileLocation" => $_,
					"FileName" => ExtractName($_),
				    }
				    );
	    $counter++;
	}
    }

    die "Cannot find any relevant file in $curr_dir. Please ensure that the files have .out extensions\n"
	if ($counter == 0);
    return %file_list;
}
	
sub ExtractName(@) {
    my ($result) = basename($_[0]);
    $result =~ s/\.out$//;
    $result =~ s/\.in$//;
    return $result;
}

sub RunApbs() {
    my ($curr_file, $isvalid, $out_cmd, $counter, $file_loc);
    my ($was_sucessful) = 0;

    for $counter (sort Numerically keys %file_list) {
	$curr_file = $file_list{$counter}{"FileName"};
	$file_loc = $file_list{$counter}{"FileLocation"};
	print "Performing calculations on $curr_file...";
	if ($is_input) {
	    $file_loc = RunJaguar($file_loc);
	}
	next if (! $file_loc);
	$out_cmd = "/ul/tpascal/scripts/apbs/jag2apbs.pl $file_loc $use_esp $curr_file";
	if (open CONVERT, "$out_cmd |") {
	    $isvalid = 0;
	    while (<CONVERT>) {
		if ($_ =~ /^\s+(\w+):\s+(\-?\d+\.\d+)\s+(.+)$/) {
		    $isvalid = 1;
		    $MOL_INFO{$curr_file}{"JAGUAR"}{$1} =  $2;
		    $MOL_INFO{"UNITS"}{$1} = $3; # units
		}
	    }
	    close CONVERT;
	} else {
	    print "Error while running $out_cmd: $!\n";
	}
	
	if (! $isvalid) {
	    delete $MOL_INFO{$curr_file};
	    print "WARNING: Failed to obtain valid information while running $out_cmd. Ignoring and continuing\n";
	} else {
	    GetMolNames($mol_name);
	    $isvalid = 0;
	    my ($molv) = ComputeACC($curr_file . ".pqr", 1); #MOLECULAR VOLUME
	    $MOL_INFO{$curr_file}{"APBS"}{"MOL_VOL"} = $molv;
#	    $molv = 2.00;
	    $out_cmd = "/ul/tpascal/scripts/apbs/doapbs.pl $curr_file" . ".pqr $use_fine $calc_epsilon $molv >& _apbstemp";
	    system ($out_cmd);
	    if (open APBS, "_apbstemp") {
		while(<APBS>) {
		    if ($_ =~ /^Electrostatic energy:\s+(\-?\d+\.\d+.+)\s+kcal\/mol$/) {
			$isvalid = 1;
			$MOL_INFO{$curr_file}{"APBS"}{"ELECTROSTATICS"} = $1;
		    }elsif ($_ =~ /Using Solute Dielectric Constant of (\d\.*\d*)/) {
			$MOL_INFO{$curr_file}{"APBS"}{"DIELECTRIC"} = $1;
		    }
		}
		close APBS;
		if (! $isvalid) {
		    print "Invalid data from $out_cmd. Ignoring\n";
		} else {
		    $MOL_INFO{$curr_file}{"APBS"}{"SASA"} = ComputeACC($curr_file . ".pqr", 0); # SASA unweighted
		    WriteOutput($was_sucessful);
		    $was_sucessful = 1;
		    print "Done\n";   
		}
		delete $MOL_INFO{$curr_file};
	    }else {
		print "Error while running $out_cmd: $!\n";
	    }
	}
    }
    system "rm -f _apbstemp";
    die "ERROR: No valid files found\n"
	if (! $was_sucessful);
}

sub ComputeACC(@) {
    my ($in_file, $is_weighted) = @_;
    my ($out_cmd, $in_text, $returnstr);

    $returnstr = 0;
    $out_cmd = "/ul/tpascal/bin/acc ";
    if ($is_weighted) {
	$out_cmd .= "-mol ";
    } else {
	$out_cmd .= "-sasa ";
    }
    $out_cmd .= "$in_file >& _sasa";

    if (! system($out_cmd)) {
	if(open(SASA, "_sasa")) {
	    while (<SASA>) {
		chomp;
		if ($is_weighted and $_ =~ /\s+Approx. volume = (\d+\.\d+)\s+A/) {
		    $returnstr = $1;
		} elsif (! $is_weighted and $_ =~ /\s+Approx. SASA = (\-?\d+\.\d+\S+)\s+A/) {
		    $returnstr = $1;
		}
	    }
	}
	close SASA;
    }
#    system "rm -f _sasa";
    return $returnstr;
}

sub WriteOutput(@) {

    my ($counter, $diff, $out_text, $counter1, $out_file);
    my ($not_first_pass) = $_[0];

    if ($use_esp) {
	$out_file = "> comparison_esp";
    } else {
	$out_file = "> comparison_mulliken";
    }

    if ($calc_epsilon) {
	$out_file .= "_epsilon";
    }

    $out_file .= ".dat";

    if (! $not_first_pass) {
	$out_text = sprintf("%40s%45s%50s%10s%10s%15s\n", " ", CenterText("JAGUAR", 45), 
			    CenterText("APBS", 50), " ", " ", "EXP");
	
	$out_text .= sprintf("%40s", " ");
	for $counter (keys %{ $MOL_INFO{"UNITS"} }) {
	    if ($counter ne "NAME" and $counter ne "EXP") {
		if ($counter eq "ELECTROSTATICS") {
		    $out_text .= sprintf("%15s", $counter);
		} else {
		    $out_text .= sprintf("%10s", $counter);
		}
	    }
	}
	$out_text .= sprintf("%20s%10s%10s%10s%10s%10s%15s\n", "ELECTROSTATICS", "MOL_VOL", "SASA", "DIELEC", "DIFF", "DIFF %", " ");

	$out_text .= sprintf("%30s", CenterText("MOLECULE", 15));
	for $counter (keys %{ $MOL_INFO{"UNITS"} }) {
	    if ($counter ne "UNITS" and $counter ne "EXP") {
		if ($counter eq "ELECTROSTATICS") {
		    $out_text .= sprintf("%15s", $MOL_INFO{"UNITS"}{$counter});
		} else {
		    $out_text .= sprintf("%10s", $MOL_INFO{"UNITS"}{$counter});
		}
	    }
	}
	$out_text .= sprintf("%20s%10s%10s%10s%10s%10s%15s\n", "kcal/mol", "[A2]", "[A2]", " ", " ", " ", "kcal/mol");
    } else {
	$out_file = ">" . $out_file;
    }	

    for $counter (keys %MOL_INFO) {
	if ($counter ne "UNITS") {
	    $diff = $MOL_INFO{$counter}{"APBS"}{"ELECTROSTATICS"} - 
		$MOL_INFO{$counter}{"JAGUAR"}{"ELECTROSTATICS"};
	    $out_text .= sprintf("%30s", $MOL_INFO{$counter}{"NAME"});
	    for $counter1 (keys %{ $MOL_INFO{$counter}{"JAGUAR"} }) {
		if ($counter1 eq "ELECTROSTATICS") {
		    $out_text .= sprintf("%15.5G", $MOL_INFO{$counter}{"JAGUAR"}{$counter1});
		} else {
                    $out_text .= sprintf("%10.5G", $MOL_INFO{$counter}{"JAGUAR"}{$counter1});
		}
	    }
            for $counter1 (keys %{ $MOL_INFO{$counter}{"APBS"} }) { 
                if ($counter1 eq "ELECTROSTATICS") { 
                    $out_text .= sprintf("%20.5G", $MOL_INFO{$counter}{"APBS"}{$counter1});
                } else {
                    $out_text .= sprintf("%10.5G", $MOL_INFO{$counter}{"APBS"}{$counter1});
                } 
            } 

	    $out_text .= sprintf("%10.3G%10.2f%15.5f\n", $diff, abs($diff * 100/$MOL_INFO{$counter}{"APBS"}{"ELECTROSTATICS"}), $MOL_INFO{$counter}{"EXP"});
	}
    }
	
    open OUTDATA, $out_file or die "Cannot write to $out_file: $!\n";
    print OUTDATA $out_text;
    close OUTDATA;

    $outfile = $out_file;
}

sub RunJaguar(@) {
    my ($in_file) = $_[0];
    my ($return_str, @infile_text, $counter, $insert_loc, $out_cmd, $do_cmd, $jagid);

    system "cp $in_file ./ >& _junk";
    $in_file = basename($in_file);
    $do_cmd = 0;
    system "mkdir -p outs";
    print "Running Jaguar...";
    $counter = 0;
    $insert_loc = 0;
    if (open READIN, $in_file) {
	while (<READIN>) {
	    $infile_text[$counter] = $_;
	    if ($_ =~/gcharge/ and $insert_loc > -1) {
		$insert_loc = $counter;
	    }elsif ($_ =~/isolv=2/) {
		$insert_loc = -1;
	    }
	    $counter++;
	}
	close READIN;

	if (! $insert_loc) {
	    print "Error in Jaguar inputfile $in_file...";
	}elsif ($insert_loc > -1) {
	    $infile_text[$insert_loc] .= "isolv=2\n";
	    system "mkdir -p newins";
	    $in_file = "./newins/" . basename($in_file);
	    if (open WRITEOUT, "> $in_file") {
		for (@infile_text) {
		    print WRITEOUT $_;
		}
		close WRITEOUT;
		$do_cmd = 1;
	    } else {
		print "Error while writing to $in_file: $!\n";
	    }
	    close WRITEOUT;
	} else {
	    $do_cmd = 1;
	}	    
    }

    if (! $do_cmd) {
	print "Error reading $in_file\n";
	return $return_str;
    } else {
	$out_cmd = "jaguar run $in_file |";
	$return_str = basename($in_file);
	$return_str =~ s/\.in//;
	if (open INFILE, $out_cmd) {
	    while (<INFILE>) {
		if ($_ =~/^JobId:\s+(\S+)/) {
		    $jagid = $1;
		}
	    }
	    close INFILE;
	    $do_cmd = 0;
	    sleep 5;
	    do {
		$out_cmd = "tail -50 $return_str" . ".out >& jagstatus";
#		print "$out_cmd";
		system ($out_cmd);
		open STATUSFILE, "jagstatus";
		while (<STATUSFILE>) {
		    if ($_ =~ /Job $return_str completed/) {
			$return_str .= ".out";
			$do_cmd = 1;
			last;
		    } elsif ($_ =~ /Jaguar cannot recover from this error and will now abort/) {
			$do_cmd = 2;
			last;
		    }
		}
		close STATUSFILE;
	    } until ($do_cmd > 0);

	    if ($do_cmd == 1) {
		system "mv $return_str ./outs/";
		$return_str = "./outs/" . $return_str;
		return $return_str;
	    } else {
		print "Error reading $return_str\n";
		return "";
	    }
	} else {
	    print "Error in Jaguar Job...\n";
	    return "";
	}
    }
}
	    

sub Numerically {
    ($a <=> $b);
}

sub GetMolNames(@) {
    my ($in_file) = $_[0];
    my ($counter);

    for $counter (keys %MOL_INFO) {
	$MOL_INFO{$counter}{"NAME"} = $counter;
	$MOL_INFO{$counter}{"EXP"} = 0.00;
    }

    if ($in_file and -e $in_file and -T $in_file) {
	if (open INFILE, $in_file) {
	    while (<INFILE>) {
		if ($_ =~ /^(\w+)\s+(\S+)\s+(\-?\d+\.\d+)/) {
		    if ($MOL_INFO{$1}) {
			$MOL_INFO{$1}{"NAME"} = $2;
			$MOL_INFO{$1}{"EXP"} = $3;
		    }
		}
	    }
	}
    }
}

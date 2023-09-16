#!/usr/bin/perl -w

BEGIN
{
    push @INC, "/ul/tpascal/scripts/";
}

use strict;
use Packages::General;
use File::Basename;
use Packages::GetParms;

sub ValidateFields();
sub GetTrajFileList();
sub GetHBondData(@);
sub Numerically;
sub GetTemp(@);
sub GetAvg(@);
sub CreateParmFile(@);
sub ParseInfoFile(@);
sub WriteOutput();

die "Usage: $0 parm_file calc_hbond tempreature_spacing trajectory_extention directories\n"
    if (! @ARGV || $#ARGV < 4);

my ($parm_file, $calc_hbond, $temp_increment, $trj_ext, @working_dirs) = @ARGV;
my (%HBonds, @Files, $curr_file, $out_data, $P_File);

print "--=== START PROGRAM ===--\n";
ValidateFields();
print "Getting list of trajectories...";
GetTrajFileList();
print "Done\n";

system "mkdir -p datfiles";
print "Calculating H-Bonds...";
for $curr_file (@Files) {
    open INFOFILE, ">> dumpfile.txt" || die "Cannot open dumpfile.txt: $!\n";
    print INFOFILE "\nUsing $curr_file\n";
    close INFOFILE;
    GetHBondData($curr_file->{"TRAJFILE"}, $curr_file->{"OUTFILE"});
    print ".";
}
print "Done\n";

print "Writing Composite HBond Data...";
WriteOutput();
print "Done\n--===END PROGRAM===--\n";
    
sub ValidateFields() {
    my ($counter, $valid_dirs);

    die "Expected boolean value for calc_hbond\n"
    if (! IsBoolean($calc_hbond));

    die "Expected integer or decimal for tempreature\n"
	if (! IsDecimal($temp_increment) && ! IsInteger($temp_increment));
    die "Cannot locate parameter file $parm_file: $!\n"
	if (! -e $parm_file);

    $P_File = Packages::GetParms->new();
    if (! $P_File->IsValidParams($parm_file)) {
        die "Error in Paramater file\n";
    }

    $P_File->{"Files"}->{"trajectory"} = "";

    for $counter (0 .. $#working_dirs) {
	if (! -e $working_dirs[$counter]) {
	    print "ERROR: " . $working_dirs[$counter] . " does not exist...Ignoring\n";
	    $working_dirs[$counter] = "";
	} else {
	    $valid_dirs = 1;
	}
    }

    die "Cannot find any valid directories\n"
	if (! $valid_dirs);
}

sub GetTrajFileList() {
    my ($curr_dir, @trj_list, $curr_trj, $curr_outFile, $rec);
    my ($ls_cmd);

    for $curr_dir(@working_dirs) {
	@trj_list = ();
	if (-d $curr_dir) {
	    $ls_cmd = "ls $curr_dir/*." . $trj_ext;
	    open OUTFILE, "$ls_cmd |" || die "Cannot execute ls command: $!\n";
	    while (<OUTFILE>) {
		chomp;
		push @trj_list, $_;
	    }
	    close OUTFILE;
	    if ($#trj_list > -1) {
		for $curr_trj (@trj_list) {
		    $curr_outFile = $curr_trj;
		    $curr_outFile =~ s/$trj_ext/out/;
		    if (-e $curr_outFile) {
			$rec = (
				{
				    "TRAJFILE" => $curr_trj,
				    "OUTFILE"  => $curr_outFile,
				}
				);
			push @Files, $rec;
		    }
		}
	    }
	}
    }

    die "Cannot locate any valid trajectories and output file\n"
	if ($#Files == -1);
}

sub GetHBondData(@) {
    my ($trj_file, $info_file)  = @_;
    my (%Temps, $info_step, $trj_step, $valid_info);
    my ($sys_cmd, $curr_temp, $valid_file, $rec);

    my $store_file = "./datfiles/" . basename($trj_file);
    $store_file =~  s/\.$trj_ext/_hbond.dat/;
     
# First Get list of tempreatures

    $valid_info = 0;
    open INFOFILE, $info_file or die "Cannot open $info_file: $!\n";
    while (<INFOFILE>) {
	chomp;
	if ($_ =~ /\s+ntpr\s*\=\s*(\d+)/) {
	    $info_step = $1;
	}
	if ($_ =~ /\s+ntwx\s*\=\s*(\d+)/) {
	    $trj_step = $1;
	}
	if ($_ =~ /R M S  F L U C T U A T I O N S/) {
	    last;
	}elsif ($_ =~ /NSTEP\s+\=\s*(\d+)\s+(.+)/) {
	    $rec = ParseInfoFile($2);
	    if ($rec) {
		$Temps{$1} = $rec;
		$valid_info = 1;
	    }
	}
    }
    close INFOFILE;

    die "ERROR: Couldn't find information about the info and/or trajectory",
    " time step in $info_file\n"
	if (! $info_step || ! $trj_step);
    die "ERROR: Unable too extract tempreature information from $info_file\n"
	if (! $valid_info);

# Create Parm File
    CreateParmFile($trj_file);

# Run HBond Anal
    $sys_cmd = "/ul/tpascal/scripts/hbond_anal.pl ./parm >> dumpfile.txt";

    if ($calc_hbond) {
	$sys_cmd = system($sys_cmd);
	if ($sys_cmd) {
	    die "Error executing hbond_anal.pl\n";
	}else {
	    system "mv hbond.dat $store_file";
	}
    }
#Get HBond Info
    $valid_file = 0;
    if (open INFILE, $store_file) { 
	while (<INFILE>) {
	    chomp;
	    if ($_ =~ /^\s+(\d+)\s+(\d+\.\d+)/) {
		$valid_file = 1;
		$curr_temp = GetTemp($1 * $trj_step, \%Temps);
		if ( $HBonds{$curr_temp} ) {
		    $HBonds{$curr_temp} = GetAvg($HBonds{$curr_temp} . " " . $2);
		}else {
		    $HBonds{$curr_temp} = $2;
		}		
	    }
	}
	close INFILE;
    }else {
	print  "Cannot open h-bond data file $store_file: $!\n";
    }

    print "ERROR: No valid data found in hbond.dat for $trj_file\n"
	if (! $valid_file);
}

sub ParseInfoFile(@) {
    my ($input_string) = $_[0];
    $input_string =~ s/\=//g;
    my (%rec, @valid_data, $counter, $index);
    my (@data_array) = split /\s+/, $input_string;

    @valid_data = ("TIME", "TEMP", "PRESS");

    if ($#data_array > -1) {
	for $counter (0 ..$#data_array) {
	    for $index (@valid_data) {
		if ($data_array[$counter] =~ /$index/) {
		    if ($data_array[$counter + 1] && $data_array[$counter + 1] =~ /(\d+\.\d+)/) {
			$rec{$index} = $1;
		    }
		}
	    }
	}
    }

    return \%rec;
}

sub CreateParmFile(@) {
    my ($trj_file) = $_[0];
    my ($outText);

    $outText = "Mol: " . $P_File->{"Molecule"}->{"major_groove"};
    $outText .= ":" . $P_File->{"Molecule"}->{"minor_groove"};
    $outText .= "\nTotal bases: " . $P_File->{"Molecule"}->{"total_bases"};
    $outText .= "\nBases at end: " . $P_File->{"Molecule"}->{"bases_at_end"};
    $outText .= "\nCrossovers: ";
    for (@{ $P_File->{"Molecule"}->{"crossovers"} }) {
	$outText .= "$_ ";
    }
    $outText .= "\nTopology file: " . $P_File->{"Files"}->{"topology"};
    $outText .="\nTrajectory file: " . $trj_file;
    $outText .="\n3 prime in: " . $P_File->{"Molecule"}->{"is3PrimeIn"} ;

    open OUTFILE, "> parm" || die "Cannot create parm: $!\n";
    print OUTFILE "$outText\n";
    close OUTFILE;

}

sub GetTemp(@) {
    my ($curr_temp, $Temp_Data) = @_;
    my ($prev_temp, $next_temp) = (0, 0);
    my ($returnval, $counter);

    if ($Temp_Data->{$curr_temp}) {
	$returnval =  $Temp_Data->{$curr_temp}->{"TEMP"};
    }else {
	for $counter (sort Numerically keys %{ $Temp_Data }) {
	    if (! $prev_temp) {
		$prev_temp = $Temp_Data->{$counter}->{"TEMP"};
	    }

	    if ($counter < $curr_temp) {
		$prev_temp = $Temp_Data->{$counter}->{"TEMP"};
	    }elsif ($counter > $curr_temp) {
		$next_temp = $Temp_Data->{$counter}->{"TEMP"};
	    }

	    if ($prev_temp && $next_temp) {
		$returnval = GetAvg($prev_temp . " " . $next_temp);
		last;
	    }
	}
    }

    if (! $returnval) {
	$returnval = $prev_temp;
    }

    return $returnval;
}		

sub GetAvg(@) {
    my (@vals) = split /\s+/, $_[0];
    my ($returnval, $counter, $sum);

    $counter = 0;
    for (@vals) {
	$sum += $_;
	$counter++;
    }


    $returnval = $sum/$counter;

    return $returnval;
}

sub Numerically {
    ($a<=>$b);
}

sub WriteOutput() {
    my ($out_data, $curr_temp, $avg_hbonds, $avg_temp);

    $curr_temp = 0;

    open OUTFILE, "> HB_Data.dat" or die "Cannot create HB_Data.dat: $!\n";
    for $out_data (sort Numerically keys %HBonds) {
	if ($curr_temp == 0) {
	    $curr_temp = $out_data;
	    $avg_hbonds = "";
	    $avg_temp = "";
	}

	if (($out_data - $curr_temp) > $temp_increment) {
	    chop($avg_hbonds);
	    chop($avg_temp);
	    $avg_hbonds = GetAvg($avg_hbonds);
	    $avg_temp = GetAvg($avg_temp);
	    printf OUTFILE "%10.3f%10.2f\n", $avg_temp, $avg_hbonds;
	    $curr_temp = 0;
	}else {
	    $avg_hbonds .= ($HBonds{$out_data} . " ");
	    $avg_temp .= ($out_data . " ");
	}
    }

    $avg_hbonds = GetAvg($avg_hbonds);
    $avg_temp = GetAvg($avg_temp);
    printf OUTFILE "%10.3f%10.2f\n", $avg_temp, $avg_hbonds;
    close OUTFILE;
}

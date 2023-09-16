#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}

use strict;
use Packages::General;
use File::Basename;

sub ValidateFields();
sub GetTempList();
sub ReadMDInfo();
sub FixInputFiles(@);
sub CreateMDCmd(@);
sub RunMD(@);
sub IsDecimal(@);


if (! @ARGV || $#ARGV < 4) {
    die "usage: $0 restart_file topology start_temp end_temp increment [mpi_cmd]\n";
}

my ($rst_file, $topo, $tempi, $tempo, $inc, @p_cmd) = @ARGV;
my ($t_counter, $new_temp, $old_temp, $curr_rst, $in_file, $my_cmd, @TData, $do_parallel);
ValidateFields();
GetTempList();

system "mkdir -p heating";
system "mkdir -p equilibration";
system "cp $rst_file ./";
system "cp $topo ./";
$old_temp = $tempi;
$curr_rst = "../" . basename($rst_file);
$topo = "../" . basename($topo);

for (@p_cmd) {
    $do_parallel .= " " . $_;
}


print "-== START DNA HEATING PROTOCOL ==-\n\n";
print "Using parallel cmd: $do_parallel\n"
    if ($do_parallel);

for $t_counter (@TData) {
    $new_temp = $t_counter;
    chdir "./heating";
    $in_file = FixInputFiles($old_temp, $new_temp, 1);
    ($curr_rst, $my_cmd) = 
	CreateMDCmd($old_temp, $new_temp, $in_file, 1, $curr_rst);
    $curr_rst = "../heating/" . $curr_rst;
    print "HEATING: $old_temp K -> $new_temp K...";
    RunMD($my_cmd);
    print "Done\n";
    $old_temp = $new_temp;
    $new_temp = ReadMDInfo();

    chdir "../equilibration";
    $in_file = FixInputFiles($old_temp, $new_temp, 0);
    ($curr_rst, $my_cmd) = 
	CreateMDCmd($old_temp, $new_temp, $in_file, 0, $curr_rst);
    $curr_rst = "../equilibration/" . $curr_rst;
    print "EQUILIBRATING AT $new_temp K...";
    RunMD($my_cmd);
    print "Done\n";
    $old_temp = ReadMDInfo();
    chdir "../";
}
print "\n-== END ==-\n";
    
sub ValidateFields() {
    die "Invalid starting tempreature: Expected deicimal, got $tempi\n"
	if (! IsDecimal($tempi));

    die "Invalid ending tempreature: Expected deicimal, got $tempo\n"
	if (! IsDecimal($tempo));

    die "Invalid tempreature interval: Expected deicimal, got $inc\n"
	if (! IsDecimal($inc));

    die "Increment must be greater than 1K\n"
	if (int($inc) < 1);

    die "Cannot locate starting restart file $rst_file: $!\n"
	if (! -e $rst_file);

    die "Cannot locate topology file $topo: $!\n"
	if (! -e $topo);

}

sub GetTempList() {
    my ($curr_time);
    
    if ($tempi < $tempo) {
	$curr_time = $tempi;
	while ($curr_time <= $tempo) {
	    $curr_time += $inc;
	    push @TData, $curr_time;
	}
    }else {
	$curr_time = $tempo;
	while ($curr_time >= $tempi) {
	    $curr_time -= $inc;
	    push @TData, $curr_time;
	}
    }	
}

sub ReadMDInfo() {
    my ($inText, $result);
    
    open MDINFO, "mdinfo" || die "Cannot open file mdinfo: $!\n";
    while (<MDINFO>) {
	chomp;
	$inText = $_;
	if ($inText =~ /^ NSTEP.+TEMP\(K\) =\s*(\d+\.\d+)/) {
	    $result = $1;
	    last;
	}
    }
    close MDINFO;

    return $result;
}

sub FixInputFiles(@) {
    my ($tmp_i, $tmp_o, $isHeat) = @_;
    my ($inFile, $inText, $outText, $saveNm);

    $inFile = "/ul/tpascal/nano/PX-molecules/AMBER/Melting/procedure/";

    if ($isHeat) {
	$inFile .= "md_heating.in";
	$saveNm = "md_heating_" . $tmp_i . "_" . $tmp_o . ".in";
    }else {
	$inFile .= "md_nvt.in";
	$saveNm = "md_heating_" . $tmp_o . ".in";
    }

    system "cp $inFile $saveNm";

    open INFILE, $saveNm || die "Cannot open $saveNm: $!\n";
    while (<INFILE>) {
	chomp;
	$inText = $_;
	$inText =~ s/temp0  = /temp0  = $tmp_o/g;
	$inText =~ s/tempi  = /tempi  = $tmp_i/g;
	$inText =~ s/value1=/value1= $tmp_i/g;
	$inText =~ s/value2=/value2= $tmp_o/g;
	$outText .= $inText . "\n";
    }
    close INFILE;

    open OUTFILE, "> $saveNm" || die "Cannot modify $saveNm: $!\n";
    print OUTFILE $outText;
    close OUTFILE;

    return $saveNm;
}

sub CreateMDCmd(@) {
    my ($tmp_i, $tmp_o, $inFile, $isMelt, $coords) = @_;
    my ($returnval, $nmPrefix, $new_rst, $sander_cmd);

    if ($isMelt) {
	$nmPrefix = "md_melt_" . $tmp_i . "_" . $tmp_o;
    }else {
	$nmPrefix = "md_nvt_" . $tmp_o;
    }

    $new_rst = "$nmPrefix" . ".restart";

    if (! $do_parallel) {
	$sander_cmd = "sander";
    } else {
	$sander_cmd = "/ul/maiti/amber7/exe/sander";
    }

    $returnval = "";
    $returnval .= "$sander_cmd -O -i $inFile ";
    $returnval .= "-o $nmPrefix" . ".out ";
    $returnval .= "-p $topo ";
    $returnval .= "-c $coords ";
    $returnval .= "-ref $coords ";
    $returnval .= "-x $nmPrefix" . ".trj ";
    $returnval .= "-r $new_rst";

    return ($new_rst, $returnval);
}

sub RunMD(@) {
    my ($md_cmd) = $_[0];

    $md_cmd = $do_parallel . " " . $md_cmd;

    if (system($md_cmd)) {
	die "An Error has occured\n";
    }
}

sub IsDecimal (@) {
    my ($inString) = $_[0];
    
    $inString =~ /^\d+\.\d+$/ ? 
	return 1 :
	return 0;
}


#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use IPC::Open2;
use Packages::General qw(FileTester Trim IsInteger TrjSelections STDev);
use Packages::AMBER qw(getTopInfo getOpts);
use Getopt::Std qw(getopt);

sub init;
sub numerically { ($a<=>$b); }
sub getSASA;
sub saveData;
sub createVMDCmds;
sub execCmd;
sub parseResults;

my ($totAtms, $SELECT, $SASA);
my ($topFile, $trjFile, $atomSelection, $saveName);

$|++;
&init;
print "Calculating SASA for \"$atomSelection\"...";
$SASA = getSASA($atomSelection, $SELECT);
print "Done\nCreating $saveName data file..";
saveData($SASA, $saveName);
print "Done\n";

sub init {
    my ($AMBEROPTS, $i, %VALID, $trjSelection, $DATA, %OPTS, $isValid);
    $isValid = 0;

    getopt('tcasf', \%OPTS);
    ($topFile, $trjFile, $atomSelection, $trjSelection, $saveName) = ($OPTS{t},$OPTS{c},$OPTS{a},$OPTS{s},$OPTS{f});
    for $i ($topFile, $trjFile, $atomSelection, $trjSelection) {
	die "usage: $0 -t amber topology file -c amber coordinate file -a solute selection\n" . 
	    "-s \"trajectory selection\" -f [saveName]\n" if (! defined($i));
    }

    print "Initializing...";
    FileTester($topFile);
    FileTester($trjFile);

    $SELECT = TrjSelections($trjSelection) if ($trjSelection ne "*");
    die "ERROR: Trajectory selection cannot be \"*\"!\n" if (! keys %{ $SELECT });
    $saveName = "residuals.dat" if (! $saveName);
    
    if ($atomSelection =~ /wat|nucleic|protein/i) {
	$isValid = 1;
    } else {
	$AMBEROPTS = &getOpts;
	print "Parsing AMBER topology file $topFile...";
	($DATA, $totAtms) = getTopInfo($topFile, $AMBEROPTS);
	for $i (keys %{ $DATA->{ATOMS} }) {
	    if ($DATA->{ATOMS}{$i}{ATMNAME} eq $atomSelection) {
		$atomSelection = "name " . $atomSelection;
		$isValid = 1;
		last;
	    } elsif ($DATA->{ATOMS}{$i}{FFTYPE} eq $atomSelection) {
		$atomSelection = "type " . $atomSelection;
		$isValid = 1;
		last;
	    }
	}
    }

    die "ERROR: $atomSelection is not a valid atom name/type!\n" if (! $isValid);
    print "Done\n";
}

sub getSASA {
    my ($atomSelect, $trjSelection) = @_;
    my ($vmdCmds, $i, $vmdCmd, $DATA);

    $vmdCmd = "/exec/VMD/bin/vmd -dispdev none -nt -e vmd_cmds.cmd >& data.dump";
    
    $vmdCmds = createVMDCmds($atomSelect, $topFile, $trjFile, $trjSelection);

    open OUTDATA, "> vmd_cmds.cmd" or die "ERROR: Cannot read from data.dump: $!\n";
    for $i (@{ $vmdCmds }) {
	print OUTDATA "$i\n";
    }
    close OUTDATA;

    if (system($vmdCmd)) {
	die "Error while executing $vmdCmd\n";
    }

    $DATA = parseResults("data.dump");
    system("rm -fr vmd_cmds.cmd data.dump");

    return $DATA;
}

sub createVMDCmds {
    my ($types, $topName, $trjName, $trjSelection) = @_;
    my ($atm, @CMDS, $i, @tmp, $sasaCommand);

    $i = 1;
    @tmp = sort numerically keys %{ $trjSelection };

    $CMDS[0] = "mol new $topName type parm7 first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all";
    $CMDS[1] = "mol addfile $trjName type crdbox first $tmp[0] last $tmp[$#tmp] step 1 filebonds 1 autobonds 1 waitfor all";
    $CMDS[2] = "set atomSelect [atomselect top \"$types\"]";

    $sasaCommand = "set atomSelect [atomselect top \"$types\" frame";
    
    for $i (@tmp) {
	push @CMDS, "$sasaCommand $i]";
	push @CMDS, "measure sasa 1.4 \$atomSelect";
    }

    push @CMDS, "quit";

    return \@CMDS;
}
    
sub execCmd {
    my ($cmdStr, $reader, $writer) = @_;
    my ($result);

    print "Command $cmdStr\n";
    print $writer "$cmdStr\n";
    chomp($result = readline<$reader>);
    while ($result && $result =~ /Info/) {
	chomp($result = <$reader>);
    }
    die "Error detected: $result\n"
	if ($result =~ /wrong|error/i);

    return $result;
}

sub parseResults {
    my ($dataFile) = $_[0];
    my ($frame, $data, %DATA, $dataStr);

    $frame = 0;
    open INDATA, $dataFile or die "ERROR: Cannot open file $dataFile:$!\n";
    while (<INDATA>) {
	chomp;
	$data = $_;
	if ($data =~ /^atomselect(\d+)/) {
	    $frame = $1;
	} elsif ($frame > 0 && $data =~ /^(\d+\.\d+)/) {
	    $DATA{VALS}{$frame} = $1;
	    $dataStr .= "$1 ";
	}
    }
    close INDATA;

    die "ERROR: No valid data found!\n" if (! keys %DATA);
    chop $dataStr;
    ($DATA{STATS}{AVG}, $DATA{STATS}{STDEV}, $DATA{STATS}{TOT}) = STDev($dataStr);
    return \%DATA;
}

sub saveData {
    my ($DATA, $save) = @_;
    my ($i);

    open OUTDATA, "> $save" or die "ERROR: Cannot create data file $save: $!\n";
    for $i (sort numerically keys %{ $DATA->{VALS} } ) {
	printf OUTDATA "%8d %8.5f\n", $i, $DATA->{VALS}->{$i};
    }
    printf OUTDATA "%-8s %8.5f\n%-8s %8.5f\n","#AVG", $DATA->{STATS}{AVG}, "#STDdev", $DATA->{STATS}{STDEV};
    close OUTDATA;
}

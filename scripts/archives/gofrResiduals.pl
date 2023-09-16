#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}


use strict;
use IPC::Open2;
use Packages::General qw(FileTester Trim IsInteger TrjSelections);
use Packages::AMBER qw(getTopInfo getOpts);

sub init;
sub numerically;
sub getResiduals;
sub saveData;
sub createVMDCmds;
sub execCmd;
sub parseResults;
sub calcResiduals;

my ($totAtms, $TYPES, $SELECT, $RESIDUALS);

die "usage: $0 amberTopFile amberTrjFile atomType1 atomType2 \"trajectory selection\" [saveName]\n"
    if (! @ARGV || $#ARGV < 4);

my ($topFile, $trjFile, $type1, $type2, $selection, $saveName) = @ARGV;

$|++;
print "Initializing...";
$TYPES = &init;
print "Done\nCalculating $type1 <-> $type2 radial distribution residuals...";
$RESIDUALS = getResiduals($TYPES, $SELECT);
print "Done\nCreating $saveName data file..";
saveData($RESIDUALS, $saveName);
print "Done\n";

sub init {
    my ($OPTS, $i, %VALID, $atm, $DATA);
    FileTester($topFile);
    FileTester($trjFile);

    $SELECT = TrjSelections($selection);    
    $saveName = "residuals.dat" if (! $saveName);
    $OPTS = &getOpts;
    print "Parsing AMBER topology file $topFile...";
    ($DATA, $totAtms) = getTopInfo($topFile, $OPTS);
    
    $type1 = Trim($type1);
    $type2 = Trim($type2);

    for $atm ($type1, $type2) {
	$VALID{$atm}{VAL} = 0;
	for $i (keys %{ $DATA->{ATOMS} }) {
	    if ($DATA->{ATOMS}{$i}{ATMNAME} eq $atm) {
		$VALID{$atm}{VAL} = 1;
		$VALID{$atm}{FIELD} = "name";
		last;
	    } elsif ($DATA->{ATOMS}{$i}{FFTYPE} eq $atm) {
		$VALID{$atm}{VAL} = 1;
		$VALID{$atm}{FIELD} = "type";
		last;
	    }
	}
	if ($atm =~ /wat|nucleic|protein/i) {
	    $VALID{$atm}{VAL} = 1;
	    $VALID{$atm}{FIELD} = "";
	}
    }

    for $atm (keys %VALID) {
	die "ERROR: $atm is not a valid atom name/type!\n"
	    if (! $VALID{$atm}{VAL});
    }
    $DATA = ();
    return \%VALID;
}

sub getResiduals {
    my ($atomTypes, $trjSelection) = @_;
    my ($vmdCmds, $i, $vmdCmd, $DATA, %RESDATA);

    $vmdCmd = "/exec/VMD/bin/vmd -dispdev none -nt -e vmd_cmds.cmd >& data.dump";
    
    $vmdCmds = createVMDCmds($atomTypes, $topFile, $trjFile, $trjSelection);

    open OUTDATA, "> vmd_cmds.cmd" or die "ERROR: Cannot read from data.dump: $!\n";
    for $i (@{ $vmdCmds }) {
	print OUTDATA "$i\n";
    }
    close OUTDATA;

    if (system($vmdCmd)) {
	die "Error while executing $vmdCmd\n";
    }

    $DATA = &parseResults("data.dump");

    for $i (1 .. $#{ $DATA }) {
	$RESDATA{$i} = calcResiduals(\%{ $DATA->[($i - 1)] }, \%{ $DATA->[$i] });
    }

    $DATA = ();
    system("rm -fr vmd_cmds.cmd data.dump");
    return \%RESDATA;
}

sub createVMDCmds {
    my ($types, $topName, $trjName, $trjSelection) = @_;
    my ($atm, @CMDS, $i, @tmp, $gofrCommand);

    $i = 1;
    @tmp = sort numerically keys %{ $trjSelection };

    $CMDS[0] = "mol new $topName type parm7 first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all";
    $CMDS[1] = "mol addfile $trjName type crdbox first $tmp[0] last $tmp[$#tmp] step 1 filebonds 1 autobonds 1 waitfor all";

    for $atm (keys %{ $types }) {
	push @CMDS, "set t$i [atomselect top \"$types->{$atm}{FIELD} $atm\"]";
	$i++;
    }

    $gofrCommand = "measure gofr \$t1 \$t2 rmax 35 usepbc 0 selupdate 0";
    
    for $i (0 .. ($#tmp - 1)) {
	push @CMDS, "$gofrCommand first " . ($tmp[$i] - 1) . " last " . ($tmp[($i + 1)] - 1);
    }

    push @CMDS, "quit";

    return \@CMDS;
}

sub numerically {
    ($a<=>$b);
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
    my (@tmp, %DATA, $i, @RET, $data, $counter);

    $counter = 0;
    open INDATA, $dataFile or die "ERROR: Cannot open file $dataFile:$!\n";
    while (<INDATA>) {
	chomp;
	$data = $_;
	if ($data =~ /\{/) {
	    $data =~ s/\{//g;
	    @tmp = split /\}/, $data;
	    die "ERROR: Invalid response from VMD: $data\n" if ($#tmp < 1);
	    %DATA = ();
	    for $i (0 .. $#tmp) {
		$tmp[$i] = Trim($tmp[$i]);
		@{ $DATA{$i} } = split /\s+/, $tmp[$i];
	    }
	    for $i (0 .. $#{ $DATA{0} }) {
		$RET[$counter]{$DATA{0}[$i]} = $DATA{1}[$i];
	    }
	    $counter++;
	}
    }
    close INDATA;

    return \@RET;
}

sub calcResiduals {
    my ($data1, $data2) = @_;
    my ($i, $sum);
    
    $sum = 0;
    for $i (keys %{ $data1 }) {
	$sum += ($data1->{$i} - $data2->{$i})**2;
    }

    $sum = sqrt($sum);
    return $sum;
}

sub saveData {
    my ($DATA, $save) = @_;
    my ($i);

    open OUTDATA, "> $save" or die "ERROR: Cannot create data file $save: $!\n";
    for $i (sort numerically keys %{ $DATA } ) {
	printf OUTDATA "%8d%8.5f\n", $i, $DATA->{$i};
    }
    close OUTDATA;
}

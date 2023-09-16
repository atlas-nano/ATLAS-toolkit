#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo);
use LAMMPS qw(ParseLAMMPSTrj ParseLAMMPSLogFile GetLammpsByteOffset);
use General qw(FileTester TrjSelections GetStats);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub saveLogData;
sub storeEng;
sub getTotAbsCharge;
sub numerically { ($a<=>$b); }

my ($bgfFile, $trjFile, $logFile, $saveFile);
my ($lSELECT, $tSELECT, $BGF, $BONDS, $pStr, $LOGDATA, $ENGDATA, $totAbsCharge);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($BGF, $BONDS) = GetBGFFileInfo($bgfFile, 0);
$totAbsCharge = getTotAbsCharge($BGF);
print "Done\nParsing LAMMPS log file $logFile...";
&ParseLAMMPSLogFile($logFile, $lSELECT, \&saveLogData, undef);
print "Done\n";
&GetLammpsByteOffset($tSELECT, $trjFile, scalar keys %{ $BGF });
$pStr = "Parsing LAMMPS trajectory $trjFile...";
ParseLAMMPSTrj($BGF, $trjFile, $tSELECT, "energy", \&storeEng, $pStr, undef);
print "Computing atom statistics...";
&calcStats($ENGDATA);
print "Done\nSaving data to $saveFile...";
&writeData($ENGDATA, $saveFile);
print "Done\n";

sub writeData {
   my ($data, $savename) = @_;
   my ($i);

   open OUTDATA, "> $savename" or die "ERROR: Cannot create $savename: $!\n";
   printf OUTDATA "%-8s%12s%12s%12s%12s%12s%12s\n",qw(Atom VDW Coul HB HB_dre Bond Total);
   #printf OUTDATA "%-8s%20s%20s%20s%20s%20s%20s\n",qw(Atom VDW Coul HB HB_dre Bond Total) if ($writeStdev);
   for $i (sort numerically keys %{ $data }) {
        printf OUTDATA "%-8d", $i;
        for (qw(vdw coul hb hb_dre bond)) {
            printf OUTDATA "%12.3f", 0;
        }
        printf OUTDATA "%12.3f\n", $data->{$i}{STATS}{AVG};
    }
    close OUTDATA;
}

sub calcStats {
    my ($data) = $_[0];
    my ($i);
    
    for $i (keys %{ $data }) {
	$data->{$i}{STATS} = GetStats($data->{$i}{DATA});
    }
}

sub storeEng {
    my ($DATA, $bgfInfo, $fileHandle) = @_;
    my ($totEng, $logEng, $tStep, $i, $pmeContrib, $atomList, $count);

    # distribute the total missing long range electrostatics based on the magnitude of the charge
    $tStep = $DATA->{TIMESTEP}[0];
    die "ERROR: Timestep $tStep not found in log file!\n" if (! exists($LOGDATA->{$tStep}));
    $logEng = $LOGDATA->{$tStep}{poteng};
    @{ $atomList } = keys %{ $DATA->{ATOMS} };
    for $i (@{ $atomList }) {
	push @{ $ENGDATA->{$i}{DATA} }, $DATA->{ATOMS}{$i}{ENERGY};
	$totEng += $DATA->{ATOMS}{$i}{ENERGY};	
    }
    $pmeContrib = ($logEng - $totEng)/$totAbsCharge;
    $totEng = 0;
    $count = $#{ $ENGDATA->{ $atomList->[0] }{DATA} };
    for $i (@{ $atomList }) {
	$ENGDATA->{$i}{DATA}[$count] += abs($bgfInfo->{$i}{CHARGE}) * $pmeContrib;
	$totEng += $ENGDATA->{$i}{DATA}[$count];
    }
    print "";
}

sub saveLogData {
    my ($data, $i, $FLEPTR) = @_;
    my ($field, $offset);

    for $field (keys %{ $data }) {
	$LOGDATA->{$i}{$field} = $data->{$field};
    }
}

sub getTotAbsCharge {
    my ($bgf) = $_[0];
    my ($i, $tot);

    for $i (keys %{ $bgf }) {
	$tot += abs($bgf->{$i}{CHARGE});
    }
    
    return $tot;
}

sub init {
    my (%OPTS, $trjSelect);
    getopt('btlso',\%OPTS);
    for ("b", "t", "l") {
	die "usage: $0 -b bgf file -t lammps trj file -l lammps log file -s (trj selection) -o (save file)\n"
	    if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($bgfFile, $trjFile, $logFile, $saveFile, $trjSelect) = ($OPTS{b}, $OPTS{t}, $OPTS{l}, $OPTS{o}, $OPTS{s});
    for ($bgfFile, $trjFile, $logFile) {
	FileTester($_);
    }
    $trjSelect = "*" if (! defined($trjSelect));
    if (! defined($saveFile)) {
	$saveFile = basename($trjFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_eng.dat";
    }
    $tSELECT = TrjSelections($trjSelect); 
    $lSELECT = TrjSelectiions($trjSelect) if ($trjSelect ne "*");
    print "Done\n";
}
	

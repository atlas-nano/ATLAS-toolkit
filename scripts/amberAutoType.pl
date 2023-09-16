#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester);

sub init;
sub typeAtoms;

my ($mol2File, $saveName, $prefix, $charge, $atomType);
my ($rFF);

$|++;
&init;
print "Typing atoms...";
$rFF = typeAtoms($mol2File, $saveName, $prefix, $charge, $atomType);
print "Done\nPlease use the following forcefields with these files: \"$rFF\"\n";

sub typeAtoms {
	my ($mol2File, $saveName, $prefix, $netcharge, $at) = @_;
	my ($antechamberCmd, $convertCmd, $parmchkCmd, $combineCmd); 
	my ($testCmd, $fs, $refFF, $fsize);

    $antechamberCmd  = "$ENV{AMBERHOME}/bin/antechamber -dr no -pf yes -c bcc -fi mol2 -nc $netcharge " .
                        "-fo mol2 -at $at -i $mol2File -o __tmp.mol2 > /dev/null";
    die "ERROR: Cannot obtain charges/fftypes for molecules:\n$antechamberCmd\n"
        if(system($antechamberCmd));

	$refFF = "$ENV{AMBERHOME}/dat/leap/parm/parm19.dat";
	$refFF = "$ENV{AMBERHOME}/dat/leap/parm/gaff2.dat" if($at ne "amber");
	$parmchkCmd = "$ENV{AMBERHOME}/bin/parmchk2 -i __tmp.mol2 -o __frcmod " . 
                    " -f mol2 -p $refFF > /dev/null";
    system($parmchkCmd);
	$testCmd = "cp __frcmod ${prefix}.frcmod";
	system($testCmd);
	$testCmd = "cp __tmp.mol2 ${prefix}.mol2";
	system($testCmd);

    $convertCmd = "$Bin/mol2bgf.pl -m __tmp.mol2 -b $saveName > /dev/null";
    die "ERROR: Cannot create $saveName"
       if (system($convertCmd));

    system("rm -fr sqm.* __tmp.mol2 __frcmod __frcmod.ff gaff.ff");
	#now test for whether the frcmod file is empty
	$fsize = -s "${prefix}.frcmod";
	if($fsize < 70) {
		$testCmd = "rm -fr ${prefix}.frcmod";
		system($testCmd);
	} else {
		$refFF .= " ${prefix}.frcmod";
	}
	return $refFF;
}

sub init {
	my (%OPTS);

	getopt('msct',\%OPTS);
	die "usage: $0 -m mol2 file -s (save_name) -c (net charge) -t (atom type=amber[default]|gaff)\n"
	    if (! exists($OPTS{m}));
	print "Initializing...";
	($mol2File, $saveName, $charge, $atomType) = ($OPTS{m}, $OPTS{s}, $OPTS{c}, $OPTS{t});
	FileTester($mol2File);
	$prefix = basename($mol2File);
    $prefix = basename($saveName) if (defined($saveName));
    $prefix =~ s/\.\w+$//;
    $ENV{AMBERHOME} = "$Bin/../codes/amber18/bin" if (!exists($ENV{AMBERHOME}));
    $charge = 0 if (! defined($charge) or $charge !~ /^\-?\d+\.?\d*/);
	$atomType = "gaff" if(defined($atomType) and $atomType =~ /gaff/i);
	$atomType = "amber" if(! defined($atomType) or $atomType =~ /amber/i);
	$atomType = lc $atomType;
	if (! defined($saveName)) {
	    $saveName = basename($mol2File);
	    $saveName =~ s/\.\w+$//;
	    $saveName .= ".${atomType}.bgf";
	}
	print "Done\n";
}

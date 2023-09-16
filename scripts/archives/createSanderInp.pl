#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts");
}

use strict;
use Packages::General;

# This program will execute a tleap and prepare the structure for simulation
# with sander
#
# usage: createSanderInp pdbfile saveName [is_glucose] [leaprc location] [solvate?]

sub Initialize;
sub CreateLeapInput;
sub RunLeap;

my ($pdbfile, $saveName, $is_glucose, $doSolvation, $leapLoc) = @ARGV;

die "usage: $0 pdbfile [saveName] [is_glucose] [solvate?=no] [leaprc location]\n"
    if (! @ARGV);

Initialize($pdbfile, $doSolvation, $leapLoc, $saveName);
CreateLeapInput($pdbfile, $doSolvation, $leapLoc, $saveName);
RunLeap;


sub Initialize {
    my ($inFile, $doSolv, $myLeap, $myName) = @_;
    my ($myCmd);

    die "Error accessing regular file $inFile: $!\n"
	if (! -e $inFile or ! -r $inFile or ! -T $inFile);
    $doSolv = 0
	if (! defined($doSolv) or ! IsInteger($doSolv));

    if (! defined($myLeap)) {
	die "ERROR: leaprc location is not specified and cannot access AMBERHOME enviromental variable\n"
	    if (! $ENV{AMBERHOME});
	$myLeap = $ENV{AMBERHOME} . "/dat/leap/cmd/leaprc.ff03";
    }
    
    if (! defined($is_glucose) or ! IsInteger($is_glucose) or $is_glucose != 1) {
	$is_glucose = 0;
    }

    die "Error accessing regular file leaprc: $!\n"
	if (! -e $myLeap or ! -r $myLeap or ! -T $myLeap);
    $myCmd = "cp $myLeap ./leaprc";
    die "ERROR while copying $myLeap: $!\n"
	if (system($myCmd));
    
    $leapLoc = "./leaprc";

    $myName = $inFile
	if (! defined($myName));

    if ($myName =~ /(.+)\.\w{3}$/) {
	$saveName = $1;
    } else {
	$saveName = $myName;
    }
}

sub CreateLeapInput {
    my ($myPDB, $doSolv, $myLeap, $modelName) = @_;

    open MYLEAPRC, ">> $myLeap" or die "Cannot write to $myLeap: $!\n";
    
    if ($is_glucose) {
	print MYLEAPRC "loadoff ../glucose.off\n";
    }

    print MYLEAPRC "model = loadmol2 $myPDB\n";
    
    if ($doSolv) {
	print MYLEAPRC "solvatebox model TIP3PBOX 10\n";
	print MYLEAPRC "addions model Na+ 0\n";
    } else {
        print MYLEAPRC "setBox model vdw\n";
    }
    print MYLEAPRC "saveamberparm model $modelName" . ".prmtop $modelName" . ".inpcrd\n";

    print MYLEAPRC "savepdb model $modelName" . "_solv.pdb\n"
	if ($doSolv);

    print MYLEAPRC "quit\n";
    
    close MYLEAPRC;

}

sub RunLeap {
    
    print "Preparing $pdbfile for sander minimization...";
    system "tleap > junk";

    print "Done\n";

    system "rm -f junk leaprc leap.log";
    print "Created: $saveName" . ".top, $saveName" . ".crd\n";
}

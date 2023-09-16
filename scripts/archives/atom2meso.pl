#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts");
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF);
use Packages::General qw(FileTester Trim);
use File::Basename;
use Packages::MESO;

#   atomistic2Meso.pl: This script will open an atomistic bgf file
#   and create the corresponding meso-scale representation

sub main;
sub init;
sub numerically;
sub CreateMesoModel;
sub getStrand;
sub getBeadInfo;
sub isMember;
sub createHBondBead;
sub updateCOM;
sub MakeBGF;
sub buildBondList;
sub getSearchList;
sub getAtmMass;

die "usage: $0 bgf_file|directory parm_file [save_name|save_dir] [strand length]\n"
    if (! @ARGV or $#ARGV < 1);

my ($fileLoc, $parmFile, $saveName, $sLen) = @ARGV;

$|++;
&main;

sub main {
    my ($FILES, $fC, $ATOMS, $BONDS, $i, $PARMS);
    my ($MESO, $bgfFile, $BGF, $CONS, $CHAIN, $HEADER, $outStr, @tmp);

    $FILES = init($parmFile, $fileLoc, $saveName);
    print "Obtaining Bead parameters ....";
    $PARMS = GetMesoParms($parmFile);
    print "Done\n";
    
    $fC = 0;
    
    for $i(@{ $FILES }) {
        print "$i->{FILE}:....READING...";
	$outStr = "$i->{FILE}:....READING...";
	($ATOMS, $BONDS) = GetBGFFileInfo($i->{FILE}, 0);
        @tmp = keys %{ $ATOMS };
	print "" . ($#tmp + 1) . " atoms";
        print "....Creating MESO";
	$outStr .= "....Creating MESO";
        $MESO = CreateMesoModel($ATOMS, $PARMS);
        $bgfFile = $i->{SAVE};
        print "....Writing Meso";
	$outStr .= "....Writing Meso";
        ($BGF, $CONS) = MakeMesoBGF($MESO, $PARMS, $ATOMS, $BONDS);
	@tmp = keys %{ $BGF };
	print " " . ($#tmp + 1) . " meso atoms...";
	$HEADER = createHeaders(my $BOX, $bgfFile);
	addHeader($BGF, $HEADER);
	createBGF($BGF, $CONS, $bgfFile);
	print "..Done\r";
	$outStr .= "..Done";
	$outStr = length($outStr) + 30;
	printf "%${outStr}s\r", " ";
	$BGF = ();
	$CONS = ();
	$HEADER = ();
	$MESO = ();
    }
    
    print "All Tasks Completed\n";

}

sub init {
    my ($pLoc, $fLoc, $save) = @_;
    my (@FILES, @tmp, $dir, $rec);

    FileTester($pLoc);

    if (! -T $fLoc) {
	if (-d $fLoc) {
	    opendir FILEDIR, $fLoc or die "ERROR: Cannot access directory $fLoc: $!\n";
	    @tmp = grep { /\.bgf$/ && -f} map { "$fLoc/$_"} readdir FILEDIR;
	    closedir FILEDIR or die "ERROR: Cannot close directory $fLoc: $!\n";
	}
    } else {
	$tmp[0] = $fLoc;
    }

    die "ERROR: No valid pdbfiles found!\n"
	if ($#tmp == -1);

    if (! defined($save)) {
	$save = $tmp[0];
	if ($#tmp == 0) {
	    $save =~ s/\.\w+/_meso\.bgf/;
	} else {
	    $save = "meso/";
	    system "mkdir -p meso";
	}
    }

    if ($#tmp == 0) {
	$rec = (
		{
		    "FILE" => $tmp[0],
		    "SAVE" => $save,
		}
		);
	push @FILES, $rec;
    } else {
	$dir = $save . "/";
	for (@tmp) {
	    $save = basename($_);
	    $save =~ s/\.\w+/_meso\.bgf/;
	    $save = $dir . $save;
	    $rec = (
		    {
			"FILE" => $_,
			"SAVE" => $save,
		    }
		    );
	    push @FILES, $rec;
	}
    }
	    
    return \@FILES;
}


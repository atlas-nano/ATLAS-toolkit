#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts");
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo);
use Packages::General qw(FileTester Trim GetBondLength GetAngle GetTorsion CenterText STDev);
use File::Basename qw(basename dirname);
use Packages::MESO qw(GetMesoParms CreateMesoModel MakeMesoBGF);
use Packages::Math::Matrix;
use Packages::BOX qw(CreateGrid GetNeighbours GetBox);
use Cwd qw(getcwd);
use Getopt::Std qw(getopt);

#   Get_Meso.pl: This script will open a bgf file(s) and will obtain the
#   Mesoscale valence and hydrogen bond parameters for the appropriate 
#   "beads" defined in the parameter file. It will then print the statistics

sub main;
sub init;
sub Numerically { ($a<=>$b); }
sub updateStats;
sub computeStats;
sub getKeyList;
sub PrintStats;
sub PrintHist;
sub getHBondStats;
sub validHB;
sub ParseBGFFile;
sub getBeadNames;
sub showUsage;

my ($bgf, $parmFile, $save_name);

$|++;
&main;

sub main {
    my ($BGF_FILES, $fC, $start, $end, $BONDS);
    my ($i, $ATOMS, $sLen, $MESO, $bgfFile, $save);
    my ($CONS, %POS_VECS, $BEADS, $PARMS, $rec);
    my (%BONDS, %ANGLES, %TORSIONS, %HBONDS, $BNAMES);
    
    ($BGF_FILES, $save) = &init;
    print "Obtaining Bead parameters ....";
    $PARMS = GetMesoParms($parmFile);
    print "Done\n\n";
    
    $fC = 0;
    $BNAMES = getBeadNames($PARMS);
  
    $start = time();
    print "START TIMER: Parsing Structures ", scalar localtime($start), "\n";
    
    for $i (@{ $BGF_FILES }) {
        print "$i:...READING...";
        ($sLen, $ATOMS, $BONDS) = ParseBGFFile($i);
        next
            if (! $sLen);
	print "Creating MESO...";
        $MESO = CreateMesoModel($ATOMS, $PARMS);
	($MESO, $CONS) = MakeMesoBGF($MESO, $PARMS, $ATOMS, $BONDS);
	print "stats...";
	$rec = (
		{
		    "ATOMS" => $MESO,
		    "BONDS" => $CONS,
		    "SLEN"  => $sLen,
		}
		);
	updateStats($rec, $PARMS, \%BONDS, \%ANGLES, \%TORSIONS, \%HBONDS);
        print "Finished\n";
        $fC++;
	$rec = ();
	$ATOMS = ();
	$CONS = ();
    }
    
    $end = time();
    print "\nEND TIMER:   ", scalar localtime($end), "\n";
    print "ELAPSED: " . ($end - $start) . " seconds\n";
    
    print "\nSTATISTICS\n==========\n";
    $start = time();
    print "START TIMER: Computing Stats ", scalar localtime($start), "\n";
    computeStats(\%BONDS, \%ANGLES, \%TORSIONS, \%HBONDS, $save);
    $end = time();
    print "END TIMER:   ", scalar localtime($end), "\n";
    print "ELAPSED: " . ($end - $start) . " seconds\n";
    
    print "\nAll Tasks Completed\n";

}

sub getBeadNames {
    my ($parms) = $_[0];
    my ($i, %BNAMES);

    for $i (keys %{ $parms->{"BEADS"} }) {
	$BNAMES{ Trim($parms->{"BEADS"}{$i}{"NAME"}) }{"ID"} = $i;
        $BNAMES{ Trim($parms->{"BEADS"}{$i}{"NAME"}) }{"RADII"} = $parms->{"BEADS"}{$i}{"RADII"};
    }

    return \%BNAMES;
}

sub init {
    my (@BGFS, @tmp, $i, $mesoName, %OPTS, $findCmd);

    getopt('bps',\%OPTS);
    ($parmFile, $bgf, $save_name) = ($OPTS{p},$OPTS{b},$OPTS{s});

    for ($parmFile, $bgf) {
	die "usage: $0 -b bgf_file|directory -p parm_file -s [save_name|directory]\n"
	    if (! defined($_));
    }

    print "Initializing...";

    FileTester($parmFile);

    if (-e $bgf and -d $bgf) {
        opendir BGFFILES, $bgf or die "ERROR: Cannot access directory $bgf: $!\n";
        @tmp = grep { /\.bgf$/ && -f} map { "$bgf/$_"} readdir BGFFILES;
        closedir BGFFILES or die "ERROR: Cannot close directory $bgf: $!\n";
    } elsif (-e $bgf) {
	$tmp[0] = $bgf;
    } else {
	$findCmd = "find $bgf -name '*.bgf' -print";
	if (open(FINDCMD, "$findCmd |")) {
	    while (<FINDCMD>) {
		chomp;
		push @tmp, $_;
	    }
	    close FINDCMD;
	}
    }
    die "ERROR: No valid BGF files found!\n"
	if ($#tmp == -1);

    for $i (0 .. $#tmp) {
	$BGFS[$i] = $tmp[$i];
    }

    if (! $save_name) {
        $save_name = "./";
    }

    $save_name = dirname($save_name);
    print "Done\n";
    return (\@BGFS, $save_name);
}

sub updateStats {
    my ($cFile, $parms, $BONDS, $ANGLES, $TORSIONS, $HBONDS) = @_;
    my ($a1, $a2, $a3, $a4, $i, $j, $k);
    my ($l, $sKey, @atmList, $sName);
    
    for $i (keys %{ $cFile->{"ATOMS"} }) {
	$a1 = $cFile->{"ATOMS"}{$i};

	for $j (@{ $cFile->{"BONDS"}{$i} }) { #BONDS
	    $a2 = $cFile->{"ATOMS"}{$j};
	    @atmList = ($a1, $a2);
	    ($sKey, $sName) = getKeyList(\@atmList, $parms);
	    $BONDS->{$sKey}{"VALS"} .= GetBondLength($a1, $a2) . " ";
	    $BONDS->{$sKey}{"NAME"} = $sName;
		
	    for $k (@{ $cFile->{"BONDS"}{$j} }) { #ANGLES
		next if ($k == $i);
		$a3 = $cFile->{"ATOMS"}{$k};
		@atmList = ($a1, $a2, $a3);
		($sKey, $sName) = getKeyList(\@atmList, $parms);
		$ANGLES->{$sKey}{"VALS"} .= GetAngle($a1, $a2, $a3, 1) .  " ";
		$ANGLES->{$sKey}{"NAME"} = $sName;

		for $l (@{ $cFile->{"BONDS"}{$k} }) { #TORSIONS
		    next if ($l == $i || $l == $j);
		    $a4 = $cFile->{"ATOMS"}{$l};
		    @atmList = ($a1, $a2, $a3, $a4);
		    ($sKey, $sName) = getKeyList(\@atmList, $parms);
		    $TORSIONS->{$sKey}{"VALS"} .= GetTorsion($a1, $a2, $a3, $a4, 1) . " ";
		    $TORSIONS->{$sKey}{"NAME"} = $sName;
                }
            }
        }
    }

    getHBondStats($cFile->{"ATOMS"}, $HBONDS);
}

sub computeStats {
    my ($BONDS, $ANGLES, $TORSIONS, $HBONDS, $save) = @_;
    my ($i);

    printf "\n\n%-20s%20s%20s%20s\n", "COMBO", "AVERAGE", "Std Dev", "TOTAL PTS";

    for $i (0 .. 79) {
	print "-";
    }
    print "\n";

    PrintStats($BONDS, "BONDS", $save);
    PrintStats($ANGLES, "ANGLES", $save);
    PrintStats($TORSIONS, "TORSIONS", $save);
    PrintStats($HBONDS, "HBONDS", $save);
}

sub getKeyList {
    my ($atmList, $parms) = @_;
    my (%indexList, $keyStr, $nameStr, $i); 
    my ($key1, $key2, @tmp, $sameBond);
    
    $key1 = $key2 = $keyStr = $nameStr = "";
    for $i (@{ $atmList }) {
	$key1 .= sprintf("%03d",$i->{"ID"});
    }

    for $i (reverse @{ $atmList }) {
	$key2 .= sprintf("%03d",$i->{"ID"});
    }

    if ($key1 > $key2) {
	@tmp = @{ $atmList };
    } else {
	@tmp = reverse @{ $atmList };
    }
    
    for $i (@tmp) {
	$keyStr .= sprintf("%03d",$i->{"ID"});
	$nameStr .= Trim($i->{"ATMNAME"}) . "-";
    }
    chop $nameStr;

    $sameBond = substr($keyStr, -6);
    if (exists($parms->{"SAME"}{"BOND"}{$sameBond})) {
	$sameBond = $parms->{"SAME"}{"BOND"}{$sameBond};
    	substr($keyStr, -6) = $sameBond;
    }

    return ($keyStr, $nameStr);
}


    
sub PrintStats(@) {
    my ($in_obj, $header, $save) = @_;
    my ($avg, $stdev, $total, $values, $counter, $cDir, $outFile, $i);

    $cDir = getcwd;
    $header = $save . "/" . $header;
    $outFile = $save . "/" . lc($header) . "_stats.dat";
    open OUTDATA, "> $outFile" or die "ERROR: Cannot create file $outFile: $!\n";
    printf OUTDATA "%-20s%20s%20s%20s\n", "COMBO", "AVERAGE", "Std Dev", "TOTAL PTS";

    for $i (0 .. 79) {
	print OUTDATA "-";
    }
    print OUTDATA "\n";

    system "mkdir -p $header";
    chdir $header or die "Cannot access directory $header: $!\n";
    printf "\n%80s\n", CenterText($header, 80);
    for $counter (sort Numerically keys %{ $in_obj } ) {
	$values = $in_obj->{$counter}{"VALS"};
	chop $values;
	
	PrintHist($values, $in_obj->{$counter}{"NAME"});
	($avg, $stdev, $total) = STDev($values);
        printf "%-20s%20.10f%20.10f%20d\n", 
	$in_obj->{$counter}{"NAME"}, $avg, $stdev, ($total/$avg);
        printf OUTDATA "%-20s%20.10f%20.10f%20d\n", 
	$in_obj->{$counter}{"NAME"}, $avg, $stdev, ($total/$avg);
	
    }

    chdir $cDir or die "ERROR: Cannot change to current directory $cDir:$!\n";
    close OUTDATA or die "ERROR: Cannot finish creation of file $outFile\n";
}

sub PrintHist(@) {
    my ($vals, $mname) = @_;
    my (@holder, $counter, $file_name);
    
    $mname =~ s/\s+/_/g;

    $file_name = $mname . ".dat";
    @holder = split /\s+/, $vals;
#    @holder = sort @holder;
    
    open HIST, ">> $file_name" or die "Cannot write to Histogram file $file_name : $!\n";
    print HIST "# $mname\n";
    for $counter (0 .. $#holder) {
	printf HIST "%10d%20.10f\n", $counter, $holder[$counter];
    }

    close HIST;
}
	
sub getHBondStats {
    my ($ATOMS, $HBOND) = @_;
    my ($GRID, $junk, $BOX, $cellC, @hList, $hb, $parent);
    my ($i, $j, $k, $cell, $aC, $atom, $CELLS, $sType);
    my ($hb_valid, $hb_dist, $hb_angle, $saveKey, $saveName);
    
    $BOX = GetBox($ATOMS, undef, undef);
    ($GRID, $junk) = CreateGrid($ATOMS, 0, $BOX, 2.5, 0); # 2.5 A grid
    
    for $i (keys %{ $GRID } ) {
        for $j (keys %{ $GRID->{$i} }) {
            for $k (keys %{ $GRID->{$i}{$j} }) {
                $cell = $GRID->{$i}{$j}{$k};
                @hList = ();
                for $atom (@{ $cell->{"ATOMS"} }) {
                    if (exists($atom->{"PARENT"}) and $atom->{"FFTYPE"} =~ /D_\d+/) { # take only donors
                        $parent = $ATOMS->{ $atom->{"PARENTID"} };
                        push @hList, (
                                      {
                                        "ATOM"   => $atom,
                                        "PARENT" => $parent,
                                      }
                                     );
                    }
                }
                next if (! @hList);
                $CELLS = GetNeighbours($GRID, $cell);
                for $hb (@hList) {
                    $hb->{"ATOM"}{"FFTYPE"} =~ /(\d+)/;
                    $sType = "H_" . $1;
                    for $cell (@{ $CELLS }) {
                        for $atom (@{ $cell->{"ATOMS"} }) {
                            if ($atom->{"FFTYPE"} eq $sType) {
                                ($hb_valid, $hb_dist, $hb_angle) = validHB($hb->{"PARENT"}, $hb->{"ATOM"}, $atom);
                                if ($hb_valid) {
                                    $saveName = $hb->{"ATOM"}{"FFTYPE"} . "-" . $sType . "_dist";
				    $saveKey = sprintf("%03d%03d", $hb->{"ATOM"}{"ID"}, $atom->{"ID"}); 
                                    $HBOND->{$saveKey}{"VALS"} .= "$hb_dist ";
                                    $HBOND->{$saveKey}{"NAME"} = $saveName;
                                    $saveName = $hb->{"ATOM"}{"FFTYPE"} . "-" . $sType . "_angl";
				    $saveKey += 1000000;
                                    $HBOND->{$saveKey}{"VALS"} .= "$hb_angle ";
                                    $HBOND->{$saveKey}{"NAME"} = $saveName;                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    @hList = ();
    $CELLS = ();
    $GRID = ();
}

sub validHB {
    my ($parent, $donor, $acceptor) = @_;
    my ($valid, $dist, $angle);    
    
    $valid = 0;
    
    $dist = GetBondLength($parent, $acceptor);
    $angle = GetAngle($parent, $donor, $acceptor, 1);
    
    if ($dist < 6 and $angle > 1.047 and $angle < 4.19) {
        $valid = 1;
    }
    
    return ($valid, $dist, $angle);
    
}


sub ParseBGFFile(@) {
    my ($bgf_file) = $_[0];
    my ($Atom_Info, $BONDS) = GetBGFFileInfo($bgf_file, 0);
    my ($counter, $curr_res_name, $curr_res_id, $strand_indicator, $strand_length, @tmp);

    $curr_res_id = $strand_length = 0;
    $strand_indicator = "";
    @tmp = sort Numerically keys %{ $Atom_Info };
    for $counter (@tmp) {
	if ($Atom_Info->{$counter}{"RESNAME"} =~ /WAT|Na\+/i) {
	    delete $Atom_Info->{$counter};
	    next;
	}
	if ($Atom_Info->{$counter}{"RESNUM"} > $curr_res_id) {
	    $curr_res_id = $Atom_Info->{$counter}{"RESNUM"};
	    $curr_res_name = $Atom_Info->{$counter}{"RESNAME"};
	    if ($curr_res_name =~ /(\d+)/) {
		if ($strand_indicator eq "") {
		    $strand_indicator = $1;
		} elsif ($strand_indicator ne $1 and ! $strand_length) {
		    $strand_length = $curr_res_id;
		}
	    }
	}
    }

    if (! $strand_length) {
	print "ERROR processing $bgf_file: CANNOT find 5' or 3' end\n";
    }
    return ($strand_length, $Atom_Info, $BONDS);

}

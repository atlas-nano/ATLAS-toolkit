#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo);
use Packages::MESO;
use Packages::Math::Matrix;
use Packages::General qw(FileTester Trim STDev);

#   getAtomVectors:  This script will open a bgf file(s), create a mesoscale 
#   model and compute the atomic displacement vectors wrt the meso atoms. This 
#   is used to reconstruct the atomistic representation from the mesocale 
#   representation

sub main;
sub init;
sub Numerically;
sub getVectors;
sub triangulate;
sub SolveMatrix(@);
sub computeAtomVectors;
sub getBead;
sub ParseBGFFile(@);

die "usage: $0 bgf_file|directory parm_file [save_name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($bgf, $parmFile, $save_name) = @ARGV;

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
    $start = time();
    print "START TIMER: Parsing Structures ", scalar localtime($start), "\n";
    
    for $i (@{ $BGF_FILES }) {
        print "$i:....READING";
        ($sLen, $ATOMS, $BONDS) = ParseBGFFile($i);
        next
            if (! $sLen);
	print "....Creating MESO";
        $MESO = CreateMesoModel($ATOMS, $PARMS);
	($MESO, $CONS) = MakeMesoBGF($MESO, $PARMS, $ATOMS, $BONDS);
	print "..atom vectors";
	getVectors($ATOMS, $MESO, $CONS, \%POS_VECS, $fC);
        print "....Finished\n";
        $fC++;
    }
    
    $end = time();
    print "\nEND TIMER:   ", scalar localtime($end), "\n";
    print "ELAPSED: " . ($end - $start) . " seconds\n";
    
    print "\nSTATISTICS\n==========\n";
    $start = time();
    computeAtomVectors(\%POS_VECS, $PARMS, $save);
    $end = time();
    print "END TIMER:   ", scalar localtime($end), "\n";
    print "ELAPSED: " . ($end - $start) . " seconds\n";
    
    print "\nAll Tasks Completed\n";

}

sub init {
    my (@BGFS, @tmp, $i, $mesoName);

    FileTester($parmFile);

    if (-e $bgf and -d $bgf) {
        opendir BGFFILES, $bgf or die "ERROR: Cannot access directory $bgf: $!\n";
        @tmp = grep { /\.bgf$/ && -f} map { "$bgf/$_"} readdir BGFFILES;
        closedir BGFFILES or die "ERROR: Cannot close directory $bgf: $!\n";
    } elsif (-e $bgf) {
	$tmp[0] = $bgf;
    }

    die "ERROR: No valid BGF files found!\n"
	if ($#tmp == -1);

    for $i (0 .. $#tmp) {
	$BGFS[$i] = $tmp[$i];
    }

    if (! $save_name) {
        $save_name = $BGFS[0];
	$save_name =~ s/\.\w+$/_atomVecs.dat/;
    }

    return (\@BGFS, $save_name);
}

sub Numerically {
    ($a<=>$b);
}

sub getVectors {
    my ($ATOMS, $MESO, $BONDS, $VECS, $file_counter) = @_;
    my ($b1, $bC, $b2, $b3, $aC, $atom, $coord, $iC, $beadName);
    
    for $bC (keys %{ $MESO }) {
        next if exists($MESO->{$bC}{"PARENT"});
        ($b1, $b2, $b3) = triangulate($MESO, $BONDS, $bC);
        next if (! $b1 || ! $b2 || ! $b3);
        $beadName = $b1->{"ID"};
        for $aC (keys %{ $b1->{"ATOMS"} }) {
            $atom = $ATOMS->{$aC};
            $iC = SolveMatrix($atom, $b1, $b2, $b3);
            for $coord (keys %{ $iC }) {
                $VECS->{$beadName}{$b2->{"ID"}}{$b3->{"ID"}}{$atom->{"ATMNAME"}}{$coord} .= 
		    $iC->{$coord} . " ";
            }
        }
    }
}

sub triangulate {
    my ($ATOMS, $CONS, $beadC) = @_;
    my ($b1, $b2, $b3, $i);
    
    $b1 = $ATOMS->{$beadC};
    
    if (@{ $CONS->{$beadC} }) {
        $b2 = $ATOMS->{ $CONS->{$beadC}[0] };
        if ($#{ $CONS->{$beadC} } > 0) {
            $b3 = $ATOMS->{ $CONS->{$beadC}[1] };
        } else {
            $beadC = $b2->{"INDEX"};
            for $i (@{ $CONS->{$beadC} }) {
                if ($i != $b1->{"INDEX"}) {
                    $b3 = $ATOMS->{ $i };
                    last;
                }
            }
        }
    }
    
    return ($b1, $b2, $b3);
}

sub SolveMatrix(@) {
# Creates a linear system Ax = b and solves it    
    my ($curr_atom, $bead_1, $bead_2, $bead_3) = @_;
    my ($matrixA, $vectorB, $vectorX, $temp, %RESULT);

    $matrixA = new Math::Matrix ([$bead_1->{"XCOORD"}, $bead_2->{"XCOORD"}, $bead_3->{"XCOORD"}],
				[$bead_1->{"YCOORD"}, $bead_2->{"YCOORD"}, $bead_3->{"YCOORD"}],
				[$bead_1->{"ZCOORD"}, $bead_2->{"ZCOORD"}, $bead_3->{"ZCOORD"}]);
    $vectorB = new Math::Matrix ([$curr_atom->{"XCOORD"}, $curr_atom->{"YCOORD"}, $curr_atom->{"ZCOORD"}]);
    $temp =  $matrixA->concat($vectorB->transpose);
    $vectorX = $temp->solve;
				 
    $RESULT{"XCOORD"} = $vectorX->[0][0];
    $RESULT{"YCOORD"} = $vectorX->[1][0];
    $RESULT{"ZCOORD"} = $vectorX->[2][0];

    return \%RESULT;
}

sub computeAtomVectors {
    my ($arrays, $parms, $save) = @_;
    my ($outString, $b1, $b2, $b3, $dim);
    my ($avg, $stdev, $total, $atom, $aC);
    
    print "\nATOM VECTORS\n==============\n\n";
    open OUTFILE, "> $save" or die "Cannot write to file $save: $!\n";
    $outString = sprintf("%-10s%10s%10s%10s%15s%10s%15s%10s%15s%10s\n", "NAME", 
    "BEAD", "BEAD 1", "BEAD 2", "P1", "+/-", "P2", "+/-", "P3", "+/-");
    print "$outString\n";
    print OUTFILE "$outString\n";
    
    for $b1 (keys %{ $arrays } ) {
        for $b2 (keys %{ $arrays->{$b1} }) {
            for $b3 (keys %{ $arrays->{$b1}{$b2} }) {
                for $aC (keys %{ $arrays->{$b1}{$b2}{$b3} }) {
                    $atom = $arrays->{$b1}{$b2}{$b3}{$aC};
                    printf "%-10s%10s%10s%10s", Trim($aC), getBead($parms, $b1, "NAME"),
                    getBead($parms, $b2, "NAME"), getBead($parms, $b3, "NAME");
                    printf OUTFILE "%-10s%10s%10s%10s", Trim($aC), getBead($parms, $b1, "ELEMENT"),
                    getBead($parms, $b2, "ELEMENT"), getBead($parms, $b3, "ELEMENT");
                    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
                        chop $atom->{$dim};
                        ($avg, $stdev, $total) = STDev($atom->{$dim});
                        printf "%15.5f%10.5f", $avg, $stdev;
                        printf OUTFILE "%15.5f%10.5f", $avg, $stdev;
                    }
                    print "\n";
                    print OUTFILE "\n";
                }
            }
        }
    }
    
    close OUTFILE;

}

sub getBead {
    my ($parms, $beadID, $hKey) = @_;
    my ($returnStr);
    
    $returnStr = $beadID;
    
    if (exists($parms->{"BEADS"}{$beadID})) {
        if (exists($parms->{"BEADS"}{$beadID}{$hKey})) {
            $returnStr = $parms->{"BEADS"}{$beadID}{$hKey};
        }
    } elsif ($beadID < 100) { #donor
        $returnStr = "D_" . ($beadID - 50);
    } else {
        $returnStr = "H_" . ($beadID - 100);
    }
    
    return Trim($returnStr);
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

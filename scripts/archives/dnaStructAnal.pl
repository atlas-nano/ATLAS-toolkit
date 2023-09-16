#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General qw(FileTester Trim TrjSelections);
use Packages::AMBER qw(getTopInfo ParseAmberTrj getOpts GetAmberByteOffset);
use Packages::GetParms;
use Packages::HelixLayout;
use Packages::FileFormats qw(sortByRes);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType);
use Getopt::Std qw(getopt);
use Packages::ManipAtoms qw(UnwrapAtoms);

sub init;
sub determineHelixLayout;
sub setLayout;
sub numerically;
sub writeHelixData;
sub savePDB;
sub getOrder;
sub fixPDB;
sub getHelicalParms;
sub showHelp;
sub moveResults;
sub convertBox;
sub makeBox;

my ($parmFile, $selection, $filebase, $field, $getByteOffset);
my ($OPTS, $SELECT, $DATA, $totAtms, $Parm, $LAYOUT, $trjType); 
my ($printStr, $RESDATA, $HC, $topFile, $trjFile, $getSnapshot);
my ($LAMMPSOPTS);

$|++;
&init;
print "Parsing AMBER topology file $topFile...";
($DATA, $totAtms) = getTopInfo($topFile, $OPTS);
$RESDATA = sortByRes($DATA->{ATOMS});
print "Done\n";
&GetLammpsTrjType($SELECT, $trjFile, "atom", \%{ $LAMMPSOPTS }) 
    if ($Parm->{Files}{trjType} == 2);
if ($Parm->{CreatePDB}) {
    if ($Parm->{Files}{trjType} == 2) {
	$field = "atom";
    } else {
	$field = $totAtms;
    }
    $printStr = "Creating PDB from $trjType trajectory $trjFile...";
    $getByteOffset->($SELECT, $trjFile, scalar keys %{ $DATA->{ATOMS} });
    $getSnapshot->($DATA->{ATOMS}, $trjFile, $SELECT, $field, \&writeHelixData, $printStr, undef);
}
print "Obtaining helical parameters...";
&getHelicalParms;
&moveResults;
print "Done\n";

sub makeBox {
    my ($box) = $_[0];
    my (%BOX, $i, @dim);

    @dim = ("X", "Y", "Z");
    for $i (0 .. 2) {
        $BOX{$dim[$i]}{lo} = 0;
        $BOX{$dim[$i]}{hi} = $box->{$i + 2}{DATA};
        $BOX{$dim[$i] . "COORD"}{lo} = 0;
        $BOX{$dim[$i] . "COORD"}{hi} = $box->{$i + 2}{DATA};
        $BOX{$dim[$i] . "COORD"}{len} = $box->{$i + 2}{DATA};
    }

    return \%BOX;
}

sub convertBox {
    my ($lammpsBox) = $_[0];
    my (%amberBox, $i);

    $amberBox{1}{DATA} = 90;
    for $i (0 .. 2) {
        $amberBox{$i + 2}{DATA} = $lammpsBox->[$i]{hi} - $lammpsBox->[$i]{lo};
    }

    return \%amberBox;
}

sub moveResults {
    system("mkdir -p results");
    my ($file, $i, $j);

    for $i (1 .. 2) {
	for $j ("Torsions", "StrandInfo", "InterStrandInfo", "Backbone") {
	    $file = "helix${i}/${filebase}_Helix${i}_lis_${j}.dat";
	    if (-e $file && -r $file) {
		system("mv $file results");
	    }
	}
    }
}

sub savePDB {
    my ($ATOMS, $helixC, $frame) = @_;
    my ($atm, $saveName, $fmt, $i);

    $fmt = "%-5s%6d%5s%4s%2s%4s%12.3f%8.3f%8.3f\n";
    $saveName = "pdbfiles/" . $filebase . "_Helix" . 
	($helixC + 1) . "_" . $HC->{$helixC}{COUNT} . ".pdb";
    open PDBFILE, "> $saveName" || die "ERROR: Cannot create PDB file $saveName: $!\n";
    for $i (sort numerically keys %{ $ATOMS }) {
	printf PDBFILE $fmt, "ATOM", $i, $ATOMS->{$i}{ATMNAME}, $ATOMS->{$i}{RESNAME},
	uc(chr(64 + $ATOMS->{$i}{CHAIN})), $ATOMS->{$i}{RESNUM}, $ATOMS->{$i}{XCOORD},
	$ATOMS->{$i}{YCOORD}, $ATOMS->{$i}{ZCOORD};
    }
    close PDBFILE;
    $HC->{$helixC}{COUNT}++;
    return $saveName;
}

sub writeHelixData {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($i, $j, $helixC, $chainC, $regionC, $atm, $pdbFile, $chain, $BBOX);  
    my ($start, $end, $DATA, $field, $resNum, $tmp, $lastAtm, $newChain);

    if ($Parm->{Files}{trjType} == 2) { # LAMMPS
	for $i (keys %{ $ATOMS->{ATOMS} }) {
	    for $j (keys %{ $BOX->{$i} }) {
		$ATOMS->{ATOMS}{$i}{$j} = $BOX->{$i}{$j};
	    }
	}
	$BOX = convertBox($ATOMS->{"BOX BOUNDS"});
	$frameNum = $ATOMS->{FRAME};
	$ATOMS = $ATOMS->{ATOMS};
        $BBOX = makeBox($BOX);
        UnwrapAtoms($ATOMS, $BBOX, $LAMMPSOPTS->{scaled});
    }
    for $helixC (0 .. $#{ $LAYOUT }) {
	$i = $resNum = 1;
	for $chainC (0 .. $#{ $LAYOUT->[$helixC] }) {
	    $newChain = 1;
	    for $regionC (0 .. $#{ $LAYOUT->[$helixC][$chainC] }) {
		$start = $LAYOUT->[$helixC][$chainC][$regionC]{StartUnit};
		$end = $LAYOUT->[$helixC][$chainC][$regionC]{EndUnit};
		$chain = $LAYOUT->[$helixC][$chainC][$regionC]{Strand};
		$tmp = getOrder($start, $end);
		for $j (@{ $tmp }) {
		    for $atm (sort numerically keys %{ $RESDATA->{$j}{ATOMS} }) {
			if ($newChain && $ATOMS->{$atm}{ATMNAME} =~ /O\dP/) {
				next;
			}
			for $field (keys %{ $ATOMS->{$atm} }) {
			    $DATA->{$i}{$field} = $ATOMS->{$atm}{$field};
			}
			$DATA->{$i}{CHAIN} = $chain;
			$DATA->{$i}{RESNUM} = $resNum;
			if ($newChain && $DATA->{$i}{ATMNAME} eq "P") {
			    $DATA->{$i}{ATMNAME} = "H5T";
			    for $field ("XCOORD", "YCOORD", "ZCOORD") {
				$DATA->{$i}{$field} -= 1;
			    }
			}
			$i++;
			$lastAtm = $atm;
		    }
		    $resNum++;
		    $newChain = 0;
		}
	    }
	    if ($ATOMS->{$lastAtm}{ATMNAME} ne "H3T") {
		for $field (keys %{ $ATOMS->{$lastAtm} }) {
		    $DATA->{$i}{$field} = $ATOMS->{$lastAtm}{$field};
		}
		$DATA->{$i}{ATMNAME} = "H3T";
		$DATA->{$i}{CHAIN} = $chain;
		$resNum--;
		$DATA->{$i}{RESNUM} = $resNum;
		$DATA->{$i}{XCOORD} += 1;
		$resNum++;
		$i++;
	    }
	}

	$pdbFile = savePDB($DATA, $helixC, $frameNum);
	fixPDB($pdbFile);
	$DATA = ();
    }

}

sub init {
    my ($i, $cmd, %OPTS);

    getopt('pft',\%OPTS);
    ($parmFile, $filebase, $selection) = ($OPTS{p},$OPTS{f},$OPTS{t});

    for ($parmFile, $selection) {
	&showHelp if (! defined($_));
    }

    print "Initializing...";

    FileTester($parmFile);

    $OPTS = &getOpts;

    $Parm = Packages::GetParms->new();
    die "Error in Paramater file\n" if (! $Parm->IsValidParams($parmFile));

    $topFile = $Parm->{Files}{topology};
    $trjFile = $Parm->{Files}{trajectory};
    FileTester($topFile);
    FileTester($trjFile);

    if (! $Parm->{isLayout} ) {
	$LAYOUT = determineHelixLayout($Parm);
    } else {
	$LAYOUT = setLayout($Parm);
    }

    system "mkdir -p pdbfiles";
    for $i (0 .. $#{ $LAYOUT }) {
	$cmd = "mkdir -p helix" . ($i + 1);
	system "$cmd\n";
	$HC->{$i}{COUNT} = 1;
    }

    $filebase = $Parm->{Molecule}{name} if (! defined($filebase));
    $SELECT = TrjSelections($selection);
    if ($Parm->{Files}{trjType} == 1) {
        $getSnapshot = \&ParseAmberTrj;
        $trjType = "AMBER";
	$getByteOffset = \&GetAmberByteOffset; 
    } else {
	$getSnapshot = \&ParseLAMMPSTrj;
	$trjType = "LAMMPS";
	$getByteOffset = \&GetLammpsByteOffset;
    }
    print "Done\n";
}

sub determineHelixLayout() {
    my ($P_File) = $_[0];
    my (@helix);

    my $hl = Packages::HelixLayout->spawn();
    $hl->DetermineHelixLayout(
                              $P_File->{"Molecule"}->{"major_groove"},
                              $P_File->{"Molecule"}->{"minor_groove"},
                              $P_File->{"Molecule"}->{"is3PrimeIn"},
                              $P_File->{"Molecule"}->{"bases_at_end"},
                              $P_File->{"Molecule"}->{"total_bases"},
                              $P_File->{"Molecule"}->{"crossovers"}
                              );

    @helix = $hl->GetHelixInfo();
    return \@helix;
}

sub setLayout {
    my ($P_File) = $_[0];
    my (@helix, $i, $j, $helixC, $chainC, $regionC, @tmp, $k, $rec, $len);

    $helixC = $chainC = $regionC = 0;

    for $i (sort numerically keys %{ $P_File->{Layout} }) {
	$chainC = 0;
	for $j (sort numerically keys %{ $P_File->{Layout}{$i} }) {
	    $regionC = 0;
	    @tmp = split /\,/, $P_File->{Layout}{$i}{$j};
	    for $k (@tmp) {
		if ($k =~ /(\d+)\-(\d+)/) {
		    $rec = (
			    {
				"StartUnit" => $1,
				"EndUnit"   => $2,
				"Strand"    => $chainC + 1,
			    }
			    );
		} else {
		    $k = Trim($k);
		    $rec = (
			    {
				"StartUnit" => $k,
				"EndUnit"   => $k,
				"Strand"    => $chainC + 1,
			    }
			    );
		}
		$helix[$helixC][$chainC][$regionC] = $rec;
		$len = abs($rec->{EndUnit} - $rec->{StartUnit}) + 1;
		$HC->{$helixC}{LEN} += $len;
		$regionC++;
		$rec = ();
	    }
	    $chainC++;
	}
	$helixC++;
    }

    return \@helix;
}

sub getOrder {
    my ($i, $j) = @_;
    my (@tmp, $k);

    if ($i < $j) {
	for $k ($i .. $j) {
	    push @tmp, $k;
	}
    } else {
	for $k (reverse $j .. $i) {
	    push @tmp, $k;
	}
    }

    return \@tmp;

}

sub fixPDB {
    my ($pdbFile) = $_[0];
    my ($cmd);

    $cmd = "/ul/tpascal/scripts/stripNaH20.pl $pdbFile $pdbFile";
    die "ERROR: Cannot execute $cmd" if (system("$cmd >& junk"));

    $cmd = "/ul/maiti/src/util/scripts/fixcurvepdb.pl $pdbFile";
    die "ERROR: Cannot execute $cmd" if (system("$cmd >& junk"));
}

sub getHelicalParms {
    my ($i, $helixC, $cmd, $strandLen);

    for $i (0 .. $#{ $LAYOUT }) {
	$helixC = $i + 1;
	$strandLen = getStrandLen($LAYOUT->[$i]);

	chdir "helix${helixC}";
	$cmd = "/ul/tpascal/scripts/do_curve_anal.pl ../pdbfiles/" . $filebase . 
	    " 1 " . ($HC->{$i}{COUNT} - 1) . " $helixC " . ($HC->{$i}{LEN}/2);
	die "ERROR: Command aborted: $cmd" if (system("$cmd >& junk"));
	
	$cmd = "/ul/tpascal/scripts/get_helical_parms_new.pl " . $filebase .
	    "_Helix${helixC}_lis 1 " . ($HC->{$i}{COUNT} - 1);
	die "ERROR: Command aborted: $cmd" if (system("$cmd >& junk"));
	chdir "../";
    }
}

sub getStrandLen {
    my ($HELIX) = $_[0];
    my ($i, $j, $k, $counter);
    
    $counter = 0;
    for $j (@{ $HELIX->[0] }) {
	if ($j->{EndUnit} < $j->{StartUnit}) {
	    $counter += $j->{StartUnit} - $j->{EndUnit} + 1;
	} else {
	    $counter += $j->{EndUnit} - $j->{StartUnit} + 1;
	}
    }
    return $counter;
}

sub numerically {
    ($a<=>$b);
}

sub showHelp {
    my ($error_str) = "usage: $0 -p parmFile -t \"trjSelect\" -f [filebase]\n";

    $error_str .= "\nArguments:\n";
    $error_str .= sprintf("%-30s%-50s\n", "-p parmfile", "The parameter file");    
    $error_str .= sprintf("%-30s%-50s\n", "-t trajSelection", "trajectory selection. ex \":It1-10:2\" frames 1 - 10 every 2");
    $error_str .= sprintf("%-30s%-50s\n", " ", "can have multiple selections within the quotes. for all use \"*\"");
    $error_str .= sprintf("%-30s%-50s\n", "-f [filebase]", "(Optional) The prefix for the saved files");

    die "$error_str";
}


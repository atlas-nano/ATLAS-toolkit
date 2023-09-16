#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use CERIUS2 qw(parseCerius2FF);
use General qw(FileTester LoadElements);
use File::Basename;
use Getopt::Std qw(getopt);

#use BIOGRAF qw (saveBiografFF);

sub init;
sub saveBiografFF;
sub matchElement;
sub stringOrder;
sub findAddedH;
sub getType;
sub fixAtomTypes;
sub addMissingParms;
sub setDefaultParms;
sub pruneParms;
sub removeEmptyParms;
sub removeInvalidParms;

my ($FFs, $saveName, $fixOpt, $bgfFile, $pruneParms);
my ($header, $ELEMENTS, $PARMS, $i);
$ELEMENTS = LoadElements;

$|++;
&init;
for $i (@{ $FFs }) {
    print "Parsing Cerius2 FF $i...";
    $PARMS = parseCerius2FF($i, 0, $PARMS);
    print "Done\n";
}
&pruneParms($PARMS, $bgfFile) if (defined($bgfFile) and $pruneParms);
&addMissingParms($PARMS, $FFs, $bgfFile) if (defined($bgfFile));
&fixAtomTypes($PARMS) if ($fixOpt);
print "Creating Biograf FF $saveName...";
&saveBiografFF($PARMS, $saveName, $header);
print "Done\n";

sub pruneParms {
    my ($parms, $bgfFile) = @_;
    my ($i, $j, $k, $ffList, $count, $pList, $tmp);

    print "Eliminating unused parameters...";
    open BGFFILE, $bgfFile or die "ERROR: Cannot read from $bgfFile: $!\n";
    while (<BGFFILE>) {
	chomp;
	if ($_ =~ /^(HETATM|ATOM)\s+(.*)$/) {
	    @{ $tmp } = split /\s+/, $2;
	    $ffList->{ $tmp->[8] } = 1;
	}
    }
    close BGFFILE;

    for $i (keys %{ $parms->{ATOMTYPES} }) {
	delete $parms->{ATOMTYPES}{$i} if (! exists($ffList->{$i}));
    }
    for $i (keys %{ $parms }) {
	next if ($i =~ /PARMS|ATOMTYPES/);
	for $j (keys %{ $parms->{$i} }) {
	    delete $parms->{$i}{$j} if (! exists($ffList->{$j}));
	    &removeInvalidParms($parms->{$i}{$j}, $ffList);
	}
	&removeEmptyParms($parms->{$i});
    }
    print "Done\n";
}

sub removeInvalidParms {
    my ($curr, $ffList) = @_;
    my ($i);

    for $i (keys %{ $curr }) {
	last if (exists($curr->{$i}{VALS}));
	if (! exists($ffList->{$i})) {
	    delete $curr->{$i};
	    next;
	} else {
	    &removeInvalidParms($curr->{$i}, $ffList);
	}
    }
}
	
sub removeEmptyParms {
    my ($parmData) = $_[0];
    my ($i);

    for $i (keys %{ $parmData }) {
        last if (exists($parmData->{$i}{VALS}));
        next if ($i eq "counter" || $i eq "TYPE" || $i eq "Lammps");
        if (! keys %{ $parmData->{$i} }) {
            delete $parmData->{$i};
        } else {
            removeEmptyParms($parmData->{$i});
            delete $parmData->{$i} if (! keys %{ $parmData->{$i} });
        }
    }
}

sub addMissingParms {
    my ($parms, $ffList, $bFile) = @_;
    my ($convert_str, @tmp, $i, $vType, $curr, $count);

    $count = 0;
    print "Adding missing valence parameters...";
    $convert_str = "$Bin/createLammpsInput.pl -b $bFile -f \"@{$ffList}\" -s tmp";
    open INCMD, "$convert_str |" or die "ERROR: Cannot execute $convert_str: $!\n";
    while (<INCMD>) {
	chop;
	if ($_ =~ /MISSING (\w+) TERMS/) {
	   $vType = $1;
	} elsif ($vType and $_ =~ /(\S+): \(occurred/) {
	   @tmp = split /\-/,$1;
	   $curr = \%{ $parms->{$vType} };
	   for $i (@tmp) {
	 	$curr = \%{ $curr->{$i} };
	   }
	   $curr->{1} = setDefaultParms($parms, $vType, $1);
	   $count++;
	}
    }
   close INCMD;
   system("rm -fr tmp.lammps.pbs in.tmp_singlepoint in.tmp data.tmp");
   print "added $count parameters...Done\n"
}

sub setDefaultParms {
    my ($parms, $vType, $cKey) = @_;
    my ($rec);

    $cKey =~ s/\-/ /g;
    $rec = (
                {
                        "INDEX"    => 1,
                        "KEY"      => $cKey,
                }
           );
    if ($vType eq "BONDS") {
	$rec->{TYPE} = "HARMONIC";
	@{ $rec->{VALS} } = (0,0);
    } elsif ($vType eq "ANGLES") {
	$rec->{TYPE} = "THETA_HARM";
	@{ $rec->{VALS} } = (0,0);
    } elsif ($vType eq "TORSIONS") {
	$rec->{TYPE} = "SHFT_DIHDR";
	@{ $rec->{VALS} } = (0,2,0);
	$rec->{"do_scale"} = 0;
	$rec->{NUM} = 1;
	$rec->{PER} = 3;
	$rec->{CTYPE} = "";
	$rec->{"1_4scale"} = 0.5;
    } else {
	$rec->{TYPE} = "IT_JIKL";
	@{ $rec->{VALS} } = (0,0,2);
   }

    return $rec;
}

sub init {
    my (%OPTS, $ffStr);

    getopt('fsnbp',\%OPTS);
    die "usage: $0 -f \"cerius2 ff(s)\" -n ('normalize' atom type = no) -b (bgf file = add missing parms = 0) -p (write only used parameters=no. requires BGF) -s (savename)\n" 
	if (! exists($OPTS{f}));
    print "Initializing...";
    ($ffStr, $saveName, $fixOpt, $bgfFile, $pruneParms) = ($OPTS{f}, $OPTS{s}, $OPTS{n}, $OPTS{b}, $OPTS{p});
    while ($ffStr =~ /(\S+)/g) {
	push @{ $FFs }, $1 if (-e $1 and -r $1 and -T $1);
    }
    die "ERROR: No valid CERIUS2 forcefields found while searching \"$ffStr\"!\n" if (! defined($FFs));
    $fixOpt = 0 if (! defined($fixOpt) or $fixOpt !~ /1|yes/i);
    $fixOpt = 1 if ($fixOpt =~ /1|yes/i);
    undef $bgfFile if (defined($bgfFile) and (! -e $bgfFile or ! -r $bgfFile or ! -T $bgfFile));
    $pruneParms = 0 if (! defined($pruneParms) or $pruneParms !~ /1|yes/i or !defined($bgfFile));
    $pruneParms = 1 if ($pruneParms =~ /1|yes/i and defined($bgfFile));

    if (! defined($saveName)) {
	$saveName = basename($FFs->[0]);
	$saveName =~ s/\.\w+$/_biograf\.par/;
    }
    open HEADER, "$Bin/dat/biograf_header.dat" or die "ERROR: Cannot open header file!\n";
    while (<HEADER>) {
	$header .= $_;
    }
    close HEADER;
    print "Done\n";
}

sub fixAtomTypes {
    my ($parms) = $_[0];
    my ($i, $ffType, $ele, $curr);

    for $i (keys %{ $parms->{ATOMTYPES} }) {
        $curr = $parms->{ATOMTYPES}{$i};
        $ele = matchElement($curr->{ATOM});
	$curr->{LABEL} = "${1}_${2}" 
	    if ($curr->{LABEL} !~ /_/ and ($curr->{LABEL} =~ /^($ele)(.*)/ or $curr->{LABEL} =~ /^(.)(.*)/));
    }
}

sub saveBiografFF {
    my ($parms, $save, $header) = @_;
    my ($i, $ele, $curr, $addedH, @tmp, $junk); 
    my ($j, $k, $l, $val, $t, $counter);

    @tmp = sort stringOrder keys %{ $parms->{ATOMTYPES} };
    open BIOGRAF, "> $save" || die "ERROR: Cannot create Biograf FF $save: $!\n";
    # FF Options
    print BIOGRAF "$header";
    print BIOGRAF "FFLABEL    ATNO      MASS CHARG HYB BND CPK \#IH \#LP\n";
    for $i (@tmp) {
	$curr = $parms->{ATOMTYPES}{$i};
	$ele = matchElement($curr->{ATOM});
	printf BIOGRAF "%-14s%-4d%-10.4f%-6.1d%-4d%-4d%-4d%-4d%-4d\n", $curr->{LABEL}, $ele, $curr->{MASS},
	$curr->{CHARGE}, $curr->{NUMBONDS}, 3, 5, $curr->{OTHER}, $curr->{LONEPAIRS};
    }
    print BIOGRAF "*\n";

    #Added H
    print BIOGRAF "ADDED H   HYDROGEN  1IMPLCTH  2IMPLCTH  3IMPLCTH  4IMPLCTH\n";
    $addedH = findAddedH($parms->{BONDS});
    print BIOGRAF "*\nLONEPAIRS\n*\nGASTEIGER          A         B         C        X+\n*\n";
    
    # VDWS
    print BIOGRAF "VDW AT ITY       RNB      DENB     SCALE\n";
    for $i (@tmp) {
	$curr = $parms->{VDW}{$i}{$i}{1};
	printf BIOGRAF "%-6s%4d", $parms->{ATOMTYPES}{$i}{LABEL}, getType($curr->{TYPE}, "VDW");
	for $val (@{ $curr->{VALS} }) {
	    printf BIOGRAF "%10.4f", $val;
	}
	print BIOGRAF "\n";
    }
    print BIOGRAF "*\nAUTOTYPE  ELEMENT   HYBRIDIZATION RING_SIZE  REQUIREMENTS  FACTOR\n*\n";
    
    # BONDS
    print BIOGRAF "BONDSTRTCH  TYPE FORC CNST  BND DIST    DE/CUB         E         F  BOND DIP\n";
    for $i (sort stringOrder keys %{ $parms->{BONDS} }) {
	for $j (sort stringOrder keys %{ $parms->{BONDS}{$i} }) {
	    @{ $junk } = keys %{ $parms->{BONDS}{$i}{$j} };
	    $curr = $parms->{BONDS}{$i}{$j}{ shift @{ $junk } };
	    printf BIOGRAF "%-5s-%-5s%5d", $parms->{ATOMTYPES}{$i}{LABEL}, 
		$parms->{ATOMTYPES}{$j}{LABEL}, getType($curr->{TYPE}, "BONDS");
	    for $val (@{ $curr->{VALS} }) {
		printf BIOGRAF "%10.4f", $val;
	    }
	    print BIOGRAF "\n";
	}
    }
    print BIOGRAF "*\n";

    # ANGLES
    print BIOGRAF "ANGLE             TYPE FORC CNST EQUIL ANG         D         E         F\n";
    for $i (sort stringOrder keys %{ $parms->{ANGLES} }) {
	for $j (sort stringOrder keys %{ $parms->{ANGLES}{$i} }) {
	    for $k (sort stringOrder keys %{ $parms->{ANGLES}{$i}{$j} }) {
		@{ $junk } = keys %{ $parms->{ANGLES}{$i}{$j}{$k} };
		$curr = $parms->{ANGLES}{$i}{$j}{$k}{ shift @{ $junk } };
		printf BIOGRAF "%-5s-%-5s-%-5s%5d", $parms->{ATOMTYPES}{$i}{LABEL}, 
			$parms->{ATOMTYPES}{$j}{LABEL}, $parms->{ATOMTYPES}{$k}{LABEL}, 
			getType($curr->{TYPE}, "ANGLES");
		for $val (@{ $curr->{VALS} }) {
		    printf BIOGRAF "%10.4f", $val;
		}
		print BIOGRAF "\n";
	    }
	}
    }
    print BIOGRAF "*\n";

    # TORSIONS
    $counter = 0;
    print BIOGRAF "TORSION                 CASE   BARRIER    PERIOD CISMIN(1)\n";
    for $i (sort stringOrder keys %{ $parms->{TORSIONS} }) {
	for $j (sort stringOrder keys %{ $parms->{TORSIONS}{$i} }) {
	    for $k (sort stringOrder keys %{ $parms->{TORSIONS}{$i}{$j} }) {
		for $l (sort stringOrder keys %{ $parms->{TORSIONS}{$i}{$j}{$k} }) {
		    $counter++;
		    @{ $junk } = keys %{ $parms->{TORSIONS}{$i}{$j}{$k}{$l} };
		    $curr = $parms->{TORSIONS}{$i}{$j}{$k}{$l}{shift @{ $junk } };
		    $t = 0;
		    while ($t <= ($#{ $curr->{VALS} } - 2)) {
			printf BIOGRAF "%-5s-%-5s-%-5s-%-5s%5d", $parms->{ATOMTYPES}{$i}{LABEL},
				$parms->{ATOMTYPES}{$j}{LABEL}, $parms->{ATOMTYPES}{$k}{LABEL},
				$parms->{ATOMTYPES}{$l}{LABEL}, $counter;
			$curr->{VALS}[$t] *= 4;
			if ($curr->{VALS}[($t + 2)] == 180) {
			    $curr->{VALS}[($t + 2)] = 1;
			} else {
			    $curr->{VALS}[($t + 2)] = -1;
			}
			for $val (0 .. 2) {
			    printf BIOGRAF "%10.4f", $curr->{VALS}[($t + $val)];
			}
			print BIOGRAF "\n";
			$t += 3;
		    }
		}
	    }
	}
    }
    print BIOGRAF "*\n";

    # INVERSIONS
    print BIOGRAF "INVERSION (CENT AT 1ST) TYPE  FRC CNST  EQU ANGL         D         E         F\n";
    for $i (sort stringOrder keys %{ $parms->{INVERSIONS} }) {
	for $j (sort stringOrder keys %{ $parms->{INVERSIONS}{$i} }) {
	    for $k (sort stringOrder keys %{ $parms->{INVERSIONS}{$i}{$j} }) {
		for $l (sort stringOrder keys %{ $parms->{INVERSIONS}{$i}{$j}{$k} }) {
		    @{ $junk } = keys %{ $parms->{INVERSIONS}{$i}{$j}{$k}{$l} };
		    $curr = $parms->{INVERSIONS}{$i}{$j}{$k}{$l}{shift @{ $junk }};
		    printf BIOGRAF "%-5s-%-5s-%-5s-%-5s%5d", $parms->{ATOMTYPES}{$i}{LABEL},
			$parms->{ATOMTYPES}{$j}{LABEL}, $parms->{ATOMTYPES}{$k}{LABEL},
			$parms->{ATOMTYPES}{$l}{LABEL}, getType($curr->{TYPE},"INVERSIONS");
		    for $val (@{ $curr->{VALS} }) {
			printf BIOGRAF "%10.4f", $val;
		    }
		    print BIOGRAF "\n";
		}
	    }
	}
    }
    print BIOGRAF "*\nHBOND       TYPE    -DE HB     RE HB\nEND OF DATA\n";
    close BIOGRAF;
		    
}

sub findAddedH {
    my ($bonds) = $_[0];
    my ($atom1, $atom2, $ele1, $ele2);
    my ($added_h);

    for $atom1 (sort stringOrder keys %{ $bonds }) {
	$ele1 = matchElement($PARMS->{ATOMTYPES}{$atom1}{ATOM});
	for $atom2 (sort stringOrder keys %{ $bonds->{$atom1} }) {
	    $ele2 = matchElement($PARMS->{ATOMTYPES}{$atom2}{ATOM});
	    if ($ele1 == 1) {
		$added_h .= sprintf("%-10s%-10s\n",$atom2, $atom1);
	    } elsif ($ele2 == 1) {
		$added_h .= sprintf("%-10s%-10s\n",$atom1, $atom2);		
	    }
	}
    }

    return $added_h;
}

sub matchElement {
    my ($ele_name) = $_[0];
    my ($ele_num, $i);

    for $i (keys %{ $ELEMENTS }) {
	if ($ELEMENTS->{$i}{SYMBOL} eq $ele_name) {
	    $ele_num = $i;
	    last;
	}
    }

    die "ERROR: No valid element found while search for $ele_name\n"
	if (! defined($ele_num));

    return $ele_num;
}

sub stringOrder {
    ($a cmp $b);
}

sub getType {
    my ($type, $i) = @_;
    my (%TYPES) = (
		   "VDW" => {
		       "LJ_6_12" => 1,
		   },
		   "BONDS" => {
		       "HARMONIC" => 1,
		   },
		   "ANGLES" => {
		       "THETA_HARM" => 21,
		   },
		   "INVERSIONS" => {
		       "IT_JIKL" => 3,
		   },
		   );

    die "ERROR: No type $type in $i found!\n" if (! exists($TYPES{$i}{$type}));
    return $TYPES{$i}{$type};
}

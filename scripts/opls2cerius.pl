#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/Packages";
use Math::Trig;

use strict;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

use General qw(FileTester LoadElements);
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetMOL2FileInfo GetPDBFileInfo createHeaders);

sub init;
sub getOplsParams;
sub writeData;
sub getRealVDW;
sub printValence;
sub getOplsHeader;
sub createMAEfile;
sub createCOMfile;
sub findEquivalence;
sub GetMAEFileInfo;
sub checkVDWParms;
sub numerically { $a<=>$b }

my ($inputFile, $ffName, $newBGF, $readFunc, $inputType);
my ($ATOMS, $BONDS, $HEADERS, $PARMS, $FILES, $EQUIV);

$|++;
&init;
if ($inputType !~ /MAE/) {
    print "Parsing $inputType file $inputFile...";
    ($ATOMS, $BONDS, $HEADERS) = $readFunc->($FILES->{STRUCTURE},1);
    @{ $HEADERS } = @{ $ATOMS->{HEADER} } if ($inputType =~ /MOL2/);
    delete $ATOMS->{HEADER} if ($inputType =~ /MOL2/);
    &createMAEfile($ATOMS, $FILES);
} else {
    print "Parsing MAE file $inputFile...";
    die "ERROR: Cannot write to current directory!\n" if (system("cp $inputFile ./$FILES->{PREFIX}.mae"));
    $HEADERS = createHeaders(undef,$FILES->{PREFIX});
    ($ATOMS,$BONDS) = GetMAEFileInfo($FILES->{PREFIX} . ".mae");
}
print "Done\nObtaining Opls parameters...";
&createCOMfile($FILES->{PREFIX});
$PARMS = getOplsParams($FILES, $ATOMS);
print "Done\nChecking VDW parms...\n";
&checkVDWParms($PARMS);
print "Creating CERIUS2 FF $ffName...";
&writeData($PARMS, $ffName, $FILES, $EQUIV);
print "Done\nCreating updated bgf $newBGF...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $newBGF);
print "Done\n";

sub checkVDWParms {
    my ($ffdata) = $_[0];
    my ($i, $j, $k, $d0, $r0, @tmp);

    for $i (keys %{ $ffdata->{VDW} }) {
        for $j (keys %{ $ffdata->{VDW}{$i} }) {
    	    $d0 = $ffdata->{VDW}{$i}{$j}{1}{D};
	        $r0 = $ffdata->{VDW}{$i}{$j}{1}{R};
	        @tmp = sort numerically keys %{ $ffdata->{VDW}{$i}{$j} };
            $d0 = $ffdata->{VDW}{$i}{$j}{$tmp[0]}{D};
            $r0 = $ffdata->{VDW}{$i}{$j}{$tmp[0]}{R};
#    	    for $k (@tmp) {
#	        	warn "\t\tWill overwrite D0 for $i $j from $d0 to $ffdata->{VDW}{$i}{$j}{$k}{D}\n"
#		            if(! $ffdata->{VDW}{$i}{$j}{$k}{D} or ($d0/$ffdata->{VDW}{$i}{$j}{$k}{D})>1.00001 or ($d0/$ffdata->{VDW}{$i}{$j}{$k}{D})<0.99999);
#                	    	warn "\t\tWill overwrite R0 for $i $j from $r0 to $ffdata->{VDW}{$i}{$j}{$k}{R}\n"
#	        	    if(! $ffdata->{VDW}{$i}{$j}{$k}{R} or ($r0/$ffdata->{VDW}{$i}{$j}{$k}{R})>1.00001 or ($r0/$ffdata->{VDW}{$i}{$j}{$k}{R})<0.99999);
#		        $d0 = $ffdata->{VDW}{$i}{$j}{$k}{D};
#    		$r0 = $ffdata->{VDW}{$i}{$j}{$k}{R};
#	        }
    	    delete($ffdata->{VDW}{$i}{$j});
	        $ffdata->{VDW}{$i}{$j}{D} = $d0;
	        $ffdata->{VDW}{$i}{$j}{R} = $r0;
	    }
    }
}

sub writeData {
    my ($ffdata, $outfile, $files, $equiv) = @_;
    my ($offDiag, $i, $j, $header);

    $header = getOplsHeader($ffdata, $files, $equiv);
    open OUTDATA, "> $outfile" || die "ERROR: Cannot write to $outfile: $!\n";
    for $i (keys %{ $ffdata->{ATOMTYPES} }) {
	next if exists($ffdata->{VDW}{$i}{$i});
	$ffdata->{VDW}{$i}{$i}{R} =  $ffdata->{VDW}{$i}{$i}{D} = 0;
    }
    print OUTDATA "${header}DIAGONAL_VDW\n";
    for $i (keys %{ $ffdata->{VDW} }) {
	for $j (keys %{ $ffdata->{VDW}{$i} }) {
	    if ($i eq $j) {
		printf OUTDATA " %-11s LJ_6_12   %21.15f   %21.15f\n",$i,$ffdata->{VDW}{$i}{$j}{R},$ffdata->{VDW}{$i}{$j}{D};
	   } else {
		$offDiag .= sprintf(" %-9s%-9s   LJ_6_12   %21.15f    %21.15f\n",$i,$j,$ffdata->{VDW}{$i}{$j}{R},$ffdata->{VDW}{$i}{$j}{D});
	   }
	}
    }
    print OUTDATA "END\n#\nATOM_TYPING_RULES\nEND\n#\nOFF_DIAGONAL_VDW\n${offDiag}END\n#\n$ffdata->{VALENCE}";
    close OUTDATA;
}

sub getOplsHeader {
    my ($ff, $files, $equiv) = @_;
    my ($header, $i, $j, $atypes, $flag, @tmp, $hyb, $ELEMENTS, $list);
    my ($errstr, $invalid);

    $invalid = 0;
    $errstr = "";
    $ELEMENTS = LoadElements();
    for $i (keys %{ $ff->{VDW} }) {
	$atypes->{$i} = 0;
	for $j (keys %{ $ff->{VDW}{$i} }) {
	    $atypes->{$j} = 0;
	}
    }

    $header = "";
    open ATOMTYPEFILE, $files->{ATOMTYPES} || die "ERROR: Cannot open atomtype file \"$files->{ATOMTYPES}\":$!\n";
    while (<ATOMTYPEFILE>) {
	chomp;
	if ($_ =~ /^\d+\s+\d+\s+\d+\.\d+\s+\S+\s+/) {
	    @tmp = split /\s+/,$_;
	    @{ $list } = ($tmp[3]);
		print "working on $tmp[3]...";
	    next if (! exists($atypes->{$tmp[3]}) && ! ($list = findEquivalence($equiv, $atypes, $tmp[3])));
		print "found!\n";
	    for $i (@{ $list }) {
		next if ($atypes->{$i});
		$atypes->{$i} = 1;
		$tmp[1] = 1 if ($tmp[1] !~ /^\d+/);
		$hyb = 1;
		if($tmp[16] eq "LIN") { $hyb = 2;
		}elsif($tmp[16] eq "TRI") { $hyb = 3;
		}elsif($tmp[16] eq "TET") { $hyb = 4;
		}elsif($tmp[16] eq "TBY") { $hyb = 6;
		}elsif($tmp[16] eq "OCT") { $hyb = 8; }
		$header .= sprintf(" %-10s %-2s %12.5f %7.4f %3d %3d %3d\n",$i,$ELEMENTS->{$tmp[1]}{SYMBOL},$tmp[2],0,$hyb,$tmp[17],0);
	    }
	}
    }
    close ATOMTYPEFILE;

    #px6
    $atypes->{PX6} = 1;
	$header .= sprintf(" %-10s %-2s %12.5f %7.4f %3d %3d %3d\n","PX6","P",19,0,0,0,0);

    #bf4
    $atypes->{BT} = 1;
	$header .= sprintf(" %-10s %-2s %12.5f %7.4f %3d %3d %3d\n","BT","B",3,0,0,0,0);

    #no3
    $atypes->{ON} = 1;
	$header .= sprintf(" %-10s %-2s %12.5f %7.4f %3d %3d %3d\n","ON","O",6,0,0,0,0);

    for $i (keys %{ $atypes }) {
	if (!$atypes->{$i}) {
	    $errstr .= "ERROR: Atomtype $i not found!\n" if (! $atypes->{$i});
	    $invalid = 1;
	}
    }
    die "$errstr" if($invalid);

    $flag = 0;
    local $/=undef;
    open HEADERFILE, $files->{OPLS_HEADER} || die "ERROR: Cannot open $files->{OPLS_HEADER}: $!\n";
    $header = <HEADERFILE> . "ATOMTYPES\n${header}END\n#\n";
    close HEADERFILE;

    return $header;
}

sub getOplsParams {
    my ($files, $atoms) = @_;
    my ($DATA, $flag, $tot, $TYPES, $count, $VALENCE, $tmp, $logFile, $prefix);
    my ($type1, $type2, $type3, $type4, $rec, $phase, $val, $execCmd, $iC);
    my ($atom1, $atom2, $atom3, $atom4);
    $prefix=$files->{PREFIX};
    $type1 = $type2 = $tot = $count = 0;
    $flag = "";
    $rec = ();

    $execCmd = "$files->{BMIN} $prefix -WAIT";
    die "ERROR: Cannot execute \"$execCmd\": $!\n" if (system("$execCmd >& /dev/null"));

    $logFile = "${prefix}.log";
    open LOGFILE, "$logFile" or die "ERROR: Cannot open log file \"$logFile\": $!\n";
    while(<LOGFILE>) {
	    chomp;
	    if ($_ =~ /transfer_params stretch/) {
            $flag = "bond";
    	} elsif ($_ =~ /transfer_params bends/) {
            $flag = "angle";
        } elsif ($_ =~ /transfer_params proper torsions/) {
            $flag = "torsion";
        } elsif ($_ =~ /transfer_params improper torsions/) {
            $flag = "improper";
        } elsif ($_ =~ /number of atoms/) {
            $flag = "charge";
	    } elsif ($flag eq "bond" && $_ =~ /atoms:\s+(\d+)\s+(\d+)/) {
            $atom1 = $1; $atom2 = $2;
            $rec = ();
        } elsif ($_ =~ /Quality of Force Field Parameters in Use/) {
            $flag = "";

        ###BONDS###
    	} elsif ($flag eq "bond" && $_ =~/(K|r0)\s+(\-?\d+\.?\d*e?\-?\+?\d*)/) {
            push @{ $rec->{VALS} }, $2;
	    } elsif ($flag eq "bond" && $_ =~ /^comment\s+\d+\s+\S+\s+\d+\s+(\S+)\-(\S+)\s+=/) {
            $type1 = $1; $type2 = $2;
	        $atoms->{$atom1}{FFTYPE} = $type1;
    	    $atoms->{$atom2}{FFTYPE} = $type2;
	        $DATA->{ATOMTYPES}{$type1} = 1;
	        $DATA->{ATOMTYPES}{$type2} = 1;
	        @{ $rec->{TYPES} } = ($type1,$type2);
    	    @{ $rec->{TYPES} } = reverse @{ $rec->{TYPES} } if ($type1 lt $type2);
	        push @{ $VALENCE->{BONDS}{DATA} }, $rec if ( ! exists($VALENCE->{BONDS}{HOLDER}{$type1}{$type2}) &&
			                            				(! exists($VALENCE->{BONDS}{HOLDER}{$type2}) || 
                                                         ! exists($VALENCE->{BONDS}{HOLDER}{$type2}{$type1})));
    	    $VALENCE->{BONDS}{HOLDER}{$type1}{$type2} = 1;
	        $VALENCE->{BONDS}{HEADER} = "BOND_STRETCH";
	        $VALENCE->{BONDS}{TYPE} = "    HARMONIC   ";

        ###ANGLES###
	    } elsif ($flag eq "angle" && $_ =~ /atoms:\s+(\d+)\s+(\d+)\s+(\d+)/) {
            $atom1 = $1; $atom2 = $2; $atom3 = $3;
            $rec = ();
    	} elsif ($flag eq "angle" && $_ =~/(K|theta0)\s+(\-?\d+\.?\d*e?\-?\+?\d*)/) {
            push @{ $rec->{VALS} }, $2;
	    } elsif ($flag eq "angle" && $_ =~ /^comment\s+\d+\s+(\S+)\-(\S+)\-(\S+)\s+=/) {
            $type1 = $1; $type2 = $2; $type3 = $3;
	        @{ $rec->{TYPES} } = ($type1,$type2,$type3);
	        $val = sprintf("%.04f",$rec->{VALS}[1]*180/pi);
	        @{ $rec->{VALS} } = ($rec->{VALS}[0],$val);
	        push @{ $VALENCE->{ANGLES}{DATA} }, $rec if ( ! exists($VALENCE->{ANGLES}{HOLDER}{$type1}{$type2}{$type3}) &&
			                            				 (! exists($VALENCE->{ANGLES}{HOLDER}{$type3}) || 
                                                          ! exists($VALENCE->{ANGLES}{HOLDER}{$type3}{$type2}) || 
                                                          ! exists($VALENCE->{ANGLES}{HOLDER}{$type3}{$type2}{$type1})));
	        $VALENCE->{ANGLES}{HOLDER}{$type1}{$type2}{$type3} = 1;
	        $VALENCE->{ANGLES}{HEADER} = "ANGLE_BEND";
	        $VALENCE->{ANGLES}{TYPE} = "    THETA_HARM ";

        ###TORSIONS###
	    } elsif ($flag eq "torsion" && $_ =~ /atoms:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
            $atom1 = $1; $atom2 = $2; $atom3 = $3; $atom4 = $4;
            $rec = ();
    	} elsif ($flag eq "torsion" && $_ =~/^\s*(\-?\d+\.?\d*e?\-?\+?\d*)/) {
            push @{ $rec->{VALS} }, $1*2;
        } elsif ($flag eq "torsion" && $_ =~ /^comment\s+\d+\s+(\S+)\-(\S+)\-(\S+)\-(\S+)\s+=/) {
            $type1 = $1; $type2 = $2; $type3 = $3; $type4 = $4;
	        next if (exists($VALENCE->{TORSIONS}{HOLDER}{$type1}{$type2}{$type3}{$type4}) ||
		            (exists($VALENCE->{TORSIONS}{HOLDER}{$type4}) && 
                     exists($VALENCE->{TORSIONS}{HOLDER}{$type4}{$type3}) &&
                     exists($VALENCE->{TORSIONS}{HOLDER}{$type4}{$type3}{$type2}) &&
                     exists($VALENCE->{TORSIONS}{HOLDER}{$type4}{$type3}{$type2}{$type1})));
	        $count = 0;
            
            $tmp = shift @{ $rec->{VALS } };
            @{ $rec->{VALS} } = (@{ $rec->{VALS} }, $tmp);   
            @{ $rec->{TYPES} } = ($type1,$type2,$type3,$type4);
            push @{ $VALENCE->{TORSIONS}{DATA} }, $rec;
	        $VALENCE->{TORSIONS}{HOLDER}{$type1}{$type2}{$type3}{$type4} = 1;
	        $VALENCE->{TORSIONS}{HEADER} = "TORSIONS";
            $VALENCE->{TORSIONS}{TYPE} = "    OPLS        ";

        ###IMPROPERS###
	    } elsif ($flag eq "improper" && $_ =~ /atoms:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
            $atom1 = $1; $atom2 = $2; $atom3 = $3; $atom4 = $4;
            $rec = ();
    	} elsif ($flag eq "improper" && $_ =~/^\s*(\-?\d+\.?\d*e?\-?\+?\d*)/) {
            push @{ $rec->{VALS} }, $1;
        } elsif ($flag eq "improper" && $_ =~ /^comment/) {
            $type1 = $atoms->{$atom3}{FFTYPE}; $type2 = $atoms->{$atom2}{FFTYPE}; $type3 = $atoms->{$atom1}{FFTYPE}; $type4 = $atoms->{$atom4}{FFTYPE};
            @{ $rec->{TYPES} } = ($type1,$type2,$type3, $type4);
	        if    (abs($rec->{VALS}[0])) { @{ $rec->{VALS} } = ($rec->{VALS}[0]*2,180,1); }
	        elsif (abs($rec->{VALS}[1])) { @{ $rec->{VALS} } = ($rec->{VALS}[1]*2,180,2); }
	        elsif (abs($rec->{VALS}[2])) { @{ $rec->{VALS} } = ($rec->{VALS}[2]*2,180,3); }
	        elsif (abs($rec->{VALS}[3])) { @{ $rec->{VALS} } = ($rec->{VALS}[3]*2,180,4); }
	        push @{ $VALENCE->{IMPROPER}{DATA} }, $rec if ( ! exists($VALENCE->{IMPROPER}{HOLDER}{$type1}{$type2}{$type3}{$type4}) &&
			                            				   (! exists($VALENCE->{IMPROPER}{HOLDER}{$type4}) || 
                                                            ! exists($VALENCE->{IMPROPER}{HOLDER}{$type4}{$type3}) ||
                                                            ! exists($VALENCE->{IMPROPER}{HOLDER}{$type4}{$type3}{$type2}) || 
                                                            ! exists($VALENCE->{IMPROPER}{HOLDER}{$type4}{$type3}{$type2}{$type1})));
	        $VALENCE->{IMPROPER}{HOLDER}{$type1}{$type2}{$type3}{$type4} = 1;
	        $VALENCE->{IMPROPER}{HEADER} = "INVERSIONS";
	        $VALENCE->{IMPROPER}{TYPE} = "    IT_JIKL    ";

        ###CHARGES###
        } elsif ($flag eq "charge" && $_ =~ /^(\d+)\s+(\-?\d+\.?\d*e?\-?\+?\d*)/) {
            $atoms->{$1}{CHARGE} = $2+0;

        ###VDWS###
        } elsif ($_ =~ /bopls_mmshare_getvdw atoms\s+(\d+)\s+(\d+)/) {
            $flag = "vdw";
            $type1 = $atoms->{$1}{FFTYPE}; $type2 = $atoms->{$2}{FFTYPE};
            ($type1,$type2) = ($type2,$type1) if ($type2 gt $type1);
            $rec = ();
        } elsif ($flag eq "vdw" && $_ =~ /params\s+(\-?\d+\.?\d*e?\-?\+?\d*)\s+(\-?\d+\.?\d*e?\-?\+?\d*)/) {
            @{ $rec->{VALS} } = ($1, $2);
        } elsif ($flag eq "vdw" && $_ =~ /comment/) {
	        $iC = scalar(keys %{ $DATA->{VDW}{$type1}{$type2} }) + 1;
	        &getRealVDW(\%{ $DATA->{VDW}{$type1}{$type2}{$iC} }, $rec->{VALS}[0], $rec->{VALS}[1]);
        }
    }
    close LOGFILE;
    die "ERROR: MacroModel logfile $logFile does not contain any valid info!\n" if (! defined($DATA));
    $DATA->{VALENCE} = printValence($VALENCE) . "COULOMBIC\n X        X           LIN-R-EPS\nEND\n";
    system("rm -fr ${prefix}-out.maegz ${prefix}.in ${prefix}.mae");
    return $DATA;
}

sub getOplsParams_old {
    my ($files, $atoms) = @_;
    my ($DATA, $flag, $tot, $TYPES, $count, $VALENCE, $tmp, $logFile, $prefix);
    my ($type1, $type2, $type3, $type4, $rec, $phase, $val, $execCmd,$iC);

    $prefix=$files->{PREFIX};
    $type1 = $type2 = $tot = $count = 0;
    $flag = "";
    $rec = ();

    $execCmd = "$files->{BMIN} $prefix -WAIT";
    die "ERROR: Cannot execute \"$execCmd\": $!\n" if (system("$execCmd >& /dev/null"));

    $logFile = "${prefix}.log";
    open LOGFILE, "$logFile" or die "ERROR: Cannot open log file \"$logFile\": $!\n";
    while(<LOGFILE>) {
	chomp;
	if ($_ =~ /opls_mmshare_getstr: (.+)$/) {
	    @{ $tmp } = split /\s+/,$1;
	    ($type1,$type2) = split /\-/,$tmp->[10];
	    $flag = "bond";
	    for ($type1, $type2) {
		if(! exists($TYPES->{NAMES}{$_})) {
		    $TYPES->{NAMES}{$_} = $count;
		    $TYPES->{COUNT}{$count} = $_;
		    $count++;
		    $TYPES->{TOT} = $count;
		}
	    }
	} elsif ($flag eq "bond" && $_ =~ /atoms\s+\=\s+(\d+)\s+(\d+)/) {
	    $atoms->{$1}{FFTYPE} = $type1;
	    $atoms->{$2}{FFTYPE} = $type2;
	    $DATA->{ATOMTYPES}{$type1} = 1;
	    $DATA->{ATOMTYPES}{$type2} = 1;
	} elsif ($flag eq "bond" && $_ =~/params\s+\=\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
	    $rec = ();
	    @{ $rec->{TYPES} } = ($type1,$type2);
	    @{ $rec->{TYPES} } = reverse @{ $rec->{TYPES} } if ($type1 lt $type2);
	    @{ $rec->{VALS} } = ($2,$1);
	    push @{ $VALENCE->{BONDS}{DATA} }, $rec if (! exists($VALENCE->{BONDS}{HOLDER}{$type1}{$type2}) &&
							! exists($VALENCE->{BONDS}{HOLDER}{$type2}{$type1}));
	    $VALENCE->{BONDS}{HOLDER}{$type1}{$type2} = 1;
	    $VALENCE->{BONDS}{HEADER} = "BOND_STRETCH";
	    $VALENCE->{BONDS}{TYPE} = "    HARMONIC   ";
	    $flag = "";
	} elsif ($_ =~ /opls_mmshare_getbnd: (.+)$/) {
	    @{ $tmp } = split /\s+/, $1;
	    ($type1,$type2,$type3) = split /\-/, $tmp->[8];
	    ($type1,$type2,$type3) = ($type3,$type2,$type1) if($type1 lt $type3);
	    $flag = "angle";
	} elsif ($flag eq "angle" && $_ =~ /params\s+\=\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
	    $rec = ();
	    @{ $rec->{TYPES} } = ($type1,$type2,$type3);
	    $val = sprintf("%.04f",$1*180/pi);
	    @{ $rec->{VALS} } = ($2,$val);
	    push @{ $VALENCE->{ANGLES}{DATA} }, $rec if (! exists($VALENCE->{ANGLES}{HOLDER}{$type1}{$type2}{$type3}) &&
							 ! exists($VALENCE->{ANGLES}{HOLDER}{$type3}{$type2}{$type1}));
	    $VALENCE->{ANGLES}{HOLDER}{$type1}{$type2}{$type3} = 1;
	    $VALENCE->{ANGLES}{HEADER} = "ANGLE_BEND";
	    $VALENCE->{ANGLES}{TYPE} = "    THETA_HARM ";
	    $flag = "";
        } elsif ($_ =~ /opls_mmshare_gettor: (.+)$/) {
	    @{ $tmp } = split /\s+/, $1;
	    ($type1,$type2,$type3,$type4) = split /\-/, $tmp->[8];
	    ($type1,$type2,$type3,$type4) = ($type4,$type3,$type2,$type1) if($type1 lt $type4);
            $flag = "torsion";
        } elsif ($flag eq "torsion" && $_ =~ /params\s+\=\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
	    next if (exists($VALENCE->{TORSIONS}{HOLDER}{$type1}{$type2}{$type3}{$type4}) ||
		     exists($VALENCE->{TORSIONS}{HOLDER}{$type4}{$type3}{$type2}{$type1}));
	    $count = 0;
	    $rec = ();
	    @{ $rec->{TYPES} } = ($type1,$type2,$type3,$type4);
	    for $_ ($2,$3,$4,$5) {
		$val = $_;
		$count++;
		next if ($val == 0);
		$phase = 180; #negative for even
		$phase = 0 if ($count %2 == 1);
		$val *= 2;
		if (! exists($rec->{VALS})) {
		    $rec->{VALS}[0] = sprintf("%9.4f %9.4f %9.4f",$val, $count, $phase);
		} else {
		    $rec->{VALS}[0] .= sprintf("\n%52s%9.4f %9.4f %9.4f","",$val, $count, $phase);
	 	}
	    }
	    next if (!exists($rec->{VALS}));
	    push @{ $VALENCE->{TORSIONS}{DATA} }, $rec;
	    $VALENCE->{TORSIONS}{HOLDER}{$type1}{$type2}{$type3}{$type4} = 1;
	    $VALENCE->{TORSIONS}{HEADER} = "TORSIONS";
	    $VALENCE->{TORSIONS}{TYPE} = "    SHFT_DIHDR  ";
	    $flag = "";
        } elsif ($_ =~ /opls_mmshare_getoop/) {
            $flag = "improper";
	} elsif ($flag eq "improper" && $_ =~ /atoms\s+\=\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
	    $type1 = $atoms->{$1}{FFTYPE};
	    $type2 = $atoms->{$2}{FFTYPE};
	    $type3 = $atoms->{$3}{FFTYPE};
	    $type4 = $atoms->{$4}{FFTYPE};
        } elsif ($flag eq "improper" && $_ =~ /params\s+\=\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
            $rec = ();
            @{ $rec->{TYPES} } = ($type1,$type2,$type3, $type4);
	    if (abs($1)) { @{ $rec->{VALS} } = ($1*2,180,1);
	    }elsif (abs($2)) { @{ $rec->{VALS} } = ($2*2,180,2);
	    }elsif (abs($3)) { @{ $rec->{VALS} } = ($3*2,180,3);
	    }elsif (abs($4)) { @{ $rec->{VALS} } = ($4*2,180,4); }
	    push @{ $VALENCE->{IMPROPER}{DATA} }, $rec if (! exists($VALENCE->{IMPROPER}{HOLDER}{$type1}{$type2}{$type3}{$type4}) &&
							   ! exists($VALENCE->{IMPROPER}{HOLDER}{$type4}{$type3}{$type2}{$type1}));
	    $VALENCE->{IMPROPER}{HOLDER}{$type1}{$type2}{$type3}{$type4} = 1;
	    $VALENCE->{IMPROPER}{HEADER} = "INVERSIONS";
	    $VALENCE->{IMPROPER}{TYPE} = "    IT_JIKL    ";
            $flag = "";
	}elsif ($_ =~/transfer_params vdw params/) {
	    $flag = "vdw";
	}elsif ($flag eq "vdw" && $_ =~ /^number of atoms (\d+)/) {
	    $tot = $1;
	} elsif ($flag eq "vdw" && $tot && $_ =~ /^\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/) {
	   $type1 = $TYPES->{COUNT}{$1};
	   $type2 = $TYPES->{COUNT}{$2};
	   next if (!$type1 || ! $type2);
	   if (! $type1 || ! $type2) {
	   	for ($type1, $type2) {
		    if(! exists($TYPES->{NAMES}{$_})) {
			$TYPES->{NAMES}{$_} = $count;
			$TYPES->{COUNT}{$count} = $_;
			$count++;
			$TYPES->{TOT} = $count;
		    }
		}
	   }
	   ($type1,$type2) = ($type2,$type1) if ($type2 gt $type1);
	   &getRealVDW(\%{ $DATA->{VDW}{$type1}{$type2}{1} }, $3, $4);
	} elsif ($_ =~/GENINT/) {
	   $flag = "";
	   $tot = 0;
	} elsif ($_ =~ /INTFFC/) {
	   $flag = "vdw";
	} elsif ($flag eq "vdw" && $_ =~ /^\s*\d+\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/) {
	   $type1 = $atoms->{$1}{FFTYPE};
	   $type2 = $atoms->{$2}{FFTYPE};
	   ($type1,$type2) = ($type2,$type1) if ($type2 gt $type1);
	   $iC = scalar(keys %{ $DATA->{VDW}{$type1}{$type2} }) + 1;
	   &getRealVDW(\%{ $DATA->{VDW}{$type1}{$type2}{$iC} }, $4, $3);
	} elsif ($_ =~ /bopls_mmshare_getvdw atoms\s+(\d+)\s+(\d+)/) {
	   $flag = "atom_vdw";
           $type1 = $atoms->{$1}{FFTYPE};
           $type2 = $atoms->{$2}{FFTYPE};
	   ($type1,$type2) = ($type2,$type1) if ($type2 gt $type1);
	} elsif ($flag eq "atom_vdw" && $_ =~ /^\s+params \s+(\S+)\s+(\S+)/) {
	   $iC = scalar(keys %{ $DATA->{VDW}{$type1}{$type2} }) + 1;
	   &getRealVDW(\%{ $DATA->{VDW}{$type1}{$type2}{$iC} }, $1, $2);
	} elsif ($_ =~ /transfer_params charages/) {
	   $flag = "charge";
	} elsif ($flag eq "charge" && $_ =~ /^(\d+)\s+(\-?\d+\.\d+e?\+?\-?\d*)/) {
	   $atoms->{$1}{CHARGE} = $2;
	}
    }
    close LOGFILE;
    die "ERROR: MacroModel logfile $logFile does not contain any valid info!\n" if (! defined($DATA));
    $DATA->{VALENCE} = printValence($VALENCE) . "COULOMBIC\n X        X           LIN-R-EPS\nEND\n";
    system("rm -fr ${prefix}-out.maegz ${prefix}.in ${prefix}.mae");
    return $DATA;
}

sub printValence {
    my ($valence) = $_[0];
    my ($retStr, $i, $j);

    for ("BONDS", "ANGLES", "TORSIONS", "IMPROPER") {
	    next if (! exists($valence->{$_}));
	    $retStr .= "$valence->{$_}{HEADER}\n";
	    for $j (@{ $valence->{$_}{DATA} }) {
	        for $i (@{ $j->{TYPES} }) {
		        $retStr .= sprintf(" %-8s", $i);
	        }
	        $retStr .= "$valence->{$_}{TYPE}";
	        if ($_ ne "TORSIONS" or $valence->{$_}{TYPE} !~ /SHFT_DIHDR/) {
		        for $i (@{ $j->{VALS} }) {
		            $retStr .= sprintf("%10.4f", $i);
		        }
	        } else {
		        $retStr .= $j->{VALS}[0];
	        }
	        $retStr .= "\n";
	    }
	    $retStr .= "END\n#\n";
    }

    return $retStr;
}

sub getRealVDW {
    my ($parm, $A, $C) = @_;

    $parm->{R} = $parm->{D} = 0;
    if($A > 0 && $C > 0) {
	$parm->{R} = (2*$A/$C)**(1/6);
	$parm->{D} = ($C/($parm->{R}**6))/2;
    }
}

sub createCOMfile {
    my ($prefix) = $_[0];
    my ($outStr);

    open OUTFILE, ">${prefix}.com" || die "ERROR: Cannot create ${prefix}.com: $!\n";
print OUTFILE <<DATA;
${prefix}.mae
${prefix}-out.maegz
 MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000
 FFLD      14      1      0      2     1.0000     0.0000     0.0000     0.0000
 EXNB       0      0      0      0    48.0000    48.0000     4.0000     0.0000
 DEBG      44      0      0      0     0.0000     0.0000     0.0000     0.0000
 BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000
 BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 READ      -1      0      0      0     0.0000     0.0000     0.0000     0.0000
 ELST      -1      0      0      0     0.0000     0.0000     0.0000     0.0000
 WRIT       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 END        0      0      0      0     0.0000     0.0000     0.0000     0.0000
DATA

    close OUTFILE;
}

sub createMAEfile {
    my ($atoms, $files) = @_;
    my ($maefile, $i, $convertCmd, $prefix);
   
    $prefix = $files->{PREFIX};
    open JAGIN, "> ${prefix}.in" || die "ERROR: Cannot create ${prefix}.in: $!\n";
    print JAGIN <<DATA;
&gen
basis=6-31g*
&
entry_name: ${prefix}
&zmat
DATA

    for $i (sort { ($a<=>$b); } keys %{ $atoms }) {
	printf JAGIN "%-5s %17.13f  %17.13f %17.13f\n",$atoms->{$i}{ATMNAME},
		$atoms->{$i}{XCOORD},$atoms->{$i}{YCOORD},$atoms->{$i}{ZCOORD};
    }
    print JAGIN "&\n";
    close JAGIN;
    
    $convertCmd = "$files->{JAGCONVERT} -omae ${prefix}.mae -ijin ${prefix}.in";
    die "ERROR: Cannot execute \"$convertCmd\"\n" if (system("$convertCmd >& /dev/null"));
}

sub findEquivalence {
    my ($equiv, $types, $str) = @_;
    my ($i, $found);

    return if (! exists($equiv->{$str}));
    for $i (keys %{ $equiv->{$str} }) {
	push @{ $found }, $i if (exists($types->{$i}) && $types->{$i} == 0);
    }
    return $found;
}

sub GetMAEFileInfo {
    my ($maeFile) = $_[0];
    my ($atoms, $bonds, $flag, @tmp, $rec, $inp);

    $flag = 0;
    open MAEFILE, $maeFile or die "ERROR: Cannot open $maeFile: $!\n";
    while (<MAEFILE>) {
	chomp;
	next if ($_ =~ /:::/);
	if ($_ =~ /i_m_template_index/) {
	    $flag = 1;
	} elsif ($_ =~ /i_m_to_rep/) {
	    $flag = 2;
	} elsif ($flag == 1 && $_ =~ /^\s*\d+\s+\d+\s+\-?\d+\.\d+/) {
	    $inp = $_;
	    $inp =~ s/\"\s+\"/X/g;
	    $inp =~ s/\"//g;
	    $inp =~ s/^\s+//;
	    $inp =~ s/\s+$//;
	    @tmp = split /\s+/, $inp;
	    $rec = ( {
			"INDEX"     => $tmp[0],
			"XCOORD"    => $tmp[2],
			"YCOORD"    => $tmp[3],
			"ZCOORD"    => $tmp[4],
			"RESNUM"    => $tmp[5],
			"RESNAME"   => $tmp[7],
			"CHAIN"     => $tmp[8],
			"CHARGE"    => $tmp[18],
			"ATMNAME"   => $tmp[17],
			"LONEPAIRS" => 0,
		      });
	    $atoms->{$tmp[0]} = $rec;
	    $atoms->{$tmp[0]}{ATMNAME} = $tmp[16] if ($tmp[17] !~ /[a-zA-Z]/);
	    $bonds->{$tmp[0]} = ();
	} elsif ($flag == 2 && $_ =~ /^\s*\d+\s+(\d+)\s+(\d+)\s+(\d)\s+\d\s+\d/) {
	    push @{ $bonds->{$1} }, $2;
	    push @{ $atoms->{$1}{ORDER} }, $3;
	} else {
	    $flag = 0;
	}
    }
    die "ERROR: No valid data found file parsing MAE file\n" if (! $atoms);

    return ($atoms, $bonds);
}

sub init {
    my (%OPTS);
    getopt('its',\%OPTS);

    die "usage: $0 -i inputfile -t (inputType = pbd|bgf|mol2|mae) -s (saveName)\n" if (! exists($OPTS{i}));

    print "Initializing...";
    ($inputFile, $inputType, $ffName) = ($OPTS{i}, $OPTS{t}, $OPTS{s});
    FileTester($inputFile);

    if (! exists($ENV{SCHRODINGER})) {
	print "setting SCHRODINGER environment variable to /exec/schrodinger_2009...";
	$ENV{SCHRODINGER} = "/exec/schrodinger_2016";
    }
    $FILES->{STRUCTURE} = $inputFile;
    $FILES->{BMIN} = "$ENV{SCHRODINGER}/bmin";
    $FILES->{JAGCONVERT} = "$ENV{SCHRODINGER}/utilities/jagconvert";
    #$FILES->{ATOMTYPES} = "$ENV{SCHRODINGER}/mmshare-v30011/data/atom.typ";
    #$FILES->{ATOMTYPES} = "$ENV{SCHRODINGER}/mmshare-v18111/data/atom.typ";
    $FILES->{ATOMTYPES} = "$ENV{SCHRODINGER}/mmshare-v3.6/data/atom.typ";
    $FILES->{ATOMTYPES} = "/exec/schrodinger_2016/mmshare-v3.6/bin/Linux-x86_64/../../data/atom.typ";
    #$FILES->{ATOMTYPES} = "/project/exec/schrodinger_2013/mmshare-v24012/data/atom.typ";
    $FILES->{OPLS_HEADER} = "/ul/tpascal/ff/OPLSFF_header.ff";
    for (keys %{ $FILES }) {
	die "ERROR: Cannot locate $_ file $FILES->{$_}. Program cannot continue\n"
	    if (! -e $FILES->{$_});
    }
    if (! defined($ffName)) {
	$ffName = basename($inputFile);
	$ffName =~ s/\.\w+$//;
	$ffName .= "_opls.ff";
    }
    $FILES->{PREFIX} = $ffName;
    $FILES->{PREFIX} =~ s/\.\w+$//;
    
    $newBGF = $ffName;
    $newBGF =~ s/\.\w+$//;
    $newBGF .= ".bgf";

    die "ERROR: Cannot determine atom type from $inputFile\n" 
	if($inputFile =~ /\.(\w+)$/ && ! defined($inputType));
    $inputType = $1 if (! defined($inputType));
    die "ERROR: Invalid input type \"$inputType\". Expected pdb|mol2|mae|bgf\n" 
	if ($inputType !~ /^(pdb|mol2|bgf|mae)$/i);
    #open EQUIVFILE, "/project/exec/schrodinger_2013/mmshare-v24012/data/f15/f15_oplsaa.type" or die "ERROR: Cannot open opls type file: $!\n";
    #open EQUIVFILE, "$ENV{SCHRODINGER}/impact-v7.3/data/opls2000/oplsaa.type" or die "ERROR: Cannot open opls type file: $!\n";
    #open EQUIVFILE, "$ENV{SCHRODINGER}/mmshare-v3.6/data/f15/f15_oplsaa.type" or die "ERROR: Cannot open opls type file: $!\n";
    #open EQUIVFILE, "$ENV{SCHRODINGER}/mmshare-v18111/data/f15/oplsaa.type" or die "ERROR: Cannot open opls type file: $!\n";
    #open EQUIVFILE, "/exec/schrodinger_2015/impact-v67011/data/opls2000/oplsaa.type" or die "ERROR: Cannot open opls type file: $!\n";
    open EQUIVFILE, "/exec/schrodinger_2016/mmshare-v3.6/bin/Linux-x86_64/../../data/f15/f15_oplsaa.type" or die "ERROR: Cannot open opls type file: $!\n";
    while (<EQUIVFILE>) {
	chomp;
	if ($_ =~ /^\d+\s+(\S+)\s+(\S+)/) {
	    $EQUIV->{$1}{$2} = 1;
	}
    }
    close EQUIVFILE;
    $EQUIV->{C2}{CZ} = 1;
    $EQUIV->{C2}{CQ} = 1;
    $EQUIV->{O2}{OZ} = 1;
    $EQUIV->{O0}{OY} = 1;
    $EQUIV->{O2}{OS} = 1;
    $EQUIV->{C2}{C5A} = 1;
    $EQUIV->{C2}{C56A} = 1;
    $EQUIV->{C2}{C56B} = 1;
    $EQUIV->{C2}{C5B} = 1;
    $EQUIV->{C2}{C5X} = 1;
    $EQUIV->{Fm}{F} = 1;
    $EQUIV->{C2}{CRA} = 1;
    $EQUIV->{C2}{CRAB} = 1;
    $EQUIV->{C2}{CRP} = 1;
    $EQUIV->{S0}{SY} = 1;
    
    $inputType =uc $inputType;
    if ($inputType =~ /BGF/) {
	$readFunc = \&GetBGFFileInfo;
    } elsif ($inputType =~ /MOL2/) {
        $readFunc = \&GetMOL2FileInfo;
    } elsif ($inputType =~ /PDB/) {
        $readFunc = \&GetPDBFileInfo;
    }

    print "Done\n";
}

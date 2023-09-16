package Packages::FileFormats;

use strict;
use Packages::ManipAtoms qw(GetMols);
use Packages::General qw(Trim);
use File::Basename qw(basename);
require Exporter;

our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(createBGF GetBGFFileInfo addElement addCon addHeader sortByRes createMOL2 GetMSIFileInfo
createHeaders addBoxToHeader GetPDBFileInfo GetMOL2FileInfo GetResAtoms GetBondList createPDB createPQR 
PrintCoord DeleteAtoms AddAtom UpdateBGF GetSystemCharge AddMass GetBGFAtoms GetMass AMBER2MOL2Types);
$VERSION = "1.00";

sub numerically { ($a<=>$b); }

sub GetMSIFileInfo {
   my ($msiName, $saveHeaders) = @_;
   my (%ATOMS, %BONDS, @HEADERS, %AtomMap, $mode, $i, $j);  
   my ($rec, $index, $atom1, $atom2, $counter, %DATA);

   $HEADERS[0] = "BGVER 322";
   $HEADERS[1] = "";
   open MSIFILE, $msiName or die "ERROR: Cannot open $msiName: $!\n";
   while (<MSIFILE>) {
	chomp;
	$_ =~ s/\#//;
	if ($_ =~ /Model/) {
	    $mode = 1; #headers
	} elsif ($_ =~ /^\s*\((\d+) Atom/) {
	    $mode = 2; #atoms
	    if ($rec) {
		$rec->{RESNAME} = "RES";
		$rec->{RESNUM} = 444;
                $rec->{CHARGE} = 0;
                $rec->{LABEL} = "ATOM";
                $rec->{LONEPAIRS} = 0;
		for $i (keys %{ $rec }) {
		    $DATA{$counter}{$i} = $rec->{$i};
		}
	    }
	    $rec = ();
            $counter = $1;
	} elsif ($_ =~ /^\s*\(\d+ Bond/) {
	    $mode = 3; #bonds
            if ($rec) {
                $rec->{RESNAME} = "RES";
                $rec->{RESNUM} = 444;
		$rec->{CHARGE} = 0;
		$rec->{LABEL} = "ATOM";
		$rec->{LONEPAIRS} = 0;
                for $i (keys %{ $rec }) {
                    $DATA{$counter}{$i} = $rec->{$i};
                }
            }
            $rec = ();
            if ($atom1 and $atom2) {
		$DATA{$atom1}{BONDS}{$atom2} = 1;
		$DATA{$atom2}{BONDS}{$atom1} = 1;
                $DATA{$atom1}{NUMBONDS}++;
                $DATA{$atom2}{NUMBONDS}++;
            }
            undef($atom1);
            undef($atom2);

	} elsif ($mode == 1 and $_=~ /Label\s+.*(\w+)/) {
	    $HEADERS[1] = "DESCRP $1";
	} elsif ($mode == 1 and $_ =~ /A C \"Save-comment \d+\" \"(.+)/) {
	    push @HEADERS, "REMARK $1";
	} elsif ($mode == 2 and $_ =~ /FFType \"(.+)\"/) {
	    $rec->{FFTYPE} = $1;
	} elsif ($mode == 2 and $_ =~ /XYZ\s*\((\-?\d+\.?\d*)\s*(\-?\d+\.?\d*)\s*(\-?\d+\.?\d*)/) {
	    $rec->{XCOORD} = $1;
	    $rec->{YCOORD} = $2;
	    $rec->{ZCOORD} = $3;
	} elsif ($mode == 2 and $_ =~ /Id\s*(\d+)/) {
	    $index = $1;
	    $AtomMap{$counter} = $1;
	    $rec->{INDEX} = $1;
	} elsif ($mode == 2 and $_ =~ /Label\s*\"(.+)\"/) {
	    $rec->{ATMNAME} = $1;
	} elsif ($mode == 2 and $_ =~ /PerParent (\d+)/) {
	    $rec = ();
	    $AtomMap{$counter} = $1;
	} elsif ($mode == 3 and $_ =~ /Atom1 (\d+)/) {
	    $atom1 = $1 if (exists($DATA{$1}));
	    $atom1 = $AtomMap{$1} if (! exists($DATA{$1}));
	} elsif ($mode == 3 and $_ =~ /Atom2 (\d+)/) {
	    $atom2 = $1 if (exists($DATA{$1}));
	    $atom2 = $AtomMap{$1} if (! exists($DATA{$1}));
	}
    }
    close MSIFILE;
    die "ERROR: $msiName does not contain any valid information!\n" if (! %DATA);

    if ($rec) {
	$rec->{RESNAME} = "RES";
	$rec->{RESNUM} = 444;
	$rec->{CHARGE} = 0;
	$rec->{LABEL} = "ATOM";
	$rec->{LONEPAIRS} = 0;
	for $i (keys %{ $rec }) {
	    $DATA{$counter}{$i} = $rec->{$i};
	}
    }
    $rec = ();
    if ($atom1 and $atom2) {
	$DATA{$atom1}{BONDS}{$atom2} = 1;
	$DATA{$atom2}{BONDS}{$atom1} = 1;
	$DATA{$atom1}{NUMBONDS}++;
	$DATA{$atom2}{NUMBONDS}++;
    }

    for $i (keys %DATA) {
	$atom1 = $DATA{$i}{INDEX};
	for $j (keys %{ $DATA{$i}{BONDS} }) {
	    $atom2 = $DATA{$j}{INDEX};
	    push @{ $BONDS{$atom1} }, $atom2;
	}
	delete $DATA{$i}{BONDS};
	for $j (keys %{ $DATA{$i} }) {
	    $ATOMS{$atom1}{$j} = $DATA{$i}{$j};
	}
    }

    if (defined($saveHeaders) and $saveHeaders =~ /^1/) {
	return (\%ATOMS, \%BONDS, \@HEADERS);
    } else {
	return (\%ATOMS, \%BONDS);
    }
}
sub createBGF {
    my ($AtmData, $CON, $file_name, $saveRadii) = @_;
    my ($atom, $Atm_no, $bond, $label, $i, $x, $y, $z);

    if (! defined($saveRadii)) {
	$saveRadii = 0;
    }

    open OUTDATA, "> $file_name" or die "Cannot create $file_name: $!\n";
    if (exists($AtmData->{"HEADER"})) {
	for (@{ $AtmData->{"HEADER"} }) {
	    print OUTDATA $_ . "\n";
	}
	delete $AtmData->{"HEADER"};
    }
    print OUTDATA "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,f10.5)\n";

    for $Atm_no (sort numerically keys %{ $AtmData }) {
	$atom = \%{ $AtmData->{$Atm_no} };
	if (! exists($atom->{"LABEL"})) {
	    $label = "ATOM";
	} else {
	    $label = $atom->{"LABEL"};
	}
        if (! exists($atom->{"NUMBONDS"})) {
	    $atom->{"NUMBONDS"} = 0;
        }
	$atom->{"NUMBONDS"} = scalar(@{ $CON->{$Atm_no} }) if (defined($CON->{$Atm_no}));
        #if (length($atom->{"ATMNAME"}) > 4) {
	    #$atom->{"ATMNAME"} = substr($atom->{"ATMNAME"},0,4);
        #}
        for $i ("RESONANCE", "OCCUPANCY", "RADII") {
	    $atom->{$i} = 0 if (! exists($atom->{$i}));
	}
	$atom->{RESNAME} = substr($atom->{RESNAME},0,3) if (length($atom->{RESNAME}) > 3);
        $atom->{CHAIN} = "X" if (! defined($atom->{CHAIN}));
	($x, $y, $z) = ($atom->{"XCOORD"}, $atom->{"YCOORD"}, $atom->{"ZCOORD"});
	$x -= $atom->{OFFSET}{XCOORD} if (exists($atom->{OFFSET}{XCOORD}));
        $y -= $atom->{OFFSET}{YCOORD} if (exists($atom->{OFFSET}{YCOORD}));
        $z -= $atom->{OFFSET}{ZCOORD} if (exists($atom->{OFFSET}{ZCOORD}));

	printf OUTDATA "%-6s %5d %-5s%4s %1s %5d%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f",
	$label, $Atm_no, $atom->{"ATMNAME"}, $atom->{"RESNAME"}, $atom->{CHAIN}, $atom->{"RESNUM"}, 
	$x, $y, $z, $atom->{"FFTYPE"}, $atom->{"NUMBONDS"}, $atom->{"LONEPAIRS"}, $atom->{"CHARGE"}; 
	if ($saveRadii) {
	    printf OUTDATA "%2d%4d%8.3f", $atom->{"RESONANCE"}, $atom->{"OCCUPANCY"}, $atom->{"RADII"};
	}
	print OUTDATA "\n";
    }

    #print OUTDATA "FORMAT CONECT (a6,12i6)\nFORMAT ORDER (a6,i6,13f6.3)\n";
    print OUTDATA "FORMAT CONECT (a6,12i6)\n";
    
    for $atom (sort numerically keys %{ $AtmData }) {
	printf OUTDATA "CONECT%6d", $atom;
	for $bond (@{ $CON->{$atom} }) {
	    printf OUTDATA "%6d", $bond;
	}
	print OUTDATA "\n";
	for ("ORDER", "DISPX", "DISPY", "DISPZ") {
	    if (exists($AtmData->{$atom}{$_})) {
		printf OUTDATA "$_ %6d", $atom;
		for $bond (@{ $AtmData->{$atom}{$_} }) {
		    $bond = 1 if ($bond eq "ar");
		    $bond = 2 if ($bond eq "am");
		    printf OUTDATA "%6d", $bond;
		}
		print OUTDATA "\n";
	    }
	}
    }
    
    print OUTDATA "END\n";
    close OUTDATA;

}

sub GetBGFFileInfo {

    my ($in_file, $saveheaders, $verbose) = ($_[0],$_[1], $_[2]);
    my (%ATOMS, %BONDS, @CON, $counter, $in_data, $atm_patern, %BBOX, @ORDER, $nb);
    my ($calcBox, @HEADERS, $chain, $atomC, $numHetAtms, $fftype, $radii, $molid);
    my (@tmp, $i, $j, $k, $conStr);

    if (defined($saveheaders) and $saveheaders =~ /^1/) {
       $saveheaders = 1;
    } else {
       $saveheaders = 0;
    }

    $numHetAtms = 0;

    $atm_patern = '^ATOM\s*(\d+)\s+(.{5})\S?\s+(\S+)\s\w?\s*(\d+)\s*' . 
	'(\-?\d+\.\d{5})\s*(\-?\d+\.\d{5})\s*(\-?\d+\.\d{5})\s*(\S+)\s+' . 
	'(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s*(\d+\s+\d+\s+\d+\.\d+)?';
    open BGFFILE, $in_file or die "Cannot open BGFfile $in_file: $!\n";
    $molid = 1;
    while (<BGFFILE>) {
	chomp;
	$in_data = $_;
	$in_data =~ s/^HETATM/ATOM  /;
	if ($in_data =~ /^XTLGRF|BIOGRF|DESCRP|REMARK|REMARK|FORCEFIELD|PERIOD|AXES|SGNAME|CELLS/) {
	    push @HEADERS, $in_data;
	} elsif ($in_data =~ /$atm_patern/) {
	    $atomC = Trim($1);
	    $fftype = $8;
	    die "ERROR: Invalid BGF format for atom $atomC. Aborting.\nLine is $_\n" if (! $fftype);
	    $ATOMS{$atomC} = (
			    {
				"INDEX"       => $atomC,
				"ATMNAME"     => $2,
				"RESNAME"     => $3,
				"RESNUM"      => $4,
				"XCOORD"      => $5,
				"YCOORD"      => $6,
				"ZCOORD"      => $7,
				"FFTYPE"      => $fftype,
				"NUMBONDS"    => $9,
				"LONEPAIRS"   => $10,
				"CHARGE"      => $11,
				"MOLECULEID"  => $molid,
                             }
                           );
	    $radii = 0;
	    $BONDS{$atomC} = ();
	    if ($12) {
		($ATOMS{$atomC}{"OCCUPANCY"},
		$ATOMS{$atomC}{"RESONANCE"},
		$ATOMS{$atomC}{"RADII"}) = split /\s+/,$12;
		$radii = $ATOMS{$atomC}{"RADII"};
	    }
            if ($in_data =~ /^ATOM\s+\d+\s+\S+\s+\S+\s+(\w)\s+/) {
                $ATOMS{$atomC}{"CHAIN"} = $1;
            } else {
                $ATOMS{$atomC}{"CHAIN"} = "A";
            }
            if ($_ =~ /^HETATM/) {
                $ATOMS{$atomC}{"LABEL"} = "HETATM";
		$numHetAtms++;
            } else {
                $ATOMS{$atomC}{"LABEL"} = "ATOM";
            }

        } elsif ($in_data =~ /^CRYSTX\s+(.+)/) {
            print "Using Box Information from file..." if (! defined($verbose) || $verbose eq "1");
	    push @HEADERS, $in_data;
        
	} elsif ($in_data =~ /^CONECT(.+)$/) {
	    $i = $1;
	    @CON = ();
	    while ($i =~ /(.{6})/g) {
		$j = $1;
		next if ($j !~ /\S+/);
		$j =~ s/\s+//g;
		push @CON, $j;
	    }
	    #for $i (0 .. $#CON) {
		#$CON[$i] =~ s/\s+//g;
	    #}
	    $atomC = shift @CON;
	    next if (! exists($ATOMS{$atomC}));
	    #@CON = sort  @CON;
	    for $counter (@CON) {
		next if (! exists($ATOMS{$counter}));
	        push @{ $BONDS{$atomC} }, $counter;
	    }
	    $ATOMS{$atomC}{NUMBONDS} = ($#CON + 1);
	} elsif ($_ =~ /^(ORDER|DISP.{1})\s+(.+)$/) {
	    @ORDER = split /\s+/, $2;
	    $atomC = shift @ORDER;
	    $nb = 0;
	    for $counter (@ORDER) {
		push @{ $ATOMS{$atomC}{$1} }, $counter;
		$nb += $counter;
	    }
	    $ATOMS{$atomC}{NUMBONDS} = $nb;
	} elsif ($_ =~ /^TER|ENDMDL/) {
	    $molid++;
	}
    }
    close BGFFILE;

    die "ERROR: $in_file does not contain any ATOM/CONNECTION information!\n"
    	if (! %ATOMS or (! %BONDS and $atomC > $numHetAtms));

    &GetMols(\%ATOMS, \%BONDS);
    if (! $saveheaders) {
	return (\%ATOMS, \%BONDS);
    } else {
        return (\%ATOMS, \%BONDS, \@HEADERS);
    }
}

sub StoreExtrema {
    my ($x, $y, $z, $radii, $par) = @_;
    my ($counter, $tester);

    @{ $tester } = keys %{ $par };
    if ($#{ $tester } == -1) {
	for $counter ("X", "Y", "Z") {
	    $par->{$counter} = (
				{
				    "lo" => 99999.999,
				    "hi" => -99999.999,
				}
				);
	}
    } else {
	$par->{"X"}{"hi"} = ($x + $radii)
	    if (($x + $radii) > $par->{"X"}{"hi"});
	$par->{"X"}{"lo"} = ($x - $radii)
	    if (($x - $radii) < $par->{"X"}{"lo"});
	$par->{"Y"}{"hi"} = ($y + $radii)
	    if (($y + $radii) > $par->{"Y"}{"hi"});
	$par->{"Y"}{"lo"} = ($y - $radii)
	    if (($y - $radii) < $par->{"Y"}{"lo"});
	$par->{"Z"}{"hi"} = ($z + $radii)
	    if (($z + $radii) > $par->{"Z"}{"hi"});
	$par->{"Z"}{"lo"} = ($z - $radii)
	    if (($z - $radii) < $par->{"Z"}{"lo"});
    }
}

sub addHeader {
    my ($ATMS, $HDATA) = @_;

    for (@{ $HDATA }) {
	push @{ $ATMS->{"HEADER"} }, $_;
    }
}

sub addCon {
    my ($ATOMS, $CONNECTIONS) = @_;
    my ($atom, $bond);
    for $atom (keys %{ $CONNECTIONS }) {
	for $bond (@{ $CONNECTIONS->{$atom} }) {
	    $ATOMS->{$atom}{"CONNECT"}{$bond} = "";
	    $ATOMS->{$bond}{"CONNECT"}{$atom} = "";
	}
    }
}
    
sub addElement {
    my ($DATA, $field) = @_;
    my ($element, $hashKey, $eleChar);

    for $hashKey (keys %{ $DATA }) {
	if ($DATA->{$hashKey}{$field} =~ /([A-Za-z]+)/) {
	    $eleChar = $1;
	} else {
	    $eleChar = $DATA->{$hashKey}{$field};
	}
	
	$DATA->{$hashKey}{"ELEMENT"} = substr($eleChar, 0, 1);
	$DATA->{$hashKey}{"INDEX"} = $hashKey;
    }
}

sub createHeaders {
    my ($BBOX, $title) = @_;
    my (@HEADER, %BOX, $muser, $mdate);

    $muser = $ENV{USER};
    $mdate = scalar(localtime);

    # get the box information
    if ($BBOX) {
	%BOX = (
		"X" => {
		    "ANGLE" => 90,
		    "VAL"   => $BBOX->{XCOORD}{hi} - $BBOX->{XCOORD}{lo},
		},
		"Y" => {
		    "ANGLE" => 90,
		    "VAL"   => $BBOX->{YCOORD}{hi} - $BBOX->{YCOORD}{lo},
		},
		"Z" => {
		    "ANGLE" => 90,
		    "VAL"   => $BBOX->{ZCOORD}{hi} - $BBOX->{ZCOORD}{lo},
		},
		);
    }

    $HEADER[0] = "BIOGRF 200";
    $HEADER[1] = "DESCRP $title";
    $HEADER[2] = "REMARK Created by $muser on $mdate"; 
    $HEADER[3] = "FORCEFIELD AMBER";
    if (%BOX) {
	$HEADER[4] = "PERIOD 111";
	$HEADER[5] = "AXES   ZYX";
	$HEADER[6] = "SGNAME P 1                  1    1";
	$HEADER[7] = sprintf("CRYSTX %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f", 
			     $BOX{"X"}{"VAL"}, $BOX{"Y"}{"VAL"}, $BOX{"Z"}{"VAL"},
			     $BOX{"X"}{"ANGLE"}, $BOX{"Y"}{"ANGLE"}, $BOX{"Z"}{"ANGLE"});
	$HEADER[8] = sprintf("CELLS %5d%5d%5d%5d%5d%5d", -1, 1, -1, 1, -1, 1);
    }

    return \@HEADER;
}

sub addBoxToHeader {
    my ($HEADER, $BOX) = @_;
    my ($rec);

    $rec  = "PERIOD 111\nAXES   ZYX\nSGNAME P 1                  1    1\nCRYSTX";
    for ("X", "Y", "Z") {
	$rec .= sprintf("%11.5f", ($BOX->{$_}{"hi"} - $BOX->{$_}{"lo"}));
    }
    $rec .=  "   90.00000   90.00000   90.00000\nCELLS    -1    1   -1    1   -1    1";

    push @{ $HEADER }, $rec;


}

sub GetPDBFileInfo {

    my ($infile, $include_solvent, $res_nm, $atom_no, $hasRadii) = @_;
    my ($inText, %FData, $rec, $isvalid, %BONDS);
    $isvalid = 0;
    open INFILE, $infile or die "Cannot open $infile: $!\n";
    while (<INFILE>) {
	chomp;
	$inText = $_;
	if ($inText =~ /^(ATOM  |HETATM)(.{5})(.{6})\s*(\S+)\s+\w*\s+(\d+)\s+(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d*\.*\d*)\s*(\d*\.*\d*)/) {
	    $res_nm = $4;
	    $atom_no = $2;
	    $rec = (
		    {
			"HEADER"   => $1,
			"LABEL"    => $3,
			"ATMNAME"  => $3,
			"RES_NAME" => $res_nm,
			"RESNAME"  => $res_nm,
			"RES_ID"   => $5,
			"RESNUM"   => $5,
			"XCOORD"   => $6,
			"YCOORD"   => $7,
			"ZCOORD"   => $8,
			"INDEX"    => $2,
			"ATMNUM"   => $2,
		    }
		    );
            $atom_no =~ s/\s+//g;
            next if ($atom_no !~ /^\d+$/);
	    
	    if ($9) {
		$rec->{"CHARGE"} = $9;
		$rec->{"RADII"} = $10;
	    }
	    if ($res_nm =~ /WAT|NA|MG/i) {
                $rec->{"SOLUTE"} = 0;
		if ($include_solvent) {
		    $FData{$atom_no} = $rec;
		}
	    } else {
                $rec->{"SOLUTE"} = 1;
		$FData{$atom_no} = $rec;
	    }

	    $isvalid = 1;
	    if (defined($hasRadii) and $hasRadii eq "1") {
		if ($inText =~ /(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*/g) {
		    $FData{$atom_no}{"CHARGE"} = $1;
		    $FData{$atom_no}{"RADII"} = $2;
		} else {
		    $isvalid = 0;
		}
	    }
	} elsif ($inText =~ /^CONECT\s+(\d+)\s+(\.+)/) {
	     $atom_no = $1;
	     $rec = $2;
	     $BONDS{$atom_no} = ();
	     while ($rec =~ /(\d+)/g) {
		push @{ $BONDS{$atom_no} }, $1;
	     }
	}
    }
    close INFILE;

    die "Invalid PDB file $infile\n"
	if (! $isvalid);

    return (\%FData, \%BONDS);
}

sub GetMOL2FileInfo {
    my ($in_file, $calcBox) = @_;
    my ($isAtm, $isBond, $atm_patern, $bond_patern, $isHeader); 
    my (%ATOMS, $description, $order, %CONN, $tmp);

    $isAtm = $isBond = $isHeader = 0;
    
    $atm_patern = '^\s+(\d+)\s+(\S+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)' . 
	'\s+(\-?\d+\.\d+)\s+(\S+)\s+(\d+)\s+(\S*)\s*(\-?\d+\.\d+)';
    $bond_patern = '^\s+\d+\s+(\d+)\s+(\d+)\s+(\w+)';

    open MOL2FILE, $in_file or die "Cannot open MOL2 file $in_file: $!\n";
    while (<MOL2FILE>) {
	chomp;
	if ($_ =~ /@<TRIPOS>MOLECULE/) {
	    $isHeader = 1;
	    $isAtm = $isBond = 0;
	} elsif ($_ =~ /^@<TRIPOS>ATOM/) {
	    $isAtm = 1;
	    $isBond = $isHeader = 0;
	} elsif ($_ =~ /^@<TRIPOS>BOND/) {
	    $isBond = 1;
	    $isAtm = $isHeader = 0;
	} elsif ($isHeader and $_ =~ /(\S+)/) {
	    $description .= "$1 ";
	} elsif ($_ =~ /$atm_patern/ and $isAtm) {
	    $ATOMS{$1} = (
			  {
			      "ATMNAME"     => $2,
			      "XCOORD"      => $3,
			      "YCOORD"      => $4,
			      "ZCOORD"      => $5,
			      "FFTYPE"      => $6,
			      "LONEPAIRS"   => 0,
			      "CHARGE"      => $9,
			      "RESNAME"     => $8,
			      "RESNUM"      => $7,
			      "INDEX"       => $1,
			  }
			  );
	    $tmp  = $1;
	    $ATOMS{$tmp}{RESNAME} =~ s/\d+$//;
	} elsif ($_ =~ /$bond_patern/ and $isBond) {
	    #if (! IsInteger($3)) {
		#$order = 1;
	    #} else {
		$order = $3;
	    #}
	    push @{ $CONN{$1} }, $2;
	    push @{ $CONN{$2} }, $1;
	    push @{ $ATOMS{$1}{"ORDER"} }, $order;
	    push @{ $ATOMS{$2}{"ORDER"} }, $order;
	
	    $ATOMS{$1}{"NUMBONDS"} = ($#{ $ATOMS{$1}{"CONNECT"} } + 1);
	    $ATOMS{$2}{"NUMBONDS"} = ($#{ $ATOMS{$2}{"CONNECT"} } + 1);
	}
    }
    close MOL2FILE;

    push @{ $ATOMS{"HEADER"} }, "BIOGRF 200";
    push @{ $ATOMS{"HEADER"} }, "DESCRP $description";
    push @{ $ATOMS{"HEADER"} }, "REMARK BGF file created by mol2bgf";
    push @{ $ATOMS{"HEADER"} }, "FORCEFIELD DREIDING";
    return (\%ATOMS, \%CONN);
}

sub createMOL2 {
    my ($ATOMS, $BONDS, $fName, $saveResInfo) = @_;
    my ($atom, $bond, @atmList, $bList, $aList, $resID); 
    my ($resName, $RESLIST, $i, $j, $order, $title);

    $title = substr(basename($fName),0,3);
    open MOL2, "> $fName" || die "ERROR: Cannot create MOL2 file $fName: $!\n";
    $i = $j = 0;
    $saveResInfo = 0 if (! defined($saveResInfo));
    @atmList = sort numerically keys %{ $ATOMS };
    $bList = "";
    for $atom (@atmList) {
	$j++;
	if (! $saveResInfo) {
	    $resID = 1;
	    $resName = "<1>";
	} else {
	    $resID = $ATOMS->{$atom}{RESNUM};
	    $resName = $ATOMS->{$atom}{RESNAME};
	}
	$aList .= sprintf("%7d %-4s %13.4f%10.4f%10.4f %-4s%7d%4s%16.5f\n",
			  $atom, $ATOMS->{$atom}{ATMNAME}, $ATOMS->{$atom}{XCOORD},
			  $ATOMS->{$atom}{YCOORD}, $ATOMS->{$atom}{ZCOORD},
			  $ATOMS->{$atom}{FFTYPE}, $resID, $resName, $ATOMS->{$atom}{CHARGE});
	$RESLIST->{ $ATOMS->{$atom}{RESNAME} } = 1;
	for $bond (0 .. $#{ $BONDS->{$atom} }) {
	    next if ($atom > $BONDS->{$atom}[$bond]);
	    $order = 1;
	    if (exists($ATOMS->{$atom}{ORDER})) {
		$order = $ATOMS->{$atom}{ORDER}[$bond];
	    }
	    $i++;
	    $bList .= sprintf("%6d %4d %4d %-5s\n", $i, $atom, $BONDS->{$atom}[$bond], $order);
	}
    }
    print MOL2 "\@<TRIPOS>MOLECULE\n$title\n";
    print MOL2 "$j $i 0 0 0\n";
    print MOL2 "SMALL\nUSER_CHARGES\n\n\@<TRIPOS>ATOM\n";
    print MOL2 "${aList}@<TRIPOS>BOND\n${bList}";
    close MOL2;    
}

sub IsInteger {
    my ($testStr) = $_[0];
    my ($returnStr);

    $returnStr = 1;
    if  ($testStr !~ /^\d+$/) {
	$returnStr = 0;
    }	

    return $returnStr;
}

sub GetResAtoms {
    my ($atom, $BONDS, $bondList) = @_;
    my ($tmp, $bond, $i);
    
    $bondList->{$atom} = 1;
    for $bond (@{ $BONDS->{$atom} }) {
	next
	    if (exists($bondList->{$bond}));
	$bondList->{$bond} = 1;
	$tmp = GetResAtoms($bond, $BONDS, $bondList);
	for $i (keys %{ $tmp }) {
	    $bondList->{$i} = 1;
	}
    }
    
    return $bondList;
}

sub PrintCoord {
    my ($atom) = $_[0];
    my ($returnStr);
    
    $returnStr = "(";
    
    for ("XCOORD", "YCOORD", "ZCOORD") {
	$returnStr .= sprintf("%.3f, ", $atom->{$_});
    }
    
    $returnStr = substr($returnStr, 0, -2) . ")";
   
    return $returnStr;
}

sub DeleteAtoms {
    my ($atmList, $atoms, $bonds) = @_;
    my ($atomC, $bondC, $i);
    
    for $atomC (keys %{ $atmList }) {
	delete $atoms->{$atomC};
	delete $bonds->{$atomC};
	#for $i (keys %{ $bonds }) {
	#    for $bondC (0 .. $#{ $bonds->{$i} }) {
	#	if ($bonds->{$i}[$bondC] and $bonds->{$i}[$bondC] == $atomC) {
	#	    delete $bonds->{$i}[$bondC];
	#	}
	#    }
	#}
    }
}

sub AddAtom {
    my ($currAtom, $atomIndex, $ATOMS, $BONDS) = @_;
    
    if (exists($ATOMS->{$atomIndex})) {
	$atomIndex = getLastAtom($ATOMS) + 1;
    }
    
    %{ $ATOMS->{$atomIndex} } = %{ $currAtom };
    @{ $BONDS->{$atomIndex} } = (); 
}

sub getLastAtom {
    my ($ATOMS) = $_[0];
    my (@atmList) = sort numerically keys %{ $ATOMS };
    
    return $atmList[$#atmList];
}

sub UpdateBGF {
    my ($ATOMS, $BONDS) = @_;
    my (%atmList, %bondList, $index, $atomC, $bondC);
    
    $index = 1;
    for $atomC (sort numerically keys %{ $ATOMS }) {
	%{ $atmList{$index} } = %{ $ATOMS->{$atomC} };
	$atmList{$index}{"INDEX"} = $atomC;
	$ATOMS->{$atomC}{"INDEX"} = $index;
	$index++;
    }
    
    for $atomC (keys %{ $BONDS }) {
	for $bondC (@{ $BONDS->{$atomC} }) {
            if (exists($ATOMS->{$bondC})) {
                push @{ $bondList{ $ATOMS->{$atomC}{"INDEX"} } }, $ATOMS->{$bondC}{"INDEX"};
            }
	}
    }
    
    $ATOMS = ();
    $BONDS = ();
    
    return (\%atmList, \%bondList);
}

sub GetSystemCharge {
    my ($ATOMS) = $_[0];
    my ($atom, $totCharge);

    $totCharge = 0;
    for $atom (keys %{ $ATOMS }) {
	$totCharge += $ATOMS->{$atom}{"CHARGE"};
    }
    
    return $totCharge;
}


sub AddMass {
    my ($atoms, $parms) = @_;
    my ($i, $fftype, $mass);
    
    for $i (keys %{ $atoms }) {
        $fftype = $atoms->{$i}{"FFTYPE"};
        die "ERROR: Atom $atoms->{$i}{ATMNAME} does not have a force field type\n"
            if (! defined($fftype));
        $mass = $parms->{"ATOMTYPES"}{$fftype}{"MASS"};
        die "ERROR: Atom $atoms->{$i}{ATMNAME} has an invalid force field type $fftype\n"
            if (! defined($mass));
        $atoms->{$i}{"MASS"} = $mass;
    }
}

sub GetMass {
    my ($parms) = $_[0];
    my ($i, $MASSES);

    for $i (keys %{ $parms->{"ATOMTYPES"} }) {
        $MASSES->{$i} = $parms->{"ATOMTYPES"}{$i}{"MASS"};
    }
    return $MASSES;
}

sub GetBGFAtoms {
    my ($atmList, $atoms, $bonds) = @_;
    my ($field, $index, $atomC, %BGF, %CONS, %TMP, $i);

    $atomC = 1;

    for $index (sort numerically keys %{ $atmList }) {
        next if (! exists($atoms->{$index}));
	$BGF{$atomC} = \%{ $atoms->{$index} };
	$atoms->{$index}{"FINDEX"} = $atomC;
	$BGF{$atomC}{"OINDEX"} = $index;
	$atomC++;
    }

    for $atomC (keys %BGF) {
	@{ $CONS{$atomC} } = ();
	for $index (@{ $bonds->{$BGF{$atomC}{"OINDEX"} } }) {
	    $i = $atoms->{$index}{"FINDEX"};
	    if ($i and exists($BGF{$i})) {
		push @{ $CONS{$atomC} }, $i;
	    }
	}
    }

    return (\%BGF, \%CONS, $atmList);

}

sub getAtmList {
    my ($atoms, $selectList, $field) = @_;
    my ($i, $atom, %sorted, %atmList);

    for $i (keys %{ $atoms }) {
	$sorted{ uc($atoms->{$i}{$field}) }{$i} = 1;
    }

    for $i (keys %{ $selectList }) {
	$i = uc($i);
	if (exists($sorted{$i})) {
	    for $atom (keys %{ $sorted{$i} }) {
		$atmList{$atom} = 1;
	    }
	}
    }

    return \%atmList;
}

sub sortByRes {
    my ($ATOMS, $select) = @_;
    my ($i, %RES, $resNum);

    %{ $select } = %{ $ATOMS } if (! defined($select));
    for $i (keys %{ $ATOMS }) {
	next if (! exists($select->{$i}));
	$resNum = $ATOMS->{$i}{RESNUM};
	$RES{$resNum}{ATOMS}{$i} = 1;
    }

    return \%RES;
}

sub GetBondList {
    my ($ATOMS, $BONDS, $SELECT) = @_;
    my ($atomC, $atom, $resAtms, $resC, %RESINFO, $i);

    $resC = 1;
    for $atomC (keys %{ $ATOMS }) {
	$atom = \%{ $ATOMS->{$atomC} };
	next if (keys %{ $SELECT } and ! exists($SELECT->{$atomC}));
	next if (exists($atom->{"bondpointer"}));
	$resAtms = GetResAtoms($atomC, $BONDS, ());
	$atom->{"bondpointer"} = $resC;
	$RESINFO{$resC}{$atomC} = 1;
	for $i (keys %{ $resAtms }) {
	    next if (keys %{ $SELECT } and ! exists($SELECT->{$i}));
	    $ATOMS->{$i}{"bondpointer"} = $resC;
	    $RESINFO{$resC}{$i} = 1;
	}
	$resC++;
    }

    return \%RESINFO;
}

sub createPDB {
    my ($atoms, $bonds, $saveName) = @_;
    my ($atmName, $resName, $i, $bondList, $j);
    my ($fmt) = '%-6s%5d %4s %3s %4d %11.3f%8.3f%8.3f';
    open PDBFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    for $i (sort numerically keys %{ $atoms }) {
        $atmName = $atoms->{$i}{"ATMNAME"};
        $resName = $atoms->{$i}{"RESNAME"};
	if (length($atmName) > 4) {
	    $atmName =~ /(\S+)/;
	    $atmName = $1;
	    $atmName =~ /^(.{4})/;
	    $atmName = $1;
	}
        $resName = substr($resName, 0,2) . substr($resName,-1,1) if (length($resName) > 3);
        printf PDBFILE "$fmt\n", "ATOM", $i, $atmName, $resName, $atoms->{$i}{RESNUM},
        $atoms->{$i}{XCOORD}, $atoms->{$i}{YCOORD}, $atoms->{$i}{ZCOORD};
	$bondList .= sprintf("%-6s%5d", "CONECT", $i);
	for $j (@{ $bonds->{$i} }) {
	    $bondList .= sprintf("%5d",$j);
	}
	$bondList .= "\n";
    }

    print PDBFILE $bondList if (keys %{ $bonds });

    close PDBFILE;
}

sub createPQR {
    my ($atoms, $bonds, $saveName) = @_;
    my ($atmName, $resName, $i, $bondList, $j);
    my ($fmt) = '%-6s%5d %4s %3s %5d %11.3f%8.3f%8.3f%8.3f%8.3f';
    open PQRFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    for $i (sort numerically keys %{ $atoms }) {
	die "ERROR: Radii for atom $i not found.. Aborting\n" if (! exists($atoms->{$i}{RADII}));
        $atmName = $atoms->{$i}{"ATMNAME"};
        $resName = $atoms->{$i}{"RESNAME"};
	$atmName =~ s/\s//g if (length($atmName) > 4);
        if (length($atmName) > 4) {
            $atmName =~ /(.{4})$/;
            $atmName = $1;
        }
        $resName = substr($resName, 0,2) . substr($resName,-1,1) if (length($resName) > 3);
        printf PQRFILE "$fmt\n", "ATOM", $i, $atmName, $resName, $atoms->{$i}{RESNUM},
        $atoms->{$i}{XCOORD}, $atoms->{$i}{YCOORD}, $atoms->{$i}{ZCOORD}, 
	$atoms->{$i}{CHARGE}, $atoms->{$i}{RADII};
        $bondList .= sprintf("%-6s%5d", "CONECT", $i);
        for $j (@{ $bonds->{$i} }) {
            $bondList .= sprintf("%5d",$j);
        }
        $bondList .= "\n";
    }

    print PQRFILE "CONECT\n$bondList" if (keys %{ $bonds });

    close PQRFILE;
}

sub AMBER2MOL2Types {
    my ($atoms) = $_[0];
    my ($i, %ffmap, $ffType, $convertFile);
    
    $convertFile = "/ul/tpascal/scripts/dat/amber2mol2Types.perldata";
    die "ERROR: Cannot read $convertFile: $!\n" 
	if (! -e $convertFile or ! -r $convertFile);
    scalar eval `cat $convertFile`;
    for $i (keys %{ $atoms }) {
	$ffType = uc($atoms->{$i}{FFTYPE});
	die "ERROR: atom $i ($atoms->{$i}{ATMNAME}) has a atomtype $ffType with no match!\n" 
	    if (! exists($ffmap{$ffType}));
	$atoms->{$i}{FFTYPE} = $ffmap{$ffType};
    }
}

1;

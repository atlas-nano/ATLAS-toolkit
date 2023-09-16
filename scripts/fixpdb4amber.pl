#! /usr/bin/perl -w
# fixpdb.pl - fixes the names of the bases in the namot2 pdb file
# to make it compatible with biograf and MSCFF4.2
# e.g changes GUA to G, THY to  T etc.

# Variable Declaration Section

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub GetPDBFileInfo;
sub updateAtoms;
sub createPDB;
sub numerically { ($a<=>$b); }

my ($pdbFile, $saveName);
my ($ATOMS, $BONDS);

$|++;
&init;
print "Parsing regular PDB file $pdbFile...";
($ATOMS, $BONDS) = GetPDBFileInfo($pdbFile);
print "Done\nUpdating atom names...";
&updateAtoms($ATOMS);
print "Create Xleap compatible PDB file $saveName...";
&createPDB($ATOMS, $BONDS, $saveName);
print "Done\n";

sub init {
    my (%OPTS);
    getopt('ps',\%OPTS);
    die "usage: $0 -p pdb file -s [save name]\n" if (! exists($OPTS{p}));
    print "Initializing...";
    ($pdbFile, $saveName) = ($OPTS{p}, $OPTS{s});
    die "ERROR: Cannot access PDB file $pdbFile: $!\n" if (! -e $pdbFile or ! -r $pdbFile or ! -T $pdbFile);

    if (! defined($saveName)) {
	$saveName = basename($pdbFile);
	$saveName =~ s/\.\w+$//;
	$saveName .= "_amber.pdb";
    }
    print "Done\n";
}

sub GetPDBFileInfo {

    my ($infile, $include_solvent, $res_nm, $atom_no, $hasRadii) = @_;
    my ($inText, %FData, $rec, $isvalid, %BONDS);
    $isvalid = 0;
    open INFILE, $infile or die "Cannot open $infile: $!\n";
    while (<INFILE>) {
        chomp;
        $inText = $_;
        if ($_ =~ /^(ATOM  |HETATM)\s*(\d+)(.{6})\s*(\S+)\s+[A-Za-z]?\s*(\d+)\s+(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)/) {
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

sub updateAtoms {
    my ($atoms) = $_[0];
    my ($atmMap, $i, $atmName, $resMap, $resName, $RES, $resID, $count);

    $count = 1;
    $atmMap = (
	    {
		"O5T"  => "  O5* ",
		"H2A*" => " 1H2* ",
		"H2B*" => " 2H2* ",
		"HN6A" => " 1H6  ",
		"HN6B" => " 2H6  ",
		"HN2A" => " 1H2  ",
		"HN2B" => " 2H2  ",
                "HN4A" => " 1H4  ",
                "HN4B" => " 2H4  ",
		"O3T"  => "  O3* ",
	     },
	   );
    $resMap = (
	     {
		"ADE"  => "DA",
		"ade"  => "DA",
                "CYT"  => "DC",
                "cyt"  => "DC",
                "GUA"  => "DG",
                "gua"  => "DG",
                "THY"  => "DT",
                "thy"  => "DT",
                "URA"  => "RU",
                "ura"  => "RU",
		"HSE"  => "HIS",
		"MG2"  => "MG2",
		"Na"   => "Na",
	      }
	    );

    for $i (values %{ $atoms }) {
	$resName = $i->{RESNAME};
	$resID = $i->{RESNUM};
	next if (exists($RES->{"${resName}_${resID}"}) and $i->{ATMNAME} !~ /HE|HB/);
	$RES->{"${resName}_${resID}"} = $resName;
	$RES->{"${resName}_${resID}"} = $resMap->{ $resName } if (exists($resMap->{ $resName }));
	if ($i->{ATMNAME} =~ /HB/) {
	    $RES->{"${resName}_${resID}"} .= "5";
	} elsif ($i->{ATMNAME} =~ /HE/) {
	    $RES->{"${resName}_${resID}"} .= "3";
	}
    }

    for $i (sort numerically keys %{ $atoms }) {
	next if ($atoms->{$i}{ATMNAME} =~ /HB/ or $atoms->{$i}{RESNAME} eq "POM");
	if ($atoms->{$i}{ATMNAME} =~ /HE/) {
	    $count++;
	    next;
	}
	$atoms->{$i}{CHAIN} = chr(64 + $count);
	$atmName = $atoms->{$i}{ATMNAME};
	$atmName =~ s/\s+//g;	
	$resID = $atoms->{$i}{RESNUM};
	$resName = $atoms->{$i}{RESNAME};
	if (exists($atmMap->{$atmName})) {
	    $atoms->{$i}{ATMNAME} = $atmMap->{$atmName};
	} elsif ($atmName !~ /^H/) {
	    $atoms->{$i}{RESNAME} =  $RES->{"${resName}_${resID}"};
	    $atoms->{$i}{ISAMBER} = 1;
	}
    }
}

sub createPDB {
    my ($atoms, $bonds, $saveName) = @_;
    my ($atmName, $resName, $i, $bondList, $j);
    my ($fmt) = '%-6s%5d%-6s%3s %1s %3d %11.3f%8.3f%8.3f';
    open PDBFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    for $i (sort numerically keys %{ $atoms }) {
	next if (! exists($atoms->{$i}{ISAMBER}));
        $atmName = $atoms->{$i}{"ATMNAME"};
        $resName = $atoms->{$i}{"RESNAME"};
        $atmName = substr($atmName, 0,5) . substr($atmName,-1,1) if (length($atmName) > 6);
        $resName = substr($resName, 0,2) . substr($resName,-1,1) if (length($resName) > 3);
        printf PDBFILE "$fmt\n", "ATOM", $i, $atmName, $resName, $atoms->{$i}{CHAIN}, $atoms->{$i}{RESNUM},
        $atoms->{$i}{XCOORD}, $atoms->{$i}{YCOORD}, $atoms->{$i}{ZCOORD};
	next if (! exists($bonds->{$i}));
        $bondList .= sprintf("%-6s%5d", "CONECT", $i);
        for $j (@{ $bonds->{$i} }) {
            $bondList .= sprintf("%5d",$j);
        }
        $bondList .= "\n";
    }

    print PDBFILE "CONECT\n$bondList" if (keys %{ $bonds });

    close PDBFILE;
}


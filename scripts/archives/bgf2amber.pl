#!/usr/bin/perl -w
BEGIN {
    push @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo);
use Packages::AMBER qw(AmberLib);
use File::Basename;
use Packages::General qw(FileTester);
use Packages::CERIUS2 qw(parseCerius2FF);

sub init;
sub getFileData;
sub createLeapFile;
sub getBondList;
sub matchElement;

die "usage: $0 bgf_file cerius2ff [savePrefix]\n"
    if (! @ARGV || $#ARGV < 1);

my ($bgf_file, $cerius2FF, $savePrefix) = @ARGV;
my ($LIBS, $ATOMS, $BONDS, $RES, $PARMS);

print "Initializing...";
&init;
print "Getting BGF Info from $bgf_file...";
($ATOMS, $BONDS) =  GetBGFFileInfo($bgf_file, 0);
print "Done\nParsing CERIUS2 forcefield $cerius2FF...";
$PARMS = parseCerius2FF($cerius2FF);
print "Done\nCreating LEAP and PDB files...";
open LEAPFILE, "> leapadd" or die "Cannot create leapadd: $!\n";
open PDBFILE, "> ${savePrefix}.pdb" or die "Cannot create ${savePrefix}.pdb: $!\n";
$RES = getFileData;
createLeapFile($RES, $LIBS);
close PDBFILE;
close LEAPFILE;

print "Done\n";

sub init {
    FileTester($bgf_file);

    if (! defined($savePrefix)) {
	$savePrefix = basename($bgf_file);
	$savePrefix =~ s/\.bgf$//;
    }
    $LIBS = &AmberLib;
}

sub getFileData {
    my ($atom, %RES, $atmname, $resname, $element); 
    my ($rescounter, $oldRes, $resIndex, $oldResIndex);
    
    $rescounter = 1;
    $oldRes = "";

    for (sort { $a <=> $b } keys %{ $ATOMS }) {
	$atom = \%{ $ATOMS->{$_} };
	if ($oldRes eq "") {
	    $oldRes = $atom->{"RESNAME"};
	    $oldResIndex = $atom->{"RESNUM"};
	}
	$atmname = $atom->{"ATMNAME"};
	$resname = $atom->{"RESNAME"};
	$resIndex = $atom->{"RESNUM"};

	if (length($atmname) > 3) {
	    $atmname = substr($atmname, 0,2) . substr($atmname,-1,1);
	}
	$atom->{"ATMNAME"} = $atmname;
	
	if (length($resname) > 3) {
	    $resname = substr($resname,0,2) . substr($resname, -1,1);
	}
	$atom->{"RESNAME"} = $resname;
	
	if ($resname ne $oldRes) {
	    $rescounter++;
	    $oldRes = $resname;
	} elsif ($resIndex != $oldResIndex) {
	    $rescounter++;
	    $oldResIndex = $resIndex;
	} elsif ($atom->{"FFTYPE"} eq "OW") {
	    $rescounter++;
	}

	$atom->{"RESNUM"} = $rescounter;

	$RES{$resname}{$rescounter}{$atmname} = $_;
	$element = $PARMS->{ATOMTYPES}{$atom->{FFTYPE}}{ATOM};
	next if (lc($element) eq "h");
	printf PDBFILE "%-6s%5d%4s%5s%6d", "ATOM", $_, $atmname, $resname,
	$atom->{"RESNUM"};
	printf PDBFILE " %11.3f%8.3f%8.3f\n", $atom->{"XCOORD"}, $atom->{"YCOORD"},
	$atom->{"ZCOORD"};
    }

    return \%RES;
}

sub createLeapFile {
    my ($resData, $amberLibs) = @_;
    my ($atomName, $atom, $resName, $atmList, $bondList, $element, $i);

    print LEAPFILE "$savePrefix = loadpdb ${savePrefix}.pdb\n";

    for $resName (keys %{ $resData} ) {
	next if exists($amberLibs->{$resName});
	for $i (keys %{ $resData->{$resName} }) {
	    for $atomName (keys %{ $resData->{$resName}{$i} }) {
		$atom = \%{ $ATOMS->{ $resData->{$resName}{$i}{$atomName} } };
		print LEAPFILE "set ${savePrefix}.${i}.${atomName} type $atom->{FFTYPE}\n";
		print LEAPFILE "set ${savePrefix}.${i}.${atomName} charge $atom->{CHARGE}\n";
		$element = $PARMS->{ATOMTYPES}{$atom->{FFTYPE}}{ATOM};
		die "ERROR: Cannot find force field type $atom->{FFTYPE} in cerius2 force field $cerius2FF\n" if (! defined($element));
		print LEAPFILE "set ${savePrefix}.${i}.${atomName} element " . uc($element) . "\n";
		
		$bondList .= getBondList($atom);
	    }
	    print LEAPFILE "\n";
	}
	print LEAPFILE "\n";
    }
    print LEAPFILE "\n$bondList\nsetbox $savePrefix vdw\n";
    print LEAPFILE "saveamberparm $savePrefix ${savePrefix}.top ${savePrefix}.crd\n";
}

sub getBondList {
    my ($atom) = $_[0];
    my ($atomName, $returnStr, $bond, @bondAtms, $line);
    
    $returnStr = "";
    $atomName = $atom->{"ATMNAME"};

    if (! exists($BONDS->{$atom->{INDEX}}) or $#{ $BONDS->{$atom->{INDEX}} } == -1) {
	return "";
    }
    @bondAtms = @{ $BONDS->{$atom->{INDEX}} };
    $line = "bond ${savePrefix}." . $atom->{RESNUM} . ".${atomName}";
    for $bond (@bondAtms) {
	next
	    if ($atom < $bond );
	$returnStr  .= "${line} ${savePrefix}." . $ATOMS->{$bond}{RESNUM} . "." . $ATOMS->{$bond}{"ATMNAME"} . "\n";
    }

    return $returnStr;
}

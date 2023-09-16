#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use MolData;
use Getopt::Std qw(getopt);
use General qw(GetSelections ShowSelectionInfo);

sub init;
sub parseAmberLib;
sub writeData;
sub findSaltBridges;

my ($fileName, $saveName, $fileType, $maxDist);
my ($molStructure, $DATA, $ATOMS, $nD);

$|++;
&init;
print "Gettting data from $fileType file $fileName...";
$molStructure->read($fileName, $fileType);
print "Done\nSelecting atoms/residues...";
$ATOMS = &parseAtomSelection;
print "Done\nFinding salt bridges (max dist: $maxDist)..."; 
($DATA, $nD) = findSaltBridges($ATOMS, $maxDist);
print "Done\nWriting data to $saveName...";
&writeData($DATA, $fileType, $saveName);
print "Done\n";
if (defined($nD)) {
    print "No salt bridge parnters found for the following res ids... consider renaming to CYS:\n";
    for (keys %{ $nD }) {
	print " $_ ";
	print "\n";
    }
}

sub writeData {
    my ($data, $fileType, $saveFile) = @_;
    my ($outData, $i, $res1, $res2);

    if ($fileType eq "pdb") {
	#write xleap data
	$outData = "#for leap\n";
	for $i (@{ $data }) {
	    $res1 = $i->{atom1}->resid;
	    $res2 = $i->{atom2}->resid;
	    $outData .= "crossLink prot.${res1} connect2 prot.${res2} connect2\n";
	}
	open OUTDATA, "> $saveFile" or die "ERROR: Cannot write to $saveFile: $!\n";
	print OUTDATA $outData;
	close OUTDATA;
    } else {
	$molStructure->write($saveFile, $fileType);
    }
}

sub findSaltBridges {
    my ($atoms, $dist) = @_;
    my (@atomList, $i, $j, $currDist, $maxDist, $jVal, $sDATA, $rec, $count, $NONE, $hydrogen);

    @atomList = keys %{ $atoms };
    $i = $j = $count = 0;
    while ($i < $#atomList) {
	$maxDist = 999999;
	undef($jVal);
	for $j (($i + 1) .. $#atomList) {
	    $currDist = $molStructure->dist($atoms->{$atomList[$i]}, $atoms->{$atomList[$j]});
	    if($currDist < $dist and $currDist <= $maxDist) {
		$maxDist = $currDist;
		$jVal = $j;
	    }
	}
	if (defined($jVal)) {
	    $atoms->{$atomList[$i]}->bondlist->create($atoms->{$atomList[$jVal]});
	    $atoms->{$atomList[$jVal]}->bondlist->create($atoms->{$atomList[$i]});
	    #find hydrogens if any and delete them
	    undef($hydrogen);
	    $hydrogen = $molStructure->findAttachedHydrogen($atomList[$i]);
	    $molStructure->deleteAtom($hydrogen) if (defined($hydrogen));
	    $hydrogen = $molStructure->findAttachedHydrogen($atomList[$jVal]);
	    $molStructure->deleteAtom($hydrogen) if (defined($hydrogen));
	    $rec = { "atom1" => $atoms->{$atomList[$i]},  "atom2" => $atoms->{$atomList[$jVal]} };
	    push @{ $sDATA }, $rec;
	    splice @atomList, $jVal, 1;
	    $molStructure->count->{bonds} += 2;
	    $count++;
	} else {
	    $NONE->{  $atoms->{$atomList[$i]}->resid } = 1;
	}
	splice @atomList, $i, 1;
    }
    die "ERROR: none found!\n" if (! defined($sDATA));
    print "${count} found...";
    return ($sDATA, $NONE);
}

sub parseAtomSelection {
    my ($RES, $i, $j, $tmo);
    for $i ("CYX", "CYS") {
	next if (! exists($molStructure->shash->{"resname"}{$i}));
	for $j (keys %{ $molStructure->shash->{"resname"}{$i} }) {
	    $RES->{$j} = $molStructure->shash->{"resname"}{$i}{$j} 
	    if ($molStructure->shash->{"resname"}{$i}{$j}->atmname =~ /SG/i);
	}
    }

    die "ERROR: No CYX or CYS residues found for salt bridges!\n" if (! defined($RES));
    return $RES;
}

sub init {
    my (%OPTS);

    getopt('fstd', \%OPTS);
    die "usage: $0 -f structure file -d (max distance = 2.5A) -t (file type) -s (save name)\n" if (! exists($OPTS{f}));

    print "Initialzing...";
    $molStructure =  MolData->new();
    ($fileName, $fileType, $saveName, $maxDist) = ($OPTS{f}, $OPTS{t}, $OPTS{s}, $OPTS{d});

    $molStructure->testFile($fileName);
    $maxDist = 2.5 if (! defined($maxDist) or $maxDist !~ /^(\d+\.?\d*)$/);
    if (! defined($fileType)) {
	$fileType = "bgf";
	if ($fileName =~ /\.(\w+)$/) {
	    $fileType = lc $1;
	}
    }
    $saveName = $molStructure->getFileName($fileName) if (! defined($saveName));
    if ($fileType eq "pdb" and ! exists($OPTS{s})) {
	$saveName =~ s/_mod\.\w+$//;
	$saveName .= "_saltbridges.leaprc";
    }
    print "Done\n";
}


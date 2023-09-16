#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use File::Basename;
use Packages::General;
use Packages::FileFormats;

sub initialize;
sub parseJaguarOutputFile;
sub updateCharges;

die "usage: $0 bgf_file jag_out_file [save_name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($bgfFile, $jagFile, $saveName) = @ARGV;

initialize;

my ($ATOMS, $CONS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
my ($JAGINFO) = parseJaguarOutputFile($jagFile);
updateCharges(\%{ $ATOMS }, \%{ $JAGINFO });
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $CONS, $saveName);

sub initialize {
    FileTester($bgfFile);
    FileTester($jagFile);
    
    if (! $saveName) {
	$saveName = basename($jagFile);
	$saveName =~ s/\.\w+$//g;
	$saveName .= "_mod.bgf";
    }
}

sub parseJaguarOutputFile {
    my ($inFile) = $_[0];

    my ($inline, $is_charge, @atom_id, @atom_charge); 
    my ($counter, $curr_id, $is_valid, %Jaguar_Data);

    $is_charge = $is_valid = 0;

    open JAGOUT, $inFile or die "Cannot open $inFile: $!\n";
    while (<JAGOUT>) {
	chomp;
	$inline = $_;
	
        if ($inline =~ /Atomic charges from Mulliken population analysis/) {
            $is_charge = 1;
        }elsif ($is_charge) {
#           print "Got Charge\n";
            if ($inline =~ /Atom\s+(.+)$/) {
                @atom_id = split /\s+/, $1;
            }elsif ($inline =~ /Charge\s+(.+)$/) {
                @atom_charge = split /\s+/, $1;
                for $counter (0 .. $#atom_charge) {
                    $curr_id = $atom_id[$counter];
		    $Jaguar_Data{$curr_id}{"CHARGE"} = $atom_charge[$counter];
		    $is_valid = 1;
                }
            }
	}
    }
    close JAGOUT;

    die "Error processing Jaguar output file: $inFile\n"
	if (! $is_valid);

    return \%Jaguar_Data;

}

sub updateCharges {
    my ($aData, $qmData) = @_;
    my ($atom, $atmName);

    for $atom (keys %{ $aData }) {
	$atmName = $aData->{$atom}{"ATMNAME"};
	$atmName =~ s/\d+//g;
	$atmName .= $atom;
	if (! exists($qmData->{$atmName}) ) {
	    print "ERROR: Atom $atom ($atmName) not found in jaguar output file. Ignoring...\n";
	    next;
	}
	$aData->{$atom}{"CHARGE"} = $qmData->{$atmName}{"CHARGE"};
    }
}

sub findAtom {
    my ($atmName, $atmList) = @_;
    my ($currAtm, $i, $returnVal, $j);

    $returnVal = $j = 0;
    for $i (0 .. $#{ $atmList }) {
	if ($ATOMS->{ $atmList->[$i] }{"ATMNAME"} eq $atmName) {
	    $returnVal = $i;
	    $j = $i;
	    last;
	}
    }

    return ($returnVal, $j);
}

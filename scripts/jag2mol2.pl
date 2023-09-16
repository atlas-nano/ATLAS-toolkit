#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use Getopt::Std qw(getopt);
use FileFormats qw(createMOL2 GetMOL2FileInfo);
use File::Basename qw(basename);

sub init;
sub ReadJagOut;
sub updateJagMOL2;
sub createHeader;

my ($jagFile, $FIELDS, $mol2File, $saveFile, $chargeType);
my ($ATOMS, $BONDS, $HEADERS, $JAGDATA, $SCFenergy, $cRMS);

$|++;
&init;
print "Parsing MOL2 file $mol2File...";
($ATOMS, $BONDS) = GetMOL2FileInfo($mol2File, 0);
@{ $HEADERS } = @{ $ATOMS->{HEADER} } if (exists($ATOMS->{HEADER}));
delete $ATOMS->{HEADER};
print "Done\nParsing Jaguar output file $jagFile...";
$JAGDATA = ReadJagOut($jagFile, $chargeType);
print "Done\nUpdating $FIELDS->{LIST}...";
&updateJagMOL2($ATOMS, $JAGDATA, $FIELDS);
print "Done\nCreating mol2 file $saveFile...";
&createMOL2($ATOMS, $BONDS, $saveFile, 1);
print "Done\n";

sub createHeader {
    my ($file, $energy, $rms) = @_;
    my ($name);
    my (@headers) = ("BIOGRF  321", "DESCRP ");
    
    $energy = 0 if (! defined($energy));
    $rms = 0 if (! defined($rms));

    $name = basename($file);
    $name =~ s/\.\w+$//;
    $headers[1] .= "$name";
    $headers[2] = "REMARK Created by $ENV{USER} @ $ENV{HOSTNAME} at " . scalar(localtime);
    $headers[3] = sprintf("REMARK : E = %12.5f hartrees RMS energy = %12.5f hartrees",$energy, $rms);
    $headers[4] = "FORCEFIELD DREIDING";
    
    return \@headers;
}

sub updateJagMOL2 {
    my ($atoms, $jagInfo, $fieldList) = @_;
    my ($i, $j);
    
    for $i (keys %{ $atoms }) {
	die "Atom $atoms->{$i}{ATMNAME} is not present in jaguar file\n"
	    if (! exists($jagInfo->{$i}));
	for $j (keys %{ $fieldList }) {
	    next if ($j eq "LIST");
	    $atoms->{$i}{$j} = $jagInfo->{$i}{$j};
	}
    }
}

sub init {
    my (%OPTS, $fieldList, %myFields);
    getopt('msjfc',\%OPTS);
    for ("m", "j") {
	die "usage: $0 -m mol2 file -j jaguar output -f [fields] -c [charge type = mul or esp] -s [save name]\n" 
	    if (! exists($OPTS{$_}));
    }
    ($mol2File, $jagFile, $fieldList, $saveFile, $chargeType) = ($OPTS{m}, $OPTS{j}, $OPTS{f}, $OPTS{s}, $OPTS{c});
    print "Initializing...";
    FileTester($mol2File);
    FileTester($jagFile);

    if (! defined($saveFile)) {
	$saveFile = basename($mol2File);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_updated.mol2";
    }

    $fieldList = "XCOORD YCOORD ZCOORD CHARGE" if (! defined($fieldList) or $fieldList eq "*");
    $myFields{XCOORD} = 1;
    $myFields{YCOORD} = 1;
    $myFields{ZCOORD} = 1;
    $myFields{CHARGE} = 1;

    while ($fieldList =~ /(\w+)/g) {
	if (exists($myFields{uc($1)})) {
	    $FIELDS->{$1} = 1;
	}
    }
    die "ERROR: No valid field found!\n" if (! defined($FIELDS));
    $FIELDS->{LIST} = $fieldList;
    $chargeType = "mul" if (! defined($chargeType) or $chargeType !~ /^mul|esp|ele/i);
    if ($chargeType =~ /^mul/i) {
	$chargeType = "Mul";
    } else {
	$chargeType = "elec";
    }
    
    print "Done\n";
}

sub ReadJagOut {
    my ($jag_out, $cType) = @_;
    my ($inline, $is_radii, $hash_key, $is_geometry, $is_charge);
    my (@atom_id, @atom_charge, $counter, $is_valid, $curr_id);
    my (%Jaguar_Data, $resCounter, $atomCounter);
    
    $is_radii = $is_geometry = $is_charge =  $is_valid =  $atomCounter = 0;
    $resCounter = 1;
    open JAGOUT, $jag_out or die "Cannot read from $jag_out: $!\n";
    while (<JAGOUT>) {
        chomp;
        $inline = $_;
        if ($inline =~ /(atom               x                 y                 z)|(final geometry)/) {
            $is_geometry = 1;
	    $atomCounter = 0;
        } elsif ($inline =~ /principal moments of inertia/) {
            $is_geometry = 0;
        }elsif ($is_geometry) {
            if ($inline =~ /([A-Z]{1}[a-z]?)(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
		$atomCounter++;
                $hash_key = $atomCounter;

                $Jaguar_Data{$hash_key}{"XCOORD"} = $3;
                $Jaguar_Data{$hash_key}{"YCOORD"} = $4;
                $Jaguar_Data{$hash_key}{"ZCOORD"} = $5;
                $Jaguar_Data{$hash_key}{"ATMNAME"} = "${1}${2}";
                $is_valid = 1;
            }
        } elsif ($inline =~ /Atomic charges from $cType/) {
            $is_charge = 1;
	    $atomCounter = 0;
        } elsif ($inline =~ /sum of atomic charges/) {
            $is_charge = 0;
        }elsif ($is_charge && $is_valid) {
            if ($inline =~ /Atom\s+(.+)$/) {
                @atom_id = split /\s+/, $1;
            }elsif ($inline =~ /Charge\s+(.+)$/) {
                @atom_charge = split /\s+/, $1;
                for $counter (0 .. $#atom_charge) {
		    $atomCounter++;
                    $atom_id[$counter] =~ /(\d+)/;
                    $curr_id = $1;

                    if ($Jaguar_Data{$atomCounter}) {
#                       print "Wrote Charge\n";
                        $Jaguar_Data{$atomCounter}{"CHARGE"} = $atom_charge[$counter];
                        $is_valid = 1;
                    } else {
                        print "Charge found for Non Existant Atom: $atom_id[$counter]\n";
                    }
                }
            }
        } elsif ($inline =~ /van der Waals radii/) {
            $is_radii = 1;
            $counter = 0;
        } elsif ($inline =~ /Number of Lewis structures found/i) {
            $is_radii = 0;
        } elsif ($is_radii and $inline =~ /^\s*[A-Z]{1}[a-z]?(\d+)\s+(\d+\.\d+)/) {
            $Jaguar_Data{$1}{"RADII"} = $2;
        } elsif ($inline =~ /SCF energy: .* (\-?\d+\.\d+) hartrees/) {
	    $SCFenergy = $1;
	} elsif ($inline =~ /RMS Error\s+(\-?\d+\.\d+\E?\-?\d*)/) {
	    $cRMS = $1;
	}
    }
    close JAGOUT;

    die "Error reading Jaguar Output $jag_out\n" if (! $is_valid);

    return (\%Jaguar_Data);
}

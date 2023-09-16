#!/usr/bin/perl -w
BEGIN {
    push @INC, "/ul/tpascal/scripts/";
}

use strict;
use Packages::General;
use Packages::FileFormats;

sub Initialize;
sub TransMol;
sub Numerically;

die "usage: $0 bgf_file num_seperations seperation_width\n"
    if (! @ARGV or $#ARGV < 2);

my ($bgf_file, $tot, $width) = @ARGV;
Initialize;

my ($MOL1, $CON1, $HEADERS) = GetBGFFileInfo($bgf_file, 1);
my ($MOL2, $CON2) = GetBGFFileInfo($bgf_file, 0);
my ($i, $TOT, $CON, $dist, $CP1, $CP2);

$bgf_file =~ s/\.bgf$//;

for $i (1 .. $tot) {
    $dist = $i * $width/2;
    $CP1 = TransMol($MOL1, ($dist * -1));
    $CP2 = TransMol($MOL2, $dist);

    ($TOT, $CON) = CombineMols($CP1, $CP2, $CON1);
    addHeader($TOT, $HEADERS);
    createBGF($TOT, $CON, $bgf_file . "_" . $i . ".bgf");
}

sub Initialize {

    FileTester($bgf_file);
    die "Error: Expected Integer for num_seperations. Got $tot\n"
	if (! IsInteger($tot));
    die "Error: Expected Decimal for seperation_width. Got $width\n"
	if (! IsDecimal($width));

}

sub TransMol {
    # Transalate the COM of the molecule by distance from origin
    my ($molData, $distance) = @_;
    my ($atom, $dimDisplace, %ATOMS, $sign);

    if ($distance < 0) {
	$sign = -1;
    } else {
	$sign = 1;
    }

    $dimDisplace = $sign * sqrt(($distance ** 2)/3);
    for $atom (keys %{ $molData }) {
	%{ $ATOMS{$atom} } = %{ $molData->{$atom} };
	
	for ("XCOORD", "YCOORD", "ZCOORD") {
	    $ATOMS{$atom}{$_} += $dimDisplace;
	}
    }

    return \%ATOMS;
}    

sub WriteBGFFile(@) {
    my ($file_nm, $ATM_DATA, $CONNS, $HEADERS, $index) = @_;
    my ($counter, $fmt, $fmt1, $fmt2, $fmt3, $atom, $bond);
    
    $file_nm .= "_$index";

    $fmt = "%-6s %5d %-5s %3s %1s %5d%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f%2d%4d%10.5f\n";
    $fmt1 = "%-6s%5d  %4s %3s  %4d    %8.3f%8.3f%8.3f %5.1f %5.1f\n";
    $fmt2 = "%-6s%5d %4s %3s %1s%4d %11.3f%8.3f%8.3f %5.1f %5.1f\n";
    $fmt3 = "%-6s%5d %4s %3s %1s%4d %11.3f%8.3f%8.3f\n";
    open OUTFILE, "> $file_nm.bgf" or die "Cannot write to $file_nm: $!\n";

    for $counter (@{ $HEADERS }) {
	print OUTFILE "$counter\n";
    }

    print OUTFILE "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,i2,i4,f10.5)\n";
    for $counter (sort Numerically keys %{ $ATM_DATA }) {
	$atom = \%{ $ATM_DATA->{$counter} };
	
	printf OUTFILE $fmt, 
	$atom->{"LABEL"}, $counter, $atom->{"ATMNAME"}, $atom->{"RESNAME"}, 
	$atom->{"CHAIN"}, $atom->{"RESNUM"}, $atom->{"XCOORD"}, $atom->{"YCOORD"}, 
	$atom->{"ZCOORD"}, $atom->{"FFTYPE"}, $atom->{"NUMBONDS"}, $atom->{"LONEPAIRS"}, 
	$atom->{"CHARGE"}, $atom->{"OCCUPANCY"}, $atom->{"RESONANCE"}, $atom->{"RADII"};
    }

    print OUTFILE "FORMAT CONECT (a6,14i6)\n";

    for $bond (sort Numerically keys %{ $CONNS }) {
	printf OUTFILE "CONECT%5d", $bond;
	for (@{ $CONNS->{$bond} }) {
	    printf OUTFILE "%6d", $_;
	}
	print OUTFILE "\n";
    }

    close OUTFILE;

    open OUTFILE, "> $file_nm.pqr" or die "Cannot write to PQR file test.pqr: $!\n";
    for $counter (sort Numerically keys %{ $ATM_DATA }) {
	printf OUTFILE $fmt1, "ATOM", $counter, $ATM_DATA->{$counter}{"ATMNAME"},
	$ATM_DATA->{$counter}{"RESNAME"}, $ATM_DATA->{$counter}{"RESNUM"},
	$ATM_DATA->{$counter}{"XCOORD"}, $ATM_DATA->{$counter}{"YCOORD"}, 
	$ATM_DATA->{$counter}{"ZCOORD"}, $ATM_DATA->{$counter}{"CHARGE"}, 
	$ATM_DATA->{$counter}{"RADII"};
    }

    for $bond (sort Numerically keys %{ $CONNS }) {
	printf OUTFILE "CONECT%5d", $bond;
	for (@{ $CONNS->{$bond} }) {
	    printf OUTFILE "%6d", $_;
	}
	print OUTFILE "\n";
    }
    close OUTFILE;
}

sub Numerically {
    ($a<=>$b);
}

#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General;
use Packages::FileFormats;

sub numerically;
sub main;

die "usage: $0 bgfFile [include_solvent=false]\n"
    if (! @ARGV);


my ($bgfFile, $incSolv) = @ARGV;
my (@tmp);

&main;

sub main {
    my ($i, @RES, $resID, $oldRes);
    FileTester($bgfFile);
    $incSolv = 0 if (! defined($incSolv) or Trim($incSolv) !~ /^1$/);
    
    print "Parsing BGF file $bgfFile...";
    my ($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
    print "Done\n";

    @tmp = sort numerically keys %{ $ATOMS };
    $oldRes = $ATOMS->{ $tmp[0] }{"RESNUM"};
    push @RES, $ATOMS->{ $tmp[0] }{"RESNAME"}; 

    for $i (@tmp) {
	$resID = $ATOMS->{$i}{"RESNUM"};
	if ($oldRes != $resID) {
	    $oldRes = $resID;
	    if ($incSolv or $ATOMS->{$i}{"RESNAME"} !~ /WAT|Na|Cl/i) {
		push @RES, $ATOMS->{$i}{"RESNAME"};
	    }
	}
    }

    for $i (@RES) {
	print "$i\n";
    }
}

sub numerically {
    ($a<=>$b);
}

#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester GetSelections);
use ManipAtoms qw(GetAtmList);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub getPBbonds;
sub createPBBonds;
sub numerically { ($a<=>$b); }

my ($ATOMS, $BONDS, $HEADERS, $SELECT, $PBbondList, $bgfFile, $saveFile, $DIMS);
my ($SELECTIONS);

$|++;
&init;
if (defined($SELECT)) {
    print "Parsing atom/residue selection...";
    $SELECTIONS = GetSelections($SELECT, 0);
    $SELECTIONS = GetAtmList($SELECTIONS, $ATOMS);
    die "ERROR: No valid atoms selected!\n" if (! keys %{ $SELECTIONS });
    print "Done\n"; 
}
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nGetting atom selection...";
$PBbondList = getPBbonds($ATOMS, $SELECT, $DIMS);
print "Done\nCreating Periodic Boundary bonds...";
createPBbonds($BONDS, $PBbondList);
print "Done\nCreate BGF file $saveFile...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub createPBbonds {
    my ($bonds, $bondList) = @_;
    my ($i, $atom1, $atom2);

    for $i (@{ $bondList }) {
	$atom1 = $i->{1};
	$atom2 = $i->{2};
	push @{ $bonds->{$atom1} }, $atom2;
	@{ $bonds->{$atom1} } = sort numerically @{ $bonds->{$atom1 } };
	push @{ $bonds->{$atom2} }, $atom1;
	@{ $bonds->{$atom2} } = sort numerically @{ $bonds->{$atom2 } };
    }
}

sub getPBbonds {
    my ($atoms, $aSelect, $dimOpt) = @_;
    my ($i, $j, $k, $l, $pointer, $rec, $tmp, $cMAP);
    my ($EXTREMA, @dlist, @BONDLIST, @maxList, @minList);

    @dlist = keys %{$dimOpt};
    for $i (keys %{ $atoms }) {
	next if(defined($aSelect) and ! exist($aSelect->{$i}));
	$pointer = \%{ $tmp };
 	for $j (@dlist) {
	    $pointer = \%{ $pointer->{ $atoms->{$i}{$j} } };
	    $EXTREMA->{$j}{min} = $atoms->{$i}{$j} if (! exists($EXTREMA->{$j}) or $EXTREMA->{$j}{min} > $atoms->{$i}{$j});
            $EXTREMA->{$j}{max} = $atoms->{$i}{$j} if (! exists($EXTREMA->{$j}{max}) or $EXTREMA->{$j}{max} < $atoms->{$i}{$j});
	}
	$pointer->{atom} = $i;
    }

    if ($#dlist == 1) { #2d so pad 1st entry with a 1
	unshift @dlist, 1;
	$cMAP->{1} = $tmp;
	$EXTREMA->{1}->{max} = 1;
	$EXTREMA->{1}->{min} = 1;
    } else { $cMAP = $tmp; }

    @maxList = keys %{ $cMAP->{ $EXTREMA->{ $dlist[0] }{max} }{ $EXTREMA->{ $dlist[1] }{max} } };
    $pointer = $cMAP->{ $EXTREMA->{ $dlist[0] }{min} }{ $EXTREMA->{ $dlist[1] }{min} };

    for $i (@maxList) {
	next if (! exists($pointer->{$i}));
	$rec = (
		{
		    "1" => $cMAP->{ $EXTREMA->{ $dlist[0] }{max} }{ $EXTREMA->{ $dlist[1] }{max} }{$i}{atom},
		    "2" => $pointer->{$i}{atom},
		}
		);
	push @BONDLIST, $rec;
    }
    die "ERROR: No valid atoms found in atom selection!\n" if (! @BONDLIST);
    return \@BONDLIST;
}

sub init {
    my (%OPTS, $atomSelect, $i, $tmp, $rec1, $rec2, $hash, $dimOpt);

    getopt('blds',\%OPTS);
    die "usage: $0 -b bgf file -l (atom selection) -d (dimension=xy plane) -s [save name]\n"
        if (! defined($OPTS{b}));
    print "Initializing...";
    ($bgfFile, $saveFile, $atomSelect, $dimOpt) = ($OPTS{b}, $OPTS{s}, $OPTS{l}, $OPTS{d});
    FileTester($bgfFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+//;
	$saveFile .= "_mod.bgf";
    }

    $dimOpt = "xy" if (! defined($dimOpt));
    $dimOpt = "xy" if ($dimOpt !~ /xy|yz|yz|zy|xz|zx/i);
    while ($dimOpt =~ /(x|y|z)/g) {
	$DIMS->{uc $1 . "COORD"} = 1;
    }
    if (defined($atomSelect) and $atomSelect =~ /\s+/) {
        @{ $SELECT } = split /\s+/, $atomSelect;
    } elsif (defined($atomSelect)) {
        $SELECT->[0] = $atomSelect;
    }
  
    print "Done\n";
}

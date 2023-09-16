#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF);
use Packages::General qw(FileTester GetSelections);
use Packages::ManipAtoms qw(GetAtmList);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub getPBbonds;
sub createPBBonds;
sub numerically { ($a<=>$b); }

my ($ATOMS, $BONDS, $HEADERS, $SELECT, $PBbondList, $bgfFile, $saveFile);

$|++;
$SELECT = &init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nGetting atom selection...";
$PBbondList = getPBbonds($ATOMS, $SELECT);
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
    my ($atoms, $bondSelect) = @_;
    my ($atoms1, $atoms2, $i, $tmp1, $tmp2, $j, @BONDLIST, $rec);

    for $i (@{ $bondSelect }) {
	$atoms1 = GetAtmList($i->{1}, $atoms);
	$atoms2 = GetAtmList($i->{2}, $atoms);
	if (scalar(keys %{ $atoms1}) == scalar(keys %{ $atoms2 })) {
	    @{ $tmp1 } = sort numerically keys(%{ $atoms1 });
	    @{ $tmp2 } = sort numerically keys(%{ $atoms2 });
	    for $j (0 .. $#{ $tmp1 }) {
		$rec = (
			{
			    "1" => $tmp1->[$j],
			    "2" => $tmp2->[$j],
			}
			);
		push @BONDLIST, $rec;
	    }
	}
    }
    die "ERROR: No valid atoms found in atom selection!\n" if (! @BONDLIST);
    return \@BONDLIST;
}

sub init {
    my ($selection, %OPTS, $atomSelect, $i, $tmp, @BONDLIST, $rec1, $rec2, $hash);

    getopt('bls',\%OPTS);
    for ("b", "l") {
	die "usage: $0 -b bgf file -l atom bond listing -s [save name]\n"
	    if (! defined($OPTS{$_}));
    }

    print "Initializing...";
    ($bgfFile, $saveFile, $atomSelect) = ($OPTS{b}, $OPTS{s}, $OPTS{l});
    FileTester($bgfFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+//;
	$saveFile .= "_mod.bgf";
    }

    if ($atomSelect =~ /\s+/) {
        @{ $selection } = split /\s+/, $atomSelect;
    } else {
        $selection->[0] = $atomSelect;
    }

    for $i (@{ $selection }) {
	@{ $tmp } = ();
	if ($i =~ /(\S+)\->(\S+)/) {
	    $tmp->[0] = $1;
	    $rec1 = GetSelections($tmp, 0);
	    $tmp->[0] = $2;
	    $rec2 = GetSelections($tmp, 0);
	    if (defined($rec1) and defined($rec2)) {
		$hash = (
			 {
			     1 => $rec1,
			     2 => $rec2,
			 }
			 );
		push @BONDLIST, $hash;
	    }
	}
    }

    die "ERROR: Invalid atom selection. Expected atom|range->atom|range. Got $atomSelect\n"
	if (! @BONDLIST);
    print "Done\n";
    return \@BONDLIST;
}

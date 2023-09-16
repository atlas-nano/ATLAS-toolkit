#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use Packages::General;
use Packages::FileFormats;
use strict;

sub numerically;
sub getCon;

die "usage: $0 bgf_file start_atom [save_name]\n"
    if (!@ARGV or @ARGV < 1);

my ($bgfFile, $start, $saveName) = @ARGV;

FileTester($bgfFile);
$start = Trim($start);

die "ERROR: Expected integer when reading start. Got $start\n"
    if (! IsInteger($start));

my ($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 0, 1, 1);

my (@tmp) = sort numerically keys %{ $ATOMS };
my ($atom, $MOL1, $MOL2, $CON1, $CON2);

for $atom (@tmp) {
    if ($atom < $start) {
	%{ $MOL1->{$atom} } = %{ $ATOMS->{$atom} };
	@{ $CON1->{$atom} } = getCon($BONDS->{$atom}, 1, $atom);
    } else {
	%{ $MOL2->{($atom - $start + 1)} } = %{ $ATOMS->{$atom} };
	@{ $CON2->{($atom - $start + 1)} } = getCon($BONDS->{$atom}, 2, $atom);
    }
}

if (! defined($saveName)) {
    $saveName = "out.bgf";
}

$saveName =~ s/\.\w+$//;

addHeader($MOL1, $HEADERS);
addHeader($MOL2, $HEADERS);
createBGF($MOL1, $CON1, $saveName . "_mol1.bgf");
createBGF($MOL2, $CON2, $saveName . "_mol2.bgf");


sub numerically {
    ($a<=>$b);
}

sub getCon {
    my ($bondList, $molid, $atomID) = @_;
    my ($bond, @BONDS);

    for $bond (@{ $bondList }) {
	if ($molid == 1 and $bond < $start) {
	    push @BONDS, $bond;
	} elsif ($molid == 2 and $bond > ($start - 1)) {
	    push @BONDS, ($bond - $start + 1);
	}
    }

    return @BONDS;
}

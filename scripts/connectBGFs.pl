#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester CombineMols);
use FileFormats qw(GetBGFFileInfo addHeader createBGF 
			     createHeaders GetBGFAtoms DeleteAtoms);

my ($FILES, $isFuse, $saveName);
my ($ATOMS, $BONDS, $HEADERS);

sub init;
sub showUsage;
sub loadBGFS;
sub connectSystems;
sub numerically { ($a<=>$b); }

$|++;
&init;
&loadBGFS($FILES);
print "Connecting systems...";
($ATOMS, $BONDS) = &connectSystems($FILES, $isFuse);
print "Done\nCreating BGF file $saveName...";
$HEADERS = createHeaders(undef, $saveName);
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub connectSystems {
    my ($bgfData, $isRemove) = @_;
    my ($BGFATOMS, $BGFBONDS, $select1, $select2, @tmp, $olap);
    my ($atoms1, $bonds1, $atoms2, $bonds2, $i1, $i2, $offset);

    @tmp = sort numerically keys %{ $bgfData->{2}{ATOMS} };
    if ($isRemove) {
	$select1->{$_} = 1 for (1 .. $bgfData->{1}{CONNECT});
	$select2->{$_} = 1 for ($bgfData->{2}{CONNECT} .. $tmp[$#tmp]);
    } else {
	$select1 = $bgfData->{1}{ATOMS};
	$select2 = $bgfData->{2}{ATOMS};
    }

    for $i1 ("XCOORD", "YCOORD", "ZCOORD") { 
	$offset = $bgfData->{1}{ATOMS}{ $bgfData->{1}{CONNECT} }{$i1} - 
	    $bgfData->{2}{ATOMS}{ $bgfData->{2}{CONNECT} }{$i1};
	for $i2 (keys %{ $bgfData->{1}{ATOMS} }) {
	    $bgfData->{1}{ATOMS}{$i2}{$i1} -= $offset;
	}
    }

    $olap->{$bgfData->{1}{CONNECT}} = ();
    DeleteAtoms($olap, $bgfData->{1}{ATOMS}, $bgfData->{1}{BONDS});
    ($atoms1, $bonds1, undef) = GetBGFAtoms($select1, $bgfData->{1}{ATOMS}, $bgfData->{1}{BONDS});
    ($atoms2, $bonds2, undef) = GetBGFAtoms($select2, $bgfData->{2}{ATOMS}, $bgfData->{2}{BONDS});
    ($BGFATOMS, $BGFBONDS) = CombineMols($atoms1, $atoms2, $bonds1, $bonds2);
    $i1 = $bgfData->{1}{ATOMS}{ $bgfData->{1}{CONNPARENT} }{INDEX};
    $i2 = $bgfData->{2}{ATOMS}{ $bgfData->{2}{CONNECT} }{INDEX};
    push @{ $BGFBONDS->{$i1} }, $i2;
    push @{ $BGFBONDS->{$i1} }, $i1;
    return ($BGFATOMS, $BGFBONDS);
}

sub loadBGFS {
    my ($BGFDATA) = $_[0];
    my ($i, $currBGF, @tmp, $j, $parent);

    for $i (1, 2) {
	$currBGF = $BGFDATA->{$i};
	print "Parsing bgffile $currBGF->{NAME}...";
	($currBGF->{ATOMS}, $currBGF->{BONDS}, undef) = GetBGFFileInfo($currBGF->{NAME});
	die "ERROR: Cannot find connect atom $currBGF->{CONNECT}!\n"
	    if (! exists($currBGF->{ATOMS}{ $currBGF->{CONNECT} }));
	@tmp = @{ $currBGF->{BONDS}{ $currBGF->{CONNECT} } };
	$parent = $currBGF->{CONNECT};
	for $j (@tmp) {
	    if ($i == 1 and $j < $currBGF->{CONNECT}) {
		$parent = $j;
		last;
	    }
	}
	$currBGF->{CONNPARENT} = $parent;
	print "Done\n";
    }
}

sub init {
    my (%OPTS, $usage, $conInfo);
    
    $usage = &showUsage;
    getopt('abfsc',\%OPTS);
    for ("a", "b", "c") {
	die "$usage" if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    FileTester($OPTS{a});
    FileTester($OPTS{b});
    ($FILES->{1}{NAME}, $FILES->{2}{NAME}, $conInfo, $saveName, $isFuse) = 
	($OPTS{a}, $OPTS{b}, $OPTS{c}, $OPTS{s}, $OPTS{f});
    if (! defined($saveName)) {
	$saveName = basename($FILES->{1}{NAME});
	$saveName =~ s/\.\w+$//;
	$saveName .= "_" . basename($FILES->{2}{NAME});
	$saveName =~ s/\.\w+$//;
	$saveName .= ".bgf";
    }
    if ($conInfo =~ /(\d+)\s*\->\s*(\d+)/) {
	$FILES->{1}{CONNECT} = $1;
	$FILES->{2}{CONNECT} = $2;
    } else {
	die "ERROR: Expected atom -> atom for connect info. Got \"$conInfo\"!\n";
    }
    $isFuse = 0 if (! defined($isFuse));
    $isFuse = 1 if ($isFuse =~ /^(1|yes)/i);
    print "Done\n";
}

sub showUsage {
    return "usage: $0 -a file1 -b file2 -c conn info -f [is fused] -s [savename]\n" .
	"Options\n" .
	"\t-a: file1. Location of bgf file 1\n\t-b: file2: Location of bgf file 2\n" .
	"\t-c: connection info: Atom numbers to connect. expected: \"atom1 -> atom2\"\n" .
	"\t-f: fused. Specifies whether atoms after file1 and before file2 connect atoms are removed.\n" .
	"\t\tDefault: no\n" .
	"\t-s: save name. The name of the connected(fused) files. Default: file1_file2.bgf\n";
}

#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF createHeaders addBoxToHeader);
use Packages::General qw(FileTester CombineMols CoM);
use Packages::BOX qw(GetBox);
use Packages::ManipAtoms qw(GetMols CenterSystem);

sub init;
sub replicateCell;
sub transMol;
sub getPBCbonds;
sub createPBCbonds;
sub deleteBond;

my ($replicate, $bgfFile, $saveFile, $doCenter, $bondImages);
my ($CENTER, $ATOMS, $BONDS, $HEADERS, $BOX);
$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$BOX = GetBox($ATOMS, undef, $HEADERS);
print "Done\nReplicating cell $replicate->{STRING}...";
($ATOMS, $BONDS, $BOX) = replicateCell($ATOMS, $BONDS, $BOX, $replicate, $doCenter, $bondImages);
&GetMols($ATOMS, $BONDS);
print "Done\nCreating BGF file $saveFile...";
$HEADERS = createHeaders(undef, $saveFile);
$CENTER = CoM($ATOMS);
for (keys %{ $BOX }) {
    $CENTER->{$_ . "COORD"} -= (($BOX->{$_}{hi} + $BOX->{$_}{lo})/2 + $BOX->{$_}{lo} - $CENTER->{$_ . "COORD"});
}
#&CenterSystem($ATOMS, $CENTER, ());
&addBoxToHeader($HEADERS, $BOX);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub transmol {
    my ($unit, $box, $dim, $disp) = @_;
    my ($i,%ATOMS, $offset);
    
    $offset = $box->{$dim}{len} * $disp;
    for $i (keys %{ $unit } ) {
	%{ $ATOMS{$i} } = %{ $unit->{$i} };
    }
    for $i (keys %{ $unit }) {
	$ATOMS{$i}{$dim . "COORD"} += $offset;
    }
    return \%ATOMS;
}

sub replicateCell {
    my ($atoms, $bonds, $box, $dims, $centerMol, $createNewbonds) = @_;
    my ($i, $j, $cellAtoms, $cellBonds, $unitAtoms, $pbcbonds, $offset);
    
    if ($centerMol) {
	for $i ("X", "Y", "Z") {
	    $j = -1 * int(($dims->{$i} -1)/2);
	    for ($j .. -1) {
		$atoms = transmol($atoms, $box, $i, $_);
	    }
	}
    }

    for $i ("X", "Y", "Z") {
	$pbcbonds = getPBCbonds($atoms, $bonds, $box) if ($createNewbonds);
	$unitAtoms = ();
	$cellBonds = ();
	$offset = 0;
	for $j (keys %{ $atoms }) {
	    %{ $unitAtoms->{$j} } = %{ $atoms->{$j} };
	    @{ $cellBonds->{$j} } = @{ $bonds->{$j} };
	    $offset = 1;
	}
	for $j (1 .. ($dims->{$i} - 1)) {
	    $cellAtoms = transmol($unitAtoms, $box, $i, $j);
	    ($atoms, $bonds) = CombineMols($atoms, $cellAtoms, $bonds, $cellBonds);
	    createPBCbonds($atoms, $bonds, $pbcbonds->{"${i}COORD"}, $offset) if(exists($pbcbonds->{"${i}COORD"}));
	}
	$box->{$i}{hi} = ($box->{$i}{len} * $dims->{$i});
	$box->{$i}{lo} = 0;
	$box->{$i}{len} = ($box->{$i}{len} * $dims->{$i});
    }
    return ($atoms, $bonds, $box);
}

sub createPBCbonds {
    my ($atoms, $bonds, $bondlist, $atomOffset) = @_;
    my ($i, $j, $atom1, $atom2, $atom3, $atom4, @list);

    for $i (keys %{ $bondlist }) {
	@list = keys %{ $bondlist->{$i} };
	for $j (0 .. $#list) {
	    $atom1 = $i;
	    $atom2 = $list[$j];
	    $atom3 = $atom1 + $atomOffset;
	    $atom4 = $atom2 + $atomOffset;
	    #for atom1 , delete bond to atom2 and form bond to atom4
	    deleteBond($bonds->{$atom1}, $atom2);
	    push @{ $bonds->{$atom1} }, $atom4;
	    #for atom2, delete bond to atom1 and form bond to atom3
	    deleteBond($bonds->{$atom2}, $atom1);
	    push @{ $bonds->{$atom2} }, $atom3;
	    #for atom3, delete bond to atom4 and form bond to atom2
            deleteBond($bonds->{$atom3}, $atom4);
            push @{ $bonds->{$atom3} }, $atom2;
            #for atom4, delete bond to atom3 and form bond to atom1
            deleteBond($bonds->{$atom4}, $atom3);
            push @{ $bonds->{$atom4} }, $atom1;
	    #now update bondlist
	    delete $bondlist->{$i}{$atom2};
	    $bondlist->{$i}{$atom4} = 1;
	}
    }
}

sub deleteBond {
    my ($atomBonds, $bondAtom) = @_;
    my ($i);

    for $i (0 .. @{ $atomBonds }) {
	if ($atomBonds->[$i] == $bondAtom) {
	     splice @{ $atomBonds }, $i, 1;
	     last;
	}
    }
}

sub getPBCbonds {
    my ($atoms, $bonds, $box) = @_;
    my ($i, $j, $k, $dist, $atom1, $atom2, $bondlist, $sign);

    for $i (keys %{ $bonds }) {
        $atom1 = \%{ $atoms->{$i} };
        $j = 0;
        while ($j <= $#{ $bonds->{$i} }) {
	    next if($i > $bonds->{$i});
            $atom2 = \%{ $atoms->{$bonds->{$i}[$j]} };
	    for $k ("XCOORD","YCOORD","ZCOORD") {
		$dist = $atom1->{$k} - $atom2->{$k};
		$sign = $dist/abs($dist);
		if ($sign < 0) {
		    $dist *= -1;
		}
		if ($dist > $box->{$k}{len}/2) {
		    if ($sign < 0) {
		    	$bondlist->{$k}{$i}{$bonds->{$i}[$j]} = 1;
		    } else {
			$bondlist->{$k}{$bonds->{$i}[$j]}{$i} = 1;
		    }
		    last;
		}
	    }
	}
    }

    return $bondlist;
}

sub init {
    my (%OPTS, $rStr);
    getopt('bdsci',\%OPTS);
    die "usage: $0 -b bgf file -d \"x y z\" replication\n" . 
	"\t-s (save name) -c (center=no) -i (make periodic bonds=no)\n" 
	if (! exists($OPTS{b}) or ! exists($OPTS{d}));
    print "Initializing...";
    ($bgfFile, $rStr, $saveFile, $doCenter, $bondImages) = ($OPTS{b}, $OPTS{d}, $OPTS{s}, $OPTS{c}, $OPTS{i});
    FileTester($bgfFile);
    if ($rStr !~ /(\d+)\s+(\d+)\s+(\d+)/) {
	die "ERROR: Expected integers for x,y and z replication. Got \"$rStr\"\n";
    } else {
	$replicate = (
		      {
			  "X"      => $1,
			  "Y"      => $2,
			  "Z"      => $3,
			  "STRING" => "${1}x${2}x${3}",
		      }
		      );
    }

    if (! defined($saveFile)) {
	$saveFile = basename ($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_replicate_" . $replicate->{STRING} . ".bgf";
    }
    $doCenter = 1;
    $doCenter = 0 if (! defined($doCenter));
    $doCenter = 0 if ($doCenter =~ /0|no/i);
    $bondImages = 1;
    $bondImages = 0 if (! defined($bondImages));
    $bondImages = 0 if ($bondImages =~ /0|no/i);
    print "Done\n";
}

#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub usage;
sub updateAtomTypes;
sub updateFF;

my ($bgfFile, $saveName, $FFs);
my ($ATOMS, $BONDS, $HEADERS, $atomTypeMap);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
print "Done\nWriting Polygraf formatted FFTYPES to $saveName...";
$atomTypeMap = updateAtomTypes($ATOMS);
&updateFFs($FFs, $atomTypeMap) if (defined($FFs));
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub updateFFs {
    my ($ffList, $aMap) = @_;
    my ($i, $j, $ldata);

    local $/ = undef;
    print "updating ffs..";
    for $i (@{ $ffList }) {
	open FF, $i or die "ERROR: Cannot open $i: $1\n";
	$ldata = <FF>;
	close FF;
	for $j (keys %{ $aMap }) {
	    $ldata =~ s/ $j / $aMap->{$j}/g;
	}
	open FF, "> $i" or die "ERROR: Cannot write too $i: $1\n";
	print FF $ldata;
	close FF;
   }
}

sub updateAtomTypes {
    my ($atoms) = $_[0];
    my ($i, $fftype, $map);

    for $i (keys %{ $atoms }) {
	$fftype = $atoms->{$i}{FFTYPE};
	$fftype =~ /(.)(.*)/;
	$atoms->{$i}{FFTYPE} = uc($1) . "_" . uc($2);
	$map->{$fftype} = $atoms->{$i}{FFTYPE};
    }
    return $map;
}

sub init {
    my (%OPTS, $ffList);

    getopt('bsf',\%OPTS);
    ($bgfFile, $saveName, $ffList)  = ($OPTS{b}, $OPTS{s}, $OPTS{f});

    &usage if (! exists($OPTS{b}));

    print "Initializing...";
    FileTester($bgfFile);
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$/\.polygraf\.bgf/;
    }
    $ffList = "" if (! defined($ffList));
    while ($ffList =~ /(\S+)/g) {
	if (-e $1 and -r $1 and -T $1) {
	    system("cp $1 ./");
	    push @{ $FFs }, basename($1);
	}
    }
    print "Done\n";
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -s (save_name) -f (forcefield)
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
DATA
die "\n";

}

#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(FileTester LoadElements);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub usage;
sub fixFFtypes;

my ($bgfFile, $saveName, $fftype);
my ($ATOMS, $BONDS, $HEADERS, $ELEMENTS);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
print "Done\nFixing $fftype FFTYPES to be based on atom names...";
&fixFFtypes($ATOMS, $fftype);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub  fixFFtypes {
    my ($atoms, $ftype) = @_;
    my ($i, $j, $test_str, $is_found, $newfftype);

    for $i (keys %{ $atoms }) {
        next if($atoms->{$i}{FFTYPE} ne $ftype);
        $j = length($atoms->{$i}{ATMNAME});
        $test_str = $atoms->{$i}{ATMNAME};
        $is_found = 0;
        while($j>1 && ! $is_found) {
            $newfftype = searchElement($test_str);
            if (! $newfftype) {
                chop $test_str;
                $j--;
            } else { $is_found = 1; last; }
        }
        die "ERROR: Cannot find element that corresponds to atom $i: $atoms->{$i}{FFTYPE}\n"
            if(!$is_found);
        $atoms->{$i}{FFTYPE} = $newfftype;
    }

}

sub searchElement {
    my ($search_str) = $_[0];
    my ($i, $element);

    for $i (keys %{ $ELEMENTS }) {
        if (uc($search_str) eq uc($ELEMENTS->{$i}{SYMBOL})) {
            $element = $ELEMENTS->{$i}{SYMBOL};
            last;
        }
    }

    return $element;
}

sub init {
    my (%OPTS);
    getopt('bsf',\%OPTS);
    ($bgfFile, $saveName, $fftype) = ($OPTS{b},$OPTS{s},$OPTS{f});
    for ($bgfFile, $fftype) {
        &usage if (! defined($_));
    }
    print "Initializing...";
    FileTester($bgfFile);
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$/_mod\.bgf/;
    }

    $ELEMENTS = LoadElements();
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -f fftype
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
  fftype: fftype to change based on atom name
DATA
die "\n";

}

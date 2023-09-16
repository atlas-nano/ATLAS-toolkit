#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF GetBGFAtoms);
use Packages::General qw(GetSelections FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub usage;
sub getCons;

my ($bgfFile, $saveName, $selection);
my ($ATOMS, $BONDS, $SELECTIONS, $BGF, $CONS, $tmp, $HEADERS, $BOX);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile,0);
print "Done\n";

print "Parsing atom/residue selection...";
$SELECTIONS = GetSelections($selection, 0);
print "Done\n";

print "Selecting relevant atoms...";
($BGF, $CONS, $tmp) = GetBGFAtoms($SELECTIONS, $ATOMS, $BONDS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $CONS });
print "Done\n";

print "Creating BGF file $saveName...";
($HEADERS) = createHeaders($BOX, $saveName);
addHeader($BGF,$HEADERS);
createBGF($BGF, $CONS, $saveName);
print "Done\n";

sub init {
    my (%OPTS, $select);
    getopt('bos',\%OPTS);
    ($bgfFile, $saveName, $select) = ($OPTS{b},$OPTS{s},$OPTS{o});
    for ($bgfFile, $select) {
        &usage if (! defined($_));
    }
    print "Initializing...";
    FileTester($bgfFile);
    if ($select =~ /\s+/) {
	@{ $selection } = split /\s+/, $select;
    } else {
	$selection->[0] = $select;
    }
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$/_mod\.bgf/;
    }
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -o options
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
  options:
    [^][:][I|T|N][a|r]
    a   - atom
    r   - residue
    IaX - atom number X
    IrX - residue index X
    TaX - atom type X
    NaX - atom name X
    NrX - residue name X
    Use ":" to specify a range, eg. :Tr1-8 :Ia3-66
    Use "^" to exclude a selection. You can use multiple combinations
    range and exclusion enclosed in quotes, eg, "^:TrIP-IM ^:Ia23-45"
    to exclude residues of type IM and IP and atoms 23-45
DATA

die "\n";

}

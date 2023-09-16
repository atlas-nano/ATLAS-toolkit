#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo addHeader createBGF createHeaders addBoxToHeader insertHeaderRemark);
use BOX qw(GetBox);

sub usage;
sub removeDimDisplay;

my ($bgfFile, $saveName, $dim);
my ($ATOMS, $BONDS, $HEADERS, $BOX);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
print "Done\nRemoving any $dim settings...";
&removeDimDisplay($ATOMS, $BONDS, $dim);
print "Done\nCreating BGF file $saveName...";
&insertHeaderRemark($HEADERS, "REMARK $bgfFile make 2d periodic removed $dim");
$BOX = GetBox($ATOMS, undef, $HEADERS);
&addBoxToHeader($HEADERS, $BOX);
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub removeDimDisplay {
	my ($atoms, $bonds, $dim) = @_;
	my ($i, @list, $j, $k);

	for $i (keys %{ $atoms }) {
		if (exists($atoms->{$i}{"DISP${dim}"})) {
			@list = @{ $atoms->{$i}{"DISP${dim}"} };
			$j = 0;
			while ($j <= $#list) {
				if ($list[$j] == 0) {
					$j++;
					next;
				}
				splice @{ $bonds->{$i} }, $j, 1;
				for $k ("DISPX", "DISPY", "DISPZ") {
					splice @{ $atoms->{$i}{$k} }, $j, 1 if(exists($atoms->{$i}{$k}));
				}
				splice @list, $j, 1;
			}
			delete $atoms->{$i}{"DISP${dim}"};
		}
			
	}
}

sub init {
	my (%OPTS, $select);
	getopt('bds',\%OPTS);
	($bgfFile, $saveName, $dim) = ($OPTS{b},$OPTS{s},$OPTS{d});
	for ($bgfFile, $dim) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	FileTester($bgfFile);
	if ($dim =~ /^(x|y|z)/i) {
		$dim = uc $1;
	} else {
		die "ERROR: Expected x|y|z for dim. Got \"$dim\"\n";
	}
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/_mod\.bgf/;
	}
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -d dimension
Arguments:
  bgf_file: name of bgf_file
  dimension: either x|y|z
  save_name: name of file to save
DATA

die "\n";

}

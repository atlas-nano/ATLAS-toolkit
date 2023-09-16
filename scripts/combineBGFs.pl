#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF);
use General qw(FileTester CombineMols);
use File::Basename qw(basename);

my ($mol1, $mol2, $tot, $saveName);
my ($mSys, $ATOMS, $BONDS, $HEADERS, $sATOMS, $sBONDS);

$|++;
&init;
($sATOMS, $sBONDS, $HEADERS) = loadBGFs($mSys);
($ATOMS, $BONDS) = mergeSystems($sATOMS, $sBONDS);
print "Creating $saveName...";
&addHeaders(\%{ $ATOMS }, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub addHeaders {
	my ($atoms, $headers) = @_;
	my ($i);

	for $i (1 .. $tot) {
		&addHeader(\%{ $atoms }, $headers->{$i});
	}
}

sub mergeSystems {
	my ($aList, $bList) = @_;
	my ($i, $atoms, $bonds);

	print "Combining systems...";
	%{ $atoms } = %{ $aList->{1} };
	%{ $bonds } = %{ $bList->{1} };
	for $i (2 .. $tot) {
		($atoms, $bonds) = CombineMols($atoms, $aList->{$i}, $bonds, $bList->{$i});
	}
	print "Done\n";
	return ($atoms, $bonds);
}

sub loadBGFs {
	my ($fList) = @_;
	my ($i, $atoms, $bonds, $headers);

	for $i (1 .. $tot) {
		print "Parsing file $i $fList->{$i}...";
		($atoms->{$i}, $bonds->{$i}, $headers->{$i}) = GetBGFFileInfo($fList->{$i}, 1);
		print "Done\n";
	}
	return ($atoms, $bonds, $headers);
}

sub init {
	my (%OPTS, $index, $sStr);
	getopt('bs',\%OPTS);

	&usage if (! exists($OPTS{b}));
	print "Initializing...";
	$tot=0;
	while ($OPTS{b} =~ /(\S+)/g) {
		$mSys->{++$tot}=$1;
		FileTester($mSys->{$tot});
		$sStr .= $mSys->{$tot};
		$sStr =~ s/\.\w+$//;
		$sStr .= "_";
	}
	die "ERROR: Need at least 2 bgf files to merge!"
		if ($tot<2);

	if (! exists($OPTS{s})) {
		$sStr =~ s/_$//;
		$saveName = "${sStr}.merged.bgf"
	} else {
		$saveName = $OPTS{s};
	}
	print "Done\n";
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b "bgf_files" -s save_name
Arguments:
  bgf_files: name and locations of bgf_files
  save_name: name of file to save
DATA
die "\n";

}

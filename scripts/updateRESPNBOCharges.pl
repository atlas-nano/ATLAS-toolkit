#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use MolData;
use Getopt::Std qw(getopt);

sub init;
sub getRESPCharges;
sub getTeraChemNBOCharges;
sub updateCharges;

my ($molStructure, $chargeFile, $saveName, $fileType, $fileName);
my ($CHARGES, $chargeFunc);

$|++;
&init;

print "Gettting data from $fileType file $fileName...";
$molStructure->read($fileName, $fileType);
print "Done\nParsing charge file $chargeFile....";
$CHARGES = $chargeFunc->($chargeFile);
print "Done\nUpdating charges...";
&updateCharges($molStructure, $CHARGES);
print "Done\nWriting to $fileType file $saveName...";
$molStructure->write($saveName, $fileType);
print "Done\n";

sub updateCharges {
	my ($molData, $chargeData) = @_;
	my ($tot, $i, $atomTot);

	$tot = scalar(keys %{ $chargeData });
	$atomTot = $molData->count->{"atoms"};
	die "ERROR: Unequal number of atoms in data file ($atomTot) vs charges ($tot)!\n"
	if ($tot != $atomTot);

	for $i (1 .. $tot) {
		$molData->atoms($i)->("charge", $chargeData->{$i});
	}
}

sub getTeraChemNBOCharges {
	my ($infile) = $_[0];
	my (%DATA, $valid);

	$valid = 0;
	open CHRG, $infile or die "ERROR: Cannot open TeraChem NBO charge file $infile: $!\n";
	while(<CHRG>) {
		chomp;
		$valid = 1 if ($_ =~ /Summary of Natural Population Analysis/);
		$DATA{$1}=  $2 
			if ($valid and $_ =~ /^\s+[a-zA-Z]+\s*(\d+)\s+(\-?\d+\.\d+)\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+/);
		last if ($_ =~ / \* Total \*/);	
	}
	close CHRG;
	die "ERROR: invalid TeraChem file\n"
		if (! $valid);
	return \%DATA;
}

sub getRESPCharges {
	my ($infile) = $_[0];
	my (%DATA, $count);
	
	$count = 0;
	open CHRG, $infile or die "ERROR: Cannot open resp charge file $infile: $!\n";
		while (<CHRG>) {
		chomp;
		if ($_ =~ /^\s*\-?\d+\.\d+/) {
			while ($_ =~ /(\-?\d+\.\d+)/g) {
				$count++;
				$DATA{$count} = $1;
			}
		}
	}
	close CHRG;

	die "ERROR: No valid data while reading $infile!\n" if (! $count);

	return \%DATA;
}

sub init {
	my (%OPTS, $chargeType);

	getopt('ftscd', \%OPTS);
	for ("f", "c") {
	die "usage: $0 -f structure_file -c charge_file -t (structure_file_type) -d (charge_file_type) -s (savename)\n" 
		if (! exists($OPTS{$_}));
	}

	print "Initialzing...";
	$molStructure =  MolData->new();
	($fileName, $fileType, $saveName, $chargeFile, $chargeType) = ($OPTS{f}, $OPTS{t}, $OPTS{s}, $OPTS{c}, $OPTS{d});
	$molStructure->testFile($chargeFile);
	if (! defined($fileType)) {
		$fileType = "bgf";
		if ($fileName =~ /\.(\w+)$/) {
			$fileType = lc $1;
		}
	}
	$molStructure->testFile($fileName);
	$chargeType = "resp" if (! defined($chargeType));
	die "ERROR: Expected resp|terachem for charge_file_type. Got \"$chargeType\"\n"
		if($chargeType !~ /(resp|terachem)/i);
	$chargeType = lc $1;
	$chargeFunc = \&getTeraChemNBOCharges;
	if ($chargeType eq "resp") {
		$chargeFunc = \&getRESPCharges;
	}

	$saveName = $molStructure->getFileName($fileName) if (! defined($saveName));
	print "Done\n";
}

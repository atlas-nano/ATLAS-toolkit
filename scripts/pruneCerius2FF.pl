#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use strict;
use FileFormats qw(GetBGFFileInfo);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub usage;
sub pruneFF;
sub saveFF;

my ($bgfFile, $saveName, $ffFile);
my ($newFF, $ATOMS, $BONDS, $HEADERS);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
print "Done\nPrunning FF ${ffFile}....";
$newFF = pruneFF($ffFile,$ATOMS);
print "Done\nSaving pruned ff ${saveName}...";
&saveFF($newFF, $saveName);
print "Done\n";

sub saveFF {
	my ($ffdata, $sfile) = @_;

	open OUTFF, "> $sfile" or die "ERROR: Cannot create $sfile: $!\n";
	print OUTFF $ffdata;
	close OUTFF;
}

sub pruneFF {
	my ($ff, $atoms) = @_;
	my ($fflist, $i, $ffstr, $header, $headerstr);
	my ($chead, $num, $shouldPrint, $tmp, $outStr, $instr, $sstr, $thead);

	$outStr = "";
	for $i (keys %{ $atoms }) {
		$fflist->{ $atoms->{$i}{FFTYPE} } = 1;
	}
	$ffstr = "";
	for $i (keys %{ $fflist }) {
		$i =~ s/\+/\\\+/;
		$i =~ s/\*/\\\*/;
		$ffstr .= "$i|";
	}
	$ffstr .= "X";

	$header = (
				{
					"ATOMTYPES" => 1,
					"DIAGONAL_VDW" => 1,
					"ATOM_TYPING_RULES" => 1,
					"OFF_DIAGONAL_VDW" => 2,
					"BOND_STRETCH" => 2,
					"ANGLE_BEND" => 3,
					"UREY_BRADLEY" => 3,
					"TORSIONS" => 4,
					"INVERSIONS" => 4,
					"GENERATOR" => 1,
				}
			);
	$headerstr = "";		
	for $i (keys %{ $header }) {
		$headerstr .= "$i|";
	}
	chop $headerstr;

	open FFSTR, $ff or die "ERROR: Cannot open $ff: $!\n";
	while (<FFSTR>) {
		chomp;
		$instr = $_;
		if ($instr  =~ /^($headerstr)/) {
			$chead = $1;
			$num = $header->{$chead};
			$outStr .= "$instr\n";
		} elsif ($instr =~ /END/) {
			$outStr .= "$instr\n";
			undef $chead;
		} elsif (defined($chead)) {
			$sstr = $instr;
			$sstr =~ s/^\s*//;
			@{ $tmp } = split /\s+/,$sstr;
			if ($chead eq "TORSIONS" and $#{ $tmp } > 4) {
				@{ $thead } = ($tmp->[0], $tmp->[1], $tmp->[2], $tmp->[3]);
			} elsif ($chead eq "TORSIONS" and $#{ $tmp } < 4) {
				unshift @{ $tmp }, @{ $thead };
			}
			$shouldPrint = 1;
			for $i (0 .. ($num-1)) {
				$shouldPrint = 0 if ($tmp->[$i] and $tmp->[$i] !~ /^($ffstr)$/);
			}
			$outStr .= "$instr\n" if ($shouldPrint);
		} else {
			$outStr .= "$instr\n";
		}
	}
	close FFSTR;
	return $outStr;
}

sub init {
	my (%OPTS, $atomSel);
	getopt('bfs',\%OPTS);
	($bgfFile, $saveName, $ffFile) = ($OPTS{b},$OPTS{s},$OPTS{f});
	for ($bgfFile, $ffFile) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	FileTester($bgfFile);
	FileTester($ffFile);
	if (! defined($saveName)) {
		$saveName = basename($ffFile);
		$saveName =~ s/\.\w+//;
		$saveName .= ".prune.bgf";
	}
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -f ff_file -s (save_ff_name) 
Arguments:
  bgf_file: name of bgf_file
  ff_file: name of cerius2 formatted forcefield file
  save_ff_name: name of file to save
DATA
die "\n";

}

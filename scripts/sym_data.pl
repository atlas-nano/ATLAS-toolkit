#!/usr/bin/perl -w

use strict;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub parseDataFile;
sub writeData;
sub numerically {($a<=>$b)}

my ($datafile, $savefile, $sym);
my ($DATA);

$|++;
&init;

print "Reading data file $datafile...";
&parseDataFile($datafile,\%{ $DATA }, $sym);
print "Done\nSymmeterizing data by reflecting $sym->{str} on $sym->{ncols} columns...";
&sym($DATA,$sym);
print "Done\nWriting data to $savefile...";
&writeData($DATA,$savefile);
print "Done\n";

sub writeData {
	my ($data, $sfile) = @_;
	my ($x, $y, $i);

	open OUTFILE, "> $sfile" or die "ERROR: Cannot write to $sfile: $!\n";
	for $x (sort numerically keys %{ $data }) {
		for $y (sort numerically keys %{ $data->{$x} }) {
			printf OUTFILE "%s %s ",$x,$y;
			for $i (@{ $data->{$x}{$y} }) {
				printf OUTFILE "%.4f ",$i;
			}
			printf OUTFILE "\n";
		}
		printf OUTFILE "\n";
	}
	close OUTFILE;
}

sub sym {
	my ($data, $sym) = @_;
	my ($x, $xs, $y, $ys, $z, $c, $tot, $mid, $i, $symdata);

	$tot = $sym->{ncols}-1;
	$mid = int($tot/2);

	for $x (keys %{ $data }) {
		$xs = -$x;
		$xs = $x if (!exists($sym->{AXIS}{X}));
		for $y (keys %{ $data->{$x} }) {
			$ys = -$y;
			$ys = $y if (!exists($sym->{AXIS}{y}));
			next if (!exists($data->{$xs}) or ! exists($data->{$xs}{$ys}));
			for $i (0 .. $mid) {
				$symdata->{$x}{$y}[$i] = ($data->{$x}{$y}[$i]+$data->{$xs}{$ys}[$tot-$i])/2;
			}
		}
	}
	%{ $data } = %{ $symdata };
}

sub parseDataFile {
	my ($infile, $data, $sym) = @_;
	my ($x,$y, @vals, $inStr);

	open DATAFILE, $infile or die "ERROR: Cannot open $infile: $!\n";
	while(<DATAFILE>) {
		chomp;
		$inStr = $_;
		if ($inStr =~ /^(\-?\d+\.?\d*)\s+(.*)$/) {
			$x = $1;
			($y,@vals) = split /\s+/,$2;
			next if ($#vals<1);
			@{ $data->{$x}{$y} } = @vals;
			$sym->{ncols} = scalar(@vals);
		}
	}
	close DATAFILE;
	die "ERROR: No valid data read from $infile!\n"
		if (! keys %{ $data });
}

sub init {
	my (%OPTS,$symStr);

	&getopt('dws',\%OPTS);
	die "usage: $0 -d datafile -s (sym_axis=x|y|xy(default)) -w (save_name)\n" 
		if (!exists($OPTS{d}));
	die "ERROR: Cannot access $OPTS{d}: $!\n" 
		if (!-e $OPTS{d} or ! -r $OPTS{d} or ! -T $OPTS{d});
	print "Initializing...";
	($datafile, $savefile, $symStr) = ($OPTS{d}, $OPTS{w}, $OPTS{s});
	if (!defined($savefile)) {
		$savefile = $datafile;
		$savefile =~ s/\.\w+$//;
		$savefile .= ".sym.dat";
	}
	$symStr = "xy" if (!defined($symStr));
	while ($symStr =~ /(x|y)/gi) {
		$sym->{AXIS}{uc $1} = 1;
	}
	die "ERROR: No valid axis value found while searching \"$symStr\"\n"
		if (! defined($sym));
	for (keys %{ $sym->{AXIS} }) {
		$sym->{str} .= "$_";
	}
	print "Done\n";
}

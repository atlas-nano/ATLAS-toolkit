#!/usr/bin/perl -w

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub parseDataFile;
sub createHistrogram;
sub writeHistrogram;

my ($datafile, $ivals, $savefile);
my ($DATA, $HIST);

$|++;
&init;
print "Parsing $datafile...";
$DATA = parseDataFile($datafile);
print "Done\nCreating histograms...";
$HIST = createHistogram($DATA, $ivals);
print "Done\nWriting histogram to $savefile...";
&writeHistogram($HIST, $DATA, $ivals, $savefile);
print "Done\n";

sub findVal {
	my ($hist, $rec, $num, $cval) = @_;
	my ($mval, $i, $val, $found);

	$val = 0;
	$found = 1;
	$mval = \%{ $hist };
	for $i (keys %{ $rec }) {
		if (!exists($mval->{$rec->{$i}})) { 
			$found = 0;
			last;
		}
		$mval = \%{ $mval->{$rec->{$i}} };
	}

	$val = $mval->{$cval} if ($found and exists($mval->{$cval}));
	return $val;

}

sub createOutStr {
	my ($rec, $data) = @_;
	my ($outStr, $i);

	for $i (sort { $a <=> $b } keys %{ $rec }) {
		$outStr .= sprintf("%.3G ",($rec->{$i}+$data->{BOUNDS}{$i}{MIN})*$data->{BOUNDS}{$i}{BIN});
	}
	return $outStr;
}

sub dumpVals {
	my ($hist,$data,$rec,$num,$index,$OUTFILE) = @_;
	my ($i, $curr, $start, $val,$outStr);

	if ($num == $index) {
		$start = $data->{BOUNDS}{$index}{MIN};
		$outStr = createOutStr($rec,$data);
		for $i (1 .. $data->{BOUNDS}{$index}{NUM}) {
			$curr = $start + $i*$data->{BOUNDS}{$index}{BIN};
			$val = findVal($hist, $rec, $num, $i);
			printf $OUTFILE "${outStr}%.3G %10d\n",$curr,$val;
		}
		print $OUTFILE "\n";
		$outStr = "";
	} else {
		$start = $data->{BOUNDS}{$index}{MIN};
		for $i (1 .. $data->{BOUNDS}{$index}{NUM}) {
			$curr = $start + $i*$data->{BOUNDS}{$index}{BIN};
			$rec->{$index} = $i;
			&dumpVals($hist,$data,$rec,$num,$index+1,$OUTFILE);
		}
	}

}

sub writeHistogram {
	my ($hist, $data, $bins, $sfile) = @_;
	my ($i, $j, $rec, $num);

	for $i (keys %{ $data->{BOUNDS} }) {
		$data->{BOUNDS}{$i}{NUM} = sprintf("%.0f", ($data->{BOUNDS}{$i}{MAX}-$data->{BOUNDS}{$i}{MIN})/$bins->[$i]);
		$data->{BOUNDS}{$i}{BIN} = $bins->[$i];
	}
	$i = 0;
	$num = scalar(@{ $bins }) - 1;
	open OUTFILE, "> $sfile" or die "ERROR: Cannot write to $sfile: $!\n";
	&dumpVals($hist,$data,\%{ $rec }, $num,$i,\*OUTFILE);
	close OUTFILE;
}

sub binData {
	my ($hist, $data, $bins, $bounds, $index) = @_;
	my ($i, $j, $ibin);

	for $i (keys %{ $data }) {
		$ibin = sprintf("%.0f",($i-$bounds->{$index-1}{MIN})/$bins->[$index-1]);
		if (keys %{ $data->{$i} }) {
			&binData(\%{ $hist->{$ibin} },$data->{$i},$bins,$bounds,$index+1);
		} else {
			$hist->{$ibin}++;
		}
	}
}

sub createHistogram {
	my ($data, $bins) = @_;
	my ($i, $j, $HIST, $ibin);

	if (scalar(keys %{ $data->{BOUNDS} }) > scalar(@{ $bins })) {
		$i = scalar(@{ $bins });
		$j = scalar(keys %{ $data->{BOUNDS} }) - 1;
		for ($i .. $j) {
			push @{ $bins }, ($data->{BOUNDS}{$_}{MAX}-$data->{BOUNDS}{$_}{MIN})/25;
		}
	}

	&binData(\%{ $HIST },$data->{VALS},$bins,$data->{BOUNDS},1);
	return $HIST;
}

sub parseDataFile {
	my ($dfile) = $_[0];
	my (%DATA, $valid, $i, $vals, $rec);
	

	$valid = $i = 0;
	open DATAFILE, $dfile or die "ERROR: Cannot open $dfile: $!\n";
	while (<DATAFILE>) {
		if ($_ =~ /^\s*(\-?\d+\.?\d*e?\-?\d*)\s*(.+)/) {
			$vals = $2; 
			$rec = \%{ $DATA{VALS}{$1} };
			$valid = 1;
			$i = 0;
			$DATA{BOUNDS}{$i}{MIN} = $1 if (!exists($DATA{BOUNDS}{$i}{MIN}) or $DATA{BOUNDS}{$i}{MIN} > $1);
			$DATA{BOUNDS}{$i}{MAX} = $1 if (!exists($DATA{BOUNDS}{$i}{MAX}) or $DATA{BOUNDS}{$i}{MAX} < $1);
			while ($vals =~ /(\-?\d+\.?\d*e?\-?\d*)/g) {
				$rec = \%{ $rec->{$1} };
				$i++;
				$DATA{BOUNDS}{$i}{MIN} = $1 if (!exists($DATA{BOUNDS}{$i}{MIN}) or $DATA{BOUNDS}{$i}{MIN} > $1);
				$DATA{BOUNDS}{$i}{MAX} = $1 if (!exists($DATA{BOUNDS}{$i}{MAX}) or $DATA{BOUNDS}{$i}{MAX} < $1);
			}
		}
	}
	die "ERROR: No valid values found while searching $dfile!\n" if (! $valid);

	return \%DATA;
}

sub init {
	my (%OPTS, $ivalString);

	getopt('dbs',\%OPTS);

	die "usage: $0 -d data_file -b [binsize=0.1] -s [save name]\n" if (! exists($OPTS{d}));
	print "Initializing...";
	die "ERROR: Cannot locate $OPTS{d}: $!\n" if (! -e $OPTS{d} or ! -r $OPTS{d} or ! -T $OPTS{d});
	($datafile, $ivalString, $savefile) = ($OPTS{d}, $OPTS{b}, $OPTS{s});
	if (! defined($savefile)) {
		$savefile = basename($datafile);
		$savefile =~ s/\.\w+$//;
		$savefile .= ".hist.dat";
	}
	@{ $ivals } = ();
	$ivalString = "" if (! defined($ivalString));
	while($ivalString =~ /(\d+\.?\d*)/g) {
		push @{ $ivals }, $1;
	}
	print "Done\n";
}

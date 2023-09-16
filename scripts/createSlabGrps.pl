#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use warnings;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester);
use ManipAtoms qw(GetAtmList GetMols);


sub init;
sub parseAtomPosFile;
sub saveSlabDataByField;
sub numerically { ($a<=>$b); }

my ($bgfFile, $atmPosFile, $saveFile, $num, $shellDist, $FIELD, $sfield, $iOff);
my ($avgPos, $ATOMS, $BONDS, $HEADERS, $MOLS);

$|++;
&init;
print "Parsing bgf file...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$MOLS = GetMols($ATOMS, $BONDS);
print "Done\nParsing atom position file $atmPosFile...";
&parseAtomPosFile($ATOMS, $atmPosFile);
print "Done\nCreating $num groups...";
&saveSlabDataByField($ATOMS, $MOLS, $num, $shellDist, $FIELD, $sfield, $avgPos, $iOff);
print "Done\nWriting BGF file $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub saveSlabDataByField {
	my ($atoms, $mols, $tot, $dists, $field, $subF, $aPos, $uoffset) = @_;
	my ($dens, $i, @list, $min, $max, $index);
	my ($j, $count, $inc, $k, $avg, $offset);

	$offset = $uoffset;
	$offset = 64 + $uoffset if ($field =~ /CHAIN/);

	for $i (keys %{ $atoms }) {
		$dens->{$atoms->{$i}{POS}{AVG}}++;
	}
	@list = sort numerically keys %{ $dens };
	$max = $list[$#list];
	$min = $list[0];
	if(defined($dists) and $list[$#list] < $dists->[$#{ $dists }]) {
		$list[$#list+1] = $dists->[$#{ $dists }]; #if the user specified more distances than we find differences
							#make the greatest entry in the list the largest specified distance
	}

	if (! defined($dists)) {
	#need to make equal intervals
		$num = $#list+1 if(($#list+1)<$num);
		$inc = (($max-$min)/$num); # increment based on distance distributions
		$dists->[0] = $min + 0.9*$inc;
		$dists->[$num-1] = $max;
		for $i (1 .. ($num-2)) {
			$dists->[$i] = $min+$inc*($i+0.9);
		}
	}

	for $i (keys %{ $mols }) {
		$avg = $count = 0;
		for $j (keys %{ $mols->{$i}{MEMBERS} }) {
			$count++;
			$avg += $atoms->{$j}{POS}{AVG};
		}
		$avg = sprintf("%.0f",$avg/$count);
		$index = $offset + $num;
		for $k (0 .. ($num-1)) {
			$index = $offset + $k+1;
			last if ($avg < $dists->[$k]);
		}
		for $j (keys %{ $mols->{$i}{MEMBERS} }) {
			if($field =~ /CHAIN|RESNAME/) {
				$atoms->{$j}{$field} = chr($index);
			} else {
				$atoms->{$j}{$field} = $index;
			} 
		}
	}
}

sub parseAtomPosFile {
	my ($atoms, $infile) = @_;
	my ($count, $tot);

	$tot = scalar(keys %{ $ATOMS });
	open ATOMPOSFILE, $infile or die "ERROR: Cannot open $infile: $!\n";
	while(<ATOMPOSFILE>) {
		chomp;
		if ($_ =~ /^\s+(\d+)\s+(\d+\.?\d*)\s+\-?(\d+\.\d+)/) {
			next if (! exists($atoms->{$1}));
			$atoms->{$1}{POS}{AVG} = $2;
			$avgPos += $2;
			$atoms->{$1}{POS}{STDEV} = $3;
			$count++;
		}
	}
	close ATOMPOSFILE;
	die "ERROR: $infile has $count entries while bgf file has $tot atoms!\n"
	if($count != $tot);
	$avgPos /= $tot;
}

sub init {
	my (%OPTS, $select, @tmp, $shellDistStr, $field);

	getopt('basndfi',\%OPTS);
	die "usage: $0 -b bgf file -a atom postion file -f (field: chain(default)|molsize|moleculeid|resnum|resname) -n [numgroups = 10] -d [shell distances] -s [save name] -g [group subfield=none] -i [index offset]\n"
	if(! defined($OPTS{b}) or ! defined($OPTS{a}));

	print "Initializing...";
	($bgfFile, $atmPosFile, $saveFile, $num, $shellDistStr, $field, $sfield, $iOff) = ($OPTS{b}, $OPTS{a}, $OPTS{s}, $OPTS{n}, $OPTS{d}, $OPTS{f}, $OPTS{g}, $OPTS{i});
	FileTester($bgfFile);
	FileTester($atmPosFile);
	die "ERROR: Expected positive integer for $num, got \"$num\"\n" if (!defined($shellDistStr)  and ($num !~ /^\d+$/ or !$num));
	if (! defined($saveFile)) {
		$saveFile = basename($bgfFile);
		$saveFile =~ s/\.\w+//;
		$saveFile .= "_slab.bgf";
	}
	undef $shellDistStr if(defined($shellDistStr) and $shellDistStr !~ /^\d+\.?\d*/);
	if (defined($shellDistStr)) {
		$num = 0;
		while ($shellDistStr =~ /(\d+\.?\d*)/g) {
			push @{ $shellDist }, $1;
			$num++;
		}
		if (exists($OPTS{n}) and $num > $OPTS{n}) {
			while ($num > $OPTS{n}) {
				splice @{ $shellDist }, ($#{ $shellDist } -1), 1;
				$num--;
				last if ($num == $OPTS{n});
				splice @{ $shellDist }, 2, 1;
				$num--;
			}
		}
	}
	if (defined($field) and $field =~ /(molsize|molecule|resname|resid|resnum|chain)/i) {
		$field = $1;
	} else {
		$field = "chain";
	}
	$field = "moleculeid" if ($field eq "molecule");
	$field = "resnum" if ($field eq "resid");
	$FIELD = uc $field;
	if(defined($sfield) and $sfield =~ /(molsize|molecule|resname|resid|resnum|chain)/i) {
		$sfield = uc $1;
	} elsif (defined($sfield)) {
		undef $sfield;
	}
	$iOff = 0 if (! defined($iOff) or $iOff !~ /^\d+$/);
	if($iOff > 0) {
		print "...offset: $iOff...";
	}
	print "Done\n";
}

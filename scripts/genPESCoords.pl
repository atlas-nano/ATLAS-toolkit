#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester GetSelections);

sub init;
sub showusage;
sub parseJagInFile;
sub genCoords;
sub createStructure;
sub writeInFile;

my ($inFile, $SPACING, $CONSTRAINTS, $prefix, $maefile, $offset);
my ($DATA, $COORDS, $tot, $prec, $isOptimize);

$|++;
&init;
print "Parsing Jaguar Input file $inFile..";
($DATA, $tot) = parseJagInFile($inFile, $CONSTRAINTS->{FROZEN});
print "Done\nGenerating coords based on options...";
$COORDS = genCoords($DATA, $SPACING, $CONSTRAINTS->{PLANE}, $prefix);
print "Done\nCreating Jaguar Input files...";
&writeInFile($COORDS);
print "Done\n";

sub writeInFile {
    my ($fdata) = $_[0];
    my ($i, $filename);
    
    for $i (@{ $fdata }) {
	$filename = $i->{NAME} . ".in";
	open OUTFILE, "> $filename" or die "ERROR: Cannot create $filename: $!\n";
	print OUTFILE $i->{DATA};
	close OUTFILE;
    }
}

sub createStructure {
    my ($data, $entryName, $dimA, $offsetA, $dimB, $offsetB, $cpFlag) = @_;
    my ($rec, $i, $coord);

    for $i (@{ $data->{HEADER} }) {
	if ($i =~ /entry_name/) {
	    $i = "entry_name: $entryName";
	} elsif ($i =~ /^MAEFILE:(.*)$/) {
	    $i = "MAEFILE: ${entryName}.mae";
	    $maefile = $1 if (! defined($maefile));
	}
	$rec .= "$i\n";
    }

    for $i (1 .. $tot) {
	$coord->{x} = sprintf("%17.13f",$data->{ATOMS}{$i}{XCOORD});
	$coord->{y} = sprintf("%17.13f",$data->{ATOMS}{$i}{YCOORD});
	$coord->{z} = sprintf("%17.13f",$data->{ATOMS}{$i}{ZCOORD});
	if (exists($data->{ATOMS}{$i}{FROZEN})) { 
	    if ($isOptimize == 1) {
		$coord->{x} .= "#";
		$coord->{y} .= "#";
		$coord->{z} .= "#";
	    }
	} else {
	    $coord->{$dimA} += $offsetA;
	    $coord->{$dimA} = sprintf("%17.13f",$coord->{$dimA});
	    $coord->{$dimA} .= "#" if ($isOptimize == 1);
	    if (defined($dimB)) {
		$coord->{$dimB} += $offsetB;
		$coord->{$dimB} = sprintf("%17.13f",$coord->{$dimB});
		$coord->{$dimB} .= "#" if ($isOptimize == 1);
	    }
	}
	if (defined($cpFlag)) {
	    if ($cpFlag == 1 and exists($data->{ATOMS}{$i}{FROZEN})) {
		$rec .= sprintf("%-6s $coord->{x} $coord->{y} $coord->{z}\n",$data->{ATOMS}{$i}{ATMNAME}. "@");
	    } elsif ($cpFlag == 2 and ! exists($data->{ATOMS}{$i}{FROZEN})) {
		$rec .= sprintf("%-6s $coord->{x} $coord->{y} $coord->{z}\n",$data->{ATOMS}{$i}{ATMNAME}. "@");
	    } else {
		$rec .= sprintf("%-6s $coord->{x} $coord->{y} $coord->{z}\n",$data->{ATOMS}{$i}{ATMNAME});	
	    }
	} else {
	    $rec .= sprintf("%-6s $coord->{x} $coord->{y} $coord->{z}\n",$data->{ATOMS}{$i}{ATMNAME});
	}
    }
	    
    for $i (@{ $data->{FOOTER} }) {
	$rec .= "$i\n";
    }
    
    return $rec;
}

   
sub genCoords {
    my ($data, $spacing, $constraints, $entryPrefix) = @_;
    my ($i, @FDATA, $rec, $a, $b, $entry, $j);

    for $i (@{ $spacing }) {
	for $a (keys %{ $constraints }) {
	    if (! keys %{ $constraints->{$a} }) {
		$rec = ();
		$entry = $rec->{NAME} = "${entryPrefix}_${a}_" . sprintf("%${prec}", $i) . "_dimer.0";
		$rec->{DATA} = &createStructure($data, $entry, $a, ($i - $offset), undef, 0);
		push @FDATA, $rec;
		if ($isOptimize == 2) {
		    $rec = ();
                    $entry = $rec->{NAME} = "${entryPrefix}_${a}_" . sprintf("%${prec}", $i) . "_a_cp.0";
                    $rec->{DATA} = &createStructure($data, $entry, $a, ($i - $offset), undef, 0, 1);
                    push @FDATA, $rec;
		    $rec = ();
                    $entry = $rec->{NAME} = "${entryPrefix}_${a}_" . sprintf("%${prec}", $i) . "_b_cp.0";
                    $rec->{DATA} = &createStructure($data, $entry, $a, ($i - $offset), undef, 0, 2);
                    push @FDATA, $rec;
		}
	    } else {
		for $j (@{ $spacing }) {
		    $rec = ();
		    for $b (keys %{ $constraints->{$a} }) {
			$entry = $rec->{NAME} = "${entryPrefix}_${a}_" . sprintf("%${prec}", $i) .
			"_${b}_" . sprintf("%${prec}", $j);
			$rec->{DATA} = &createStructure($data, $entry, $a, ($i - $offset), $b, $j);
			push @FDATA, $rec;
		    }
		}
	    }
	}
    }
    return \@FDATA;
}

sub parseJagInFile {
    my ($fileName, $frozen) = @_;
    my (%JAGDATA, $start, $count);

    $start = $count = 0;
    open JAGFILE, $fileName or die "ERROR: Cannot open $fileName: $!\n";
    while (<JAGFILE>) {
	chomp;
	if ($_ =~ /^\s*\&zmat/o) {
	    $start = 1;
	    push @{ $JAGDATA{HEADER} }, $_;
	} elsif ($start and $_ =~ /^\s*\&/) {
	    $start = 0;
	    push @{ $JAGDATA{FOOTER} }, $_;
	} elsif ($start and $_ =~ /^\s*(\w+)\s+(\-?\d+\.\d*)\#?\s+(\-?\d+\.\d*)\#?\s+(\-?\d+\.\d*)\#?(.+)/) {
	    $count++;
	    $JAGDATA{ATOMS}{$count} = (
				       {
					   "ATMNAME" => $1,
					   "XCOORD"  => $2,
					   "YCOORD"  => $3,
					   "ZCOORD"  => $4,
				       }
				       );
	    if (defined($5)) {
		$JAGDATA{ATOMS}{$count}{OTHER} = $5;
	    }
	    $JAGDATA{ATOMS}{$count}{FROZEN} = 1 if exists($frozen->{$count});
	} else {
	    if (exists($JAGDATA{ATOMS})) { #comes after atoms
		push @{ $JAGDATA{FOOTER} }, $_;
	    } else { # comes before
		push @{ $JAGDATA{HEADER} }, $_;
	    }
	}
    }
    die "ERROR: $fileName is not a valid Jaguar Input file!\n" if (! %JAGDATA);
    return (\%JAGDATA, $count);
}

sub init {
    my (%OPTS, $usage, $spacingInfo, $constraintInfo);
    my ($start, $end, $interval, $tmp, $i, $frozen);

    $usage = &showusage;
    getopt('icsfoa',\%OPTS);
    for ("i", "s", "f") {
	die "$usage" if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($inFile, $spacingInfo, $constraintInfo, $frozen, $isOptimize, $offset) = 
	($OPTS{i}, $OPTS{s}, $OPTS{c}, $OPTS{f}, $OPTS{o}, $OPTS{a});
    FileTester($inFile);
    $prefix = basename($inFile);
    $prefix =~ s/\.\w+$//;
    
    $isOptimize = 0 if (! defined($isOptimize));
    if ($isOptimize =~ /^(\d+)/) {
	if ($1 == 1 or $1 == 2) {
	    $isOptimize = $1;
	} else {
	    $isOptimize = 0;
	}
    } else {
	$isOptimize = 0;
    }

    while ($frozen =~ /(\S+)/g) {
	$tmp = $1;
	if ($tmp =~ /(\d+)\-(\d+)/) {
	    $start = $1;
	    $end = $2;
	    ($start, $end) = ($end, $start) if ($start > $end);
	} else {
	    $start = $end = $1;
	}
	for $i ($start .. $end) {
	    $CONSTRAINTS->{FROZEN}{$i} = 1;
	}
    }
    die "ERROR: Expected a-b for frozen atoms. Got \"$frozen\"\n" if (! $CONSTRAINTS);

    $prec = 0;
    while ($spacingInfo =~ /(\-?\d+\.?\d*)\-(\-?\d+\.?\d*):(\d+\.?\d*)/g) {
	$start = $1;
	$end = $2;
	$interval = $3;
	if ($interval =~ /\.(\d+)/) {
	    $prec = length($1) if (length($1) > $prec);
	}
	($start, $end) = ($end, $start) if ($start > $end);
	if ($interval > 0) {
	    while ($start <= $end) {
		push @{ $SPACING }, $start;
		$start += $interval;
	    }
	}
    }
    die "ERROR: Invalid spacing info got while searching \"$spacingInfo\"!\n" if (! $SPACING);
    $prec = ".${prec}f";

    $constraintInfo = "xy" if (! defined($constraintInfo));
    $tmp = ();
    $tmp = \%{ $CONSTRAINTS->{PLANE} };
    while ($constraintInfo =~ /(x|y|z)/ig) {
	$tmp->{lc($1)} = ();
	$tmp = \%{ $tmp->{lc($1)} };
    }

    $offset = 0 if (! defined($offset) or $offset !~ /^\d+\.?\d*/);

    print "Done\n";
}

sub showusage {
    return "usage: $0 -i input template file -f frozen atoms -s spacing info -c [constrainted plane] -o [optimize=yes]\n" . 
	"options:\n\t-i input template file: required. location of Jaguar input template file\n" . 
	"\t-f frozen atoms: reequied. atoms to keep constant. Expected 'a-b' for atoms a to b\n" .
	"\t-s spacing info: required. 'a-b:c' for a to b every c\n" . 
	"\t-c constrained plane: optional. '[x|y|z|xy|yz]'plane to constain unfrozen atoms a. Default xy\n" .
	"\t-o optimize: optional. (SinglePoint = 0 (default), GOpt = 1, CP = 2)\n" .
	"\t-a plane-plane offset: optional. Current distance seperating planes. Default 0\n";
}

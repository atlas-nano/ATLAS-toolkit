#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo);
use General qw(GetSelections FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(SplitAtomsByMol GetMols);

sub init;
sub getJagAtomData;
sub getJagInOpts;
sub createCPFiles;
sub getMols;
sub numerically { ($a<=>$b); }
sub writeFile;
sub updateAtomNames;

my ($bgfFile, $jagOutFile, $jagInFile, $sname);
my ($ATOMS, $BONDS, $HEADER, $mode, $MOLS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, undef) = GetBGFFileInfo($bgfFile, 1);
$MOLS = GetMols($ATOMS, $BONDS);
print "Done\n";
if ($mode == 1) {
    &getJagAtomData($jagOutFile, $ATOMS);
} else {
    &updateAtomNames($ATOMS);
}
print "Getting Jaguar input options from $jagInFile...";
$HEADER = getJagInOpts($jagInFile);
print "Done\nCreating $sname CP input and mae files...";
&createCPFiles($HEADER, $ATOMS, $MOLS, $sname);
print "Done\n";

sub writeFile {
    my ($filePrefix, $fileData) = @_;
    my ($maeCmd) = "jaguar babel -ijagin ${filePrefix}.in -omacmod ${filePrefix}.mae";

    open OUTDATA, "> ${filePrefix}.in" or die "ERROR: Cannot create ${filePrefix}.in: $!\n";
    print OUTDATA $fileData;
    close OUTDATA;

    die "ERROR while executing $maeCmd\n" if (system($maeCmd));
}

sub createCPFiles {
    my ($optHeader, $atomData, $molData, $prefix) = @_;
    my ($molName, $i, $outFile, $j, $outData, $atomName);

    for $i (1..2) {
	$molName = "${prefix}_" . lc(chr(64 + $i));

	# write the cp .in file
	$outData = "MAEFILE: " . $ENV{PWD} . "/${molName}_cp.mae\n\&gen\n${optHeader}" . 
	    "\&\nentry_name: $molName\n\&zmat\n";
	for $j (sort numerically keys %{ $atomData }) {
	    $atomName = $atomData->{$j}{ATMNAME};
	    $atomName =~ s/\s//;
	    $atomName .= "@" if (! exists($molData->{$i}{$j}));
	    $outData .= sprintf("%-5s%17.13f%17.13f%17.13f\n", $atomName,
				$atomData->{$j}{XCOORD}, $atomData->{$j}{YCOORD},
				$atomData->{$j}{ZCOORD});
	}
	$outData .= "\&\n";
	writeFile("${molName}_cp", $outData);

	#write the no_cp .in file
	$outData = "MAEFILE: " . $ENV{PWD} . "/${molName}_no_cp.mae\n\&gen\n${optHeader}" . 
	    "\&\nentry_name: $molName\n\&zmat\n";
	for $j (sort numerically keys %{ $atomData }) {
	    next if (! exists($molData->{$i}{$j}));
	    $atomName = $atomData->{$j}{ATMNAME};
	    $atomName =~ s/\s//;
	    $outData .= sprintf("%-5s%17.13f%17.13f%17.13f\n", $atomName,
				$atomData->{$j}{XCOORD}, $atomData->{$j}{YCOORD},
				$atomData->{$j}{ZCOORD});
	}
	$outData .= "\&\n";
	writeFile("${molName}_no_cp", $outData);
    }
}

sub updateAtomNames {
    my ($atomData) = $_[0];
    my ($i);

    for $i (keys %{ $atomData }) {
	$atomData->{$i}{ATMNAME} = substr($atomData->{$i}{FFTYPE}, 0, 1) . "$i";
    }
}

sub getJagInOpts {
    my ($inFile) = $_[0];
    my ($start, $optData);

    $start = 0;
    open INFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
  LINE: while (<INFILE>) {
      chomp;
      if ($_ =~ /^\&gen/) {
	  $start = 1;
      } elsif ($_ =~ /^\&/) {
	  last LINE;
      } elsif ($start) {
	  next if ($_ =~ /igeopt/);
	  $optData .= "$_\n";
      }
  }
    close INFILE;
    die "ERROR: No valid data read!\n" if (! $optData);
    return $optData;
}

sub getMols {
    my ($atomData) = $_[0];
    my ($molData, $count);

    $molData = SplitAtomsByMol($atomData,undef);
    for (keys %{ $molData }) {
	$count++;
    }

    die "ERROR: detected $count molecules. Expected 2. Aborting\n"
	if ($count != 2);
    return $molData;
}

sub getJagAtomData {
    my ($jagOut, $atomData) = @_;
    my ($inline, $is_geometry, $atomCounter, $hash_key, %Jaguar_Data, $is_valid, $i, $j);

    if (open(JAGOUT, $jagOut)) {
	print "Updating atom data from Jaguar output file $jagOut...";
	while (<JAGOUT>) {
	    chomp;
	    $inline = $_;
	    if ($inline =~ /(atom               x                 y                 z)|(final geometry)/) {
		$is_geometry = 1;
		$atomCounter = 0;
	    } elsif ($inline =~ /principal moments of inertia/) {
		$is_geometry = 0;
	    }elsif ($is_geometry) {
		if ($inline =~ /([A-Z]{1}[a-z]?)(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
		    $atomCounter++;
		    $hash_key = $atomCounter;

		    $Jaguar_Data{$hash_key}{"XCOORD"} = $3;
		    $Jaguar_Data{$hash_key}{"YCOORD"} = $4;
		    $Jaguar_Data{$hash_key}{"ZCOORD"} = $5;
		    $Jaguar_Data{$hash_key}{"ATMNAME"} = "${1}${2}";
		    $is_valid = 1;
		}
	    }
	}
	close JAGOUT;
	if ($is_valid) {
	    for $i (keys %{ $atomData }) {
		if (! exists($Jaguar_Data{$i})) {
		    $is_valid = 0;
		    last;
		}
	    }
	    if ($is_valid) {
		for $i (keys %{ $atomData }) {
		    for $j ("ATMNAME", "XCOORD", "YCOORD", "ZCOORD") {
			$atomData->{$i}{$j} = $Jaguar_Data{$i}{$j};
		    }
		}
	    } else {
		print "incompatible files";
	    }
	} else {
	    print "invalid file";
	}
    } else {
	print "unable to open file";
    }
    print "Done\n";
}

sub init {
    my (%OPTS);
    getopt('boi',\%OPTS);
    die "usage: $0 -i jaguar input file -b bgf file -o (jag output file) -s (save prefix)\n" 
	if (! exists($OPTS{i}) or ! exists($OPTS{b}));
    ($bgfFile, $jagOutFile, $jagInFile, $sname) = ($OPTS{b}, $OPTS{o}, $OPTS{i}, $OPTS{s});
    print "Initializing...";
    FileTester($jagInFile);
    FileTester($bgfFile);
    $mode = 0;
    $sname = $jagInFile if (! defined($sname));
    if (defined($jagOutFile)) {
	FileTester($jagOutFile);
	$mode = 1;
	$sname = $jagOutFile if ($sname ne $jagInFile);
    }
    print "Done\n";
    $sname = basename($sname);
    $sname =~ s/\.\w+$//;
}


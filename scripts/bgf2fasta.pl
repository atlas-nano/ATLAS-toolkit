#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo sortByRes);
use General qw(FileTester);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub getFileData;
sub printFastaSequence;

my ($ATOMS, $BONDS, $fastaFile, $FILES, $RES, $outHandle, $printOpt);

$|++;
&init;
&getFileData($FILES);
&printFastaSequence($ATOMS, $RES, $outHandle, $printOpt);
if (defined($fastaFile)) {
    print "Done\n";
    close $outHandle;
}

sub getFileData {
    my ($files) = $_[0];
    my ($i);

    for $i (0 .. $#{ $files }) {
	print "Parsing BGF file $files->[$i]...";
	($ATOMS->{$i}, $BONDS->{$i}, undef) = GetBGFFileInfo($files->[$i]);
	$RES->{$i} = sortByRes($ATOMS->{$i});
	print "Done\n";
    }

    if (defined($fastaFile)) {
	print "Writing fasta information to $fastaFile..";
	open FASTAFILE, "> $fastaFile" or die "ERROR: Cannot create $fastaFile: $!\n";
	$outHandle = \*FASTAFILE;
    } else {
	$outHandle = \*STDOUT;
    }
}

sub printFastaSequence {
    my ($atomsInfo, $resInfo, $PRINTHANDLE, $printLevel) = @_;
    my (%short, $i, $atom, @tmp, $count, $headers); 
    my ($k, $maxLen, $name, $j, $outStr, $maxK);

    %short = (
	      "ALA" => "A",
	      "CYS" => "C",
	      "CYX" => "C",
	      "ASP" => "D",
	      "GLU" => "E",
	      "PHE" => "F",
	      "GLY" => "G",
	      "HIS" => "H",
	      "HIE" => "H",
	      "HSD" => "H",
	      "HSE" => "H",
	      "HSP" => "H",
	      "ILE" => "I",
	      "LYS" => "K",
	      "LEU" => "L",
	      "MET" => "M",
	      "ASN" => "N",
	      "PRO" => "P",
	      "GLN" => "Q",
	      "ARG" => "R",
	      "SER" => "S",
	      "THR" => "T",
	      "VAL" => "V",
	      "TRP" => "W",
	      "TYR" => "Y" );
    $count = $maxLen = $maxK = 0;
    for $i (0 .. $#{ $FILES }) {
	$name = basename($FILES->[$i]);
	$name =~ s/\.\w+$//;
	$maxLen = length($name) if (length($name) > $maxLen);
	$headers->{$i} = $name;
    }
    $maxLen += 2;
    for $i (sort {$a cmp $b} keys %{ $headers }) {
	$k = $count = 0;
	$outStr->{$i}[$k] = sprintf("%-${maxLen}s", $headers->{$i});
	for $j (sort {$a<=>$b} keys %{ $resInfo->{$i} }) {
	    @tmp = keys %{ $resInfo->{$i}{$j}{ATOMS} };
	    $atom = $tmp[0];
	    if (exists($short{ $atomsInfo->{$i}{$atom}{RESNAME} })) {
		$outStr->{$i}[$k] .= $short{ $atomsInfo->{$i}{$atom}{RESNAME} };
	    } else {
		$outStr->{$i}[$k] .= "-";
	    }
	    $count++;
 	    $outStr->{$i}[$k] .=  " " if (($count % 10) == 0 and $printLevel > 0);
	    $k++ if (($count % 50) == 0);
	    $maxK = $k if ($k > $maxK);
	    $outStr->{$i}[$k] = sprintf("%-${maxLen}s",$headers->{$i}) if (($count % 50) == 0);
	}
    }
    for $i (0 .. $maxK) {
	if ($printLevel > 0) {
	    printf $PRINTHANDLE "%${maxLen}s", "";
	    for $j (1 .. 5) {
		printf $PRINTHANDLE "%10s ", ($j * 10 + ($i * 50));
	    }
	    print $PRINTHANDLE "\n";
	}
	for $j (keys %{ $headers }) {
	    next if ($#{ $outStr->{$j} } < $i);
	    print $PRINTHANDLE "$outStr->{$j}[$i]\n";
	}
    }
}

sub init {
    my (%OPTS, $fileStr);
    getopt('bfv',\%OPTS);
    die "usage: $0 -b bgfFile(s) -f [fasta file (optional)] -v [verbose = no]\n" if (! exists($OPTS{b}));
    ($fileStr, $fastaFile, $printOpt) = ($OPTS{b}, $OPTS{f}, $OPTS{v});
    print "Initializing...";
    while ($fileStr =~ /(\S+)/g) {
	if (-e $1 and -r $1 and -T $1) {
	    push @{ $FILES }, $1;
	}
    }
    die "ERROR: No valid bgf files found which seaching \"$fileStr\"!\n" if (! $FILES);
    $printOpt = 0 if (! defined($printOpt) or $printOpt !~ /^1|yes/i);
    $printOpt = 1 if ($printOpt =~ /^1|yes/i);
    print "Done\n";
}
    

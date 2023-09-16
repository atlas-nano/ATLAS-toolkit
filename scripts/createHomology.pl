#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use FileFormats qw(GetBGFFileInfo sortByRes);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub parseFastaFile;
sub testFileData;
sub mutateRes;
sub printFastas;

my ($bgfFile, $fastaFile, $saveFile);
my ($BGF, $FASTA, $RES, $CONVERTER);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($BGF, undef, undef) = GetBGFFileInfo($bgfFile, 1);
$RES = sortByRes($BGF);
print "Done\nParsing fasta file $fastaFile...";
$FASTA = parseFastaFile($fastaFile, $CONVERTER);
&testFileData($BGF, $RES, $FASTA);
print "Done\n";
&mutateRes($FASTA, $saveFile, $bgfFile);
&printFastas($FASTA);

sub printFastas {
    my ($convertData, $startBGF, $endBGF) = $_[0];
    my ($newFasta, $name1, $name2, $out1, $out2, $execStr); 
    my ($i, $count, $tot, $oldFasta, $start);

    $name1 = basename($startBGF);
    $name2 = basename($endBGF);
    $name1 =~ s/\.\w+$//;
    $name2 =~ s/\.\w+$//;
    $tot = length($name1);
    $tot = length($name2) if (length($name2) > $tot);
    @{ $oldFasta } = split /\w/, $convertData->{TEXT};

    $execStr = "/ul/tpascal/scripts/bgf2fasta.pl -b $endBGF";
    open FASTAFILE, "$execStr |" or die "ERROR: Cannot exeecute $execStr:$!\n";
    while (<FASTAFILE>) {
	chomp;
	if ($_ =~ /^Parsing BGF file/) {
	    $start = 1;
	} elsif ($start) {
	    while ($_ =~ /(A-Z)/g) {
		push @{ $newFasta }, $1;
	    }
	}
    }
    close FASTAFILE;
    die "ERROR: Invalid bgf file $endBGF\n" if (! defined($newFasta));
    $i = 0;
    $out1 = sprintf("%-${tot}s: ", $name1);
    $out2 = sprintf("%-${tot}s: ", $name2);
    while ($i <= $#{ $newFasta }) {
	$out1 .= $oldFasta->[$i];
	$out2 .= $newFasta->[$i];
	$i++;
	if (($i % 50) == 0) {
	    $out1 .= sprintf("\n%-${tot}s: ", $name1);
	    $out2 .= sprintf("\n%-${tot}s: ", $name2);
	} elsif (($i % 10) == 0) {
	    $out1 .= " ";
	    $out2 .= " ";
	}
    }
    print "FASTA Comparison\n${out1}\n${out2}\n";
}

sub mutateRes {
    my ($convertData, $outFile, $inFile) = @_;
    my ($i, $j, $printStr, $execStr);

    system("cp $inFile $outFile");
    
    for $i (keys %{ $convertData }) {
	next if ($i =~ /COUNT|TEXT/);
	$printStr = "";
	for $j (keys %{ $convertData->{$i} }) {
	    $execStr = "/ul/tpascal/scripts/amberProteinMutate.pl -f $outFile -s $outFile -t bgf -m $i -r \"Ir${j}\"";
	    print "Converting ${printStr} res $j (" . $convertData->{$i}{$j} . ") to $i\r";
	    if (system("${execStr} >& _homology.out")) {
		die "ERROR while executing \"$execStr\"\n";
	    }
	}
    }
    system("rm -fr _out");
}

sub testFileData {
    my ($bgfData, $resData, $fastaData) = @_;
    my ($i, $j, @resIDs, @atoms, $resName);

    die "ERROR: More residues specified in fasta than exists ib BGF!\n"
	if (scalar(keys %{ $RES }) < $fastaData->{COUNT});
    for $i (keys %{ $fastaData }) {
	next if ($i =~ /COUNT|TEXT/);
	@resIDs = (keys %{ $fastaData->{$i} });
	for $j (@resIDs) {
	    @atoms = (keys %{ $resData->{$j}{ATOMS} });
	    $resName = $bgfData->{$atoms[0]}{RESNAME};
	    if ($i eq $resName) {
		delete $fastaData->{$i}{$j};
	    } else {
		$fastaData->{$i}{$j} = $resName;
	    }
	}
    }
}

sub parseFastaFile {
    my ($inFile, $convertHash) = @_;
    my ($DATA, $resCount, $resName);

    open FASTAFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
    $resCount = 0;
    while (<FASTAFILE>) {
	chomp;
	next if ($_ !~ /^[A-Z]/);
	while ($_ =~ /([A-Z])/g) {
	    $resCount++;
	    $DATA->{COUNT} = $resCount;
	    $DATA->{TEXT} .= "$1";
	    next if (! exists($convertHash->{FASTA}{$1}));
	    $resName = $convertHash->{FASTA}{$1};
	    $DATA->{$resName}{$resCount} = 1;
	}
    }
    close FASTAFILE;
    die "ERROR: No valid info read!\n" if (! $DATA);
    return $DATA;
}

sub init {
    my (%OPTS, $resMap, $resInfo, $fastaName);
    getopt('bfsr',\%OPTS);
    for ("b", "f") {
	die "usage: $0 -b bgf file -f fasta file -s (save file) -r (residue name map)\n"
	    if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($bgfFile, $fastaFile, $resInfo, $saveFile) = ($OPTS{b}, $OPTS{f}, $OPTS{r}, $OPTS{s});
    FileTester($bgfFile);
    FileTester($fastaFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_mutated.bgf";
    }

    $CONVERTER->{BGF} = {
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
	"TYR" => "Y",
    };

    for (keys %{ $CONVERTER->{BGF} }) {
	$CONVERTER->{FASTA}{$CONVERTER->{BGF}{$_}} = $_;
    }

    if (defined($resInfo)) {
	while ($resInfo =~ /(\w+) \-\> (\w+)/g) {
	    $fastaName = $CONVERTER->{BGF}{$2};
	    $CONVERTER->{FASTA}{$fastaName} = uc $2 if (exists($CONVERTER->{FASTA}{$fastaName}));
	}
    }

    print "Done\n";
}

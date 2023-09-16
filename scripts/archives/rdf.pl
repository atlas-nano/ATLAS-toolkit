#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use warnings;
no warnings "recursion";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetBGFFileInfo GetBondList AddMass);
use Packages::General qw(STDev FileTester CoM GetSelections IsInteger 
			 IsDecimal TrjSelections CenterOnMol GetBondLength);
use Packages::AMBER qw(ParseAmberTrj GetAmberByteOffset);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType);
use Packages::BOX qw(GetBox ConvertBox MakeBox);
use Packages::ManipAtoms qw(UnwrapAtoms ImageAtoms GetAtmList);
use Packages::CERIUS2 qw(LoadFFs);
use constant PI => atan2(1,1) * 4;

sub init;
sub calcRDF;
sub saveFile;
sub getAtoms;
sub boxConvert;
sub saveFile;
sub numerically { ($a<=>$b); }
sub showUsage;

my ($SOLUTE, $SOLVENT, $trjFile, $BGF, $BONDS, $rMax, $printStr);
my ($getSnapshot, $getByteOffset, $trjType, $LAMMPSOPTS, $delR);
my ($BBOX, $reImage, $SELECT, $saveFile, $field, $DATA, $PARMS);

$|++;
&init;

if (! $trjType) {
    print "Calculating RDF...";
    &calcRDF($BGF, $BBOX, 1, undef);
    print "Done\n";
} else {
    $field = scalar keys %{ $BGF };
    $getByteOffset->($SELECT, $trjFile, $field);
    if ($trjType == 2) {
	&GetLammpsTrjType($SELECT, $trjFile, "atom", \%{ $LAMMPSOPTS });
	$field = "atom";
    }
    $printStr = "Calculating RDF from $trjFile...";
    $getSnapshot->($BGF, $trjFile, $SELECT, $field, \&calcRDF, $printStr, undef);
}

print "Saving RDF data to $saveFile...";
saveFile($DATA, $saveFile);
print "Done\n";

sub saveFile {
    my ($data, $saveName) = @_;
    my ($i, $j, $tot, $STATS, $norm); 
    my ($avg, $stdev, $rl, $ru, $avgFile);

    $tot = scalar keys %{ $SOLVENT };
    open OUTFILE, "> $saveName" or die "ERROR: Cannot create file $saveName: $!\n";
    for $i (sort numerically keys %{ $data }) {
	print OUTFILE "\#Snapshot \# $i\n";
	for $j (0 .. $#{ $data->{$i}{VALS} }) {
	    $rl = $j*$delR;
	    $ru = $rl + $delR;
	    $norm = $data->{$i}{CONST} * $tot * $tot * ($ru*$ru*$ru - $rl*$rl*$rl);
	    $data->{$i}{VALS}[$j] /= $norm;
	    print OUTFILE "$rl $data->{$i}{VALS}[$j]\n";
	    $STATS->[$j] .= "$data->{$i}{VALS}[$j] ";
	}
	print OUTFILE "\n";
    }

    close OUTFILE;
    if ($saveName =~ /(.*)\.(.+)$/) {
	$avgFile = "${1}_avg.${2}";
    } else {
	$avgFile = "${saveName}_avg.rdf";
    }

    open OUTFILE, "> $avgFile" or die "ERROR: Cannot create $avgFile:$!\n";
    print OUTFILE "\#AVERAGE\n";
    for $j (0 .. $#{ $STATS }) {
	chomp $STATS->[$j];
	($avg, $stdev, $i) = STDev($STATS->[$j]);
	print OUTFILE ($j*$delR) . " $avg $stdev\n";
    }
    close OUTFILE;
}

sub getAtoms {
    my ($allAtoms, $atomList) = @_;
    my (%ATOMS, $i);

    for $i (keys %{ $atomList }) {
        $ATOMS{$i} = $allAtoms->{$i};
    }

    return \%ATOMS;
}

sub calcRDF {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($MOL, $CENTER, @tmp, $i, $j, $SOLCENTER); 
    my ($bin, $dist, $nBins, $count, $prec, $rounding);

    $nBins = $rMax/$delR; # total number of bins
    while ($delR =~ /(\d|.)/g) {
	if ($1 ne ".") {
	    if ($1 > 0) {
		$prec .= "1";
	    } else {
		$prec .= "0";
	    }
	} else {
	    $prec .= ".";
	}
    }

    $rounding = $prec;
    chop $rounding;
    $rounding .= "0";

    if ($trjType == 2) { #LAMMPS
	$frameNum = $ATOMS->{TIMESTEP}[0];
        $BOX = ConvertBox($ATOMS->{"BOX BOUNDS"});
        $BBOX = MakeBox($BOX);
        $ATOMS = $ATOMS->{ATOMS};
        UnwrapAtoms($ATOMS, $BBOX, $LAMMPSOPTS->{scaled});
        @tmp = grep {!/COORD/i} keys %{ $BGF->{1} };
        for $i (keys %{ $ATOMS }) {
            for $j (@tmp) {
                next if ($j =~ /COORD/i);
                $ATOMS->{$i}{$j} = $BGF->{$i}{$j};
            }
        }	
    } elsif ($trjType == 1) {
	$BBOX = MakeBox($BOX);
    } else {
	$BBOX = $BOX;
    }
    
    for $i (0 .. $nBins) {
        $DATA->{$frameNum}{VALS}[$i] = 0; #initialize all elements of array
    }

    $MOL = getAtoms($ATOMS, $SOLUTE);
    $SOLCENTER = CoM($MOL);

    if ($reImage) {
        CenterOnMol($ATOMS, $SOLCENTER);

        for $i (keys %{ $SOLVENT }) {
            $MOL = getAtoms($ATOMS, $SOLVENT->{$i});
            $CENTER = CoM($MOL);
            ImageAtoms($ATOMS, $SOLVENT->{$i}, $CENTER, $BBOX);
        }

    }

    $MOL = getAtoms($ATOMS, $SOLUTE);
    $SOLCENTER = CoM($MOL);

    $count = 0;
    for $i (keys %{ $SOLVENT }) {
	$MOL = getAtoms($ATOMS, $SOLVENT->{$i});
	$CENTER = CoM($MOL);
	$dist = sprintf("%${prec}f",GetBondLength($SOLCENTER, $CENTER));
	next if ($dist > $rMax);
	$count++;
	$bin = sprintf("%${rounding}f",($dist/$delR) + 1);
	$DATA->{$frameNum}{VALS}[$bin]++;
    }
    $DATA->{$frameNum}{TOTAL} = $count;
    $DATA->{$frameNum}{CONST} =  4.0*PI / 
	(3*$BBOX->{XCOORD}{len}*$BBOX->{YCOORD}{len}*$BBOX->{ZCOORD}{len});
}

sub init {
    my (%OPTS, @tmp, $i, $j, $list, $solvTmp, $HEADERS);
    my ($solSel, $solvSel, $bgfFile, $tSel, $usage, $FF);
    getopt('mvbtdrlswif',\%OPTS);

    $usage = &showUsage;

    for ($OPTS{m},$OPTS{v},$OPTS{b}) {
	die "$usage" if (! defined($_));
    }

    print "Initializing...";
    ($solSel,$solvSel,$trjFile,$bgfFile,$rMax,$delR,$trjType,$tSel,$saveFile,$reImage,$FF) = 
	($OPTS{m},$OPTS{v},$OPTS{t},$OPTS{b},$OPTS{r},$OPTS{d},$OPTS{l},$OPTS{s},$OPTS{w},$OPTS{i},$OPTS{f});

    if (defined($trjFile) and -e $trjFile and -r $trjFile and -T $trjFile) {
	if (! defined($trjType) or lc($trjType) ne "lammps") {
	    $trjType = 1;
	    $getSnapshot = \&ParseAmberTrj;
	    $getByteOffset = \&GetAmberByteOffset;
	} else {
	    $trjType = 2;
	    $getSnapshot = \&ParseLAMMPSTrj;
	    $getByteOffset = \&GetLammpsByteOffset;
	}
	if (! defined($tSel)) {
	    $tSel = "*";
	}
	$list = TrjSelections($tSel);
	for $j (keys %{ $list }) {
	    $SELECT->{$j} = $list->{$j};
	}
	
	die "ERROR: No valid frames selected with selection $tSel!\n"
	    if (! keys %{ $SELECT } and $tSel ne "*");
    } else {
	$trjType = 0;
    }
	
    FileTester($bgfFile);

    if ($solSel =~ /\s+/) {
	@tmp = split /\s+/, $solSel;
    } else {
	$tmp[0] = $solSel;
    }
    $list = GetSelections(\@tmp, 0);
    for $j (keys %{ $list }) {
	$SOLUTE->{$j} = $list->{$j};
    }
    
    die "ERROR: No valid atoms selected for SOLUTE with selection $solSel!\n"
	if (! keys %{ $SOLUTE });

    
    if ($solvSel =~ /\s+/) {
	@tmp = split /\s+/, $solvSel;
    } else {
	$tmp[0] = $solvSel;
    }
    $list = GetSelections(\@tmp, 0);
    for $j (keys %{ $list }) {
	$solvTmp->{$j} = $list->{$j};
    }

    die "ERROR: No valid atoms selected for SOLVENT with selection $solvSel!\n"
	if (! keys %{ $solvTmp });

    $delR = 0.1 if (! defined($delR) or (! IsInteger($delR) and ! IsDecimal($delR)));

    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$/\.rdf/;
    }
    $reImage = 0 if (! defined($reImage) or $reImage !~ /^1/);
    print "Done\nParsing BGF file $bgfFile...";
    ($BGF, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
    $SOLUTE = GetAtmList($SOLUTE, $BGF);
    $solvTmp = GetAtmList($solvTmp, $BGF);
    $SOLVENT = GetBondList($BGF, $BONDS, $solvTmp);
    die "ERROR: Solvent atoms does not correspond to any in BGF file!\n"
	if (! keys %{ $SOLVENT });
    $BBOX =GetBox($BGF, undef, $HEADERS);
    &boxConvert($BBOX);
    if (! defined($rMax) or (! IsInteger($rMax) and ! IsDecimal($rMax))) {
	$rMax = 0;
	for $i ("X", "Y", "Z") {
	    $rMax = $BBOX->{$i}{len} if ($BBOX->{$i}{len} > $rMax);
	}
    }
	
    print "Done\n";
    if (exists($OPTS{f}) and -e $FF) {
	print "Parsing CERIUS2 forcefield $FF...";
	if ($FF =~ /\s/) {
	    @tmp = split /\s+/, $FF;
	} else {
	    $tmp[0] = $FF;
	}
	
	$PARMS = LoadFFs(\@tmp);
	AddMass($BGF, $PARMS);
    }
}

sub boxConvert {
    my ($box) = $_[0];
    
    for ("X", "Y", "Z") {
	%{ $box->{$_ . "COORD"} } = %{ $box->{$_} };
    }
}

sub showUsage {
    my ($usage) = "usage: $0 -m solute atoms -v solvent atoms -b bgf file " . 
	"\nOptional parameters:" . 
	"\n\t-w save name -r max distance -d del r -f cerius2 forcefield" . 
	"\n\t-t trajectory -s traj selection -l traj type (amber(default)|lammps)\n";
    return $usage;
}

#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General;
use Packages::FileFormats;
use Packages::LAMMPS;
use Packages::AMBER;

sub numerically;
sub parseLammpsData;
sub createAmberTrj;
sub createLammpsTrj;
sub unwrapAtoms;
sub determinePos;
sub init;
sub sortByRes;
sub determineIfScaled;

die "usage: $0 lammps_traj starting_bgf (snapshot: start stop end) [save_name] [out_type=lammps]\n"
    if (! @ARGV or @ARGV < 3);

my ($lammpsTrj, $bgfFile, $selection, $saveName, $outType) = @ARGV;
my ($printAtoms, $SELECT, $scaled);

init;
$|++;

print "Parsing BGF file $bgfFile...";
my ($BGF, $BONDS) = GetBGFFileInfo($bgfFile, 0);
my ($atmBondList) = sortByRes($BGF, $BONDS);
print "Done\n";

$scaled = determineIfScaled($lammpsTrj);

my $pStr = "Creating $outType trajectory $saveName...";
open OUTDATA, "> $saveName" or die "Cannot write to trajectory file $saveName: $!\n";
print OUTDATA "TITLE: LAMMPS Trajectory file create by lammpsTrj2amber on " . scalar(localtime(time())) . "\n";
parseLammpsData($BGF, $lammpsTrj, $atmBondList, \*OUTDATA, $pStr, $SELECT);
close OUTDATA or die "ERROR: Cannot close $saveName: $!\n";

sub determineIfScaled {
    my ($lammpsFile) = $_[0];
    my ($scaled, $counter, $inStr);

    $counter = $scaled = 1;

    open LAMMPS, $lammpsFile or die "ERROR: Cannot open LAMMPS trajectory file $lammpsFile: $!\n";
    while (<LAMMPS>) {
	chomp;
	$inStr = $_;
	if ($inStr =~ /^\d+\s\d+\s(\S+)\s(\S+)\s(\S+)/) {
	    if ($1 > 2 || $2 > 2 || $3 > 2) {
		$scaled = 0;
		last;
	    }
	    $counter++;
	}
	last if ($counter > 5);
    }

    close LAMMPS or die "ERROR: Cannot close LAMMPS trajectory file $lammpsTrj: $!\n";

    return $scaled;
}

sub sortByRes {
    my ($ATOMS, $BONDS) = @_;
    my ($atomC, $atom, $resAtms, $resC, %RESINFO, $i);

    $resC = 1;
    for $atomC (keys %{ $ATOMS }) {
	$atom = \%{ $ATOMS->{$atomC} };
	next
	    if (exists($atom->{"bondpointer"}));
	$resAtms = GetResAtoms($atomC, $BONDS, ());
	$atom->{"bondpointer"} = $resC;
	$RESINFO{$resC}{$atomC} = 1;
	for $i (keys %{ $resAtms }) {
	    $ATOMS->{$i}{"bondpointer"} = $resC;
	    $RESINFO{$resC}{$i} = 1;
	}
	$resC++;
    }

    return \%RESINFO;
}

sub parseLammpsData {
    my ($bgfInfo, $lammpsFile, $resinfo, $outFile, $printStr) = @_;
    my ($inStr, $BOX, $counter, $ATOMS, @dim, $validTStep, $index, $snapC);
    my ($atomC, $HEADERS, $currPos, $progress, $start, $strLen, $tot, $endFrame);
    my ($filesize) = -s $lammpsFile;

    
    $ATOMS = $BOX = $HEADERS = ();
    $counter = $validTStep = $index = $snapC = 0;
    $start = 0;

    $strLen = 0;
    if (defined($SELECT)) {
	@dim = sort numerically keys %{ $SELECT };
	$tot = $#dim + 1;
	$endFrame = pop @dim;
    } else {
	$tot = 0;
    }

    @dim =  ("X", "Y", "Z", "Xi", "Yi", "Zi");
    open LAMMPS, $lammpsFile or die "ERROR: Cannot open LAMMPS trajectory file $lammpsFile: $!\n";
    while (<LAMMPS>) {
	chomp;
	$inStr = $_;
	if ($inStr =~ /^ITEM: TIMESTEP$/) { # either timestep of #atoms
	    if ($ATOMS and $BOX and $HEADERS) {
		#print "SAVING\n";
		$printAtoms->($ATOMS, $BOX, $resinfo, $HEADERS, $outFile);
		$counter = 0;
		$index++;
		if (! $tot) {
		    $currPos = tell(LAMMPS);
		    $strLen = PrintProgress($currPos, $filesize, $start, $printStr);
		} else {
		    $strLen = PrintProgress($index, $tot, $start, $printStr);
		    if ($snapC == $endFrame) {
			$BOX = $ATOMS = $HEADERS = ();
			$validTStep = 0;
			last;
		    }
		}
	    }
	    $BOX = $ATOMS = $HEADERS = ();
	    $validTStep = 0;
	    $counter = 0;
	    #print "$inStr\n";
	    push @{ $HEADERS }, $inStr;
	    $snapC++;
	    if (! $tot || exists($SELECT->{$snapC})) {
		$validTStep = 1;
		if (! $start) {
		    $start = time();
		}
	    }
	} elsif ($validTStep) {
	    if ($inStr =~ /^(\-?\d+\.?\d*) (\-?\d+\.?\d*)$/) { #box info
		push @{ $HEADERS }, $inStr;
		$BOX->{$dim[$counter]} = (
					  {
					      "lo" => $1,
					      "hi" => $2,
					  }
					  );
		$counter++;
	    } elsif ($inStr =~ /^(\d+)\s(\d+)(\s\-?\d+.+)/) {
		$atomC = $1;
		$ATOMS->{$atomC}{"bondpointer"} = $bgfInfo->{$atomC}{"bondpointer"};
		$ATOMS->{$atomC}{"TYPE"} = $2;
		$counter = 0;
		$inStr = $3;
		while ($inStr =~ /\s(\-?\d+\.?\d*e?\-?\d*)/g and $counter < 6) {
		    $ATOMS->{$atomC}{$dim[$counter]} = $1;
		    $counter++;
		}
	    } else {
		#print "$inStr\n";
		push @{ $HEADERS }, $inStr;
	    }
	}
    }
    close LAMMPS or die "ERROR: Cannot close $lammpsFile: $!\n";

    if ($ATOMS and $BOX) {
	$printAtoms->($ATOMS, $BOX, $resinfo, $HEADERS, $outFile);
    }
    printf "$printStr%-${strLen}s\n", "Done";
}

sub createLammpsTrj {
    my ($ATOMS, $BOX, $bondList, $HEADERS, $OUTFILE) = @_;
    my ($i, $atomC, $atom, $dim, $index, $j, $bLen);

    for $dim ("X", "Y", "Z") {
	$bLen->{$dim} = $BOX->{$dim}{"hi"} - $BOX->{$dim}{"lo"};
    }

    for $i (@{ $HEADERS }) {
	print $OUTFILE "$i\n";
    }

    unwrapAtoms($ATOMS, $BOX, $bLen);

    for $atomC (keys %{ $ATOMS }) {
	$atom = \%{ $ATOMS->{$atomC} };
	print $OUTFILE "$atomC " . $atom->{"TYPE"} . " ";
	if (! exists($atom->{"imaged"})) {
	    $j = $atom->{"bondpointer"};
	    determinePos($ATOMS, $atomC, $BOX, \%{ $bondList->{$j} }, $bLen);
	}
	$index = "";
	for $dim ("X", "Y", "Z") {
	    $index .= $atom->{$dim . "i"} . " ";
	    printf $OUTFILE "%.5f ", $atom->{$dim};
	}
	chop $index;
	printf $OUTFILE "$index\n";
    }
    
}    

sub unwrapAtoms {
    my ($ATOMS, $BOX, $bLen) = @_;
    my ($atomC, $atom, $dim);


    for $atomC (keys %{ $ATOMS }) {
	$atom = \%{ $ATOMS->{$atomC} };
	for $dim ("X", "Y", "Z") {
	    $atom->{$dim} +=  $atom->{$dim . "i"};
	    if ($scaled) {
		#$atom->{$dim} = ($atom->{$dim} * $bLen->{$dim}) + $BOX->{$dim}{"lo"};
	    }
	}
    }
}

sub determinePos {
    my ($ATOMS, $i, $BOX, $RES, $bLen) = @_;
    my ($dim, $resIndex, $atom, $atomC, $index);


    for $dim ("X", "Y", "Z") {
	$resIndex->{$dim} = 0; # set the dimension index of the residue to 0
	                       # now search through all the atoms in the residue
                               # and determine the smallest index required to be
                               # in the box
	for $atomC (keys %{ $RES }) {
	    $atom = \%{ $ATOMS->{$atomC} };
	    $index = $atom->{$dim . "i"}; # get the integer offset
	    if ($index != 0) {
		if ($index > 0) { # if the atom coordinate lies to right of box
		    $resIndex->{$dim} = $index
			if ($resIndex->{$dim} > $index || $resIndex->{$dim} == 0); # make the residue index the atom index 
                                                                                   # if it is the smallest
		} else {
		    $resIndex->{$dim} = $index
			if ($resIndex->{$dim} < $index || $resIndex->{$dim} == 0);
		}
	    }
	    $atom->{"imaged"} = 1;

	}
	for $atomC (keys %{ $RES }) {
	    $index = $atom->{$dim . "i"};	    
	    $index -= $resIndex->{$dim};
	    $atom = \%{ $ATOMS->{$atomC} };
	    #$atom->{$dim} -= $index * $bLen->{$dim} - $BOX->{$dim}{"lo"};
	    #$atom->{$dim} /= $bLen->{$dim};
	}
    }
}
	    	
sub createAmberTrj {
    my ($ATOMS, $BBOX, $bondList, $HEADERS, $OUTFILE) = @_;
    my ($atomC, $counter, @tmp, $dim, $BOX);

    $counter = 0;
    @tmp = sort numerically keys %{ $ATOMS };
    if (defined($BBOX)) {
	for $dim ("X", "Y", "Z") {
	    $BOX->{$dim} = $BBOX->{$dim}{"hi"} - $BBOX->{$dim}{"lo"};
	}
    }

    for $atomC (@tmp) {
	for $dim ("X", "Y", "Z") {
	    $counter++;
	    if ($scaled and defined($BBOX)) {
		$ATOMS->{$atomC}{$dim} *= $BOX->{$dim}
	    }
	    if (exists($ATOMS->{$atomC}{$dim . "i"}) and defined($BBOX)) {
		printf $OUTFILE "%8.3f", $ATOMS->{$atomC}{$dim} + ($BOX->{$dim} * $ATOMS->{$atomC}{$dim . "i"});
	    } else {
		printf $OUTFILE "%8.3f", $ATOMS->{$atomC}{$dim};
	    }
	    if ($counter == 10) {
		print $OUTFILE "\n";
		$counter = 0;
	    }
	}
    }
    if (defined($BBOX)) {
	for $dim ("X", "Y", "Z") {
	    printf $OUTFILE "%8.3f", $BOX->{$dim};
	}
	print $OUTFILE "\n";
    }
}

sub numerically {
    ($a<=>$b);
}

sub init {
    my ($i);

    FileTester($lammpsTrj);
    FileTester($bgfFile);
    
    if (! defined($scaled)) {
	$scaled = 0;
    }

    if (! defined($outType)) {
	$outType = "lammps";
    }
    

    $scaled = Trim($scaled);
    if (! IsInteger($scaled) or $scaled !~ /^1|0$/) {
	print "WARNING: Expected boolean for scaled. Got $scaled.. \nSetting to default value of 0\n";
    }

    if ($outType !~ /(lammps|amber)/i) {
	print "WARNING: Expected lammps or amber for out_type. Got $outType...\n" . 
	    "Setting to default value: lammps\n";
	$outType = "lammps";
    } else {
	$outType = $1;
    }

    if (! defined($saveName)) {
	$saveName = $lammpsTrj;
    }

    $outType = Trim(lc($outType));

    if ($outType eq "lammps") {
	$printAtoms = \&createLammpsTrj;
	$saveName =~ s/\.\w+$/\.lammpstrj/;
    } else {
	$printAtoms = \&createAmberTrj;
	$saveName =~ s/\.\w+$/\.crd/;
    }
    $outType = uc($outType);
    
    print "SNAPSHOT SELECTION: ";
    if (Trim($selection) =~ /^(\d+)\s+(\d+)\s+(\d+)/) {
	print "start: $1 end: $2 every: $3...";
	$i = $1;
	while ($i < $2) {
	    $SELECT->{$i} = 1;
	    $i += $3;
	}
    } else {
	$SELECT = ();
	print "using all";
    }
    print "\n";
}


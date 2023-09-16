#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use warnings;
no warnings "recursion";
use Packages::FileFormats qw(GetBGFFileInfo sortByRes addHeader createBGF AddMass GetBondList GetBGFAtoms);
use Packages::General qw(FileTester IsInteger GetBondLength CoM GetSoluteAtoms LoadFFs);
use Packages::BOX qw(GetBox CreateGrid GetSurface GetRadii GetNeighbours);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use Packages::ManipAtoms qw(GetMols);

sub init;
sub updateBGF;
sub numerically { ($a<=>$b); };
sub findWatShells;
sub createFile;
sub addBulk;
sub createShellBGFs;
sub getOffset;
sub showUsage;
sub getShakeDOF;

my ($bgfFile, $cerius2FF, $shellStr, $saveName);
my ($ATOMS, $BONDS, $HEADERS, $PARMS, $BOX, $GRID, $BBOX, $atm_counter, $SOLUTE);
my ($SURFACE, $molCharge, $watShells, $increment, $RES, $WSHELLS, $bList, $OFFSET);

$|++;
$watShells = &init;
print "Getting atom data from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$bList = GetMols($ATOMS, $BONDS);
$SOLUTE = GetSoluteAtoms($ATOMS, $bList);
$RES = sortByRes($ATOMS);
print "Done\nGetting forcefield parameters from cerius2 file $cerius2FF...";
$PARMS = LoadFFs($cerius2FF);
AddMass($ATOMS, $PARMS);
print "Done\nGetting system dimensions...";
$BOX = GetBox($ATOMS, $PARMS, $HEADERS);
$OFFSET = getOffset($ATOMS, $BOX);
updateBGF($ATOMS, $PARMS, $OFFSET, 1);
print "Done\nCreating grid...\n";
($GRID, $BBOX, $atm_counter) = CreateGrid($ATOMS, 0, $BOX, $increment, 1);
print "\n";
$WSHELLS = findWatShells($SOLUTE, $GRID, $watShells, $RES, $ATOMS);
addBulk($WSHELLS, $ATOMS);
print "Creating $saveName...";
createFile($WSHELLS, $saveName);
updateBGF($ATOMS, $PARMS, $OFFSET, -1);
print "Done\n";
createShellBGFs($ATOMS, $WSHELLS, $saveName);

sub createShellBGFs {
    my ($atoms, $shells, $savePrefix) = @_;
    my ($MOLECULE, $CONS, $atmList, $i, $bgfName, @tmp, $atom, $printStr);

    $savePrefix =~ s/\.\w+$//;
    @tmp = sort numerically keys %{ $shells };
    $printStr = "Creating BGF file";

    for $i (0 .. $#tmp) {
	$MOLECULE = ();
	$CONS = ();
	$atmList = ();
	if ($i == 0) {
	    $bgfName = "${savePrefix}_solute.bgf";
	} elsif ($i == 1) {
	    $bgfName = "${savePrefix}_ions.bgf";
	}elsif ($i == $#tmp) {
	    $bgfName = "${savePrefix}_bulk.bgf";
	} else {
	    $bgfName = "${savePrefix}_shell" . ($i - 1) .".bgf";
	}
	print "${printStr} ${bgfName}...\r";
	for $atom (keys %{ $shells->{$tmp[$i]} }) {
	    $atmList->{$atom} = 1;
	}
	($MOLECULE, $CONS, undef) = GetBGFAtoms($atmList, $ATOMS, $BONDS);
	addHeader($MOLECULE, $HEADERS);
	createBGF($MOLECULE, $CONS, $bgfName);
    }
    printf "${printStr}s...Done%50s\n", "";

}

sub addBulk {
    my ($shells,$atoms) = @_;
    my ($i, $bulkC, @tmp, $found, $j);

    if (keys %{ $shells }) {
	@tmp = sort numerically keys %{ $shells };
	$bulkC = $tmp[$#tmp] + 1;
    } else {
	$bulkC = 2;
    }

    for $i (keys %{ $atoms }) {
	next if ($atoms->{$i}{VISITED});
	if ($atoms->{$i}{RESNAME} =~ /WAT/) {
	    $shells->{$bulkC}{$i} = $atoms->{$i}{RESNUM}; #bulk waters
	} else {
	    $shells->{1}{$i} = $atoms->{$i}{RESNUM}; #ions
	}
    }
}
	
sub createFile {
    my ($wshells, $save) = @_;
    my (@tmp, @tmp2, $i, $j, $start, $prev, $groups, $counter, $dofStr);
    
    @tmp = keys %{ $wshells };
    open OUTDATA, "> $save" || die "ERROR: Cannot create $save: $!\n";
    print OUTDATA "Total Groups: " . ($#tmp + 1) . "\n";

    @tmp = sort numerically keys %{ $wshells };
    $counter = 0;
    for $i (0 .. $#tmp) {
	@tmp2 = sort numerically keys %{ $wshells->{$tmp[$i]} };
        $dofStr .= getShakeDOF($wshells->{$tmp[$i]}, $ATOMS) . " ";
	print OUTDATA "Group " . ($i + 1) . " Atoms " . ($#tmp2 + 1) . "\n";
	$start = $prev = -1;
	$groups = "";
	$counter = 0;
	for $j (@tmp2) {
	    if ($start == -1) {
		$start = $j;
	    } elsif (($j - $prev) > 1) {
		$groups .= "${start} - ${prev} ";
		$start = $j;
		$counter++;
	    }
	    $prev = $j;
	    if ($counter == 10) {
		$counter = 0;
		$groups = substr($groups, 0, -2);
		$groups .= "\n";
	    }
	}
	$groups .= "${start} - ${prev} ";
	$groups = substr($groups, 0, -1);
	print OUTDATA "$groups\n";
    }
    print OUTDATA "Constraints\n$dofStr\n";
    close OUTDATA;
}

sub findWatShells {
    my ($soluteAtms, $grid, $shells, $resData, $atoms) = @_;
    my ($i, $start, $counter, $currCell, $NEIGHS, $x, $y, $z); 
    my ($water, $neighCell, $dist, %WATERS, $atom, $k, @tmp);

    $start = 0;
    for $i (0 .. $#{ $shells }) {
	print "Getting water molecules in from $start - $shells->[$i]...";
	for $counter (keys %{ $soluteAtms }) {
	    $WATERS{0}{$counter} = $atoms->{$counter}{RESNUM};
	    $atoms->{$counter}{VISITED} = 1;
	    $x = $atoms->{$counter}{CELL}{XINDEX};
	    $y = $atoms->{$counter}{CELL}{YINDEX};
	    $z = $atoms->{$counter}{CELL}{ZINDEX};
	    $currCell = $grid->{$x}{$y}{$z};
	    $NEIGHS = GetNeighbours($grid, $currCell, ($i + 1));
	    for $neighCell (@{ $NEIGHS }) {
		next if (! exists($neighCell->{WATERS}));
		for $water (@{ $neighCell->{WATERS} }) {
		    @tmp = keys %{ $resData->{ $water->{RESNUM} }{ATOMS} };
		    next if (exists($atoms->{ $tmp[0] }{VISITED}));
		    for $atom (@{ $currCell->{ATOMS} }) {
			for $k (@tmp) {
			    $dist = GetBondLength($atom, $atoms->{$k});
			    if ($dist >= $start && $dist <= $shells->[$i]) {
				for $k (@tmp) {
				    $WATERS{($i + 2)}{$k} = $water->{RESNUM};
				    $atoms->{$k}{VISITED} = 1;
				}
				last;
			    }
			}
		    }
		}
	    }
	}
	print "found " . scalar(keys %{ $WATERS{($i + 2)} })/3 . " molecules\n";
	$start = $shells->[$i];
    }
    
    return (\%WATERS);
}
    
sub getOffset {
    my ($atoms, $box) = @_;
    my (%OSET, $CENTER, $dim);

    $CENTER = CoM($atoms);
    for $dim ("X", "Y", "Z") {
	$OSET{$dim} = (($box->{$dim}{lo} + $box->{$dim}{hi})/2) - $CENTER->{$dim . "COORD"};
    }
    
    return \%OSET;
}

sub updateBGF {
    my ($atoms, $parms, $OFFSET, $direction) = @_;
    my ($counter, $dim);

    for $counter (keys %{ $atoms }) {
	$atoms->{$counter}{RADII} = GetRadii($ATOMS->{$counter}, $PARMS);
	for $dim ("X", "Y", "Z") {
	    $atoms->{$counter}{$dim . "COORD"} += ($direction * $OFFSET->{$dim});
	}
    }
}

sub init {
    my (@tmp, $i, $j, %SHELL, $counter, %OPTS, $usage, $FF);

    getopt('bfsw',\%OPTS);

    ($bgfFile, $FF, $shellStr, $saveName) = ($OPTS{b},$OPTS{f},$OPTS{w},$OPTS{s});
    $usage = &showUsage;
    for ("b", "f", "w") {
	die "$usage\n" if (! $OPTS{$_});
    }
    
    print "Initializing...";
    FileTester($bgfFile);
    
    if ($shellStr =~ /\s+/) {
	@tmp = split /\s+/, $shellStr;
    } else {
	$tmp[0] = $shellStr;
    }

    $counter = 1;
    $increment = 0;
    for $j (@tmp) {
	if ($j =~ /(\d+\.*\d+)/) {
	    $SHELL{$counter} = $1;
	    if ($counter == 1) {
		$increment = $1;
	    } elsif (abs($SHELL{$counter} - $SHELL{($counter - 1)}) > $increment) {
		$increment = abs($SHELL{$counter} - $SHELL{($counter - 1)});
	    }
	    $counter++;
	}
    }
    die "Expected at least one decimal(integer) value for the shell limts. Got \"$shellStr\"\n" if (! keys %SHELL);
    
    if (! defined($saveName)) {
	$saveName = basename($bgfFile);
	$saveName =~ s/\.\w+$/_watshell\.dat/;
    }

    if ($FF =~ /\s/) {
        @tmp = split /\s+/, $FF;
    } else {
        $tmp[0] = $FF;
    }
    $cerius2FF = ();
    for $i (@tmp) {
        push @{ $cerius2FF }, $i if (-e $i);
    }

    @tmp = sort numerically values %SHELL;
    print "grid length: $increment, total watershells: " . ($#tmp + 1) . "...Done\n";
    return \@tmp;
}

sub getShakeDOF {
    my ($select, $atoms) = @_;
    my ($i, $dof);
    $dof = 0;
    for $i (keys %{ $select }) {
        if ($atoms->{$i}{FFTYPE} =~ /^H/) { #if this is a hydrogen
            $dof++;
            if (${ $atoms->{$i}{MOLSIZE} } == 3 ) {
                $dof += 0.5;
            }
        }
    }
    return $dof;
}

sub showUsage {
    my ($usage) = "This script will break up the bgf file into:  solute, ions, water shell 1...n, and bulk\n\n" . 
	"usage: $0 -b bgf file -f force field file -w water shell distance(s) -s [save name]\n" .
	"\toptions:\n" .
	"\t-b bgf file: the location of the bgf file\n" .
	"\t-f force field file: the location of the CERIUS2 formatted force field file.\n" .
	"\t    the bgf file above must correctly typed with this force field\n" . 
	"\t-w water shell distances(s): the distance from the surface of the 1st water shell.\n" .
	"\t    you can specify multiple water shell distances by enclosing them in \"\" qoutes\n" .
	"\t-s [save name]: (optional) the name to save the group information generated by the script\n" .
	"\t    if not specified, will be bgf_file.grps. Each water shell will be save_name_shellx.bgf\n";
    return $usage;
}

#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
    unshift @INC, "/ul/tpascal/scripts/Packages";
}

use Packages::General qw(CoP FileTester TrjSelections IsInteger GetTorsion Trim);
use Packages::AMBER qw(getTopInfo ParseAmberTrj getOpts GetAmberByteOffset);
use Packages::MESO qw(GetMesoParms);
use Packages::FileFormats qw(sortByRes);
use Math::Polynomial::Solve qw(GetLeastSquaresPlane);
use Math::Trig;
use File::Basename;
use strict;

sub init;
sub doAnal;
sub sysInit;

die "usage: $0 amberTop amberTrj \"trajectory selection\" parm_file residue1 residue2 [saveName]\n"
    if (! @ARGV || $#ARGV < 5);

my ($topFile, $trjFile, $selection, $parmFile, $res1, $res2, $saveName) = @ARGV;
my ($SELECT, $DATA, $totAtms, $RESDATA, $printStr); 
my ($MOLOPTS, $PARMS, $OPTS);

print "Initializing...";
&init;
print "Done\nParsing AMBER topology file $topFile...";
($DATA, $totAtms) = getTopInfo($topFile, $OPTS);
$RESDATA = sortByRes($DATA->{ATOMS});
print "Done\nSetting up system...";
$MOLOPTS = &sysInit($DATA->{ATOMS}, $RESDATA, $PARMS, $res1, $res2);
print "Done\n";

$printStr = "Parsing AMBER trajectory $trjFile...";
open OUTDATA, "> $saveName" || die "ERROR: Cannot create $saveName: $!\n";
printf OUTDATA "%-8s%8s%8s\n","#TSTEP", "STACK ", "ROTATE ";
&GetAmberByteOffset($SELECT, $trjFile, $totAtms);
ParseAmberTrj(undef, $trjFile, $SELECT, $totAtms, \&doAnal, $printStr, \*OUTDATA);
close OUTDATA;
print "Created $saveName\n";

sub doAnal {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($i, $j, $k, $MOL, $dihdr);
    my ($a1, $a2, $b1, $b2, $c1, $c2);

    for $i ($res1, $res2) {
	for $j ("ATOMS", "PHO") {
	    for $k (keys %{ $MOLOPTS->{$i}{$j} }) {
		$MOL->{$i}{$j}{$k} = $ATOMS->{$k};
	    }
	    $MOL->{$i}{"${j}_COP"} = CoP($MOL->{$i}{$j});
	}
    }
    
    #stacking distance
    $dihdr = GetTorsion($MOL->{$res1}{ATOMS_COP}, $MOL->{$res1}{PHO_COP},
			$MOL->{$res2}{PHO_COP}, $MOL->{$res2}{ATOMS_COP}, 1);
    printf $fileHandle "%8d%8.3f", $frameNum, ($dihdr * 180/3.142);

    # rotational orientation
    for $i ($res1, $res2) {
	GetLeastSquaresPlane($MOL->{$i}, $MOL->{$i}{ATOMS_COP});
    }
	
    $a1 = $MOL->{$res1}{NORMAL}{a};
    $b1 = $MOL->{$res1}{NORMAL}{b};
    $c1 = $MOL->{$res1}{NORMAL}{c};
   
    $a2 = $MOL->{$res2}{NORMAL}{a};
    $b2 = $MOL->{$res2}{NORMAL}{b};
    $c2 = $MOL->{$res2}{NORMAL}{c};

    $dihdr = ($a1*$a2 + $b1*$b2 + $c1*$c2);
    $dihdr /= (sqrt($a1**2 + $b1**2 + $c1**2) * sqrt($a2**2 + $b2**2 + $c2**2));
    $dihdr = acos($dihdr);
    printf $fileHandle "%8.3f\n", ($dihdr*180/3.142);
    
}

sub init {
    $|++;
    FileTester($topFile);
    FileTester($trjFile);
    FileTester($parmFile);
    $SELECT = TrjSelections($selection);
    $OPTS = &getOpts;
    $PARMS = GetMesoParms($parmFile);

    for ($res1, $res2) {
	$_ = Trim($_);
	die "ERROR: Expected integer for residue. Got \"$_\"\n" if (! IsInteger($_));
    }

    if (! defined($saveName)) {
	$saveName = basename($topFile);
	$saveName =~ s/\.\w+$//;
	$saveName .= ".dat";
    }
}

sub sysInit {
    my ($atoms, $res, $parms, $id1, $id2) = @_;
    my ($i, $j, @atms, $resName, @tmp, @tmp2, $k, $l, %MOL);
    
    for $i ($id1, $id2) {
	die "ERROR: Residue $i does not exists!\n" if (! exists($res->{$i}));
	@tmp = keys %{ $res->{$i}{ATOMS} };
	$resName = $atoms->{$tmp[0]}{RESNAME};
	for $j (keys %{ $parms->{BEADS} }) {
	    next if ($parms->{BEADS}{$j}{NAME} =~ /SUO|SUP/ || ! $parms->{BEADS}{$j}{MEMBERS});
	    if ($parms->{BEADS}{$j}{RES} =~ /$resName/i || $parms->{BEADS}{$j}{NAME} eq "PHO") { # found residue
		@atms = split /\,/, $parms->{BEADS}{$j}{MEMBERS};
		for $k (@tmp) {
		    for $l (0 .. $#atms) {
			if (lc($atoms->{$k}{ATMNAME}) eq lc($atms[$l])) {
			    if ($parms->{BEADS}{$j}{NAME} eq "PHO") {
				$MOL{$i}{PHO}{$k} = 1;
			    } else {
				#print "res $i is $parms->{BEADS}{$j}{NAME}...";
				$MOL{$i}{ATOMS}{$k} = 1;
				$MOL{$i}{NAME} = $parms->{BEADS}{$j}{NAME};
			    }
			    splice(@atms, $l, 1);
			    last;
			}
		    }
		}
	    }
	}
    }

    print "res $id1 is $MOL{$res1}{NAME}, res $id2 is $MOL{$res2}{NAME}...";
    die "ERROR: No valid atoms found!\n" if (! keys %{ $MOL{$id1} } || ! keys %{ $MOL{$id2} });

    return \%MOL;
}

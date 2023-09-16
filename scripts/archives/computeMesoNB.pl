#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::AMBER qw(ParseAmberTrj getOpts getTopInfo GetAmberByteOffset);
use Packages::General qw(FileTester TrjSelections CoM GetSelections GetBondLength STDev);
use Packages::ManipAtoms qw(GetAtmList);

sub init;
sub showUsage;
sub parseMPSimEngFile;
sub calcCom;
sub numerically { ($a<=>$b); }

my ($trjFile, $topFile, $atomSelection, $saveFile, $createAvg);
my ($DATA, $totAtms, $printStr, $engFile, $ENERGY, $COMDATA, $SELECT, $AVGOPTS);

$|++;

&init;
print "Obtaining energies from $engFile->[0]...";
$ENERGY = parseMPSimEngFile($engFile, $SELECT);
print "Done\n";

&GetAmberByteOffset($SELECT, $trjFile, $totAtms);
$printStr = "Calculating COM from AMBER trj $trjFile...";
ParseAmberTrj($DATA->{ATOMS}, $trjFile, $SELECT, $totAtms, \&calcCom, $printStr, undef);

die "ERROR: No overlap between trjfile and engfile!\n" if (! keys %{ $COMDATA });
print "Creating COM vs Energy file ${saveFile}.dat...";
&writeData($COMDATA, $saveFile);
print "Done\n";

sub writeData {
    my ($data, $outFile) = @_;
    my ($i, %AVG, $stats);
    
    open OUTFILE, "> ${outFile}.dat" or die "ERROR: Cannot create ${outFile}.dat: $!\n";
    open COMFILE, "> ${outFile}_profile.dat" or die "ERROR: Cannot create ${outFile}_profile.dat: $!\n";
    for $i (sort numerically keys %{ $data }) {
	if (exists($AVGOPTS->{$i})) {
	    $AVG{ $AVGOPTS->{$i} }{COM} .= "$data->{$i}{COM} "; 
	    $AVG{ $AVGOPTS->{$i} }{ENG} .= "$data->{$i}{ENG} ";
	}
	printf OUTFILE "%8.3f %12.3f\n", $data->{$i}{COM}, $data->{$i}{ENG};
	printf COMFILE "%8d %8.3f\n", $i, $data->{$i}{COM};
    }

    if ($createAvg) {
	open AVGFILE, "> ${outFile}_avg.dat" or die "ERROR: Cannot creat ${outFile}_avg.dat: $!\n";
	for $i (sort numerically keys %AVG ) {
	    chop $AVG{$i}{COM};
	    chop $AVG{$i}{ENG};
	    $stats = ();
	    for ("COM", "ENG") {
		($stats->{$_}{AVG}, $stats->{$_}{STDEV}, $stats->{$_}{TOT}) = STDev($AVG{$i}{$_});
	    }
	    printf AVGFILE "%8.3f %12.3f # %8.3f %8.3f\n", $stats->{COM}{AVG}, $stats->{ENG}{AVG},
		$stats->{COM}{STDEV}, $stats->{ENG}{STDEV};
	}
	close AVGFILE;
    }
   
    close COMFILE;
    close OUTFILE;
}

sub calcCom {
    my ($ATOMS, $BOX, $frameNum, $junk) = @_;
    my ($GROUP, $i, $j, $comDist);
    
    next if (! exists($ENERGY->{$frameNum}));
    for $i (1, 2) {
	for $j (keys %{ $atomSelection->{$i} }) {
	    die "ERROR: Atom $j is not valid!\n" if (! exists($ATOMS->{$j} ));
	    %{ $GROUP->{$i}{$j} } = %{ $ATOMS->{$j} };
	}
	$GROUP->{$i}{COM} = CoM($GROUP->{$i});
    }

    $comDist = GetBondLength($GROUP->{1}{COM}, $GROUP->{2}{COM});
    $COMDATA->{$frameNum} = (
			     {
				 "COM" => $comDist,
				 "ENG" => $ENERGY->{$frameNum},
			     }
			     );
}
sub parseMPSimEngFile {
    my ($ENGFILES, $trjSelect) = @_;
    my (%ENG, $inFile);

    for $inFile (@{ $ENGFILES }) {
	open MPSIMFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
	while (<MPSIMFILE>) {
	    chomp;
	    if ($_ =~ /^\s+(\d+)\s+(\-?\d+\.\d+)/) {
		if (! keys %{ $trjSelect } or exists($trjSelect->{$1})) {
		    $ENG{$1} += $2;
		}
	    }
	}
    }
    close MPSIMFILE;

    die "ERROR: $inFile does not contain any valid data!\n" if (! %ENG);
    
    return \%ENG;
}

sub showUsage {
    my ($usage) = "usage: $0 -e eng file -p top file -f amber trajectory -t trajectory selection " . 
	"-a atom selection -o [trj avg options] -s [saveFile]\n" .
	"Options:\n\t-f amber trajectory: location of amber trajectory file (with box information)\n" . 
	"\t-p top file: location of amber topology file\n" .
	"\t-e eng file: location of the energy per timestep file (created with /ul/tpascal/scripts/dnaEnergyAnal.pl)\n" .
	"\t-t trajectory selection: For all frames in the trajectory: \"*\". A integer specifies as specific frame\n" .
	"\t\tTo specify a range :Ita-b:c will select frames a to b, every c frames\n" .
	"\t-a atom selection: the atoms to be in group 1 and group 2. Expects \"group1 group2\".\n" . 
	"\t\tExpects [:][I|N|T][a|r]x[-y:z].\n" . 
	"\t\twill select a Integer|Name|Type of Atom|Residue of x (to y step z if : is specified)\n" .
	"\t-o [trj avg options]: (Optional) Specifies the frames to average over. Expected x-y:z\n" .
	"\t\twhere x: start, y: end and z is total frames per block. So 20-30:45 will average every 20th\n" .
	"\t\tto 30th frames,  with 45 frames per interval\n" .
	"\t-s [saveFile]: (Optional) Name of the file to save\n";
    die $usage . "\n";
}

sub init {
    my (%FOPTS, $tSelection, $aSelection, @tmp, $OPTS, $avg, $i, $j, $interval, $EFILES);

    getopt('ftsaepo',\%FOPTS);
    ($trjFile, $tSelection, $aSelection, $saveFile, $EFILES, $topFile, $avg) = 
	($FOPTS{f}, $FOPTS{t}, $FOPTS{a}, $FOPTS{s}, $FOPTS{e}, $FOPTS{p}, $FOPTS{o});
    
    for ($trjFile, $tSelection, $aSelection, $EFILES, $topFile) {
	&showUsage if (! defined($_));
    }

    print "Initializing...";
    $createAvg = 0;
    FileTester($trjFile);
    while ($EFILES =~ /(\S+)/g) {
	if (-e $1 and -r $1 and -T $1) {
	    push @{ $engFile }, $1;
	}
    }
    die "ERROR: No valid energy files found while search $EFILES\n" if (! @{ $engFile });
    FileTester($topFile);
    $SELECT = TrjSelections($tSelection);

    $OPTS = &getOpts;
    if (! defined($saveFile)) {
	$saveFile = basename($trjFile);
	$saveFile =~ s/\.\w+$/_com\.dat/;
    }
    print "Parsing AMBER topology file $topFile...";
    ($DATA, $totAtms) = getTopInfo($topFile, $OPTS);
    if ($aSelection =~ /(.+)\s+(.+)/) {
	@tmp = ($1);
	$atomSelection->{1} = GetAtmList(GetSelections(\@tmp, 0), $DATA->{ATOMS});
	die "ERROR: Atom selection $1 is invalid!\n" if (! keys %{ $atomSelection->{1} });
	@tmp = ($2);
	$atomSelection->{2} = GetAtmList(GetSelections(\@tmp, 0), $DATA->{ATOMS});
 	die "ERROR: Atom selection $2 is invalid!\n" if (! keys %{ $atomSelection->{2} });
   } else {
	print "ERROR: Invalid atom selection!\n";
	&showUsage;
    }
    
    $saveFile =~ s/\.\w+//;
    
    if ($avg =~ /^(\d+)\-(\d+):(\d+)/) {
	$interval = $3;
	$interval = $2 if ($2 > $3);
	for $i (1 .. $interval) {
	    for $j ($1 .. $2) {
		$AVGOPTS->{ ($j + (($i - 1) * $3)) } = $i;
	    }
	}
    }

    $createAvg = 1 if (keys %{ $AVGOPTS });

    print "Done\n";
}
  

#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use Packages::General qw(FileTester);
use strict;
use warnings;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::LAMMPS qw(CreateLAMMPSTrj);
use Packages::AMBER qw(CreateAmberTrj);

sub init;
sub parseSeqquestNEBout;
sub writeData;
sub numerically { ($a<=>$b) };
sub showHelp;
sub createLammpsHeaders;

my ($nebFile, $trjFile, $DATA, $trjType, $savePrefix, $BOX, $atmCount);

$|++;
&init;
print "Processing Seqquest NEB output $nebFile...";
($DATA, $BOX, $atmCount) = parseSeqquestNEBout($nebFile);
print "Done\nCreating $trjType trajectory file ${trjFile}...";
writeData($DATA, $trjFile, $trjType, $savePrefix, $BOX, $atmCount);
print "Done\nCreated ${savePrefix}.energies ${savePrefix}.forces\n";

sub createLammpsHeaders {
    my ($box, $numAtms) = @_;
    my (%DATA); 
		  
    $DATA{"TIMESTEP"}[0] = 0;
    $DATA{"NUMBER OF ATOMS"}[0] = $numAtms;
    $DATA{"BOX BOUNDS"}[0]{lo} = $box->{X}{lo};
    $DATA{"BOX BOUNDS"}[0]{hi} = $box->{X}{hi};
    $DATA{"BOX BOUNDS"}[1]{lo} = $box->{Y}{lo};
    $DATA{"BOX BOUNDS"}[1]{hi} = $box->{Y}{hi};
    $DATA{"BOX BOUNDS"}[2]{lo} = $box->{Z}{lo};
    $DATA{"BOX BOUNDS"}[2]{hi} = $box->{Z}{hi};

    return \%DATA;
}

sub writeData {
    my ($data, $trjFile, $trjType, $prefix, $box, $numAtms) = @_;
    my ($trjWriter, $INFO, $i);

    open COORDFILE, "> $trjFile" or die "ERROR: Cannot create $trjType trajectory file $trjFile: $!\n";
    if (lc($trjType) eq "amber") {
	print COORDFILE "Title: ${prefix} NEB trajectory\n";
	$trjWriter = \&CreateAmberTrj;
	$INFO = $box;
    } else {
	$trjWriter = \&CreateLAMMPSTrj;
	ScaleAtoms($DATA->{COORDS}, $BOX);
	$INFO = createLammpsHeaders($BOX, $numAtms);
    }

    for $i (sort numerically keys %{ $DATA->{COORDS} }) {
	$INFO->{TIMESTEP}[0] = $i if ($trjType eq "lammps");
	$trjWriter->($DATA->{COORDS}{$i}, $INFO, \*COORDFILE);
    }
    close COORDFILE;

    open FORCEFILE, "> ${prefix}.forces" or die "ERROR: Cannot create ${prefix}.forces: $!\n";
    for $i (sort numerically keys %{ $DATA->{FORCES} }) {
	$INFO->{TIMESTEP}[0] = $i if ($trjType eq "lammps");
	$trjWriter->($DATA->{FORCES}{$i}, $INFO, \*FORCEFILE);
    }
    close FORCEFILE;

    open ENGFILE, "> ${prefix}.energies" or die "ERROR: Cannot create ${prefix}.energies: $!\n";
    for $i (sort numerically keys %{ $DATA->{ENERGY} }) {
	printf ENGFILE "%10d %12.6f\n", $i, $DATA->{ENERGY}{$i};
    }
    close ENGFILE;
}


sub parseSeqquestNEBout {
    my ($inFile) = $_[0];
    my (%DATA, $byteOffset, $inStr, $field, $imageNum, $lastFrame, $rec, $i, %BOX);
    my ($x, $y, $z);

    my ($grepCmd) = "grep -b 'NEB modified forces' $inFile";
    
    $byteOffset = -1;
    open GREPCMD, "$grepCmd |" || die "ERROR: Cannot execute $grepCmd: $!\n";
    while (<GREPCMD>) {
        chomp;
	if ($_ =~ /^(\d+)/) {
	    $byteOffset = $1;
	}
    }
    close GREPCMD;
    
    die "ERROR: File $inFile does not contain any valid data!\n" if ($byteOffset == -1);

    $imageNum = $lastFrame = -1;
    open NEBFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
    die "ERROR: $!\n" if (! seek(NEBFILE,$byteOffset,0));
    print "" if (<NEBFILE>);
    while (<NEBFILE>) {
	chomp;
	$inStr = $_;
	if ($inStr =~ /^(imagef - image number|image number in NEB)/) {
	    $field = "frame_num";
	} elsif ($inStr =~ /^eimage - energy of image/) {
	    $field = "frame_energy";
	} elsif ($inStr =~ /^forcei - forces/) {
	    $field = "forces";
	    $i = 1;
	} elsif ($inStr =~ /^einitial - energy for start-point/) {
	    $field = "frame_energy";
	    $imageNum = 0;
	} elsif ($inStr =~ /^efinal - energy for end-product/) {
	    $field = "frame_energy";
	    $imageNum = $lastFrame + 1;
	} elsif ($inStr =~ /^coordinates for/) {
	    $field = "coords";
	    $i = 1;
	} elsif ($field eq "frame_num" && $inStr =~ /\s+(\d+)/) {
	    $imageNum = $1;
	    $lastFrame = $1 if ($1 > $lastFrame);
	} elsif ($imageNum > -1) {
	    if ($field eq "frame_energy" && $inStr =~ /^\s+(\-?\d+\.\d+)/) {
		$DATA{ENERGY}{$imageNum} = $1;
	    } elsif ($field eq "forces" || $field eq "coords") {
		if ($inStr =~ /^\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)/) {
		    $x = $1;
		    $y = $2;
		    $z = $3;

		    if ($field eq "coords") {
			$x *= 0.529177;
			$y *= 0.529177;
			$z *= 0.529177;
			$BOX{X}{hi} = $x if (! exists($BOX{X}{hi}) || $BOX{X}{hi} < $x);
			$BOX{Y}{hi} = $y if (! exists($BOX{Y}{hi}) || $BOX{Y}{hi} < $y);
			$BOX{Z}{hi} = $z if (! exists($BOX{Z}{hi}) || $BOX{Z}{hi} < $z);
			$BOX{X}{lo} = $x if (! exists($BOX{X}{lo}) || $BOX{X}{lo} > $x);
			$BOX{Y}{lo} = $y if (! exists($BOX{Y}{lo}) || $BOX{Y}{lo} > $y);
			$BOX{Z}{lo} = $z if (! exists($BOX{Z}{lo}) || $BOX{Z}{lo} > $z);
		    }

		    $rec = (
			    {
				"XCOORD" => $x,
				"YCOORD" => $y,
				"ZCOORD" => $z,
			    }
			    );

		    $DATA{uc($field)}{$imageNum}{$i} = $rec;
		    $i++;
		}
	    }
	}
    }
    close NEBFILE;
    
    die "ERROR: No valid data found in $inFile!\n" if (! keys %DATA);
    return (\%DATA, \%BOX, ($i - 1));
}
    
sub init {
    my (%OPTS, $usage);

    getopt('ost',\%OPTS);
    
    ($nebFile, $savePrefix, $trjType) = ($OPTS{o},$OPTS{s},$OPTS{t});

    $usage = &showHelp;
    die "$usage\n" if (! defined($nebFile));
    
    if (! defined($savePrefix)) {
	$savePrefix = basename($nebFile);
	$savePrefix =~ s/\.\w+$/_seqquestNEB/;
    }

    print "Initializing...";
    FileTester($nebFile);
    $trjType = "amber" if (! defined($trjType));
    $trjType = lc($trjType);
    
    if ($trjType =~ /(lammps|amber)/) {
	$trjType = $1;
	if ($trjType eq "lammps") {
	    $trjFile = $savePrefix . ".lammpstrj";
	} else {
	    $trjFile = $savePrefix . ".crd";
	}
    }

    $trjType = uc($trjType);

    die "ERROR: Invalid trajectory type specified!\n$usage\n" if (! defined($trjFile));
}

sub showHelp {
    my ($usage) = "usage: $0 -o seqquest neb out file -s saved file prefix -t [trajectory type (amber is default)]\n";
    $usage .= "options:\n";
    $usage .= "\t-o seqquest neb file: location of Seqquest NEB output file. Will fail if line " . 
	"\"NEB modified forces\" not present\n";
    $usage .= "\t-s saved file prefix: the prefix for the save trajectory. Each image will be a seperate frame.\n" .
	"\t\tThis prefix can also include the path information: i.e. /ul/tpascal/seq_results/enthanol is acceptable\n";
    $usage .= "\t-t [trajectory type]: (optional) Sepecifies the type of the trajectory to create. The default is AMBER.\n" . 
	"\t\tTher other choices is Lammps.";
    return $usage;
}
    

#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester TrjSelections);
use LAMMPS qw(ParseLAMMPSLogFile);
use FileFormats qw(GetBGFFileInfo);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub numerically { ($a<=>$b); }
sub runLAMMPS;
sub execCmd;
sub getBestStructures;
sub getLammpsEnergy;
sub GetRange;
sub saveFile;

my ($monomer, $DIMS, $datfile, $lammps, $prefix, $dirname); 
my ($monomerEng, $axis, $forcefield, $inputTemplate);
my ($DATA, $numAtoms, $tmp, @junk, $EXTREMA);

$|++;
&init;
print "Parsing bgf file $monomer...";
($tmp, undef, undef) = GetBGFFileInfo($monomer, 0);
@junk  = sort numerically(keys %{ $tmp });
$numAtoms = pop @junk;
print "Done\n";
($DATA, $dirname) = &runLAMMPS($monomer, $forcefield, $inputTemplate, $DIMS, $lammps, $axis, $prefix, $numAtoms);
print "\nGetting lowest energy structures...";
$EXTREMA = getBestStructures($DATA, $axis);
print "Done\nSaving datafile $datfile...";
&writeData($DATA, $EXTREMA, $datfile, $DIMS);
system("rm -fr $dirname");
print "Done\n";

sub writeData {
    my ($data, $best, $file, $dimOpts) = @_;
    my ($prefix, @dim, $i, $filename, $j, %dimRange);
    
    for $i ("x", "y", "z") {
	if (! exists($dimOpts->{$i})) {
	    $dimRange{$i}[0] = 0;
	} else {
	    @{ $dimRange{$i} } = sort numerically keys %{ $dimOpts->{$i} };
	}
    }

    @dim = grep {!/$axis/} ("x", "y", "z");
    $prefix = basename($file);
    $prefix =~ s/\.\w+//;
    system("mkdir -p ${prefix}_structures ${prefix}_structures/$dim[0]_best ${prefix}_structures/$dim[1]_best");
    for $i ("MIN", "MAX") {
	$filename = "${prefix}_structures/${i}.bgf";
	system("cp $best->{$i}{FILE} $filename");
    }
    for $i (@dim) {
	for $j (keys %{ $best->{$i} }) {
	    $filename = "${prefix}_structures/${i}_best/${j}.bgf";
	    system("cp $best->{$i}{$j}{FILE} $filename");
	}
    }
    &saveFile($data, $file, \%dimRange, \@dim);
    $dimRange{$dim[0]} = ([$dim[0], $dim[1]]);
    &saveFile($best, "${prefix}_structures/minima", \%dimRange, \@dim);
}

sub saveFile {
    my ($data, $file, $dimRange, $dim) = @_;
    my ($i, $j);

    $file =~ s/\.\w+\/?$//;

    open ENG, "> ${file}_energy.csv" or die "ERROR: Cannot create ${file}_energy.dat: $!\n";
    open DIST, "> ${file}_dist.csv" or die "ERROR: Cannot create ${file}_dist.dat: $!\n";
    print ENG ",";
    print DIST ",";
    for $i (@{ $dimRange->{$dim->[1]} }) {
	print ENG "$i,";
	print DIST "$i,";
    }
    print ENG "\n";
    print DIST "\n";
    for $i (@{ $dimRange->{$dim->[0]} }) {
	print ENG "$i,";
	print DIST "$i,";
	for $j (@{ $dimRange->{$dim->[1]} }) {
	    if (! exists($data->{$i}{$j})) {
		print ENG ",";
		print DIST ",";
	    } else {
		if (exists( $data->{$i}{$j}{BEST} )) {
		    printf ENG "%12.3f,", $data->{$i}{$j}{BEST}{ENERGY};
		    printf DIST "%12.3f,", $data->{$i}{$j}{BEST}{DIST};
		} else {
		    printf ENG "%12.3f,", $data->{$i}{$j}{ENERGY};
		    printf DIST "%12.3f,", $data->{$i}{$j}{DIST};
		}
	    }
	}
	print ENG "\n";
	print DIST "\n";
    }
    close ENG;
    close DIST;
}

sub getBestStructures {
    my ($data, $axis) = @_;
    my (@dim, $i, $j, $k, %BEST, $min);

    @dim = grep {!/$axis/} ("x", "y", "z");

    for $i (keys %{ $data }) {
	for $j (keys %{ $data->{$i} }) {
	    $min = ();
	    for $k (keys %{ $data->{$i}{$j}{ENERGIES} }) {
		if (! $min or $data->{$i}{$j}{ENERGIES}{$k} < $min->{ENERGY}) {
		    $min->{ENERGY} = $data->{$i}{$j}{ENERGIES}{$k};
		    $min->{DIST} = $k;
		    $min->{FILE} =  $data->{$i}{$j}{FILES}{$k};
		    $min->{J} = $j;
		}
	    }
	    $data->{$i}{$j}{BEST} = \%{ $min };
	    $BEST{$dim[0]}{$i} = \%{ $min } if (! exists($BEST{$dim[0]}{$i}) or 
						$min->{ENERGY} < $BEST{$dim[0]}{$i}{ENERGY} or
						($min->{ENERGY} == $BEST{$dim[0]}{$i}{ENERGY} and
						 $min->{J} > $BEST{$dim[0]}{$i}{J}));
	    $BEST{$dim[1]}{$j} = \%{ $min } if (! exists($BEST{$dim[1]}{$j}) or 
						$min->{ENERGY} < $BEST{$dim[1]}{$j}{ENERGY} or
						($min->{ENERGY} == $BEST{$dim[1]}{$j}{ENERGY} and
						 $min->{J} > $BEST{$dim[1]}{$j}{J}));
	    $BEST{MIN} = \%{ $min } if (! exists($BEST{MIN}) or 
					$min->{ENERGY} < $BEST{MIN}{ENERGY} or
					($min->{ENERGY} == $BEST{MIN}{ENERGY} and
					 $min->{J} > $BEST{MIN}{J}));
	    $BEST{MAX} = \%{ $min } if (! exists($BEST{MAX}) or 
					$min->{ENERGY} > $BEST{MAX}{ENERGY} or
					($min->{ENERGY} == $BEST{MAX}{ENERGY} and
					 $min->{J} > $BEST{MAX}{J}));
	}
    }
    return \%BEST;
}

sub GetRange {
    my ($range) = $_[0];
    my (%RANGE, $start, $stop, $interval, $curr, $dist);
    my ($prec);

    if ($range =~ /:(-?\d+\.?\d*)-(-?\d+\.?\d*):(\d+\.?\d*)/) {
	($start, $stop, $interval) = ($1, $2, $3);
    } elsif ($range =~ /(-?\d+\.?\d*)/) {
	($start, $stop, $interval) = ($1, $1, $1);
    }

    $prec = 0;
    if ($interval =~ /\./g) {
        while ($interval =~ /(\d)/g) {
            $prec++
	}
	$prec = "%0.${prec}";
    }

    $curr = $start;
    while ($curr <= $stop) {
	if ($interval =~ /\./) {
	    $dist = sprintf("${prec}f", $curr);
	    $RANGE{$dist} = "${prec}f";
	} else {
	    $RANGE{$curr} = 1;
	}
	$curr += $interval;
    }
    
    die "ERROR: Invalid range specified:  \"$range\"\n" if (! %RANGE);
    return \%RANGE;
}

sub getLammpsEnergy {
    my ($logFile) = $_[0];
    my ($potEng);

    open LOGFILE, $logFile or die "ERROR: Cannot open $logFile: $!\n";
    while (<LOGFILE>) {
	chomp;
	if ($_ =~ /Poteng\s+\=\s+(\-?\d+\.\d+)/i) {
	    $potEng = $1;
	    last;
	}
    }
    close LOGFILE;
    die "ERROR: $logFile does not contain any valid info!\n" if (! defined($potEng));
    return ($potEng - (2*$monomerEng));
}

sub runLAMMPS {
    my ($monomerFile, $forcefield, $template, $dimOpts, $lammpsBinary, $axis, $prefix, $tot) = @_;
    my (%dimRange, $i, $j, $k, @dim, %DATA, $LOGDATA);
    my ($execScript, $bgfFile, $printStr);

    my ($duplicateMolScript) = "/ul/tpascal/scripts/transmol.pl";
    my ($transMolScript) = "/ul/tpascal/scripts/modifyAtomData.pl";
    my ($createLammpsInputScript) = "/ul/tpascal/scripts/createLammpsInput.pl";

    $tot++;
    for $i ("x", "y", "z") {
	if (! exists($dimOpts->{$i})) {
	    $dimRange{$i}[0] = 0;
	} else {
	    @{ $dimRange{$i} } = sort numerically keys %{ $dimOpts->{$i} };
	}
    }

    @dim = grep {!/$axis/} ("x", "y", "z");

    $prefix .= "_" . basename($forcefield) . "_pes";
    system("rm -fr getPES.log");
    for $i (@{ $dimRange{$dim[0]} }) {
	for $j (@{ $dimRange{$dim[1]} }) {
	    $printStr = "Running PES calculation for $dim[0] $i $dim[1] $j ...";
	    print "$printStr\r";
	    $execScript = "mkdir -p ${prefix}/$dim[0]_${i}/$dim[1]_${j}";
	    &execCmd($execScript);
	    for $k (@{ $dimRange{$axis} }) {
		$bgfFile = "${prefix}/$dim[0]_${i}/$dim[1]_${j}/${axis}_${k}.bgf";
		$execScript = "$duplicateMolScript $monomerFile $k $axis $bgfFile";
		&execCmd($execScript);
		$execScript = "$transMolScript -a \"Ia<${tot}\" -s $bgfFile -w $bgfFile -f \"" .
		    uc($dim[0]) . "COORD:+${i} " . uc($dim[1]) . "COORD:+${j}\"";
		&execCmd($execScript);
		$execScript = "$createLammpsInputScript -b $bgfFile -f $forcefield -s ${prefix}_${axis}_${k}";
		&execCmd($execScript);
		$execScript = "$lammpsBinary -in $template -var datafile data.${prefix}_${axis}_${k} " . 
		    "-log ${prefix}_${axis}_${k}_log.lammps";
		&execCmd($execScript);
		$DATA{$i}{$j}{ENERGIES}{$k} = getLammpsEnergy("${prefix}_${axis}_${k}_log.lammps");
		$DATA{$i}{$j}{FILES}{$k} = $bgfFile;
		system("rm -fr in.${prefix}_${axis}_${k}* data.${prefix}_${axis}_${k} " . 
		       "${prefix}_${axis}_${k}_lammps.script ${prefix}_${axis}_${k}_log.lammps");
	    }
	}
    }

    print "%-" . (length($printStr)+10) . "s\n", "Running PES calculations...Done";
    return (\%DATA, $prefix);
}

sub execCmd {
    my ($cmdStr) = $_[0];
    $cmdStr .= " >> getPES.log";
    if (system($cmdStr)) {
        die "ERROR: Cmd: $_[0]\nCheck getPES.log\n";
    }
}

sub init {
    my (%OPTS, $usage, $tmp);
    
    $usage =  "usage: $0 -b bgf monomer -l lammps binary -a axis -t lammps input template file -f forcefield\n" . 
	"Optional:\n\t-s (save file) -e lammps monomer energy and any combination of\n" . 
	"\t-x xrange (either a or :a-b:c), -y yrange, -z zrange\n";
    
    getopt('ebxyzlsaft',\%OPTS);
    for ("a", "b", "l", "f", "t") {
	die "$usage" if (! exists($OPTS{$_}));
    }
    die "$usage" if (! exists($OPTS{x}) and ! exists($OPTS{y}) and ! exists($OPTS{z}));
    print "Initializing...";
    FileTester($OPTS{b});
    die "ERROR: Cannot access $OPTS{l}: $!\n" if (! -e $OPTS{l});
    FileTester($OPTS{f});
    FileTester($OPTS{t});

    die "ERROR: Expected \"x|y|z\" for axis, got \"$axis\"\n" if ($OPTS{a} !~ /^(x|y|z)/i);
    $axis = lc($1);

    $monomer = $OPTS{b};
    $lammps = $OPTS{l};
    $forcefield = $OPTS{f};
    $inputTemplate = $OPTS{t};
    $monomerEng = 0;
    $monomerEng = $OPTS{e} if (exists($OPTS{e}) and $OPTS{e} =~ /^\-?\d+\.?\d*/);

    for ("x", "y", "z") {
	next if (! exists($OPTS{$_}));
	$DIMS->{$_} = GetRange($OPTS{$_});
    }
    die "ERROR: No valid selections found!\n" if (! $DIMS);
    if (! exists($OPTS{s})) {
	$datfile = $monomer;
	$datfile =~ s/\.\w+//;
	$datfile .= "_pes.dat";
    } else {
	$datfile = $OPTS{s};
    }
    $prefix = basename($datfile);
    $prefix =~ s/\.\w+$//;

    print "Done\n";
}
	

#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/.libs/Math-ematica-1.108";
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo);
use Packages::MESO;
use Packages::General qw(FileTester Trim STDev CenterText);
use Math::ematica qw(:PACKET :TYPE :FUNC);
use File::Basename;
use Packages::CERIUS2 qw(saveCeriusFF);
use constant PI => atan2(1,1) * 4;

#  genMesoFF.pl: This will take the files in the specified directory
#  and run Mathematica to obtain a non-linear curve fit to obtain valence
#  parameters. It will then print the results in a Cerius2 format

sub main;
sub init;
sub getVType;
sub getFFTypes;
sub getValenceCoeffs;
sub mathLinkCmd;
sub readDatFile;
sub numerically;
sub getFitFunc;
sub extractVals;
sub createCeriusFF;
sub printResults;
sub addHBondAtoms;

die "usage: $0 dat_file|directory meso_parm_file [ff_save_name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($datLoc, $parmFile, $saveName) = @ARGV;

$|++;
&main;

sub main {
    my ($DATFILES, $i, $PARMS, $data, $table, $pFF); 
    my ($dumpFile, $rec, $FFPARMS, $FF, $ELEMENTS);
    
    $|++;

    $dumpFile = "dump.dat";
    $DATFILES = &init;
    $pFF = 0;

    open DUMPFILE, "> $dumpFile" or die "Cannot write to file $dumpFile: $!\n";    
    for $i (@{ $DATFILES }) {
	$i->{"RESULTS"} = ();
	$rec = getValenceCoeffs($i, \*DUMPFILE);
	if ($rec) {
	    $i->{"RESULTS"} = $rec;
	    print "$i->{FILE}:...Done!\n";
	    $pFF = 1;
	} else {
	    print "$i->{FILE}:...Failed!\n";
	}
    }
    close DUMPFILE or die "ERROR: Cannot finish creating file $dumpFile: $!\n";

    if ($pFF) {
	print "Creating Cerius2 FF...";
	($FF, $ELEMENTS) = createCeriusFF($DATFILES);
	print "Done\n"; 
	print "Saving Cerius2 FF...";
	saveCeriusFF($FF, $saveName, $ELEMENTS);
	print "Done\n";
	printResults($DATFILES);
    } else {
	print "Cannot save FF, No valid data\n";
    }
}

sub printResults {
    my ($data) = $_[0];
    my (%STATS, $i, $parms, @NAMES, $type, @headers);

    @NAMES = ("BONDS", "ANGLES", "TORSIONS");
    for $i (@{ $data }) {
	next if (! exists($i->{"RESULTS"}) or ! keys %{ $i->{"RESULTS"} });
	$type = $i->{VTYPE};
	$type--;
	$type = $NAMES[$type];
	$i->{TITLE} = Trim($i->{TITLE});
	for $parms (keys %{ $i->{"RESULTS"}{"PARMS"} }) {
	    next if ($parms eq "a0");
	    $STATS{$type}{VALS}{ $i->{TITLE} }{$parms} = $i->{"RESULTS"}{"PARMS"}{$parms};
	    $STATS{$type}{HEADERS}{$parms} = 1;
	}
	$STATS{$type}{VALS}{ $i->{TITLE} }{"variance"} = $i->{"RESULTS"}{"VARIANCE"};
	$STATS{$type}{HEADERS}{"variance"} = 1;
    }

    for $type (@NAMES) {
	next if (! exists($STATS{$type}));
	print "" . CenterText("$type STATS", 40) . "\n";
	printf "%20s", " ";
	@headers = ();
	for $parms (keys %{ $STATS{$type}{HEADERS} }) {
	    printf "%-10s", $parms;
	    push @headers, $parms;
	}
	print "\n";
	
	for $i (keys %{ $STATS{$type}{VALS} }) {
	    printf "%-20s", $i;
	    for $parms (@headers) {
		printf "%-10.3f", $STATS{$type}{VALS}{$i}{$parms};
	    }
	    print "\n";
	}
	printf "\n";
    }
}

sub createCeriusFF {
    my ($data) = $_[0];
    my (%FF, $i, $atom, $j, @atomList, $fftype, %ELEMENTS, $counter);
    
    $counter = 0;
    for $i (@{ $data }) {
	next if (! exists($i->{"RESULTS"}) or ! keys %{ $i->{"RESULTS"} });
	@atomList = sort numerically keys %{ $i->{"ATOMS"} };
	for $j (@atomList) {
	    $atom = $i->{"ATOMS"}{$j};
	    if (! exists($FF{"atoms"}{ $atom->{"FFTYPE"} })) {
		$FF{"atoms"}{ $atom->{"FFTYPE"} }{"VALS"}[0] = (
								{
								    "element" => $atom->{"ELEMENT"},
								    "mass"    => $atom->{"MASS"},
								    "hybrid"  => $atom->{"NUMBONDS"},
								    "r"       => $atom->{"RADII"},
								    "e"       => 0,
								    "name"    => $atom->{"NAME"},
								}
								);
		$ELEMENTS{ $atom->{"ELEMENT"} }{"SYMBOL"} = $atom->{"ELEMENT"};
	    }
	}
	$fftype =  $i->{"ATOMS"}{ shift @atomList }{"FFTYPE"};
	#$FF{ $i->{"RESULTS"}{"TYPE"} }{$fftype} = ();
	$atom = \%{ $FF{ $i->{"RESULTS"}{"TYPE"} }{$fftype} };
	for $j (@atomList) {
	    $fftype =  $i->{"ATOMS"}{$j}{"FFTYPE"};
	    #$atom->{$fftype} = ();
	    $atom = \%{ $atom->{$fftype} };
	}
        if ($i->{VTYPE} == 3) {
            $counter++;
            $atom->{counter} = $counter;
            $FF{"torsionOrders"}{$counter} = "0 1 2 3";
	    if ($i->{RESULTS}{PARMS}{kp} < 0) { #negative torsion fix
		$i->{RESULTS}{PARMS}{kp} *= -1;
		#if ($i->{RESULTS}{PARMS}{p0} > PI) { #> 180deg so angle = angle - 180
		    $i->{RESULTS}{PARMS}{p0} -= PI;
		#} else { # else angle = angle + 180
		#    $i->{RESULTS}{PARMS}{p0} += PI;
		#}
	    }
        }
	for $j (keys %{ $i->{"RESULTS"}{"PARMS"} }) {
	    $atom->{"VALS"}[0]{$j} = $i->{"RESULTS"}{"PARMS"}{$j};
	}
    }

    return (\%FF, \%ELEMENTS);
}


sub getValenceCoeffs {
    my ($data, $dumpfile) = @_;
    my ($link, $error, $rawData, $binData, $rec, $histCat, $histData);

    $link = new Math::ematica '-linklaunch', '-linkname', 'math -mathlink';

    print "$data->{FILE}:...reading\r";    
    $rawData = readDatFile($data->{FILE}, $data->{VTYPE});
    printf "$data->{FILE}:%25s\r", " ";

    print "$data->{FILE}:...creating bins\r";
    ($binData, $histData, $histCat) = createBins($rawData, $data->{VTYPE});
    printf "$data->{FILE}:%25s\r", " ";

    print "$data->{FILE}:...creating histogram\r";
    $error = createHist($binData, $histData, $rawData, $histCat, $link, $dumpfile);
    printf "$data->{FILE}:%25s\r", " ";
    return 0 if ($error);
    
    print "$data->{FILE}:..fitting\r";
    $rec = getFit($data, $rawData, $link, $dumpfile);
    printf "$data->{FILE}:%25s\r", " ";
    return 0 if (! $rec->{VARIANCE});
   
    print "$data->{FILE}:...plotting\r";
    $error = createPlots($data, $rawData, $rec, $link, $dumpfile);
    if ($error) {
	print "$data->{FILE}:$error...\r";
	return 0;
    }
    
    return $rec;
}

sub createHist {
    my ($binData, $histData, $rawData, $histCat, $link, $dumpfile) = @_;
    my ($cmd, $error, $result);

    $cmd = "<<Graphics`Graphics`";
    mathLinkCmd($cmd, 1, $link, $dumpfile);
    $cmd = "<< Graphics`Legend`";
    mathLinkCmd($cmd, 1, $link, $dumpfile);
    $cmd = "mdata = $binData;";
    ($error, $result) = mathLinkCmd($cmd, 0, $link, $dumpfile);
    if (! $error) {
	$cmd = "countdata = {$histData};";
	mathLinkCmd($cmd, 1, $link, $dumpfile);
	$cmd = "hist = Histogram[countdata, FrequencyData -> True, " .
	    "HistogramScale -> $rawData->{AREA}, HistogramCategories -> {$histCat}];";
	mathLinkCmd($cmd, 1, $link, $dumpfile);
	return "";
    } else {
	return $error;
    }
}

sub getFit {
    my ($data, $rawData, $link, $dumpfile) = @_;
    my ($cmd, $fitFunction, $vars, $tmp1, $tmp, $rec, $i, $error, $result, $j, $tmp2, $p0Val);

    $cmd = "Needs[\"Statistics`NonlinearFit`\"];";
    mathLinkCmd($cmd, 1, $link, $dumpfile);
    ($fitFunction, $vars) = getFitFunc($data, $rawData);
    for $i (1 .. 6) {
	for $j (1 .. 6) {
	    $tmp1 = $fitFunction;
	    $tmp2 = $vars;
	    $tmp1 =~ s/periodicity/$i/;
	    $p0Val = ((PI * $j)/6);
	    $tmp2 =~ s/p0Val/$p0Val/;
	    $cmd = "fitFunc = $tmp1;";
	    mathLinkCmd($cmd, 1, $link, $dumpfile);
	    $cmd = "InputForm[{BestFitParameters, EstimatedVariance} /. " . 
		"NonlinearRegress[mdata, fitFunc, {x}, {$tmp2}]]";
	    ($error, $result) = mathLinkCmd($cmd, 0, $link, $dumpfile);
	    if ($error) {
		showError($error, $link);
		return 0;
	    }
	    $tmp = extractVals($result, $data->{VTYPE});
	    $tmp->{FITFUNC} = $tmp1;
	    $tmp->{VARS} = $tmp2;
	    if ($data->{VTYPE} != 3) {
		%{ $rec } = %{ $tmp };
		last;
	    } else {
		if (! exists($rec->{VARIANCE}) || $tmp->{VARIANCE} < $rec->{VARIANCE}) {
		    %{ $rec } = %{ $tmp };
		    $rec->{PARMS}{n} = $i;
		}
	    }
	}
	last if ($data->{VTYPE} != 3);
    }

    return $rec;
}

sub createPlots {
    my ($data, $rawData, $rec, $link, $dumpfile) = @_;
    my ($cmd, $error, $result);

    $cmd = "fitFunc = $rec->{FITFUNC};";
    mathLinkCmd($cmd, 1, $link, $dumpfile);    
    $cmd = "datafit = NonlinearFit[mdata, fitFunc, {x}, {" . $rec->{VARS} . "}]";
    mathLinkCmd($cmd, 1, $link, $dumpfile);

    $cmd = "rGraph = ListPlot[mdata, PlotStyle -> {{PointSize[0.02]}, Black}];";
    mathLinkCmd($cmd, 1, $link, $dumpfile);
    $cmd = "fitGraph = Plot[datafit, {x,$rawData->{MIN},$rawData->{MAX}}, PlotStyle -> " . 
	"{Blue}];";
    mathLinkCmd($cmd, 1, $link, $dumpfile);
    $cmd = "g = Show[{hist, rGraph, fitGraph}, Frame -> True, RotateLabel -> True, " . 
	"PlotLabel->StyleForm[\"Plot of $data->{TITLE} $rec->{TITLE}\", " . 
	"FontWeight->\"Bold\"], Ticks -> IntervalCenters, " . 
	"FrameLabel->StandardForm/\@{\"$rec->{XLABEL}\", \"ycneuqerF\"}, " . 
	"TextStyle -> {FontSize -> 8}, PlotRegion -> {{0, 1}, {0.02, 0.98}}];";
    mathLinkCmd($cmd, 1, $link, $dumpfile);
    $cmd = "Export[\"$data->{PIC}\",g, ImageResolution -> 300];";
    ($error, $result) = mathLinkCmd($cmd, 0, $link, $dumpfile);
    $link->PutFunction("Exit", 0);
    $link->EndPacket;
    $link = ();
    printf "$data->{FILE}:%25s\r", " ";
    return $error;
}

sub extractVals {
    my ($inStr, $vtype) = @_;
    my (%DATA, $i, $counter);

    if ($vtype == 1) {
	$DATA{"TYPE"} = "bonds";
	$DATA{"TITLE"} = "BOND length Distribution (A)";
	$DATA{"XLABEL"} = "BOND Length (A)";
    } elsif ($vtype == 2) {
	$DATA{"TYPE"} = "angles";
	$DATA{"TITLE"} = "ANGLE Distribution (radians)";
	$DATA{"XLABEL"} = "ANGLE (radians)";
    } else {
	$DATA{"TYPE"} = "torsions";
	$DATA{"TITLE"} = "DIHEDRAL ANGLE Distribution (radians)";
	$DATA{"XLABEL"} = "DIHEDRAL ANGLE (radians)";
    }

    $vtype += 2;
    $counter = 0;
    while ($inStr =~ /(\w+)\s+\->\s+(\-?\d+\.*\d*)/g) {
	$DATA{"PARMS"}{$1} = $2;
	$counter++;
    }

    if ($inStr =~ /(\d+\.\d+)\}$/) {
	$DATA{"VARIANCE"} = $1;
	$counter++;
    }
    
    if ($counter < 4) {
	%DATA = ();
    }

    return \%DATA;
}

sub getFitFunc {
    my ($parms, $data) = @_;
    my ($fitFunc, $vars, $avg);

    $avg = $data->{"MAX"} - $data->{"MIN"};
    $avg /= 2;
    $avg += $data->{"MIN"};

    if ($parms->{VTYPE} == 1) {
	$fitFunc = "a0 * Exp[-kb *(x-r0)^2/(2*0.0019872*298)]";
	$vars = "{a0,0.1}, {kb,250},{r0,$avg}";
    } elsif ($parms->{VTYPE} == 2) {
	$fitFunc = "a0 * Exp[-kt *(Cos[x] - Cos[t0])^2/(2*0.0019872*298)]";
	$vars = "{a0,0.1}, {kt,5000},{t0,$avg}";
    } elsif ($parms->{VTYPE} == 3) {
	$fitFunc = "a0 * Exp[-kp *(1 + Cos[periodicity*x - p0])/(2*0.0019872*298)]";
	$vars = "{a0, 0.1}, {kp,75}, {p0,p0Val}";
    }

    return ($fitFunc, $vars);
}

sub createBins {
    my ($data, $vtype) = @_;
    my ($nbins, $range, $bin, $binFile, $i, $max, $histData, $l); 
    my (%BINS, @tmp, $min, $tot, $binData, $binSize, $histCat);

    if ($vtype < 3) {
	$nbins = 30;
	$min = $data->{"MIN"};
	$range = $data->{"MAX"} - $min;
    } else {
	$nbins = 100;
	$min = 0;
	$range = 2 * PI;
    }
    
    $binSize = $range/$nbins;
    $tot = 0;
    $max = -999999999;

    return () if ($range == 0);
    
    for $i (keys %{ $data->{"VALS"} }) {
	$bin = int($nbins * ($i - $min)/$range)  + 1;
	$BINS{$bin}++;
    }

    $tot = 0;
    for $i (keys %BINS) {
	$tot += $BINS{$i};
    }

    $tot /= 100;
    $binData = "";
    $histData = "";
    $histCat = "";
    $data->{"AREA"} = 0;
    #$l = pop @tmp;

    for $i (1 .. $nbins) {
	$BINS{$i} = 0 if (! exists($BINS{$i}));
	$binData .= sprintf("{%.3f,%.3f},", ($min + (($i - 1) * $binSize)), ($BINS{$i}/$tot));
	$max = ($BINS{$i}/$tot) if (($BINS{$i}/$tot) > $max);
	$histData .= sprintf("%.3f,", ($BINS{$i}/$tot));
	$histCat .= sprintf("%.3f,", ($min + (($i - 1) * $binSize) - $binSize/2));
	$data->{"AREA"} += ($BINS{$i}/$tot) * $binSize;
    }
    $histCat .= "" . ($min + (($nbins * $binSize) - $binSize/2));
    chop $binData;
    chop $histData;
    $binData = "{$binData}";
    $data->{"BINMAX"} = $max + 5 * $binSize;
    $data->{"TOT"} = $tot;
    return ($binData, $histData, $histCat);
}
    
sub readDatFile {
    my ($inFile, $vtype) = @_;
    my (%DATA, $max, $min, $total, $val, $factor);

    $min = 99999999;
    $max = -99999999;
    $total = 0;

    open INDATA, $inFile or die "ERROR: Cannot open file $inFile:$!\n";
    while (<INDATA>) {
	chomp;
	if ($_ =~ /^\s+\d+\s+(-?\d+\.\d+)/) {
	    $val = $1;
	    if ($val < 0 and $vtype > 1) {
		$factor = abs(int($1/(2 * PI)));
		$factor++;
		$val += $factor * (2 * PI);
		$val = $val % (2 * PI);
	    }

	    $DATA{"VALS"}{$val} = 1;
	    $max = $val if ($val > $max);
	    $min = $val if ($val < $min);
	    $total++;
	}
    }
    
    if ($total > 30) {
	$DATA{"MAX"} = $max;
	$DATA{"MIN"} = $min;
	$DATA{"TOT"} = $total;
    }

    close INDATA or die "ERROR: Cannot close file $inFile: $!\n";
    return \%DATA;
}

sub init {
    my (@tmp, @DATFILE, @fName, $rec, $i, $PARMS);
    FileTester($parmFile);
    print "Obtaining Bead parameters ....";
    $PARMS = GetMesoParms($parmFile);
    print "Done\n\n";

    if (-d $datLoc) {
	opendir DATDIR, $datLoc or die "ERROR: Cannot access directory $datLoc: $!\n";
	@tmp = grep { /\.dat$/ && -f} map { "$datLoc/$_"} readdir DATDIR;
	closedir DATDIR or die "ERROR: Cannot close directory $datLoc:$!\n";
    } elsif (-e $datLoc and -r $datLoc and -T $datLoc) {
	$tmp[0] = $datLoc;
    } else {
	die "ERROR: No valid data files found: $datLoc:$!\n";
    }

    if (! defined($saveName)) {
	$saveName = "meso.par";
    }

    for $i (@tmp) {
	$rec = ();
	$rec->{FILE} = $i;
	$i = basename($i);
	$i =~ s/\.\w+$//;
	$rec->{PIC} = $i . ".png";
	@fName = split /\-/, $i;
	next if ($#fName < 1);
	($rec->{"ATOMS"}, $rec->{"TITLE"}) = getFFTypes(\@fName, $PARMS);
	next if (! $rec->{"ATOMS"} );
	$rec->{VTYPE} = getVType(\@fName);
	push @DATFILE, $rec;
    }

    die "ERROR: $datLoc is not a valid file(s)\n" if (! @DATFILE);
    return \@DATFILE;
}

sub getFFTypes {
    my ($atomList, $parms) = @_;
    my ($i, $j, %ATOMDATA, $isValid, $rec, $counter, $title);

    $counter = 1;
    for $i (@{ $atomList }) {
	$rec = ();
	if ($i =~ /^(D|H)(\d+)/) {
	    $isValid = 1;
	    $rec = (
		    {
			"NAME"     => $i,
			"ELEMENT"  => $1,
			"FFTYPE"   => $i,
			"MASS"     => 1.0001,
			"NUMBONDS" => 1,
			"RADII"    => 1.00001,
		    }
		    );
	} else {
	    $isValid = 0;
	    for $j (keys %{ $parms->{"BEADS"} }) {
		if ($parms->{"BEADS"}{$j}{"NAME"} eq $i) {
		    %{ $rec } = %{ $parms->{"BEADS"}{$j} };
		    $rec->{"FFTYPE"} = $rec->{"ELEMENT"};
		    $isValid = 1;
		    last;
		}
	    }
	}

	if ($isValid) {
	    $ATOMDATA{$counter} = $rec;
	    $title .= $rec->{"NAME"} . " ";
	    $counter++;
	}
    }

    return (\%ATOMDATA, $title);
}

sub getVType {
    my ($atomList) = $_[0];
    my ($returnVal);

    $returnVal = 1;

    if ($#{ $atomList } == 2) {
	$returnVal = 2;
    } elsif ($#{ $atomList } > 2) {
	$returnVal = 3;
    }

    return $returnVal;
}

sub mathLinkCmd {
    my ($inExpression, $isTerminal, $link, $dumpfile) = @_;
    my ($tmp, $err, @holder);

    print $dumpfile "$inExpression\n";
    $link->PutFunction("EvaluatePacket",1);
    if (! $isTerminal) {
        $link->PutFunction("ToString", 1);
    }
    $link->PutFunction("ToExpression", 1);
    $link->PutString("$inExpression");
    $link->EndPacket;

    while ($link->NextPacket != RETURNPKT) {
        $link->NewPacket;
    }
    if ($isTerminal == 1) {
        $link->NewPacket;
        if (! $link) {
            print "ERROR: " . $link->ErrorMessage . "\n";
        }
    } else {
        $err = $link->ErrorMessage;
        if ($isTerminal == 2) {
            @holder = $link->GetRealList;
        } else  {
            $tmp = $link->GetString;
        }
        if ($err =~ /ok so far/) {
            $err = "";
        }
        if ($isTerminal == 2) {
            return ($err, @holder);
        } else {
            return ($err, $tmp);
        }
    }
}

sub numerically {
    ($a<=>$b);
}

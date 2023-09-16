#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename dirname);
use Getopt::Std qw(getopt);
use General qw(FileTester STDev);

sub init;
sub parseThermoFile;
sub collectStats;
sub writeStats;
sub getThermoData;
sub getFieldName;
sub numerically { ($a<=>$b); }

my ($FILES, $saveFile, $FIELDS, $DATA, $STATS, $molOpts);

$|++;
&init;
$DATA = &parseThermoFiles($FILES, $molOpts);
print "Collecting stats...";
$STATS = collectStats($DATA, $FIELDS);
print "Done\nWriting stats to $saveFile...";
&writeStats($STATS, $saveFile);
print "Done\n";

sub writeStats {
    my ($stats, $outFile) = @_;
    my ($i, @grps, $j, @headers);

    @grps = sort numerically keys %{ $stats };
    @headers =  sort { ($a cmp $b) } keys %{ $stats->{ $grps[0] } };
    open OUTFILE, "> $outFile" or die "ERROR: Cannot write to $outFile: $!\n";
    printf OUTFILE "%-8s", "Group";
    for $i (@headers) {
        printf OUTFILE "%21s", $i;
    }
    print OUTFILE "\n";
    for $i (@grps) {
        printf OUTFILE "%-8d", $i;
        for $j (@headers) {
            printf OUTFILE "%13.2f %7.2f", $stats->{$i}{$j}{AVG}, $stats->{$i}{$j}{STDEV};
        }
        print OUTFILE "\n";
    }
    close OUTFILE;

}

sub collectStats {
    my ($data, $fields) = @_;
    my (@fieldList, @groups, $i, $j, $k, $wFields, $STATS, $dat);

    $wFields = join "|", keys(%{ $fields });
    @groups = sort numerically keys(%{ $data->{0} });
    @fieldList =grep {/^($wFields)_\d+$/i} keys(%{ $data->{0}{ $groups[0] } });
    die "ERROR: No valid field found while searching $wFields!\n"
	if (! @fieldList);
    for $i (@groups) {
	for $j (@fieldList) {
	    $dat = "";
	    for $k (keys %{ $data }) {
		$dat .= $data->{$k}{$i}{$j} . " ";
	    }
	    chop $dat;
	    ($STATS->{$i}{$j}{AVG},$STATS->{$i}{$j}{STDEV},$STATS->{$i}{$j}{TOT}) = STDev($dat);
	}
    }
    
    return $STATS;
}

sub parseThermoFiles {
    my ($files, $molOpts) = @_;
    my ($i, $STATS, $data);

    for $i (0 .. $#{ $files }) {
	print "Parsing thermo file $files->[$i]...\r";
	$data->{$i} = getThermoData($files->[$i], $molOpts);
    }
    print "Parsing thermo files...Parsed " . scalar(@{ $files }) . " files\n";
    return $data;
}

sub getThermoData {
    my ($inFile, $molOpts) = @_;
    my ($group, $ENG, $data, $field, @vals, $i, $j, $isEng);

    open INFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^\s*(\S+)\=\s+(.+)$/) {
	    $group = $i = $j = 0;
	    ($field, $isEng) = getFieldName($1);
	    @vals = split /\s+/, $2;
	    while ($i <= (scalar(@vals) - $molOpts->{fieldlength})) {
		for $j (keys %{ $molOpts->{fields} }) {
		    $vals[$i + $j] /= 4.182 if ($vals[$i + $j] =~ /\-?\d+\.\d+/ and $isEng);
		    next if ($field eq "mols" and $j > 2);
		    $ENG->{$group}{"${field}_${j}"} = $vals[$i + $j];
		    $ENG->{$group}{"${field}_${j}"} = 0 if ($vals[$i + $j] eq "nan");
		}
		$group++;
		$i += $molOpts->{fieldlength};
	    }
	}
    }
    close INFILE;
    die "ERROR: No valid data read from $inFile!\n" if (! defined($ENG));
    return $ENG;
}

sub getFieldName {
    my ($infield) = $_[0];
    my ($isEng) = 0;

    $infield = "S0" if ($infield eq "S(0)(cm/mol/SimBox)");
    $infield = "atoms" if ($infield eq "#_of_atoms_in_group");
    $infield = "mols" if ($infield eq "#_of_molec_in_group");
    $infield = "d0" if ($infield eq "constant__________D");
    $infield = "dof" if ($infield eq "degree__of__freedom");
    $infield = "f" if ($infield eq "fluidicity___factor");
    $infield = $1 if ($infield =~ /^(\w+)/);
    $infield =~ s/_*$//;
    $infield =~ s/temperature/temp/g;
    $infield = "V0" if ($infield eq "Vo");

    $isEng = 1 if ($infield =~ /(aq|sq|eq|v0|zpe|ec|sc|ac|sqd|cvq|cvc|emd)/i);
    return ($infield, $isEng);
}

sub init {
    my (%OPTS, $fileList, $fieldList, $findCmd, $dir, $mOpts);
    getopt('tfsm',\%OPTS);
    die "usage: $0 -t thermo files -f [fields] -s [save file] -m [mol options]\n"
	if (! exists($OPTS{t}));
    print "Initializing...";
    ($fileList, $fieldList, $saveFile, $mOpts) = ($OPTS{t}, $OPTS{f}, $OPTS{s}, $OPTS{m});
    $findCmd = "find " . dirname($fileList) . " -name '" . basename($fileList) . "' -print" if (! -e $fileList);
    $findCmd = "ls $fileList" if (-e $fileList);
    die "ERROR: No valid file found while searching $fileList!\n"
        if (! open(DATFILES, "$findCmd |"));
    while (<DATFILES>) {
        chomp;
        if (-e $_ and -r $_ and -T $_) {
	    $dir = dirname($_);
	    #next if ($dir ne ".");
            push @{ $FILES }, $_;
        }
    }
    close DATFILES;
    die "ERROR: No valid files found in path \"$fileList\"\n"
        if (! defined($FILES));
    #die "ERROR: Need more than 1 data file to compute averages!\n"
        #if (scalar(@{ $FILES } == 1));
    $fieldList = "sq zpe v0 eq s0 aq f d0 cvq emd temp dof"
	if (! defined($fieldList));
    while ($fieldList =~ /(\w+)/g) {
	$FIELDS->{$1} = 1;
    }
    if (! defined($saveFile)) {
	$saveFile = basename($FILES->[0]);
	$saveFile =~ s/\.\w+$//;
	#$saveFile =~ s/_mol$//;
	#$saveFile =~ s/_grps$//;
	$saveFile =~ s/_\d+//;
	$saveFile .= "_stats.dat";
    }
    $mOpts = "3 every 3" if (! defined($mOpts) or $mOpts !~ /^(.+) every (\d+)/i);
    $mOpts =~ /^(.+) every (\d+)/i;
    ($molOpts->{input}, $molOpts->{fieldlength}) = ($1, $2);
    while ($molOpts->{input} =~ /(\d+)/g) {
	$molOpts->{fields}{($1 - 1)} = 1;
    }
    print "Done\n";
}

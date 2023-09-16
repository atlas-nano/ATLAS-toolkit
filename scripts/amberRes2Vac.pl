#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use Packages::FileFormats qw(GetBGFFileInfo sortByRes);
use Packages::General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use strict;

sub init;
sub readAmberResFile;
sub createVACGrpFile;
sub numerically { ($a<=>$b); }

my ($ATOMS, $BONDS, $bgfFile, $RESFILES, $vacFile, $SHELLS, $RES);

$|++;
$RESFILES = &init;
print "Parsing bgf file $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
$RES = sortByRes($ATOMS);
print "Done\n";
$SHELLS = readAmberResFile($RESFILES, $RES);
print "Creating VAC group file $vacFile...";
createVACGrpFile($SHELLS, $vacFile);
print "Done\n";

sub createVACGrpFile {
    my ($shellData, $saveName) = @_;
    my ($i, $count, $j);

    open VACFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    print VACFILE "Total groups " . (scalar(keys %{ $shellData })) . "\n";
    for $i (sort numerically keys %{ $shellData }) {
	$count++;
	print VACFILE "Group $count Atoms $shellData->{$i}{COUNT} # $shellData->{$i}{FILE}\n";
	chop $shellData->{$i}{DATA};
	chop $shellData->{$i}{DATA};
	print VACFILE "$shellData->{$i}{DATA}\n";
    }
    close VACFILE;
}

sub readAmberResFile {
    my ($resFiles, $resData) = @_;
    my ($prnStr, $strLen, $inStr, $atom1, $atom2, %shellData, @tmp, $i);
    my ($count);
    
    for $i (0 .. $#{ $resFiles }) {
	$count = 0;
	$prnStr = "Reading AMBER RES file $resFiles->[$i]...";
	$strLen = length($prnStr);
	print "${prnStr}\r";
	open RESFILE, "$resFiles->[$i]" or die "ERROR: Cannot open $resFiles->[$i]: $!\n";
	while (<RESFILE>) {
	    chomp;
	    while ($_ =~ /\:(\d+\-?\d*)/g) {
		undef($atom1);
		undef($atom2);
		$inStr = $1;
		if ($inStr =~ /(\d+)\-(\d+)/) { #range
		    @tmp = sort numerically keys %{ $resData->{$1}{ATOMS} };
		    $atom1 = $tmp[0];
		    @tmp = sort numerically keys %{ $resData->{$2}{ATOMS} };
		    $atom2 = $tmp[$#tmp];
		} else {
		    @tmp = sort numerically keys %{ $resData->{$inStr}{ATOMS} };
		    $atom1 = $tmp[0];
		    $atom2 = $tmp[$#tmp];
		}
		$prnStr = "${atom1}-${atom2}";
		$count += length($prnStr);
		if ($count > 75) {
		    $prnStr .= "\n";
		    $count = 0;
		} else {
		    $prnStr .= ", ";
		}
		$shellData{$i}{DATA} .= "$prnStr";
		$shellData{$i}{COUNT} += ($atom2 - $atom1) +1;
		$shellData{$i}{FILE} = $resFiles->[$i];
	    }
	}
	close RESFILE;
	printf "%${strLen}s\r", " ";
    }

    die "ERROR: No valid Amber res data found in any of the files\n" if (! keys %shellData);
    print "Reading AMBER Res file(s)...Done\n";
    
    return \%shellData;
}

sub init {
    my (%OPTS, $resLoc, @FILES);
    
    getopt('brs',\%OPTS);
    for ("b", "r") {
	die "usage: $0 -b bgf file -r residue file(s) -s (save name)\n"
	    if (! defined($OPTS{$_}));
    }

    print "Initializing...";
    ($bgfFile, $vacFile, $resLoc) = ($OPTS{b}, $OPTS{s}, $OPTS{r});
    FileTester($bgfFile);

    open FINDCMD, "ls $resLoc |" or die "ERROR: Cannot find $resLoc: $!\n";
    while (<FINDCMD>) {
	chomp;
	push @FILES, $_;
    }
    close FINDCMD;
    
    die "ERROR: No valid res data files found!\n" if (! @FILES);

    if (! defined($vacFile)) {
	$vacFile = basename($bgfFile);
	$vacFile =~ s/\.\w+$//;
	$vacFile .= "_vac.grp";
    }
    print "Done\n";
    
    return \@FILES;
}

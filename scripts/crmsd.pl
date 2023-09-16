#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use General qw( ReadFFs FileTester GetBondLength LoadFFs);
use FileFormats qw(GetBGFFileInfo AddMass);
use Superimpose;
#use Kabash qw(SuperimposeAtomsPDL SuperimposeAtomsMR);
#use Math::MatrixReal;

my ($ffFiles, $FFdata, $ref, $i, $rms, $FILES, @RMSDATA, $saveFile, $useMass);

$|++;
$FILES = &init;
$FFdata = LoadFFs($ffFiles, 2);
print "Parsing reference bgf file $ref->{NAME}...";
($ref->{ATOMS}, $ref->{BONDS}) = GetBGFFileInfo($ref->{NAME}, 0);
&AddMass($ref->{ATOMS}, $FFdata);
print "Done\n";
for $i (@{ $FILES }) {
    print "$i->{NAME}: ";
    ($i->{ATOMS}, $i->{BONDS}) = GetBGFFileInfo($i->{NAME}, 0);
    if (scalar(keys %{ $ref->{ATOMS} }) != scalar(keys %{ $i->{ATOMS} })) {
		  print "different total atoms! Skipping...\n";
		  next;
    }
    &AddMass($i->{ATOMS}, $FFdata);
    &Kabash::SuperimposeAtoms($ref->{ATOMS}, $i->{ATOMS});
    #SuperimposeAtomsPDL($ref->{ATOMS}, $i->{ATOMS});
    $rms = &calcRMS($ref->{ATOMS}, $i->{ATOMS}, $useMass);
    print "rms: $rms...done\n";
    $i->{RMS} = $rms;
}

if (defined($saveFile)) {
    print "Saving CRMSD data to $saveFile...";
    open PRINTHANDLE, "> $saveFile" or die "ERROR: Cannot write to $saveFile: $!\n";
    &writeData($FILES, \*PRINTHANDLE);
    close PRINTHANDLE;
    print "Done\n";
} else {
    &writeData($FILES, \*STDOUT);
}

sub writeData {
    my ($fileData, $OUTDATA) = @_;
    my ($i);

    printf $OUTDATA "%-8s %8s Name\n", "#Number", "CRMSD";
    
    for $i (0 .. $#{ $fileData }) {
		  next if (! exists($fileData->[$i]{RMS}));
		  printf $OUTDATA "%8d %8.3f $fileData->[$i]{NAME}\n", ($i + 1), $fileData->[$i]{RMS};
    }
}

sub calcRMS {
    my ($ref, $curr, $massOpt) = @_;
    my ($rms, $totMass, $i, $dist, $mass);

    $rms = $totMass = $mass = 0;

    for $i (keys %{ $ref }) {
		  $dist = GetBondLength($ref->{$i}, $curr->{$i});
		  $mass = 1;
		  if (exists($ref->{$i}{MASS}) and $massOpt) {
		      $mass = $ref->{$i}{MASS};
		  }
		  $rms += $mass * $dist**2;
		  $totMass += $mass;
    }
    $rms = sqrt($rms/$totMass);
    return $rms;
}

sub init {
    my (%OPTS, $tmp, @bgfFiles);

    getopt('rfscm',\%OPTS);
    for ("r", "f", "c") {
	die "usage: $0 -r reference bgf -c cerius2 force field -f " . 
	    "\"bgf1 bgf2...\" -s [save file] -m [mass weighted crms (default 0)]\n" 
	    if (! defined($OPTS{$_}));
    }
    print "Initializing...";
    ($ref->{NAME}, $saveFile, $useMass) = ($OPTS{r}, $OPTS{s}, $OPTS{m});
    
    $useMass = 0 if (! defined($useMass)); 
    if ($useMass =~ /^1|yes|true/i) {
		  $useMass = 1;
    }
    FileTester($ref->{NAME});
    ($ffFiles, undef) = ReadFFs($OPTS{c}, undef);
    #$ffFiles = findFiles($OPTS{c});
    die "ERROR: No valid forcefield files found while searching $OPTS{c}!\n" if (! @{ $ffFiles });
    $tmp = findFiles($OPTS{f});
    die "ERROR: No valid bgf files found while searching $OPTS{f}!\n", if (! @{ $tmp });
    for (@{ $tmp }) {
		push @bgfFiles, ( { "NAME" => $_, });
    }
    print "Done\n";
    return \@bgfFiles;
}

sub findFiles {
    my ($fileList) = $_[0];
    my (@tmp, $i, $findCmd, @FILES, $file);

    @tmp = split /\s+/, $fileList;
    
    for $i (@tmp) {
		$file = basename($i);
		$i =~ s/$file$//;
		$findCmd = "find $i -name '$file' -print";
		if (open(FINDCMD, "$findCmd |")) {
			while (<FINDCMD>) {
			  chomp;
				if (-e $_ and -r $_ and -T $_) {
					push @FILES, $_;
				}
			}
	    close FINDCMD;
	  	}
    }
    return \@FILES;
}

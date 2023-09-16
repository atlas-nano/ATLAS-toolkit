#!/usr/bin/perl -w

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub getParms;
sub writeMacro;

my ($bgfFile, $saveName, $rotSym);
my ($name);
$|++;
&init;
print "Getting params...";
$name = getParms($bgfFile);
print "Done\nWriting Polygraf macro $saveName...";
&writeMacro($bgfFile, $name, $rotSym, $saveName);
print "Done\n";

sub writeMacro {
    my ($bFile, $descrp, $rSymm, $outFile) = @_;

    open OUTDATA, "> $outFile" or die "ERROR: Cannot create $outFile: $!\n";
    print OUTDATA <<DATA;
beginmacro
%
% version      : 3.21
% version date : 21:49:53 4/30/93
% link date    : 11:43:25 8/14/96
%
% Macro created on  2/17/11    0:36:46
%
Top menu/in-out
   In-Out/read
   File types/BioDesign
     "$bFile"
   In-Out/return
Top menu/simulate
   Simulate/defaults
      Defaults/energy var
         Energy var/use all invers
         Energy var/return
      Defaults/return
   Simulate/setup eex
     $descrp
     return
     $descrp
     return
     $descrp
     return
   Simulate/mechanics
      Mechanics/# of steps/atms
        "5000"
      Mechanics/minimize
      Mechanics/return
   Simulate/vibrate
      not sep trn rot
      select proprty
      thermochemistry
      units
      init
        "0.15"
      incr
        "1"
      # incr
        "400"
      rotn symm numb
        "$rSymm"
      return
      return
      calculat modes
      thermochemstry
      return
   Simulate/return
Top menu/exit
  "OK"
%
endmacro
DATA

    close OUTDATA;
}

sub getParms {
    my ($infile) = $_[0];
    my ($descrp);

    open BGFFILE, $infile or die "ERROR: Cannot open BGF file $infile: $!\n";
    while(<BGFFILE>) {
	chomp;
	if ($_ =~ /DESCRP (.*)/) {
	    $descrp = $1;
	    $descrp = substr($descrp, 0, 8) if(length($1) > 8);
	}
    }
    close BGFFILE;
    die "ERROR: Cannot find DESCRP file in BGF file!\n" if (! defined($descrp));
    return $descrp;
}
		
sub init {
    my (%OPTS);

    getopt('brs',\%OPTS);

    for ("b") {
	die "usage: $0 -b bgf file -r (rot symm number) -s (savename)\n" 
   	    if (! exists($OPTS{$_}));
    }

    print "Initializing...";
    ($bgfFile, $rotSym, $saveName) = ($OPTS{b}, $OPTS{r}, $OPTS{s});
    for ($bgfFile) {
	die "ERROR: Cannot find $_: $!\n"
	    if (! -e $_ or ! -r $_ or ! -T $_);
    }
    $rotSym = 1 if (! defined($rotSym));
    if (! defined($saveName)) {
	$saveName = basename($bgfFile);
	$saveName =~ s/\.\w+$//;
	$saveName .= ".polygraf.macro";
    }
    print "Done\n";
}

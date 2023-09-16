#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester LoadElements);
use AMBER qw(parseAmberFF);
use CERIUS2 qw(saveCeriusFF);
use File::Basename;

sub findElement;
sub addFields;

die "usage: $0 amber_ff1 amber_ff2 ...\n"
    if (! @ARGV);

my (@amberfiles) = @ARGV;

my ($amberFF, $ff, $i);

my ($ELEMENTS) = LoadElements();
for $amberFF (@amberfiles) {
    FileTester($amberFF);
    $ff = parseAmberFF($amberFF, $ff);
}
addFields($ff->{atoms}, $ELEMENTS);
my ($save) = basename($amberfiles[0]);
$save =~ s/\.\w{3}$//;
$save .= ".ff";
print "Saving $save...";
saveCeriusFF($ff, $save, $ELEMENTS);
print "Done\n";

sub addFields {
   my ($atoms, $elements) = @_;
   my ($i, $element);


   for $i (keys %{ $atoms }) {
	if ($i =~ /IP|IB/) {
	    $element = 11;
    } elsif ($i =~ /^Mg/i) {
        $element = 12;
	} elsif ($i eq "IM") {
	    $element = 19;
	} elsif ($i =~ /EP|X|Lp/i) {
		$element = 1;
    } else {
	    $element = findElement($atoms->{$i}{VALS}[0]{mass}, $elements);
	    die "ERROR: Cannot find element for type $i"
			if(!defined($element));
	}
	$atoms->{$i}{VALS}[0]{element} = $element;
	$atoms->{$i}{VALS}[0]{hybrid} = 0;
    }
}
	   
sub findElement {
   my ($sMass, $elements) = @_;
   my ($i, $elementNum, $eMass);

   for $i (keys %{ $elements }) {
	$eMass->[0] = $elements->{$i}{MASS}-.1;
	$eMass->[1] = $elements->{$i}{MASS}+.1;
	if ($sMass>$eMass->[0] and $sMass<$eMass->[1]) {
            $elementNum = $i;
            last;
        }
   }
   return $elementNum;
}


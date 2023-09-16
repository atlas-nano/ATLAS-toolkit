#!/usr/bin/perl -w
BEGIN {
    push @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General qw(FileTester LoadElements);
use Packages::AMBER qw(parseAmberFF);
use Packages::CERIUS2 qw(saveCeriusFF);
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
my ($save) = basename($amberfiles[$#amberfiles]);
$save =~ s/\.\w{3}$//;
$save .= ".ff";
saveCeriusFF($ff, $save, $ELEMENTS);

sub addFields {
   my ($atoms, $elements) = @_;
   my ($i, $element);


   for $i (keys %{ $atoms }) {
	if ($i eq "IP") {
	    $element = 11;
	} elsif ($i eq "IM") {
	    $element = 19;
        } else {
	    $element = findElement($i, $elements);
    	    if (! $element && length($i) > 1) {
	    	$element = findElement(substr($i,0,1),$elements);
	    	$element = 1 if (! $element);
	    }
	}
	$atoms->{$i}{VALS}[0]{element} = $element;
	$atoms->{$i}{VALS}[0]{hybrid} = 0;
    }
}
	   
sub findElement {
   my ($elementName, $elements) = @_;
   my ($i, $elementNum);

   for $i (keys %{ $elements }) {
	if (uc($elements->{$i}{SYMBOL}) eq uc($elementName)) {
            $elementNum = $i;
            last;
        }
   }
   return $elementNum;
}


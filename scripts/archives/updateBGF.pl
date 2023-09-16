#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats;
use Packages::General;
use Getopt::Std qw(getopt);

sub main;
sub init;
sub updateBGF;

my ($refBGF, $modBGF, $saveName, $field, $newBox);

$|++;
&main;

sub main {
    my ($refATOMS, $refBONDS, $HEADERS, $modATOMS, $tmp, $FIELDS, $HEADERS1);

    $FIELDS = &init;
    print "Obtaining reference BGF Info from $refBGF...";
    ($refATOMS, $refBONDS, $HEADERS) = GetBGFFileInfo($refBGF, 1);
    print "Done\nObtaining coordiantes/charges from $modBGF...";
    ($modATOMS, $tmp, $HEADERS1) = GetBGFFileInfo($modBGF, 1);
    print "Done\nUpdating \"$field\"...";
    updateBGF($refATOMS, $modATOMS, $FIELDS);
    print "Done\nCreating BGF file $saveName...";
    if ($newBox) {
	addHeader($refATOMS, $HEADERS1);
    } else {
	addHeader($refATOMS, $HEADERS);
    }
    createBGF($refATOMS, $refBONDS, $saveName, 1);
    print "Done\n";
}

sub updateBGF {
    my ($ref, $mod, $fieldList) = @_;
    my ($atomC, $modName, $refName, $atom, $field);

    for $atomC (keys %{ $ref }) {
	die "ERROR: Atom # $atomC ($ref->{$atomC}{ATMNAME}) is not present in reference file\n"
	    if (! exists($mod->{$atomC}));
	$modName = Trim($mod->{$atomC}{"ATMNAME"});
	$refName = Trim($ref->{$atomC}{"ATMNAME"});

	#print "WARNING: Atoms # $atomC have different names! REF: $refName, MOD: $modName\n"
	#    if ($modName ne $refName);
	$atom = $ref->{$atomC};
	for $field (@{ $fieldList }) {
	    next if (! exists($atom->{$field}));
	    $atom->{$field} = $mod->{$atomC}{$field};
	}
    }
}
	
sub init {
    my (%OPTS, @FIELDS, $i);

    getopt('brts',\%OPTS);
    ($refBGF, $modBGF, $field, $saveName) = ($OPTS{b},$OPTS{r},$OPTS{t},$OPTS{s});
    for ($refBGF, $modBGF, $field) {
	die "usage: $0 -b updateBGF -r referenceBGF -t \"fields\" -s [saveBGF]\n"
	    if (! defined($_));
    }
    
    print "Initializing...";
    FileTester($refBGF);
    FileTester($modBGF);

    if (! $saveName) {
	$saveName = $refBGF;
	$saveName =~ s/\.\w+$/_mod\.bgf/;
    }

    if ($field =~ /\s+/) {
	@FIELDS = split /\s+/, $field;
    } else {
	$FIELDS[0] = $field;
    }

    $newBox = 0;
    for $i (@FIELDS) {
	if (lc($i) eq "box") {
	     $newBox =1;
	     last;
	}
    }
    print "Done\n";
    return \@FIELDS;
}
    

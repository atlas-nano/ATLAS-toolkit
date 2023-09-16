#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester Trim);
use Getopt::Std qw(getopt);
use ManipAtoms qw(GetMols AddMolsToSelection SelectAtoms BuildAtomSelectionString);

sub numerically { ($a<=> $b); }

my ($refBGF, $modBGF, $saveName, $field, $newBox, $selection);
my ($refATOMS, $refBONDS, $HEADERS, $modATOMS, $modBONDS, $tmp, $tmp1, $FIELDS, $HEADERS1, $SELECTIONS);

$|++;
$FIELDS = &init;
print "Obtaining reference BGF Info from $refBGF...";
($refATOMS, $refBONDS, $HEADERS) = GetBGFFileInfo($refBGF, 1);
$tmp1 = &GetMols($refATOMS, $refBONDS);
$SELECTIONS = SelectAtoms($selection, $refATOMS);
&AddMolsToSelection($SELECTIONS, $refATOMS);
print "Done\nObtaining @{ $FIELDS } from $modBGF...";
($modATOMS, $modBONDS, $HEADERS1) = GetBGFFileInfo($modBGF, 1);
$tmp = &GetMols($modATOMS, $modBONDS);
die "ERROR: Expected " . scalar(keys %{ $tmp1 }) . " molecule(s) in $modBGF. Got " . scalar(keys %{ $tmp }) . "\n" 
    if (scalar(keys %{ $tmp }) != scalar(keys %{ $tmp1 }));
print "Done\nUpdating \"$field\"...";
&updateBGF($refATOMS, $modATOMS, $FIELDS, $SELECTIONS);
print "Done\nCreating BGF file $saveName...";
if ($newBox) {
    &addHeader($refATOMS, $HEADERS1);
} else {
    &addHeader($refATOMS, $HEADERS);
}
&createBGF($refATOMS, $refBONDS, $saveName);
print "Done\n";

sub updateBGF {
    my ($ref, $mod, $fieldList, $atomSel) = @_;
    my ($i, $j, $k, $f, $molid, $refMolSize, $modMolSize, $molList, $modMolNatoms);

    $modMolNatoms = scalar(keys%{ $mod });
    for $i (keys %{ $atomSel }) {
        $molid = ${ $ref->{$i}{MOLECULEID} };
        next if (exists($molList->{$molid}));
        $molList->{$molid} = 1;
        $modMolSize = ${ $mod->{$i}{MOLSIZE} };
        $refMolSize = ${ $ref->{$i}{MOLSIZE} };
        die "ERROR: Molecule ${ $ref->{$i}{MOLECULEID} } has $refMolSize atoms while the reference has $modMolSize!\n"
            if ($refMolSize != $modMolSize);
        for $j (sort numerically keys %{ $ref->{$i}{MOLECULE}{MEMBERS} }) {
           for $f (@{ $fieldList }) {
                $ref->{$j}{$f} = $mod->{$j}{$f};
            }
        } 
    }
}
	
sub init {
    my (%OPTS, @FIELDS, $i, $atomSel);

    getopt('brtsa',\%OPTS);
    ($refBGF, $modBGF, $field, $saveName, $atomSel) = ($OPTS{b},$OPTS{r},$OPTS{t},$OPTS{s},$OPTS{a});
    for ($refBGF, $modBGF, $field) {
	die "usage: $0 -b updateBGF -r referenceBGF -t \"fields\" -a (atom_selection) -s [saveBGF]\n"
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
    $atomSel = "index>0" if (! defined($atomSel));
    $selection = BuildAtomSelectionString($atomSel);
    print "Done\n";
    return \@FIELDS;
}

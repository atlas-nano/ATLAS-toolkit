#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::General qw(FileTester);

sub init;
sub execCMD;
my ($amberFile, $saveFile);
my ($scripts);
$|++;
$scripts = &init;
print "Step 1. Converting AMBER bgf $amberFile -> MOL2...";
&execCmd("$scripts->[0] -b $amberFile -m _tmp.mol2");
print "Done\nStep 2. Converting AMBER types to MOl2...";
&execCmd("$scripts->[1] -m _tmp.mol2 -s _tmp_dreiding.mol2");
print "Done\nStep 3. Converting MOL2 to BGF...";
&execCmd("$scripts->[2] -imol2 _tmp_dreiding.mol2 -obgf _tmp_dreding.bgf");
print "Done\nStep 4. Creating $saveFile...";
if (system("$scripts->[3] _tmp_dreiding.mol2 _tmp_dreding.bgf > _bgf_dreiding_typed.bgf")) {
    die "error\n";
}
&execCmd("$scripts->[4] -r _bgf_dreiding_typed.bgf -b $amberFile -t FFTYPE -s ${amberFile}_almost");
&execCmd("$scripts->[5] ${amberFile}_almost tmp-tmp-juju33.bgf");
&execCmd("$scripts->[4] -b tmp-tmp-juju22.bgf -r $amberFile -t CHARGE -s $saveFile");
my $prefix = $saveFile;
$prefix =~ s/\.\w+$//;
system ("mv _tmp_dreiding.mol2 ${prefix}.mol2");
system("rm -fr _tmp.mol2 _tmp_dreding.bgf _bgf_dreiding_typed.bgf ${amberFile}_almost tmp-tmp-juju*");
print "Done\n";

sub execCmd {
    my ($cmdStr) = $_[0];
    $cmdStr .= " >> amber2dreiding_convert.log";
    if (system($cmdStr)) {
	die "ERROR: Check amber2dreidingconvert.log\n";
    }
}

sub init {
    my (%OPTS, @scripts, $i);
    getopt('bs',\%OPTS);
    die "usage: $0 -b amber bgf file -s [save name]\n"
	if (! exists($OPTS{b}));
    print "Initializing...";
    ($amberFile, $saveFile) = ($OPTS{b}, $OPTS{s});
    FileTester($amberFile);
    if (! defined($saveFile)) {
	$saveFile = basename($amberFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_dreiding.bgf";
    }
    @scripts = ("/net/hulk/home3/tpascal/scripts/bgf2mol2.pl", 
		"/net/hulk/home3/tpascal/scripts/amber2mol2Types.pl",
		"/net/hulk/home4/maiti/babel-1.6/babel",
		"/ul/tpascal/scripts//clean_bgf.pl", 
		"/ul/tpascal/scripts/updateBGF.pl",
		"/ul/victor/utilities/bgf_tools/Lingraf_post_process.pl");
    for $i (0 .. $#scripts) {
	if ($i != 2) {
	    FileTester($scripts[$i]);
	} else {
	    die "ERROR: $scripts[$i] is not valid!\n" if (! -e $scripts[$i]);
	}
    }
    print "Done\n";
    return \@scripts;
}

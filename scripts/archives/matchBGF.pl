#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename;
use Packages::FileFormats qw(GetBGFFileInfo createBGF addHeader);
use Packages::General qw(FileTester);

sub init;
sub updateAtoms;
sub showUsage;
sub numerically { ($a<=>$b); }

my ($refATOMS, $refBONDS, $refBGF,$saveBGF); 
my ($modBGF, $modATOMS, $modBONDS, $HEADERS);
my ($updatedBGF, $updatedBONDS);

$|++;
&init;
print "Parsing reference BGF $refBGF...";
($refATOMS, $refBONDS) = GetBGFFileInfo($refBGF, 0);
print "Done\nParsing modified BGF $modBGF...";
($modATOMS, $modBONDS, $HEADERS) = GetBGFFileInfo($modBGF, 1);
print "Done\nUpdating atom numbers in modBGF to match those in refBGF...";
($updatedBGF, $updatedBONDS) = updateAtoms($modATOMS, $modBONDS, $refATOMS);
print "Done\nCreating $saveBGF...";
addHeader($updatedBGF, $HEADERS);
createBGF($updatedBGF, $updatedBONDS, $saveBGF);
print "Done\n";

sub updateAtoms {
    my ($mod, $mBonds, $ref) = @_;
    my (%MAP, $i, $atmName, %UATOMS, %UBONDS, $j);

    $j = 0;
    for $i (keys %{ $ref }) {
	$atmName = $ref->{$i}{ATMNAME};
	$MAP{$atmName} = $i;    
	$j = $i if ($i > $j);
    }

    for $i (sort numerically keys %{ $mod }) {
	$atmName = $mod->{$i}{ATMNAME};
	if (exists($MAP{$atmName})) {
	    %{ $UATOMS{$MAP{$atmName}} } = %{ $mod->{$i} };
	    $UATOMS{$MAP{$atmName}}{OINDEX} = $i;
	    $mod->{$i}{UINDEX} = $MAP{$atmName};
	} else {
	    $j++;
	    %{ $UATOMS{$j} } = %{ $mod->{$i} };
	    $UATOMS{$j}{OINDEX} = $i;
	    $mod->{$i}{UINDEX} = $j;
	}
    }
    
    for $i (sort numerically keys %UATOMS) {
	for $j (@{ $mBonds->{ $UATOMS{$i}{OINDEX} } }) {
	    push @{ $UBONDS{$i} }, $mod->{$j}{UINDEX};
	}
    }

    return (\%UATOMS, \%UBONDS);
}

sub init {
    my (%OPTS, $usage);
    getopt('rms',\%OPTS);
    ($refBGF, $modBGF, $saveBGF) = ($OPTS{r},$OPTS{m},$OPTS{s});
    $usage = &showUsage;
    for ($refBGF, $modBGF) {
	die "$usage\n" if (! defined($_));
    }

    print "Initializing...";
    FileTester($modBGF);
    FileTester($refBGF);
    if (! defined($saveBGF)) {
	$saveBGF = basename($modBGF);
	$saveBGF =~ s/\.\w+$/_mod\.bgf/;
    }
    print "Done\n";
}

sub showUsage {
    my ($usage) = "usage: $0 -r reference bgf -m modified bgf -s [savename]\n" . 
	"Options:\n\t-r reference bgf: The location of the bgf template\n\t" .
	"-m modified bgf: The location of the bgf to match to the reference\n\t" . 
	"-s [savename]: (Optional) The name to save the updated (matched) bgf as\n";
    return $usage;
}

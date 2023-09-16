#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::Math::MatrixReal;
use Packages::FileFormats;
use Packages::CERIUS2;
use Packages::General;
use Packages::Superimpose;

sub init;
sub main;
sub usage;
sub GetAvgCoords;

&usage if (! @ARGV or $#ARGV < 2);

$|++;

my ($bgf1, $bgf2, $cerius2ff, $save, $selection) = @ARGV;
my (@movFiles, @COORDS);

&main;

sub main {
    my ($refBGF, $refBONDS, $matchBGF, $HEADER, $PARMS, $file);
    my ($tmpAtoms, $tmpBonds, $saveBGF, $matchBONDS, $BONDS2, $MASSES, $tmp);

    print "Initializing...";
    &init;
    print "Done\nGetting CERIUS2 forcefield information...";
    $PARMS = parseCerius2FF($cerius2ff);
    $MASSES = GetMass($PARMS);
    print "Done\nGetting reference BGF File Information from $bgf1...";
    ($tmpAtoms, $tmpBonds, $HEADER) = GetBGFFileInfo($bgf1, 1);
    ($refBGF, $refBONDS, $tmp) = GetBGFAtoms($selection, $tmpAtoms, $tmpBonds);
    die "ERROR: No valid atoms found in selection\n" if (! keys %{ $refBGF });
    print "Done\n";
    push @COORDS, $refBGF;
    for $file (@movFiles) {
	print "Superimposing $file onto ref..\r";
	($tmpAtoms, $tmpBonds) = GetBGFFileInfo($file, 0);
	($matchBGF, $matchBONDS, $tmp) =  GetBGFAtoms($selection, $tmpAtoms, $tmpBonds);
	die "ERROR: No valid atoms found in selection\n" if (! keys %{ $matchBGF });
	SuperimposeAtoms($refBGF, $matchBGF, $MASSES);
	push @COORDS, $matchBGF;
    }
    print "\nCreating average BGF file $save...";
    $saveBGF = GetAvgCoords(\@COORDS);
    addHeader($matchBGF, $HEADER);
    createBGF($matchBGF, $BONDS2, $save);
    print "Done\n";
}

sub GetAvgCoords {
    my ($FILES) = $_[0];
    my ($BGF, $atomC, $total, $i, $dim);

    %{ $BGF } = %{ $FILES->[0] };

    $total = $#{ $FILES };
    for $atomC (keys %{ $FILES->[0] }) {
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    for $i (1 .. $total) {
		$BGF->{$atomC}{$dim} += $FILES->[$i]{$atomC}{$dim};
	    }
	    $BGF->{$atomC}{$dim} /= ($total + 1);
	}
    }

    return $BGF;
}
	    
sub init {
    my (@tmp, $i);
    FileTester($bgf1);
    FileTester($cerius2ff);
    
    if (! defined($save)) {
        $save = $bgf1;
        $save =~ s/\.\w+$/_avg.bgf/;
    }

    if (! defined($selection)) {
	$selection = "*";
    } else {
	print "Parsing atom/residue selection...";
	if ($selection =~ /\s/) {
	    @tmp = split /\s+/, $selection;
	} else {
	    @tmp = ($selection);
	}
	$selection = GetSelections(\@tmp);
    }

    if (! -d $bgf2) {
	FileTester($bgf2);
	push @movFiles, $bgf2;
    } else {
	opendir BGFFILES, $bgf2 or die "ERROR: Cannot access directory $bgf2: $!\n";
	@tmp = grep { /\.bgf$/ && -f} readdir BGFFILES;
	close BGFFILES;
	for $i (@tmp) {
	    next if ($i eq $bgf1);
	    push @movFiles, $i;
	}
	@movFiles = map { "$bgf2/$_" } @movFiles;
	die "ERROR: No valid match BGF files found in $bgf2\n" if (! @movFiles);
    }

}


sub usage {
    print STDOUT <<DATA;
usage: ref_bgf match_bgf cerius2ff [save_name] [options]
  ref_bgf: name of reference bgf file
  match_bgf: name of matching bgf file| directory of files
  cerius2ff: name of cerius2 force field for bgf files
  save_name: name of file to save (optional)
  options: (default is all)
    IaX - atom number X
    IrX - residue index X
    TaX - atom type X
    NaX - atom name X
    NrX - residue name X
   : - range, eg. :Tr1-8 :Ia3-66 
DATA

die "\n";

}

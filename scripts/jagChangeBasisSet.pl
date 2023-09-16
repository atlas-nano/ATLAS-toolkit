#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo);
use General qw(GetSelections FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub parseJagInFile;
sub writeFile;

my ($jagInFile, $basisSet, $sname);
my ($DATA);

$|++;
&init;
print "Getting Jaguar input options from $jagInFile...";
$DATA = parseJagInFile($jagInFile, $basisSet);
print "Done\nCreating Jaguar input file $sname...";
&writeFile($DATA, $sname);
print "Done\n";

sub writeFile {
    my ($data, $outfile) = @_;

    open JAGFILE, "> $outfile" or die "ERROR: Cannot create $outfile: $!\n";
    print JAGFILE $data;
    close JAGFILE;
}

sub parseJagInFile {
    my ($inFile, $basis) = @_;
    my ($start, $optData, $inStr);

    $start = 0;
    open INFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
  LINE: while (<INFILE>) {
      chomp;
      $inStr = $_;
      if ($inStr =~ /^\&zmat/) {
	  $start = 1;
	  $optData .= "$inStr\n";
      } elsif ($inStr =~ /^\&/ and $start) {
	  $optData .= "$inStr\n";
	  last LINE;
      } else {
	  next if ($inStr =~ /igeopt/);
	  next if ($inStr =~ /MAEFILE/);
	  next if ($inStr =~ /entry/);
	  $inStr = "basis=${basis}" if ($inStr =~ /basis/);
	  $optData .= "$inStr\n";
      }
  }
    close INFILE;
    die "ERROR: No valid data read!\n" if (! $optData);
    return $optData;
}

sub init {
    my (%OPTS);
    getopt('bis',\%OPTS);
    for ("i", "b") {
	die "usage: $0 -i jaguar input file -b basis set -s (save name)\n" if (! exists($OPTS{$_}));
    }
    ($jagInFile, $basisSet, $sname) = ($OPTS{i}, $OPTS{b}, $OPTS{s});
    print "Initializing...";
    FileTester($jagInFile);
    if (! defined($sname)) {
	$sname = basename($jagInFile);
	$sname =~ s/\.\w+$//;
	$sname .= "_${basisSet}.in";
	$sname =~ s/\*/s/g;
	$sname =~ s/\+/p/g;
	$sname =~ s/\(//g;
	$sname =~ s/\)//g;
    }
    print "Done\n";
}


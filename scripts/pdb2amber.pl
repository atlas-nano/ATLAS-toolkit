#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

my ($leaprc, $pdbFile, $prefix, $RMAP);

$|++;
&init;
print "Creating leaprc file using $leaprc...";
&createLeapFile($leaprc, $prefix, $pdbFile);
print "Done\nCreating AMBER $prefix files from $pdbFile...";
&runLeap($prefix);
system("rm -fr leaprc __test.pdb");
print "Done\n";

sub runLeap {
    my ($mol) = $_[0];
    if (system "$ENV{AMBERHOME}/bin/tleap > ${mol}.out") {
	print "Error Occurred. Email ${mol}.out file to tpascal\@ucsd.edu. Here is what it contains:\n";
	#system "cat ${mol}.out";
    }
}

sub createLeapFile {
    my ($leapFile, $molname, $pdbfile) = @_;

    die "ERROR: Cannot copy $leapFile to current directory!\n" if (system("cp $leapFile ./leaprc"));
	#die "ERROR: Cannot file leapfiles!" if(system("cp $ENV{AMBERHOME}/dat/leap/cmd/leaprc.protein.ff03.r1 ./leaprc"));
	#die "ERROR: Cannot file leapfiles!" if(system("cat $ENV{AMBERHOME}/dat/leap/cmd/leaprc.DNA.bsc1 > ./leaprc"));
    open MYLEAPRC, ">> leaprc" or die "Cannot write to leaprc: $!\n";
    print MYLEAPRC "gaff = loadamberparams gaff.dat\n";
	print MYLEAPRC "clearPdbResMap\n";	
    print MYLEAPRC "prot = loadpdb $pdbfile\n";
    print MYLEAPRC "saveamberparm prot $molname" . ".prmtop $molname" . ".inpcrd\n";
    print MYLEAPRC "quit\n";
    close MYLEAPRC;
}
    
sub init {
    my (%OPTS);

    getopt('lps',\%OPTS);
    die "usage: $0 -p pdb file -s (save prefix) -l (leaprc location)\n" if (! exists($OPTS{p}));
    print "Initializing...";
    ($pdbFile, $leaprc, $prefix) = ($OPTS{p}, $OPTS{l}, $OPTS{s});
    FileTester($pdbFile);
    $ENV{AMBERHOME} = "/home/tpascal/codes/amber18/" if (!exists($ENV{AMBERHOME}));
    $leaprc = "$ENV{AMBERHOME}/dat/leap/cmd/leaprc.DNA.bsc1" if (! defined($leaprc));
    if (! -e $leaprc or ! -r $leaprc or ! -T $leaprc) {
		print "invalid amber leaprc file... using default...";
        $leaprc = "$ENV{AMBERHOME}/dat/leap/cmd/leaprc.DNA.bsc1";
    }
    $prefix = basename($pdbFile) if (! defined($prefix));
    $prefix =~ s/\.\w+$//;
    print "Done\n";
	print "Using Leaprc: $leaprc\n";
}


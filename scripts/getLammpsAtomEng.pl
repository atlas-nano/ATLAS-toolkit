#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub createLammpsInput;
sub getLammpsEng;

my ($bgffile, $forcefields, $savename);

$|++;
&init;
print "Creating LAMMPS files using $bgffile and \"$forcefields\"...";
&createLammpsInput($bgffile, $forcefields, $savename);
print "Done\nSaving LAMMPS atom energy to $savename...";
&getLammpsEng($savename);
print "Done\n";

sub getLammpsEng {
    my ($sname) = $_[0];
    my ($fname, $lmpcmd, $tot, $savecmd, $outStr);

    $outStr = "";
    for ("data","in") {
	$fname = $_ . ".${sname}";
	die "ERROR: Cannot access $fname: $!\n" if (! -e $fname ||! -r $fname || ! -T $fname);
    }
    $fname = "in.${sname}_singlepoint";
    open INFILE, $fname || die "ERROR: Cannot read from $fname: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^pair_style/) {
	    $outStr .= "pair_style lj/cut 10\n";
	} elsif ($_ =~/^kspace_style/) {
	    $outStr .= "kspace_style none\n";
	} else {
	    $outStr .= "$_\n";
	}
    }
    close INFILE;
    $outStr .= <<DATA;
compute         atomEng all pe/atom
dump            4 all custom 10000 ${sname}.atom.eng id c_atomEng
run             0
DATA

    die "ERROR: Cannot access $fname: $\n" if (! -e $fname ||! -r $fname || ! -T $fname);
    open OUTFILE, "> $fname" || die "ERROR: Cannot write to $fname: $!\n";
    print OUTFILE $outStr;
    close OUTFILE;

    $lmpcmd = "/ul/tpascal/openmpi/1.4.2/64/bin/mpirun -np 1 /ul/tpascal/programs/bin/lmp_openmpi64 -in $fname -screen none";
    die "ERROR while executing \"$lmpcmd\"\n" if(system($lmpcmd));
    $tot = `cat ${sname}.atom.eng | wc | awk '{print \$1}'`;
    $tot-=9;
    $savecmd = "tail -${tot} ${sname}.atom.eng  > ${sname}";
    die "ERROR while executing \"$savecmd\"\n" if (system($savecmd));
    system("rm -fr data.${sname} in.${sname} in.${sname}_singlepoint ${sname}_lammps.script ${sname}.atom.eng");
}

sub createLammpsInput {
    my ($bgf, $ff, $sname) = @_;
    my ($scriptfile);

    $scriptfile = "/ul/tpascal/scripts/createLammpsInput.pl -b $bgf -f \"$ff\" -s $sname";
    die "ERROR: Cannot execute \"$scriptfile\"\n" if(system("${scriptfile} >& /dev/null"));
}

sub init {
    my (%OPTS);

    getopt('bfs',\%OPTS);
    for ("b", "f") {
	die "usage: $0 -b bgffile -f \"forcefield(s)\" -s (savename)"
	    if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($bgffile, $forcefields, $savename) = ($OPTS{b}, $OPTS{f}, $OPTS{s});
    FileTester($bgffile);
    if (! defined($savename)) {
	$savename = basename($bgffile);
	$savename =~ s/\.\w+$//;
	$savename .= "_atom_eng.dat";
    }
    $ENV{LD_LIBRARY_PATH} = "/ul/tpascal/programs/lib/ifcore/64";
    print "Done\n";
}

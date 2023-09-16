#!/usr/bin/perl -w

# This program will execute a tleap and prepare the structure for simulation
# with sander
#
# usage: prep_sander.pl pdbfile molname home_dir

($pdbfile, $molname, $home_dir) = @ARGV;

if (! $pdbfile or ! $molname or ! $home_dir) {
    die "usage: prep_sander.pl pdbfile molname home_dir\n";
}

-e $pdbfile or die "Cannot find $pdbfile: $!\n";
-d $home_dir or die "Cannot find $home_dir: $!\n";
-e $home_dir . "/myleaprc" or die "Cannot find leaprc file\n";

$my_cmd = "cp $home_dir/myleaprc leaprc";
system $my_cmd;

open MYLEAPRC, ">> leaprc" or die "Cannot write to leaprc: $!\n";
print MYLEAPRC "gaff = loadamberparams gaff.dat\n";
print MYLEAPRC "prot = loadpdb $pdbfile\n";
#print MYLEAPRC "chi = loadmol2 /ul/tpascal/OmpA/chitobiose_md/chitobiose_region1_amber.mol2\n";
#print MYLEAPRC "dna = combine {prot chi}\n";
#print MYLEAPRC "solvatebox dna TIP3PBOX 10\n";
#print MYLEAPRC "addions dna Na+ 0\n";
#print MYLEAPRC "setbox dna vdw\n";
print MYLEAPRC "saveamberparm prot $molname" . ".prmtop $molname" . ".rst7\n";
#print MYLEAPRC "savepdb dna $molname" . "_solv.pdb\n";
print MYLEAPRC "quit\n";

close MYLEAPRC;

print "Preparing $pdbfile for sander minimization...";
$ENV{AMBERHOME} = "/home/tpascal/codes/amber12" if (! exists($ENV{AMBERHOME}));
system "$ENV{AMBERHOME}/bin/tleap > junk";

print "Done\n";

#system "rm -f junk leaprc leap.log";
print "Created: $molname" . ".prmtop, $molname" . ".inpcrd\n";

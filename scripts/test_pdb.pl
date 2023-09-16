#!/usr/bin/perl -w
#
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use PDB;

$a=PDB->new();
$a->read(fname => "1AUA.pdb");
$b=PDB->new();
$b->read(fname => "3B7N.pdb");

@lista=$a->select(resnum => [8..221,243..271,276..299], atnam => ["CA"], chain => ["A"]);
@listb=$b->select(resnum => [4..91,98..223,245..273,278..301], atnam => ["CA"], chain => ["A"]);

$rmsd=$a->fit(moving_pdb => $b, moving_atoms => \@listb, ref_atoms => \@lista);

print "RMSD=$rmsd\n";

$b->write(fname => "3B7N_fit.pdb");

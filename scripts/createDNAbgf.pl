#!/usr/bin/perl -w

use strict;

sub ValidateInput();
sub CreateHelix(@);
sub DoCmd(@);

die "usage: $0 save_name.bgf sequence\n"
    if (! @ARGV or $#ARGV < 1);

my ($save_name, $sequence) = @ARGV;
my ($pdb_file);

ValidateInput();
print "Creating DNA Double Helix...";
$pdb_file = CreateHelix($sequence);
print "Done\nConverting PDB->BGF...";
DoCmd("/ul/tpascal/scripts/fixpdb4amber.pl $pdb_file", 1);
DoCmd("/ul/tpascal/scripts/pdbfixforxleap.pl $pdb_file", 1);
DoCmd("/ul/tpascal/scripts/amber2bgf.pl $pdb_file /ul/tpascal/ff/AMBER95.cnv $save_name 0.0 0.0 0.0", 1);
DoCmd("rm -fr $pdb_file _junk leaprc leap.log junk output", 0);
print "Done\nCreated $save_name. Enjoy!\n";

sub ValidateInput() {
    if ($ARGV[0] !~ /(\w+\.bgf)/) {
	die "ERROR: Expected save_name.bgf got $ARGV[0]\n";
    }

    while ($ARGV[1] =~ /(\w)/g) {
	if (lc($1) !~ /[a|t|c|g]/) {
	    die "ERROR: in DNA sequence. Expected [A|T|G|C] got $1\n";
	}
    }
}

sub DoCmd(@) {
    my ($inCmd, $evalCmd) = $_[0];

    if (system("$inCmd >& output") and $evalCmd) {
	die "ERROR executing cmd $inCmd\n";
    }
}

sub CreateHelix(@) {
    my ($inSequence) = $_[0];

    open OUTFILE, "> _junk";
    print OUTFILE "generate d d b $inSequence\n";
    print OUTFILE "write amber _tmp.pdb";
    close OUTFILE;

    DoCmd("/ul/tpascal/scripts/namot_cmd.pl _junk", 0);
    return "_tmp.pdb";
}

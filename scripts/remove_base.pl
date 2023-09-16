#!/usr/bin/perl -w
use strict;
use warnings;

# remove_base.pl - chops of a base from the teminal group of a crossover molecule.
# usage: remove_base.pl pdbfile total_bases number_bases_to_remove which_end

sub ProcessFile();
sub RemoveBases();
sub WriteFile();
sub DeleteBase(@);
sub numerically { ($a<=>$b); }

if (! @ARGV or $#ARGV < 2) {
    die "usage: $0 pdbfile number_bases_to_remove which_end\n";
}

my ($pdbfile, $no_base, $which_end) = @ARGV;
my ($base_file_name, %AtomArray, $no_chains, $tot_base);

-e $pdbfile or die "Cannot locate $pdbfile: $!\n";

if (! $no_base =~ /^\d+$/) {
    die "Invalid number specified for number_bases_to_remove. Expected integer, got $no_base\n";
}

if (! $which_end =~ /^[1|2|3]$/) {
    die "Invalid number specified for which_end. Expected 1, 2 or 3, got $which_end\n";
}

if ($pdbfile =~ /(.+)\.pdb/) {
    $base_file_name = $1;
} else {
    die "Invalid pdb file name: $pdbfile\n";
}

my ($new_file) = $base_file_name . "_nowatna.pdb";

system "cp $pdbfile $new_file";
system "/ul/tpascal/scripts/stripNaH20.pl $new_file";

ProcessFile();
RemoveBases();
WriteFile();

sub ProcessFile() {

    my ($rec, $counter, $base_spec, $total_bases);

    $counter = 0;
    $no_chains = 0;
    $total_bases = 0;

    open PDBFILE, $new_file or die "Cannot open $new_file: $!\n";
    while (<PDBFILE>) {
	chomp;
	if ($_ =~ /^ATOM\s+\d+\s(.+)\s+(\w+[\d+]?)\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
	    $total_bases = $3;
	    $rec = (
		    {
			"ATOMTYPE" => $1,
			"BASETYPE" => sprintf("%-3s", $2),
			"XPOS" => $4,
			"YPOS" => $5,
			"ZPOS" => $6
		    }
		);
	    $base_spec = sprintf("%03d", $3);
	    push @{ $AtomArray{$base_spec} }, $rec;
	    $counter++; 
	    if ($2 =~ /3/ and $#{ $AtomArray{$base_spec} } == 0 ) {
		$no_chains++;
	    }
	} else {
	    print "ERROR: $_\n";
	}
    }
    close PDBFILE;

    $tot_base = $total_bases;

    if ($total_bases == 0) {
	die "Error while reading pdb file $new_file: No valid atoms found.\n";
    }

    if ($no_chains == 0) {
	die "Error while reading pdb file $new_file: No chains found. Found $counter atoms\n";
    }

    print "Found $counter atoms in $no_chains chains and $total_bases bases\n";
			 
}

sub RemoveBases() {

    my ($i, $total_atoms_removed, $curr_base);

    if ($tot_base % $no_chains > 0) {
	die "Invalid total bases specified: $tot_base is not a multiple of $no_chains\n";
    }

    my ($strand_len) = $tot_base/$no_chains;

    $total_atoms_removed = 0;

    if ($which_end == 1 or $which_end == 3) {
	print "Removed $no_base bases from bottom of molecule ";
	for $i (1 .. $no_chains) {
	    if ($i % 2 == 1) {
		$curr_base = ($i - 1) * $strand_len + 1;
		$total_atoms_removed += DeleteBase($curr_base, 0, $no_base);
	    } else {
		$curr_base = $i * $strand_len;
		$total_atoms_removed += DeleteBase($curr_base, 1, $no_base);
	    }
	}
	print "($total_atoms_removed atoms)\n";
    }

    $total_atoms_removed = 0;

    if ($which_end == 2 or $which_end == 3) {
	print "Removed $no_base bases from top of molecule ";
	for $i (1 .. $no_chains) {
	    if ($i % 2 == 0) {
		$curr_base = ($i - 1) * $strand_len + 1;
		$total_atoms_removed += DeleteBase($curr_base, 0, $no_base);
	    } else {
		$curr_base = $i * $strand_len;
		$total_atoms_removed += DeleteBase($curr_base, 1, $no_base);
	    }
	}
	print " ($total_atoms_removed atoms)\n";
    }

}

sub DeleteBase(@) {

    my ($c_base, $is3prime, $no_base) = @_;
    my ($next_base, $return_val, $i);

    $return_val = 0;
    for $i (0 .. $no_base - 1) {
	if ($is3prime) {
	    $c_base = sprintf("%0" . length("$tot_base") . "d", ($_[0] - $i));
	    $return_val += ($#{ $AtomArray{$c_base} } + 1);
	    delete $AtomArray{$c_base}; 
	} else {
	    $c_base = sprintf("%0" . length("$tot_base") . "d", ($_[0] + $i));
	    $return_val += ($#{ $AtomArray{$c_base} } + 1);
	    delete $AtomArray{$c_base}; 
	}
    }

    if ( ! $is3prime) {
	$c_base = $_[0] + $no_base;
	$c_base = sprintf("%0" . length("$tot_base") . "d", $c_base);
	splice @{ $AtomArray{$c_base} }, 0, 3;
	$return_val += 3;
	for $i (0 .. $#{ $AtomArray{$c_base} }) {
	    substr($AtomArray{$c_base}->[$i]->{"BASETYPE"}, -1)  = "5";
	}
    } else {
	$c_base = $_[0] - $no_base;
	$c_base = sprintf("%0" . length("$tot_base") . "d", $c_base);
	for $i (0 .. $#{ $AtomArray{$c_base} }) {
	    substr($AtomArray{$c_base}->[$i]->{"BASETYPE"}, -1)  = "3";
	}
    }

    return $return_val;
}

sub WriteFile() {
    my ($atom_counter, $base_counter, $old_base_counter, $curr_atom, $out_file);
    $atom_counter = $base_counter = 1;

    $out_file = $base_file_name . "_w" . $no_base . "baseremoved.pdb";
    open OUTPDBFILE, "> $out_file" or die "Cannot create $out_file: $!\n";

    for $old_base_counter (sort numerically keys %AtomArray) {
	for $curr_atom (@{ $AtomArray{$old_base_counter} }) {
	    printf OUTPDBFILE "ATOM%7d %-4s%4s%6d%12.3f%8.3f%8.3f\n", $atom_counter, $curr_atom->{"ATOMTYPE"}, $curr_atom->{"BASETYPE"}, $base_counter, $curr_atom->{"XPOS"}, $curr_atom->{"YPOS"}, $curr_atom->{"ZPOS"};
	    $atom_counter++;
	}
	$base_counter++;
    }
    close OUTPDBFILE;
    print "Created file with " . ($atom_counter -1) . " atoms in $no_chains chains\n";
}

#!/usr/bin/perl -w
use strict;
use warnings;

sub GetBaseEnergies();

my (@base_energies, @base_total, $i);

if (! @ARGV or $#ARGV < 1) {
    die "usage: $0 totaleneneygy base_energy_file\n";
}

my ($total_energy, $base_e_file) = @ARGV;

-e $base_e_file or die "Could not locate $base_e_file: $!\n";

if (! $total_energy =~ /^\-?\d+\.\d+$/) {
    die "Invalid total energy, expected decmial\n";
}

GetBaseEnergies();

my ($result);

for $i (0 .. $#base_energies) {
    $result -= ($base_energies[$i] * $base_total[$i]);
}

print "$result\n";

sub GetBaseEnergies() {

    open INFILE, $base_e_file or die "Cannot open $base_e_file: $!\n";
    while (<INFILE>) {
	if ($_ =~ /^\w\s(\-\d+\.\d_)\s(\d+)/) {
	    push @base_energies, $1;
	    push @base_total, $2;
	}
    }
    close INFILE;

    if ($#base_energies < 3) {
	die "$base_e_file is invalid\n";
    }
}

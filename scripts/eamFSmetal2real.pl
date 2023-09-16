#!/usr/bin/perl -w

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub usage;
sub parseEAMFSfile;
sub writeEAMFSfile;

my ($eamfsFile, $saveFile);
my ($DATA);

$|++;
&init;
print "Parsing EAM/FS file ${eamfsFile}....";
$DATA = parseEAMFSfile($eamfsFile);
print "Done\nWriting updated EAM/FS file ${saveFile} with real units...";
&writeEAMFSfile($DATA,$saveFile);
print "Done\n";

sub writeEAMFSfile {
	my ($data, $outFile) = @_;
	my ($i, $j, $k);

	open OUTFILE, "> $outFile" or die "ERROR: Cannot write to $outFile: $!\n";
	for $i (0 .. 4) {
		printf OUTFILE "%s\n", $data->{HEADER}[$i];
	}
	for $i (1 .. $data->{Nelements}) {
		for $j ("eleNum", "mass", "latConst", "latType") { #line 1 = atomic number, mass, lattice constant, lattice type (e.g. FCC)
			printf OUTFILE "%10s ", $data->{$i}{$j};
		}
		printf OUTFILE "\n";
		#embedding function F(rho) (Nrho values)
		for $j (1 .. $data->{Nrho}) {
			printf OUTFILE "%25.16f ", $data->{Frho}{$i}[$j-1]*23.06;
			printf OUTFILE "\n" if(($j%5==0) or ($j==$data->{Nrho}));
		}
		#density function rho(r) for element I at element 1 (Nr values)
		for $j (1 .. $data->{Nelements}) {
			for $k (1 .. $data->{Nr}) {
				printf OUTFILE "%25.16f ", $data->{rho}{$i}{$j}[$k-1];
				printf OUTFILE "\n" if(($k%5==0) or ($k==$data->{Nr}));
			}
		}
	}
	#Following the Nelements sections, Nr values for each pair potential phi(r) array are listed in the same manner 
	for $i (1 .. $data->{Nelements}) {
		for $j ($i .. $data->{Nelements}) {
			for $k (1 .. $data->{Nr}) {
				printf OUTFILE "%25.16f ", $data->{phi}{$i}{$j}[$k-1]*23.06;
				printf OUTFILE "\n" if(($k%5==0) or ($k==$data->{Nr}));
			}
		}
	}
	close OUTFILE;

}

sub parseEAMFSfile {
	my ($infile) = $_[0];
	my ($data, $i, $j, $k, $tmp, $allData);

	open EAMFSFILE, $infile or die "ERROR: Cannot open $infile: $!\n";
	chomp(@{ $tmp } = <EAMFSFILE>);
	close EAMFSFILE;
	#slurm the first 5 lines into the header
	$i=0;
	while($i<5) {
		push @{ $data->{HEADER} }, shift @{ $tmp };
		$data->{HEADER}[ $#{ $data->{HEADER} } ] =~ s/^\s+//;
		$data->{HEADER}[ $#{ $data->{HEADER} } ] =~ s/\s+$//;
		$i++;
	}
	$allData = "@{$tmp}"; $allData =~ s/^\s+//; $allData =~ s/\s+$//;
	@{ $tmp } = split /\s+/, $allData;

	#line 4: Nelements Element1 Element2 … ElementN
	($data->{Nelements}, @{ $data->{elements} } )  = split /\s+/, $data->{HEADER}[3];
	#line 5: Nrho, drho, Nr, dr, cutoff
	($data->{Nrho}, $data->{drho}, $data->{Nr}, $data->{dr}, $data->{cutoff} ) = split /\s+/, $data->{HEADER}[4];

	#Following the header are Nelements sections, one for each element I, each with the following format:
	for $i (1 .. $data->{Nelements}) {
		for $j ("eleNum", "mass", "latConst", "latType") { #line 1 = atomic number, mass, lattice constant, lattice type (e.g. FCC)
			$data->{$i}{$j} = shift @{ $tmp };
		}
		#embedding function F(rho) (Nrho values)
		for $j (1 .. $data->{Nrho}) {
			push @{ $data->{Frho}{$i} }, shift @{ $tmp };
		}
		#density function rho(r) for element I at element 1 (Nr values)
		for $j (1 .. $data->{Nelements}) {
			for $k (1 .. $data->{Nr}) {
				push @{ $data->{rho}{$i}{$j} }, shift @{ $tmp };
			}
		}
	}
	#Following the Nelements sections, Nr values for each pair potential phi(r) array are listed in the same manner 
	for $i (1 .. $data->{Nelements}) {
		for $j ($i .. $data->{Nelements}) {
			for $k (1 .. $data->{Nr}) {
				push @{ $data->{phi}{$i}{$j} }, shift @{ $tmp };
			}
		}
	}

	return $data;
}

sub init {
	my (%OPTS);
	getopt('io',\%OPTS);
	($eamfsFile, $saveFile) = ($OPTS{i},$OPTS{o});
	die "usage: $0 -i EAM/FS file -o [output file]\n"
		if(!defined($eamfsFile));
	print "Initializing...";
	die "ERROR: Cannot access $eamfsFile: $!\n"
		if (! -e $eamfsFile or ! -r $eamfsFile or ! -T $eamfsFile);
	if (! defined($saveFile)) {
		$saveFile = basename($eamfsFile);
		$saveFile =~ s/\.\w+$/.realUnits\.eam.fs/;
	}
}

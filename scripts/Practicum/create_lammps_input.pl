#!/usr/bin/perl -w

use strict;

sub Initialize();
sub ValidateTemplate(@);
sub ParseParam(@);
sub CreateFiles();

die "usage: $0 template parameter_file [save_prefix]\n"
    if (! @ARGV or $#ARGV < 1);

my ($template, $p_file, $save_prefix) = @ARGV;
my (@OUTFILES, @out_data);

Initialize();
print "Validating template file $template...";
ValidateTemplate($template);
print "Done\n";

print "Parsing parameter file $p_file..."
ParseParam($p_file);
print "Done\n";

print "Writing Lammps input files....";
CreateFiles();
print "Done\n";

sub Initialize() {

$save_prefix = "inp_lammps"
    if (! $save_prefix);

die "Cannot locate parameter file $p_file: $!\n"
    if (! -e $p_file);

die "Cannot locate template file $template: $!\n"
    if (! -e $template);

die "Cannot access parameter file $p_file: $!\n"
    if (! -T $p_file or ! -r $p_file);

die "Cannot access template file $template: $!\n"
    if (! -T $template or ! -r $template);

$OUTFILE[0] = (
	       {
		   "suffix"        => "_01",
		   "timestep"      => 0.5,
		   "simlength"     => 1000,
		   "read_data"     => 1,
		   "read_restart"  => 0,
		   "fixstyle1"     => 0,
		   "assignfix"     => [],
		   "tflag"         => 100,
		   "dump"          => 100,
		   "restart"       => 100,
	       }
	       );

$OUTFILE[1] = (
	       {
		   "suffix"        => "_02",
		   "timestep"      => 1.0,
		   "simlength"     => 2000,
		   "read_data"     => 0,
		   "read_restart"  => 1,
		   "fixstyle1"     => 0,
		   "assignfix"     => [],
		   "tflag"         => 100,
		   "dump"          => 100,
		   "restart"       => 100,
	       }
	       );

$OUTFILE[2] = (
	       {
		   "suffix"        => "_04",
		   "timestep"      => 2.5,
		   "simlength"     => 2000,
		   "read_data"     => 0,
		   "read_restart"  => 1,
		   "fixstyle1"     => 0,
		   "assignfix"     => [],
		   "tflag"         => 100,
		   "dump"          => 100,
		   "restart"       => 100,
	       }
	       );

$OUTFILE[3] = (
	       {
		   "suffix"        => "_04",
		   "timestep"      => 5.0,
		   "simlength"     => 45000,
		   "read_data"     => 0,
		   "read_restart"  => 1,
		   "fixstyle1"     => 0.5E-3,
		   "assignfix"     => [1,2,3],
		   "tflag"         => 1000,
		   "dump"          => 1000,
		   "restart"       => 1000,
	       }
	       );

$OUTFILE[0] = (
	       {
		   "suffix"        => "_05",
		   "timestep"      => 5.0,
		   "simlength"     => 6000000,
		   "read_data"     => 1,
		   "read_restart"  => 0,
		   "fixstyle1"     => 0,
		   "assignfix"     => [1,2,3],
		   "tflag"         => 250,
		   "dump"          => 250,
		   "restart"       => 250,
	       }
	       );


}

sub ValidateTemplate(@) {

    my ($in_file) = $_[0];
    my ($rec);

    open INFILE, $in_file or die "Cannot access $in_file: $!\n";
    while (<INFILE>) {
	chomp;

	if ($_ =~ /^(.+)_(.+)$/) {
	    $rec = (
		    {
			"header"    => $1,
			"key"       => $2,
		    }
		    );
	} else {
	    $rec = (
		    {
			"header"    => $_;
			"key"       => "";
		    }
		    );
	}
	push @out_data, $rec;
    }
    close INFILE;
}

sub ParseParm(@) {

    my ($in_file) = $_[0];

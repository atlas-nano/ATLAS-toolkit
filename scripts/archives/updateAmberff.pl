#!/usr/bin/perl -w
# This program will take Cerius2 AMBER cnv file and update it with an AMBER .lib file
use strict;
use File::Basename;

sub initialize;
sub readCnvFile;
sub readLibFile;
sub writeOutput;
sub Trim;


die "usage: $0 conversion_file amber_lib_file [save_name]\n"
    if (! @ARGV or $#ARGV == 0);

my ($cnv_file, $lib_file, $save_file) = @ARGV;
initialize;

print "Reading Amber Conversion file $cnv_file...";
my ($LIB) = readLibFile($lib_file);
print "Done\nReading Amber library file $lib_file...";
my ($outdata) = readCnvFile($cnv_file, $LIB);
print "Done\nCreating update Conversion file $save_file...";
writeOutput($outdata);
print "Done\n";

sub initialize {
    for ($ARGV[0], $ARGV[1]) {
	die "Error accessing file $_: $!\n"
	    if (! -e $_ or ! -r $_ or ! -T $_);
    }
    if (! $save_file) {
	$save_file = basename($cnv_file);
	$save_file =~ s/\.\w{3}$/\.cnv/g;
    }
}

sub readCnvFile {
    my ($inFile, $LibData) = @_;
    my ($in_data, $ln_counter, $is_valid, $outdata, $indata); 
    my ($res, $atom, $curr_charge, $new_charge);
    $is_valid = 0;

    open INFILE, $inFile || die "Cannot load conversion file $inFile: $!\n";
    while (<INFILE>) {
	chomp;
	$indata = $_;
	if ($indata =~ /^\s+(\S+)\s(\s*\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\-?\d+\.\d+)/) {
	    $is_valid = 1;
	    $curr_charge = $6;
	    $res = Trim($1);
	    $atom = Trim($2);
	    if (exists($LibData->{$res})) {
		if (exists($LibData->{$res}{$atom})) {
		    $new_charge = $LibData->{$res}{$atom};
		    if ($curr_charge != $new_charge) {
			print "Updated RES: $res ATOM: $atom. Charge: $curr_charge -> $new_charge\n";
			$indata =~ s/$curr_charge/$new_charge/;
		    }
		}
	    }
	}
	$outdata .= "$indata\n";

    }

    close INFILE;

    die "Unable too parse conversion file $cnv_file\n"
	if (! $is_valid);
    return ($outdata);
}

sub readLibFile {
    my ($inFile) = $_[0];
    my (%LIBDATA, $isValid, $curr_res, $stop, @dat, $pattern);

    $pattern = '^\s+"(\w+)"\s+"\w+"\s+(.+)';
    $isValid = 0;

    open INFILE, $inFile or die "Cannot open AMBER lib file $inFile: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^\!entry\.(\w+).unit.atoms/) {
	    $curr_res = $1;
	    $stop = 0;
	} elsif ($_ =~ /^\!entry\.\w+.unit.atomspertinfo/) {
	    $stop = 1;
	} elsif (! $stop and $_ =~ /$pattern/) {
	    @dat = split /\s+/, $2;
	    if ($dat[5]) {
		$LIBDATA{$curr_res}{$1} = $dat[5];
		$isValid = 1;
	    }
	}

    }

    close INFILE;

    die "Error: $inFile is not a valid AMBER .lib file\n"
	if (! $isValid);

    return \%LIBDATA;
}

sub writeOutput {
    my ($filedata) = $_[0];
    open OUTFILE, "> $save_file" or die "Cannot write AMBER lib file $save_file: $!\n";
    print OUTFILE $filedata;
    close OUTFILE;

}

sub Trim {
    my ($inString) = $_[0];

    for ($inString) {
        s/^\s+//;
        s/\s+$//;
    }

    return $inString;
}

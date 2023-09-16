#!/usr/bin/perl -w
use strict;

sub ProcessFile(@);
sub PrintStats(@);
sub PrintGraphs(@);
sub STDev(@);
sub Numerically;

#    This script will process the file vel_xx.dat and create graphs/averages

die "usage: $0: input_file file_type\n"
    if (! @ARGV or $#ARGV < 1);

my ($input_file, $input_type) = @ARGV;

die "Cannot access file $input_file: $!\n"
    if (! -e $input_file or ! -T $input_file or ! -r $input_file);

my ($DATA) = ProcessFile($input_file, $input_type);

PrintStats($DATA);
#PrintGraphs($DATA);

sub ProcessFile(@) {
    my ($in_file, $input_type) = @_;
    my ($in_type, $in_block, $counter, %VELS, @holder, @headers, $index);

    $in_type = $in_block = 0;
    open INFILE, $in_file or die "Cannot open input file $in_file: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^\s+\#\s+(.+)Velocity Profile for Type\s+(\d+)\s+Block\s+\#\s+(\d+)/) {
	    $in_type = $2;
	    $in_block = $3;
	    @headers = split /\,*\s+/, $1;
	} elsif ($_ =~ /^\s+End of set/) {
	    $in_type = $in_block = 0;
	} elsif ($in_type > 0 and $in_block > 0 and $_ =~ /^\s+(\-?\d+\.\d+E[\+|\-]\d+)(.+)/) {
	    $VELS{$in_type}{$in_block}{"TIMESTEP"} = $1;
	    @holder = split /\s+/, $2;
	    $index = 0;
	    for $counter (@holder) {
		if ($counter =~ /^\-?\d+/) {
		    $VELS{$in_type}{$in_block}{"VALS"}{$headers[$index]} = $counter;
		    $index++;
		}
	    }
	}
    }
    close INFILE;

    die "Error: No Valid data found in $in_file\n"
	if (! %VELS);

    return (\%VELS);
}
	    
sub PrintStats(@) {
    my ($data) = $_[0];
    my ($counter, $index, $key_nm, %STATS);
    my ($average, $stdev, $total);

    printf "%-8s%13s%13s%5s%13s%13s%5s%13s%13s%5s\n", "Type", " ", "X", " "," ", "Y", " "," ", "Z", " ";
    printf "%-8s%13s%13s%5s%13s%13s%5s%13s%13s%5s\n", " ", "Avg", "Std.Dev", "Tot", 
    "Avg", "Std.Dev", "Tot", "Avg", "Std.Dev", "Tot";
    for $counter (sort Numerically keys %{ $data }) {
	%STATS = ();
	printf "%-8d", $counter;
	for $index (keys %{ $data->{$counter} }) {
	    for $key_nm (keys %{ $data->{$counter}{$index}{"VALS"} }) {
		$STATS{$key_nm} .= $data->{$counter}{$index}{"VALS"}{$key_nm} . " ";
	    }
	}
	for $key_nm ("X", "Y", "Z") {
	    chomp $STATS{$key_nm};
	    ($average, $stdev, $total) = STDev($STATS{$key_nm});
	    printf " %12.6G %12.6G %4d", $average, $stdev, $total;
	}
	print "\n";
    }
}

sub STDev(@) {
    my ($in_data) = $_[0];
    my (@datavalues, $n_total, $avg, $result, $i);

    $avg = $result = 0;
    @datavalues = split / /, $in_data;
    $n_total = $#datavalues;
    if ($n_total > -1) {
        $avg = 0.0;
        $result = 0.0;

        foreach $i (@datavalues) {
            $avg += $i;
        }

        $avg = $avg/($n_total + 1);

        foreach (@datavalues) {
            $result += ($_ - $avg) **2;
        }

        if ($n_total ==0) {
            $n_total = 1;
        }


        $result = sqrt($result/$n_total);
    }
    return ($avg, $result, ($n_total + 1));
}

sub Numerically {
    ($a<=>$b);
}

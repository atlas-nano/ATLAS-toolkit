#!/usr/bin/perl -w
#   This script will open a template file, and compare the nearest neighbor
#   energy to that of the files provided, providing a base-base analysis of
#   the strain
use strict;
sub CheckInput();
sub GetFiles(@);
sub GetEnergy(@);
sub CalcDifference(@);
sub Numerically;
sub PrintStats(@);
sub DetermineOptimalSequence(@);
sub GetMax(@);

die "usage: $0 template file|directory [save_name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($template, $loc, $save_name) = @ARGV;
my ($FILES, $REF_Energy, $counter, $curr_Energy, $out_string, $STATS);

CheckInput();
$STATS = ();

print "Getting file(s)...";
$FILES = GetFiles($loc);
print "Done\nObtaining Reference Energies...";
$REF_Energy = GetEnergy($template);
print "Done\nWriting to $save_name...";
open OUTFILE, ">> $save_name" or die "Cannot write to $save_name: $!\n";
printf OUTFILE "%-8s%12s%12s%12s%12s%12s\n","#", "SEQ_NAME", "REF_NAME", "SEQ_ENG", "REF_ENG", "DIFF"; 
for $counter (@{ $FILES } ) {
    $curr_Energy = GetEnergy($counter);
    next
	if (! defined($curr_Energy->{"2"}));
    print OUTFILE "#$counter\n";
    ($STATS, $out_string) = CalcDifference($curr_Energy, $REF_Energy, $STATS);
    print OUTFILE "$out_string\n";
}
close OUTFILE;
print "Done\n";
print "--===STATISTICS===--\n";
$out_string = PrintStats($STATS);

open OUTFILE, "> optimal_sequence.txt" or die "Cannot write to optimal_sequence.txt: $!\n";
print OUTFILE $out_string;
close OUTFILE;

sub CheckInput() {
    die "Error accesing template file $template: $!\n"
	if (! -e $template or ! -r $template or ! -T $template);
    die "Error acessing file/directory $loc: $!\n"
	if (! -e $loc);
    $save_name = "comparison.txt"
	if (! $save_name);

    open OUTFILE, "> $save_name" or die "Cannot write to $save_name: $!\n";
    print OUTFILE "\n";
    close OUTFILE;
}

sub GetFiles(@) {
    my ($location) = $_[0];
    my (@FILES);

    if (! -d $location) {
	$FILES[0] = $location;
    } else {
	opendir (CURRDIR, $loc) or die "Cannot access directory $loc: $!\n";
	@FILES = grep { /^\w+\.anal/ } readdir CURRDIR;
	closedir CURRDIR;
    }

    die "Error: No .anal files found in directory $loc\n"
	if ($#FILES == -1);

    return \@FILES;
}

sub GetEnergy(@) {
    my ($in_file) = $_[0];
    my (%ENG, $rec);
    
    open INFILE, $in_file or die "Cannot open $in_file: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^\s+(\d+)\s+(\S+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
	    $rec = (
		    {
			"SEQ_NAME" => $2,
			"SEQ_ENG"  => $3,
			"REF_ENG"  => $4,
			"DIFF"     => $5,
		    }
		    );
	    $ENG{$1} = $rec;
	}
    }


    return \%ENG;
}

sub CalcDifference(@) {
    my ($curr_eng, $ref_eng, $Stats) = @_;
    my ($counter, $returnStr, $eng_diff, $seq_name);
    my ($ref_total, $seq_total, $diff_total, $ref_name);

    $ref_total = $seq_total = $diff_total = 0;
    $returnStr = "";
    for $counter (sort Numerically keys %{ $curr_eng}) {
	next 
	    if (! defined ($ref_eng->{$counter}));
	$seq_name = $curr_eng->{$counter}{"SEQ_NAME"};
	$ref_name = $ref_eng->{$counter}{"SEQ_NAME"};
	$eng_diff = ($curr_eng->{$counter}{"SEQ_ENG"} - $curr_eng->{$counter}{"REF_ENG"}) - 
	    ($ref_eng->{$counter}{"SEQ_ENG"} - $ref_eng->{$counter}{"REF_ENG"});

	$Stats->{$counter}{"SEQ_NAME"} = $ref_eng->{$counter}{"SEQ_NAME"};
	if ($eng_diff < 0 and ($seq_name ne $ref_name)) {
	    $Stats->{$counter}{"NEGATIVE"}{"COUNTER"}++;
	    if (! defined($Stats->{$counter}{"NEGATIVE"}{$seq_name}) or
		$Stats->{$counter}{"NEGATIVE"}{$seq_name} < $eng_diff) {
		$Stats->{$counter}{"NEGATIVE"}{$seq_name} = $eng_diff;
	    }
	} else {
	    $Stats->{$counter}{"POSITIVE"}{"COUNTER"}++;
	    if (! defined($Stats->{$counter}{"POSITIVE"}{$seq_name}) or
		$Stats->{$counter}{"POSITIVE"}{$seq_name} < $eng_diff) {
		$Stats->{$counter}{"POSITIVE"}{$seq_name} = $eng_diff;
	    }
	}
	$ref_total += $ref_eng->{$counter}{"SEQ_ENG"};
	$seq_total += $curr_eng->{$counter}{"SEQ_ENG"};
	$diff_total += $eng_diff;
	$returnStr .= sprintf("%-8d%12s%12s%12.3f%12.3f%12.3f\n", $counter, 
			      $seq_name, $ref_eng->{$counter}{"SEQ_NAME"}, $curr_eng->{$counter}{"SEQ_ENG"}, 
			      $ref_eng->{$counter}{"SEQ_ENG"}, $eng_diff);
    }
    $returnStr .= sprintf("%-8s%12s%12s%12.3f%12.3f%12.3f\n", "#","","TOTAL", $seq_total, $ref_total, $diff_total);
    return ($Stats, $returnStr);
			  
			      
}

sub Numerically {
    ($a<=>$b);
}

sub PrintStats(@) {
    my ($stats) = $_[0];
    my ($counter, $pos, $neg, $subval, $index);
    my ($out_string);

    $out_string = sprintf("%-10s%10s\n", "OPTIMAL", "ORIGINAL");
    printf "%12s%12s%12s  %-40s\n", "SEQUENCE", "POSITIVE", "NEGATIVE", "NEG LIST";
    for $counter (sort Numerically keys %{ $stats }) {
	$pos = $neg = 0;
	$subval = "";
	if (defined $stats->{$counter}{"POSITIVE"}{"COUNTER"}) {
	    $pos = $stats->{$counter}{"POSITIVE"}{"COUNTER"};
	}
	if (defined $stats->{$counter}{"NEGATIVE"}{"COUNTER"}) {
	    $neg = $stats->{$counter}{"NEGATIVE"}{"COUNTER"};
	}

	for $index (keys %{ $stats->{$counter}{"NEGATIVE"} }) {
	    next
		if ($index eq "COUNTER");
	    $subval .= "$index:" . sprintf("%6.4f", $stats->{$counter}{"NEGATIVE"}{$index} ) . " ";
	}
	printf "%12s%12d%12d  %-40s\n", $stats->{$counter}{"SEQ_NAME"}, $pos, $neg, $subval;
	$out_string .= DetermineOptimalSequence($stats->{$counter}{"SEQ_NAME"}, $subval);
    }
    return $out_string;
}

sub DetermineOptimalSequence(@) {
    my ($curr_seq, $stats) = @_;
    my ($counter, $return_string, @Seq);

    if ($curr_seq =~ /(\d+)\'\/\w/) {
	@Seq = GetMax($curr_seq, $stats, 1);
	if ($1 == "5") {
	    $return_string = "#Helix\n";
	}
    } else {
	@Seq = GetMax($curr_seq, $stats, 0);
    }

    for $counter (@Seq) {
	$return_string .= $counter;
    }

    return $return_string;
}

sub GetMax(@) {
    my ($original_name, $optimal_list, $end_residue) = @_;
    my (@return_vals, %holder, $min_val, $min_name, @tmp, $h1, $index, @tmp2, $h2);

    $holder{"C"} = "CG";
    $holder{"G"} = "GC";
    $holder{"T"} = "TA";
    $holder{"A"} = "AT";

    $min_val = 0;
    $min_name = "";
    if ($optimal_list ne "") {
	@tmp = split /:|\/|\s+/, $optimal_list;
	$index = 0;
	while ($index < $#tmp) {
	    if ($tmp[$index + 2] < $min_val) {
		$min_val = $tmp[$index + 2];
		$min_name = $tmp[$index + 1];
	    }
	    $index += 3;
	}
    } else {
	@tmp = split /:|\/|\s+/, $original_name;
	if ($end_residue) {
	    $min_name = $tmp[1];
	} else {
	    $min_name = $tmp[0];
	}
    }

    @tmp = split /:|\/|\s+/, $original_name;
    if ($end_residue) {
	$original_name = $tmp[1];
    } else {
	$original_name = $tmp[0];
    }

    $h1 = substr($min_name, 0, 1);
    $h2 = substr($original_name, 0, 1);
    $return_vals[0] = sprintf("%-10s%10s\n", $holder{$h1}, $holder{$h2});
    if (! $end_residue) {
	$h1 = substr($min_name, 1, 1);
	$h2 = substr($original_name, 1, 1);
	$return_vals[1] = sprintf("%-10s%10s\n", $holder{$h1}, $holder{$h2});
	
    }

    return (@return_vals);
}

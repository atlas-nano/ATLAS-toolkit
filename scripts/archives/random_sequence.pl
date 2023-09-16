#!/usr/bin/perl -w
BEGIN {
    push @INC, "/ul/tpascal/scripts/";
}
use strict;
use File::Basename;
use Packages::General;
use Packages::GetParms;
                                                                                                                      
                                                                                                                      

# This script will take an input sequence, maintain the sequence at the crossovers and 
# create a random sequence subject to the %GC specified.

sub Initialize();
sub DetSequence(@);
sub GenerateSequence(@);
sub CreateFiles(@);

die "usage: $0 parameter_file template_1 template_2 percent_gc\n"
    if (! @ARGV or$#ARGV < 3);

my ($parm, $helix1, $helix2, $per_gc) = @ARGV;
my ($P_File, $h1_seq, $h2_seq, $Crossovers, $file_1, $file_2);

print "Initializing...";
Initialize();
print "Done\nDeterming Sequence....";
($h1_seq, $h2_seq, $Crossovers) = DetSequence($P_File);
print "Done\nGenerating Random Sequence...";
$file_1 = GenerateSequence($h1_seq, $Crossovers, $per_gc);
$file_2 = GenerateSequence($h2_seq, $Crossovers, $per_gc);
print "Done\nCreate Sequence Files...";
CreateFiles("TOPSEQUENCE.TXT", $file_1, "BOTTOMSEQUENCE.TXT", $file_2);
print "Done\n";
sub Initialize() {
    my ($counter);

    for $counter ($helix1, $helix2, $parm) {
	die "Error: Cannot access $counter: $!\n"
	    if (! -e $counter or ! -r $counter or ! -T $counter);
    }

    die "Error: Expected decimal for percent_gc, got $per_gc"
	if ($per_gc !~ /^\d+\.\d+/);

    $P_File = Packages::GetParms->new();
    die "Error in Paramater file\n"
	if (! $P_File->IsValidParams($parm));
    
                                                                                                                      
}

sub DetSequence(@) {
    my ($Parms) = $_[0];
    my (@Crossovers, @helix1, @helix2);
    @Crossovers = @{ $Parms->{"Molecule"}->{"crossovers"} };

    open INFILE, $helix1 or die "Cannot open helix sequence file $helix1: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^([a-zA-Z]{2})/) {
	    push @helix1, $1;
	}
    }
    close INFILE;

    open INFILE, $helix2 or die "Cannot open helix sequence file $helix2: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^([a-zA-Z]{2})/) {
	    push @helix2, $1;
	}
    }
    close INFILE;

    die "Error: $helix1 doesn't contain any relevant data\n"
	if ($#helix1 == -1);

    die "Error: $helix2 doesn't contain any relevant data\n"
	if ($#helix2 == -1);

    return (\@helix1, \@helix2, \@Crossovers);

}

sub GenerateSequence(@) {
    my ($seq1, $Cross, $percent) = @_;
    my (%CrossSpec, $counter, @tmp, $total_gc, $total_bps, @returnval, $rand_nm, $index);

    for $counter (@{ $Cross }) {
	$index = $counter - 1;
	$CrossSpec{"$index"} = 1;
    }

    $total_bps = ($#{ $seq1 } + 1) - (3 * ($#{ $Cross } + 1));
    $total_gc = $total_bps * $percent;

    if ($total_gc > int($total_gc)) {
	$total_gc++;
    }

    for $counter (0 .. ($total_gc -1)) {
	$tmp[$counter] = "GC";
    }

    for $counter ($total_gc .. ($total_bps - 1)) {
	$tmp[$counter] = "AT";
    }

    $total_bps += 3 * ($#{ $Cross } + 1);
    for $counter (0 .. ($total_bps  - 1)) {
	if (defined $CrossSpec{$counter + 1}) {
	    $returnval[$counter] = $seq1->[$counter];
	} elsif (defined $CrossSpec{$counter}) {
	    $returnval[$counter] = $seq1->[$counter];
	} elsif (defined $CrossSpec{$counter - 1}) {
	    $returnval[$counter] = $seq1->[$counter];
	} else {
	    $rand_nm = int(rand $#tmp);
	    $returnval[$counter] = $tmp[$rand_nm];
	    splice @tmp, $rand_nm, 1;
	}
    }

    return \@returnval;
	    
}

sub CreateFiles(@) {
    my (@file_data) = @_;
    my ($counter, $index, $file_nm, $data);

    $counter = 0;

    while ($counter < $#file_data) {
	$file_nm = $file_data[$counter];
	$data  = $file_data[$counter + 1];

	open OUTFILE, "> $file_nm" or die "Cannot create $file_nm: $!\n";
	print OUTFILE $P_File->{"Molecule"}->{"major_groove"} . ":" . $P_File->{"Molecule"}->{"minor_groove"} . "\n";
	for $index (@{ $data }) {
	    print OUTFILE "$index\n";
	}
	close OUTFILE;
	$counter += 2;
    }
}

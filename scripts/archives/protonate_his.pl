#!/usr/bin/perl -w
use strict;

sub GetResInfo();
sub GetCeriusFile(@);
sub ExecuteCmd();

if (!@ARGV) {
    die "usage: $0 bgffile start_res\n";
}

my ($infile, $start_res) = @ARGV;
my ($template) = "/ul/tpascal/scripts/protonate_his.template";

my (@cerius_array, $i, $tmp_file, @created_files);
-e $template or die "Cannot locate $template: $!\n";

$created_files[0] = $infile;

while ($#created_files > -1) {
    $infile = $created_files[0];
    pop @created_files;

    @cerius_array = ();
    -e $infile or die "Cannot locate $infile: $!\n";
    
    GetResInfo();
    
    for $i (0 .. $#cerius_array) {
	$tmp_file = "tmp_file." . $i;
	
	open TMPFILE, "> $tmp_file" or die "Cannot write to $tmp_file: $!\n";
	print TMPFILE "$cerius_array[$i]\n";
	close TMPFILE;
	ExecuteCmd();
	system "rm -f $tmp_file";
    }
}

sub ExecuteCmd() {
    open CERIUS2, "cerius2 -n tmp_file |" or die "Cannot execute cerius2: $!\n";
    while (<CERIUS2>) {
	chomp;
    }
    close CERIUS2;
    
}

sub GetResInfo() {
    my ($tmp_info, $inline, @his_array, $current_residue);
    my ($nd1, $nd1h, $end_atm) = (0, 0, 0);

    open INFILE, $template or die "Cannot locate $template: $!\n";
    while (<INFILE>) {
	chomp;
	$tmp_info .= $_ . "\n";
    }
    close INFILE;

    $tmp_info =~ s/his_file/$infile/;

    $current_residue = 0;
    open INFILE, $infile or die "Cannot open $infile: $!\n";
    while (<INFILE>) {
	if ($_ =~ /^ATOM\s+(\d+)\s+(\w+)\s+\HIS\s+(\d+)/) {
	    if ($2 eq "ND1") {
		print "Got HERE: $3 $current_residue $start_res\n";
		if ($current_residue < $3 and $3 > $start_res) {
		    if ($3 > 0 && $nd1h) {
			push @his_array, "$nd1 $nd1h $current_residue";
		    }
		    $nd1h = 0;
		    $current_residue = $3;
		    $nd1 = $1;
		}
	    } elsif ($2 eq "HND1") {
		$nd1h = $1;
	    } else {
		$end_atm = $1;
	    }
	} elsif ($_ =~ /^ATOM\s+(\d+)\s+\w+\s+\w+\s+\d+/) {
	    $end_atm = $1;
	}
    }
    close INFILE;
    print "HERE: $nd1 $nd1h $end_atm\n";

    if ($current_residue && $end_atm && $nd1h && $nd1) {
	push @his_array, "$nd1 $nd1h $current_residue";
    }

    if ($#his_array > -1) {
	for $i (0 .. $#his_array) {
	    $inline = GetCeriusFile((split / /, $his_array[$i]), $end_atm, $tmp_info);
	    push @cerius_array, $inline;
	}
    } else {
	die "Error: Unable too find any HIS information in $infile\n";
    }
    
}

sub GetCeriusFile(@) {
    my ($returnstr, $pattern, $i, $test_string);
    my ($hsd_fname, $hsp_fname);
    my ($nd1, $nd1h, $curr_res, $endatm, $template_info) =  @_;

    $endatm++;

    $template_info =~ s/his_nd1_atm/ATOM($nd1)/;
    $template_info =~ s/new_hydrogen/ATOM($endatm)/;
    $template_info =~ s/old_hydrogen/ATOM($nd1h)/;

# change filename
    $pattern = $infile;
    $pattern =~ s/\.bgf//;
    my (@fdata) = split /_/, $pattern;

    if ($#fdata > -1) {
	for $i (0 .. ($curr_res -2)) {
	    $test_string .= $fdata[$i] . "_";
	}
	$hsd_fname = $test_string . "HSD_";
	$hsp_fname = $test_string . "HSP_";
	for $i ($curr_res .. $#fdata) {
	    $hsd_fname .= $fdata[$i] . "_";
	    $hsp_fname .= $fdata[$i] . "_";
	}
	chop $hsd_fname;
	chop $hsp_fname;
	$hsd_fname .= ".bgf";
	$hsp_fname .= ".bgf";
    } else {
	$hsd_fname = $pattern . "_HSD.bgf";
	$hsp_fname = $pattern . "_HSP.bgf";
    }

    $template_info =~ s/hsd_file/$hsd_fname/;
    $template_info =~ s/hsp_file/$hsp_fname/;

    push @created_files, $hsd_fname;
    push @created_files, $hsp_fname;

    print "$hsd_fname $hsp_fname $#fdata\n";
    return $template_info;
}

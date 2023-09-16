#!/usr/bin/perl -w
#    This script will input a biograf format file and output a pdb file with
#    Connectivities in rasmol
use strict;
use File::Basename;

sub CheckInput();
sub GetFiles(@);
sub ParseBGFFile(@);
sub CreatePDB(@);
sub SavePDB(@);
sub Numerically;
sub PrintRasmol(@);

die "usage: $0 bgf_file|directory rasmol_cmd [pdb_file]\n"
    if (! @ARGV or $#ARGV < 1);

my ($loc, $rasmol_cmd, $save_name) = @ARGV;
my ($counter, $F_Data, $out_string, $CON, $SPACE, $FILES);
CheckInput();

$FILES = GetFiles($loc);

print "Converting files...";
for $counter (0 .. $#{ $FILES }) {
    print ".";
    ($F_Data, $CON) = ParseBGFFile($FILES->[$counter]{"FILENAME"});
     next
        if (! defined($F_Data) or ! defined($CON));
    ($out_string, $SPACE) = CreatePDB($F_Data, $CON);
    SavePDB($out_string, $FILES->[$counter]{"SAVENAME"});
    PrintRasmol($SPACE, $FILES->[$counter]{"SAVENAME"}, $FILES->[$counter]{"FILENAME"});
    system ("$rasmol_cmd < rasmol_cmd >& junk");
}

print "Done\n";


sub CheckInput() {
    die "Error: Cannot access $loc: $!\n"
	if (! -e $loc and ! -d $loc);
    system("mkdir -p pics");
    system("mkdir -p pdbs");

    die "Cannot find rasmol at $rasmol_cmd\n"
	if (! -e $rasmol_cmd);
}

sub GetFiles(@) {
    my ($curr_loc) = $_[0];
    my (@FILES, @holder, $savenm, $counter, $rec);

    if (! -d $curr_loc) {
	if (! $save_name) {
	    $savenm = basename($curr_loc);
	    $savenm =~ s/\.bgf/\.rasmol/g;
	} else {
	    $savenm = $save_name;
	}
	$rec = (
		{
		    "FILENAME" => $curr_loc,
		    "SAVENAME" => $savenm,
		}
		);
	$FILES[0] = $rec;
    } else {
	opendir CURRDIR, $curr_loc or die "Cannot access directory $curr_loc: $!\n";
	@holder = grep { /^\S+\.bgf$/ } map { "$curr_loc/$_"} readdir CURRDIR;
	closedir CURRDIR;

	for $counter (@holder) {
	    $savenm = basename($counter);
	    $savenm =~ s/\.bgf/\.rasmol/g;
	    $rec = (
		    {
			"FILENAME" => $counter,
			"SAVENAME" => $savenm,
		    }
		    );
	    push @FILES, $rec;
	}
    }

    die "Error: No bgf files found in $curr_loc\n"
	if ($#FILES == -1);

    return \@FILES;
}

sub ParseBGFFile(@) {
    my ($in_file) = @_;
    my (%ATOMS, @CON, $counter, $in_data, $atm_patern, @holder, $constr);
    
    $atm_patern = '^ATOM\s+(\d+)\s+(\w+)\s+(\w+)\s\w?\s+(\d+)\s+' .
        '(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\w+)\s+' .
        '(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s+(\d)+\s+(\d)+\s+(\d+\.\d+)';
    open BGFFILE, $in_file or die "Cannot open BGFfile $in_file: $!\n";
    while (<BGFFILE>) {
        chomp;
        $in_data = $_;
        $in_data =~ s/HETATM/ATOM/;
        if ($in_data =~ /$atm_patern/) {
            $ATOMS{$1} = (
                            {
                                "ATMNAME"     => $2,
                                "RESNAME"     => $3,
                                "RESNUM"      => $4,
                                "XCOORD"      => $5,
                                "YCOORD"      => $6,
                                "ZCOORD"      => $7,
                                "REALATM"     => $8,
                                "NUMBONDS"    => $9,
                                "LONEPAIRS"   => $10,
                                "CHARGE"      => $11,
                                "OCCUPANCY"   => $12,
                                "RESONANCE"   => $13,
                                "RADII"       => $14,
                            }
                            );
        } elsif ($in_data =~ /^CONECT\s+(.+)$/) {
	    @holder = split /\s+/, $1;
	    $constr = "";
	    for $counter (@holder) {
		$constr .= sprintf("%5d", $counter);
	    }
	    push @CON, "CONECT" . $constr . "\n";
	}
    }
    close BGFFILE;
    
    die "ERROR: $in_file does not contain any ATOM/CONNECTION information!\n"
        if ($#CON == -1 or ! %ATOMS == -1);
    
    return (\%ATOMS, \@CON);
}

sub CreatePDB(@) {
    my ($AtmInfo, $Conn) = @_;
    my ($counter, $outstring, %SPACE, $atm_letter);
    
    $outstring = "";
    for $counter (sort Numerically keys %{ $AtmInfo }) {
	$atm_letter = substr($AtmInfo->{$counter}{"ATMNAME"}, 0, 1);
	$outstring .= sprintf("%-5s%6d%2s%7s%2s%4s%12.3f%8.3f%8.3f\n", "ATOM", $counter,
			      $atm_letter, $AtmInfo->{$counter}{"RESNAME"},
			      "A", $AtmInfo->{$counter}{"RESNUM"}, $AtmInfo->{$counter}{"XCOORD"},
			      $AtmInfo->{$counter}{"YCOORD"}, $AtmInfo->{$counter}{"ZCOORD"});
	$SPACE{$atm_letter} = int(100 *  $AtmInfo->{$counter}{"RADII"});
    }
    $outstring .= "TER\n";
    
    for $counter (@{ $Conn }) {
	$outstring .= $counter;
    }
    
    
    return ($outstring, \%SPACE);
}

sub SavePDB(@) {
    my ($outtext, $save_name) = @_;

    open OUTFILE, "> pdbs/$save_name" or die "Cannot write to pdbs/$save_name: $!\n";
    print OUTFILE $outtext;
    close OUTFILE;
}

sub Numerically {
    ($a<=>$b);
}

sub PrintRasmol(@) {
    my ($cpk, $save_name, $fle_nm) = @_;
    my ($outstring, $counter);

    $outstring = "load pdb pdbs/$save_name\n";
    $outstring .= "rotate x 270\n";
    $outstring .= "rotate y 90\n";
    $outstring .= "background white\n";
	
    if ($fle_nm =~ /(\d+)\.bgf$/) {
	$counter = $1;
    } else {
	$counter = 0;
    }

    for $counter (keys %{ $cpk }) {
	$outstring .= "select *." . $counter . "\n";
	$outstring .= "cpk " . $cpk->{$counter} . "\n";
    }

    $outstring .= "write pics/pic_" . $counter . ".gif\n";
    $outstring .= "quit\n";
    open OUTFILE, "> rasmol_cmd" or die "Cannot write to rasmol_cmd: $!\n";
    print OUTFILE $outstring;
    close OUTFILE;
}

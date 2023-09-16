#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use Packages::General qw(FileTester Trim LoadElements);
use Packages::FileFormats qw(GetBGFFileInfo);
use Packages::BOX qw(GetBox);
use Packages::CERIUS2 qw(parseCerius2FF);

use strict;

# This program will open a bgf file and will write an msi file
sub addACL;
sub findElement;
sub init;
sub CreateMSI;

die "usage: $0 bgf_file cerius2FF [save_name]\n"
    if (! @ARGV || $#ARGV < 1);

my ($bgf_file, $ceriusFF, $save_name) = @ARGV;
my ($ATOMS, $BBOX, $MSI, $BONDS, $HEADERS, $PARMS, $ELEMENTS);

$|++;
print "Initializing...";
&init;
print "Parsing CERIUS2 forcefield $ceriusFF...";
$PARMS = parseCerius2FF($ceriusFF);
print "Done\nParsing BGF file $bgf_file...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgf_file, 1);
$BBOX = GetBox($ATOMS, $PARMS, $HEADERS);
addACL($ATOMS, $PARMS, $ELEMENTS);
print "Done\nWriting MSI file $save_name...";
CreateMSI($ATOMS, $BONDS, $HEADERS, $BBOX, $save_name);
print "Done\n";

sub addACL {
    my ($atoms, $parms, $elements) = @_;
    my ($i, $j, $fftype, $elementName, $elementNum);

    for $i (keys %{ $atoms }) {
	$fftype = $atoms->{$i}{FFTYPE};
	($elementName, $elementNum) = findElement($fftype, $parms, $elements);
	if (length($elementName) > 1) {
	    $elementName = uc(substr($elementName, 0, 1)) . lc(substr($elementName, 1, length($elementName)));
	}
	$atoms->{$i}{ACL} = "$elementNum $elementName";
	$atoms->{$i}{LABEL} = $atoms->{$i}{ATMNAME};
    }
}

sub findElement {
    my ($fftype, $parms, $elements) = @_;
    my ($elementName, $elementNum, $i);

    die "ERROR: $fftype is not a valid force field type!\n" if (! exists($parms->{ATOMTYPES}{$fftype}));
    $elementName = $parms->{ATOMTYPES}{$fftype}{ATOM};
    $elementNum = 0;

    for $i (keys %{ $elements }) {
	if (uc($elements->{$i}{SYMBOL}) eq uc($elementName)) {
	    $elementNum = $i;
	    last;
	}
    }

    die "ERROR: Element $elementName is not a valid element!\n" if (! $elementNum);

    return ($elementName, $elementNum);
}

sub init {
    FileTester($bgf_file);
    FileTester($ceriusFF);
    if (! $save_name) {
	$save_name = $bgf_file;
	$save_name =~ s/\.bgf//g;
	$save_name .= ".msi";
    }
    $ELEMENTS = &LoadElements;
}


sub ReadBGF(@) {
    my ($in_file) = $_[0];
    my (%MSI_File, $rec, $index, $inStr, $pTableInfo, $fftype, $BOX, $atm_patern);

    $atm_patern = '^ATOM\s+\d+\s(\s*\w+)\s+\w+\s+\w?\s+\d+\s+(\-?\d+\.\d+)\s+(\-?\d+\.' .
	'\d+)\s+(\-?\d+\.\d+)\s+(\S+)\s+\d+\s+\d+\s+(\-?\d+\.\d+)';
    open BGFFILE, $in_file or die "Cannot open BGF file $in_file: $!\n";
    while (<BGFFILE>) {
	chomp;
	$inStr = $_;
	$inStr =~ s/HETATM/ATOM/g;
	$rec = ();

	if ($inStr =~ /DESCRP\s+(\S+)/) {
	    $MSI_File{"Label"} = $1;
	} elsif ($inStr =~ /REMARK\s+(.+)$/) {
	    push @{ $MSI_File{"Comment"} }, $1;
	} elsif ($inStr =~ /ATOM\s+(\d+)/) {
	    $index = $1;
	    if ($inStr =~ /$atm_patern/) {
                $pTableInfo = "";
		if ($5 eq "Po") {
		    $pTableInfo = "84 Po";
		} elsif ($5 eq "Cs") {
		    $pTableInfo = "55 Cs";
		} elsif ($5 eq "O" or $5 eq "N" or $5 eq "S") {
		    $pTableInfo = "1 H";
		} elsif ($5 eq "Ra") {
		    $pTableInfo = "88 Ra";
		} elsif ($5 eq "Ac") {
		    $pTableInfo = "89 Ac";
		} elsif ($5 eq "La") {
		    $pTableInfo = "57 La";
		} elsif ($5 eq "Sr") {
		    $pTableInfo = "38 Sr";
		}
		$fftype = Trim($5);
		$rec = (
			{
			    "Label"   => $1,
			    "XYZ"     => "$2 $3 $4",
			    "FFType"  => $fftype,
			    "Charge"  => $6,
			    "ACL"     => $pTableInfo,
			    "Id"      => $index,
			}
			);
		push @{ $MSI_File{"Atoms"} }, $rec;
		$BOX = StoreExtrema($2, $3, $4, $6, $BOX);
#		if (! defined($MSI_File{"XMAX"}) or ($2 > $MSI_File{"XMAX"})) {
#		    $MSI_File{"XMAX"} = $2;
#		} 
#		if (! defined($MSI_File{"YMAX"}) or ($3 > $MSI_File{"YMAX"})) {
#		    $MSI_File{"YMAX"} = $3;
#		} 
#		if (! defined($MSI_File{"ZMAX"}) or ($4 > $MSI_File{"ZMAX"})) {
#		    $MSI_File{"ZMAX"} = $4;
#		} 
	    }
	} elsif ($inStr =~ /CONECT\s+(\d+)\s+(.+)$/) {
	    @{ $rec } = split /\s+/, $2;
	    $MSI_File{"Bonds"}{$1} = $rec;
	}
    }
    close BGFFILE;
    
    die "ERROR: Invalid BGF file $in_file\n"
	if (! %MSI_File);

    return (\%MSI_File, $BOX);
}

sub CreateMSI(@) {
    my ($ATOMS, $BONDS, $HEADERS, $BOX, $savename) = @_;
    my ($out_string, $counter, $index, $h_keys, $data_type);

    $out_string = "# MSI CERIUS2 DataModel File Version 4 0\n";
    $out_string .= "(1 Model\n (A I Id 1)\n";
    $out_string .= " (A I PeriodicType 100)\n";
    $out_string .= " (A D A3 (" . ($BOX->{"X"}{"hi"} - $BOX->{"X"}{"lo"}) . " 0 0))\n";
    $out_string .= " (A D B3 (0 " . ($BOX->{"Y"}{"hi"} - $BOX->{"Y"}{"lo"}) .  " 0))\n";
    $out_string .= " (A D C3 (0 0 " . ($BOX->{"Z"}{"hi"} - $BOX->{"Z"}{"lo"}) . "))\n";
    $out_string .= ' (A C SpaceGroup "1 1")' . "\n";
    
    $counter = 1;
    for $index (@{ $HEADERS }) {
	if ($index =~ /^DESC/) {
	    $index =~ s/^DESC//;
	    $out_string .= ' (A C Label "' . $index . '" )' . "\n";
	} elsif ($index =~ /^REMARK/) {
	    $index =~ s/^REMARK//;
	    $out_string .= ' (A C "Save-comment ' . $counter . '" "' . $index . '")' . "\n";
	    $counter++;
	}
    }

    $counter = 2;
    for $index (sort Numerically keys %{ $ATOMS }) {
	$out_string .= " ($counter Atom\n";
	for $h_keys ("ACL", "Charge", "FFType", "XYZ", "Id", "Label") {
	    if ($h_keys eq "XYZ") {
		$data_type = "D";
		$out_string .= "  (A $data_type $h_keys ($ATOMS->{$index}{XCOORD} $ATOMS->{$index}{YCOORD} $ATOMS->{$index}{ZCOORD}))\n";
	    } elsif ($h_keys eq "Id") {
		$data_type = "I";
		$out_string .= "  (A $data_type $h_keys $index)\n";
	    } elsif ($h_keys eq "Charge") {
		$data_type = "I";
		$out_string .= "  (A $data_type $h_keys $ATOMS->{$index}{CHARGE})\n";
	    } else {
		$data_type = "C";
		$out_string .= "  (A $data_type $h_keys " . '"' . $ATOMS->{$index}{uc($h_keys)} . '")' . "\n";
	    }
		   
	}
	$out_string .= " )\n";
	$counter++;
    }
    
    for $index (sort Numerically keys %{ $BONDS }) {
	for $h_keys (@{ $BONDS->{$index} }) {
	    if ($index < $h_keys) {
		$out_string .= " ($counter Bond\n";
		$out_string .= "  (A O Atom1 " . ($index + 1). " )\n";
		$out_string .= "  (A O Atom2 " . ($h_keys + 1) . " )\n";
		$out_string .= " )\n";
		$counter++;
	    }
	}
    }
    $out_string .= ")\n";

    open OUTFILE, "> $savename" or die "Cannot write to $savename: $!\n";
    print OUTFILE $out_string;
    close OUTFILE;
}

sub Numerically {
    ($a<=>$b);
}

sub STrim(@) {
    my ($inText) = $_[0];

    if ($inText =~ /(\w+)/) {
	$inText = $1;
    }

    return $inText;
}

sub StoreExtrema(@) {
    my ($x, $y, $z, $radii, $par) = @_;
    my ($counter);

    if (! defined($par)) {
	for $counter ("X", "Y", "Z") {
	    $par->{$counter} = (
				{
				    "lo" => 99999.999,
				    "hi" => -99999.999,
				}
				);
	}
    }
	
    $par->{"X"}{"hi"} = ($x + $radii)
	if (($x + $radii) > $par->{"X"}{"hi"});
    $par->{"X"}{"lo"} = ($x - $radii)
	if (($x - $radii) < $par->{"X"}{"lo"});

    $par->{"Y"}{"hi"} = ($y + $radii)
	if (($y + $radii) > $par->{"Y"}{"hi"});
    $par->{"Y"}{"lo"} = ($y - $radii)
	if (($y - $radii) < $par->{"Y"}{"lo"});

    $par->{"Z"}{"hi"} = ($z + $radii)
	if (($z + $radii) > $par->{"Z"}{"hi"});
    $par->{"Z"}{"lo"} = ($z - $radii)
	if (($z - $radii) < $par->{"Z"}{"lo"});

    return $par;
}

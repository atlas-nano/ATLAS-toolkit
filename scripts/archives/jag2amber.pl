#!/usr/bin/perl -w
# This program will take a Jaguar output file, read the final geometry
# and charges and create an compatible AMBER pdb type file
BEGIN {
    push (@INC, "/ul/tpascal/scripts");
}
use strict;
use Packages::FileFormats;
use File::Basename;
use Packages::General;

sub ReadJagOut;
sub ReadConvFile;
sub WriteAmberInput;
sub Initialize;
sub Numerically;
sub GetResCounter;
sub FindCnvHashKey;
sub ReadTemplateFile;

die "usage: $0 jaguar_output conversion_file #atoms [res_name] [save_name]\n"
    if (! @ARGV or $#ARGV < 2);

my ($jagOut, $convFile, $resNum, $resName, $saveName) = @ARGV;
my ($JagData, $ConvData);

Initialize();
print "Parsing Jaguar Output file $jagOut...";
$JagData = ReadJagOut($jagOut);
print "Done\nParsing Conversion file $convFile...";
$ConvData = ReadConvFile($convFile);

print "Done\nCreating $saveName" . ".pdb and myleaprc files...";
WriteAmberInput($saveName);
print "Done\n\nDo the following steps: \n1. Copy the file /ul/tpascal/amber8/leaprc to this directory\n";
print "2. Append the myleaprc file just created to the bottom of leaprc file\n";
print "3. Run xleap\n";


sub Initialize {
    $resName = "UNK"
	if (! defined($resName));

    $saveName = $jagOut
	if (! defined($saveName));

    if ($saveName =~ /(.+)\.\w{3}$/) {
	$saveName = $1;
    }

    die "ERROR: Expected Integer for atoms_in_residue. Got $resNum\n"
	if (! IsInteger($resNum));

    $resNum++;
    FileTester($jagOut);
    FileTester($convFile);
}

sub ReadJagOut(@) {
    my ($jag_out) = $_[0];
    my ($inline, $is_radii, $hash_key, $is_geometry, $is_charge);
    my (@atom_id, @atom_charge, $counter, $is_valid, $curr_id);
    my (%Jaguar_Data, $resCounter);

    $is_radii = $is_geometry = $is_charge =  $is_valid = 0;
    $resCounter = 1;
    open JAGOUT, $jag_out or die "Cannot read from $jag_out: $!\n";
    while (<JAGOUT>) {
	chomp;
	$inline = $_;
	if ($inline =~ /final geometry/) {
	    $is_geometry = 1;
	} elsif ($inline =~ /bond lengths \(angstroms\)/) {
	    $is_geometry = 0;
	}elsif ($is_geometry) {
	    if ($inline =~ /([A-Z]{1}[a-z]?)(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
		$hash_key = $1 . $2;
		$resCounter = int($2/$resNum) + 1;

		$Jaguar_Data{$resCounter}{$hash_key}{"XCOORD"} = $3;
		$Jaguar_Data{$resCounter}{$hash_key}{"YCOORD"} = $4;
		$Jaguar_Data{$resCounter}{$hash_key}{"ZCOORD"} = $5;
		
		$is_valid = 1;
	    }
	} elsif ($inline =~ /Atomic charges from/) {
	    $is_charge = 1;
	} elsif ($inline =~ /sum of atomic charges/) {
	    $is_charge = 0;
	}elsif ($is_charge && $is_valid) {
	    if ($inline =~ /Atom\s+(.+)$/) {
		@atom_id = split /\s+/, $1;
	    }elsif ($inline =~ /Charge\s+(.+)$/) {
		@atom_charge = split /\s+/, $1;
		for $counter (0 .. $#atom_charge) {
		    $atom_id[$counter] =~ /(\d+)/;
		    $curr_id = $1;
		    $resCounter = GetResCounter(\%Jaguar_Data, $atom_id[$counter]);

		    if ($Jaguar_Data{$resCounter}{$atom_id[$counter]}) {
#			print "Wrote Charge\n";
			$Jaguar_Data{$resCounter}{$atom_id[$counter]}{"CHARGE"} = $atom_charge[$counter];
			$is_valid = 1;
		    } else {
			print "Charge found for Non Existant Atom: $atom_id[$counter]\n";
		    }
		}
	    }
	}
    }
    close JAGOUT;

    die "Error reading Jaguar Output $jag_out\n"
	if (! $is_valid);

    return (\%Jaguar_Data);
}

sub ReadConvFile {
    my ($cnvFile) = $_[0];
    my ($isStart, %CNV);

    open CNVFILE, "$cnvFile" or die "Error reading regular file $cnvFile: $!\n";

    $isStart = 0;

    while (<CNVFILE>) {
	chomp;
	if ($_ =~ /^ATOMNAME	ATOMTYPE	ELEMENT/) {
	    $isStart = 1;
	} elsif ($isStart and $_ =~ /^(\w+)\s+(\S+)\s+(\w)/) {
	    $CNV{$1} = (
			{
			    "type"    => $2,
			    "element" => $3,
			}
			);
	}
    }
    close CNVFILE;

    die "ERROR: $cnvFile is not valid\n"
	if (! $isStart);

    return \%CNV;

}

sub WriteAmberInput {
    my ($savePrefix) = $_[0];
    my ($counter, $hashKey, $pdbfile, $resIndex, $cnvHash);
    my ($leapData, $pdbData);

    $counter = 0;
    $pdbfile = $savePrefix . ".pdb";

    
    $leapData = "model = loadpdb " . $pdbfile . "\n";
    
    for $resIndex (sort Numerically keys %{ $JagData }) {
    #for $resIndex (1 .. 1) {
	for $hashKey (keys %{ $JagData->{$resIndex} }) {
	    $counter++;

	    if (! $ConvData->{$hashKey}) {
		$cnvHash = FindCnvHashKey($resIndex, $hashKey);
	    } else {
		$cnvHash = $hashKey;
	    }

	    die "ERROR: Didn't find atom name in conversion file: $hashKey\n"
		if (! $ConvData->{$cnvHash});
	
	    $pdbData .= sprintf("%-5s%6d%2s%-4s%3s%6d%12.3f%8.3f%8.3f\n",
				"ATOM", $counter, " ", $hashKey, $resName, 1, $JagData->{$resIndex}{$hashKey}{"XCOORD"}, 
				$JagData->{$resIndex}{$hashKey}{"YCOORD"}, $JagData->{$resIndex}{$hashKey}{"ZCOORD"});	
	
	    $leapData .= "set model." . $resName . "." . $hashKey . " type \"" . $ConvData->{$cnvHash}{"type"} . "\"\n";
	    $leapData .= "set model." . $resName . "." . $hashKey . " charge " . $JagData->{$resIndex}{$hashKey}{"CHARGE"} . "\n";
	    $leapData .= "set model." . $resName . "." . $hashKey . " element " . $ConvData->{$cnvHash}{"element"} . "\n";
	    
	}
	$pdbData .= "TER\n";

    }
    
    $leapData .= "bondbydistance model\n";


    open OUTDATA, "> $pdbfile" or die "Cannot write to $pdbfile: $!\n";
    open LEAPDATA, "> myleaprc" or die "Cannot write to myleaprc: $!\n";
    print OUTDATA $pdbData;
    print LEAPDATA $leapData;
    close LEAPDATA;
    close OUTDATA;
}

sub FindCnvHashKey {
    my ($currRes, $currKey) = @_;
    my ($offset, $currNum);
    
    if ($currKey =~ /(\d+)/) {
	$currNum = $1;
	$offset = $currNum - (($currRes - 1) * ($resNum -1));
	$currKey =~ s/$currNum/$offset/;
    }

    return $currKey;
    
    
}

sub Numerically {
    ($a<=>$b);

}

sub GetResCounter {
    my ($myHash, $searchKey) = @_;
    my ($returnVal, $counter, $hashKey);

    $returnVal = -1;
    
  MAINLOOP: for $counter (keys %{ $myHash }) {
      for $hashKey (keys %{ $myHash->{$counter}} ) {
	  if ($hashKey eq $searchKey) {
	      $returnVal = $counter;
	      last MAINLOOP;
	  }
      }
  }

    return $returnVal;

}

sub ReadTemplateFile {
    my ($pdbFile) = $_[0];
    my ($hashKey, $isValid);
    my ($ATOMS) = GetPDBFileInfo($pdbFile, 0, " ", " ", 1);
    my ($jagName, $realName, $element, %CNV);

    $isValid = 0;
    for $hashKey (sort Numerically keys %{ $ATOMS }) {
	$realName = $ATOMS->{$hashKey}{"LABEL"};
	if ($realName =~ /([A-Z]{1}[a-z]?)/) {
	    $isValid = 1;
	    $element = $1;
	    $jagName = $1 . $hashKey;
	    $CNV{$jagName} = (
			      {
				  "type"    => $realName,
				  "element" => $element,
			      }
			      );
	    print "$jagName\t\t$realName\n";
	}
    }

    die "ERROR: $pdbFile contains invalid data\n"
	if (! $isValid);

    return \%CNV;
    
}

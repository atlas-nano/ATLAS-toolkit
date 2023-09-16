#!/usr/bin/perl -w

# Displays the list of directories and the total size
# usage dirsize.pl [whichdir]

use strict;
sub fileSize($);

my $curr_dir = ".";

if ($ARGV[0]) {
    -d $ARGV[0] or die "Cannot find $ARGV[0]: $!\n";
    $curr_dir = $ARGV[0];
}

my ($counter, $mytotal, $largest_file);
my ($instr, $lfile_path, $lfile_name);

my $dirlist_cmd = "ls -alR " . $curr_dir . "| awk '{printf\"%12d %s\\n\", \$5, \$9}' |";
open DIRLIST, "$dirlist_cmd" or die "Cannot get directory listing\n";

$largest_file = 0.0;
$lfile_name = "";
while (<DIRLIST>) {
    $instr = $_;
    chomp($instr);
    if ($instr =~ /(\d+)\s(\w+.*)$/) {
	if ($1 > 0) {
	    $mytotal += $1;
	    $counter +=1;
	    if ($1>$largest_file) {
		$largest_file = $1;
		$lfile_name =  $2;
#		print "$instr\n";
	    }
	}
    }
}
close DIRLIST;

if ($lfile_name) {
my $myflepath = "find -name $lfile_name";
open FLEPATH, "$myflepath |" or die "Cannot get file location: $!\n";
while (<FLEPATH>) {
    chomp;
    $lfile_path = $_;    
}

close FLEPATH;

print "Largest file: $lfile_name ($lfile_path): ";
$largest_file = fileSize($largest_file);
print "\n";

print "Total: ";
$mytotal = fileSize($mytotal);
print "in $counter files\n";
}
sub fileSize($) {
    my $flesz = $_[0];

    if ($flesz/1024 < 1) {
	printf "%f Bytes ", $flesz;
    } else {
	$flesz = $flesz/1024;
	if ($flesz/1024 < 1) {
	    printf "%f KBytes ", $flesz;
	} else {
	    $flesz = $flesz/1024;
	    if ($flesz/1024 < 1) {
		printf "%f MBytes ", $flesz;
	    } else {
		$flesz = $flesz/1024;
		if ($flesz/1024 < 1) {
		    printf "%f GBytes ", $flesz;
		} else {
		    $flesz = $flesz/1024;
		    printf "%f TBytes ", $flesz;
		}
	    }
	}
    }
}






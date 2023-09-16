#!/usr/bin/perl -w

use File::Basename;
use strict;

sub GetCurrFile($);

my ($in_template, $start, $end, $save_name) = @ARGV;
my ($old_name, $new_name);

if (! @ARGV || $#ARGV < 2) {
	die "usage: $0 in_template start end\n";
}

if (! defined($save_name)) {
    $save_name = $in_template;
}

system "mkdir -p newfiles";

for ($start .. $end) {
	($old_name, $new_name) = GetCurrFile($_);
	if ($new_name) {
		print "Renaming $old_name -> $new_name ...";
		system("cp $old_name $new_name") ?
			print "Failed!\n" :
			print "Sucess!\n";
	}
}

sub GetCurrFile($) {
	my ($temp) = $in_template;
	$temp =~ s/\.lis$//g;
	$temp =~ s/_\d+$//g;

	my ($curr_file) = $temp . "_" . $_[0] . ".lis";

	if (-e $curr_file) {
		my ($response) = "./newfiles/" . $save_name . "_lis" . "." . $_[0];
		return ($curr_file, $response);
	}
}


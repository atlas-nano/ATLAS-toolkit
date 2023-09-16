#!/usr/bin/perl -w
use warnings;
use strict;


BEGIN {
    push (@INC, "/ul/tpascal/.libs/blib/arch/auto/p5namot");
    push (@INC, "/ul/tpascal/.libs/PerlMagick/arch/auto/Image/Magick/");
    push (@INC, "/ul/tpascal/.libs/Net/Ssh/lib/perl5/site_perl/5.6.1");
}

use p5namot;
use Magick;
use Net::SSH qw(ssh);

p5namot::Cmd("set hush ERROR off");
p5namot::Cmd("set hush INFO off");
p5namot::Cmd("set hush REQUESTED off");
p5namot::Cmd("set hush WARNING off");
p5namot::Cmd("set WCad on");
p5namot::Cmd("set WCadFirst on");

if (!@ARGV or $#ARGV < 1) {
    die "usage: $0 scriptfile moviename [hostname]\n";
}

my ($scriptfile, $movienm, $host_name) = @ARGV;
my (@s_array, $array_ele, $filebase, $currpic, $pic_files, $last_frame);
my ($counter, $totalpics, $movie_cmd, @curr_cmd, $user_name);

sub executeSSHCmd(@);
sub ConvertPic($);

{
    
    -e $scriptfile or die "Cannot locate $scriptfile: $!\n";
    if (! $host_name) {
	$host_name = "octane2";
    }
    print "INFO: Will create movie on $host_name\n";

    print "\n1. Reading $scriptfile...";
    open SCRIPTFILE, $scriptfile or die "Cannot open $scriptfile: $!\n";
    
    while (<SCRIPTFILE>) {
	chomp;
	if ($_ =~ /render/ or $_ =~ /\#/) { } else {
	    push @s_array, $_;
	}
    }
    
    close SCRIPTFILE;
   
    $totalpics = $#s_array +1;

    print "Done\n";
    if ($totalpics < 5) {
	die "ERROR: $scriptfile is invalid\n";
    }

    $movienm =~ /(\w+)/;
    $filebase = $1;
    
    $counter = 0;
    
    if (! -d "picfiles") {
	mkdir "picfiles";
    }
    $user_name = $ENV{LOGNAME};
    
    chdir "picfiles";

    print "\n2. Converting $totalpics picture files: png -> gif...";
    for $array_ele (@s_array) {
	$currpic = sprintf("%0" . length($totalpics) . "d", $counter);
	$currpic = $filebase . "_" . $currpic; 
	p5namot::Cmd("set background black");
	p5namot::Cmd("render");
	p5namot::Cmd($array_ele);
	p5namot::Cmd("render");
	p5namot::Cmd("write png $currpic" . ".png");
	if ( ConvertPic($currpic) ) {
	    $pic_files .= " $currpic" . ".gif";
	    $movie_cmd .= " /temp1/$user_name/movie/$currpic" . ".gif";
	    $last_frame = " /temp1/$user_name/movie/$currpic" . ".gif";
	    $counter++;
	}
    }
    print "Done\n";
    
    $movie_cmd = $last_frame . $movie_cmd;

    if (! -d "/ul/$user_name/tmp") {
	mkdir "/ul/$user_name/tmp";
    }
    print "\n3. Creating movie $movienm in avi format\n";
    print "-----------------------------------------\n";
    print "NOTE: Movie will contain " . ($counter +1) . " frames at 5 frames/sec. Total length: ";
    printf("%-5.2f secs\n", ( ($counter + 1) / 5) );
    executeSSHCmd("mkdir -p /temp1/$user_name/movie", "Creating tmp directory");
    print "Copying pictures to tmp directory..";
    if (system "rcp $pic_files $host_name:/temp1/$user_name/movie/") {
	die "$!.. Failure\n";
    }
    print "Sucess\n";
    executeSSHCmd("makemovie -o /temp1/$user_name/movie/$movienm -f avi -r 5" . $movie_cmd, "Running makemovie");

    print "\n4. Moving movie file..";
    chdir "../";
    if (! system "rcp $host_name:/temp1/$user_name/movie/$movienm .") {
	print "Done\n";
    } else {
	print "There was an error: $!\n";
    }
    executeSSHCmd("rm -f /temp1/$user_name/movie/$filebase" . "*.gif", "5. Cleaning up");

}

sub executeSSHCmd(@) {
    my ($curr_cmd, $display_string) = @_;
    my ($result, $kid);

    print "$display_string...";
    $result = ssh("$user_name\@$host_name", "$curr_cmd");
#    do {
#	$kid = waitpid(-1, WNOHANG);
#    } until $kid > 0;

    if ($result) {
	die "Failure\n$!\n";
    } else { 
	print "Sucess\n"; 
        $result = 1;
    }

    return $result;

}

sub ConvertPic($) {
    my ($infile) = $_[0];
    my ($ImgHolder, $result);

    my ($outfile) = $infile . ".gif";
    $infile .= ".png";
#    print "Converting $infile -> $outfile...";
    if (-e $infile) {
	$ImgHolder = Image::Magick->new();
	$result = $ImgHolder->ReadImage($infile);
	if (! $result) {
	    $result = $ImgHolder->Write($outfile);
	    if (! $result) {
		system "rm -f $infile";
#		print "Sucess\n";
		return 1;
	    } else {
		print "WARNING: Couldn't create $outfile: $result\n"; 
		return 0; 
	    }
	} else {
	    print "WARNING: Couldn't read $infile: $result\n";
	    return 0;
	}
    } else {
	print "ERROR: Cannot locate $infile: $!... Frame ignored\n";
    }
}

#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}

# sander_movie.pl - This scripts will create a movie from a sander trajectory file
# inputs:
#       trajectory: location of trajectory file
#       start: starting frame
#       end: ending frame
#       movie_host: the (SGI) host computer on which too run moviemaker
#       has_Mg & has_Na: specifies whether that ion is present
#       topology : location of the topology file
#
# optional inputs:
#       interval: frame interval (default = 1)
#       keep_water: keep the waters in the picture (default = false)
#       keep_ions: keep the sodium and/or magnesium ions (default = false)

use strict;
use File::Basename;
use Packages::General;
use Packages::Namot;
use Packages::MovieGen;
use Packages::ManipImages;

# subroutine declaration
sub GetParms();
sub CreatePdbFiles();
sub CreatePictures();
sub ConvertPictures();
sub GetTime(@);
sub CreateSanderMovie();

# start

die "usage: $0 parameter_file movie_name movie_lenght\n"
    if (!@ARGV or $#ARGV < 2);

my ($parm_file, $movie_name, $movie_length) = @ARGV;
my ($trajectory, $start_frame, $end_frame, $topology, $total_bases, $info_file); 
my ($host, $interval, $keep_water, $keep_ions, $has_Mg, $has_Na, $create_pdb, $create_pics);
my (%PDB_Files, @Pictures, $Sander_Info, $trj_step, $movie_size);

die "Cannot find $parm_file: $!\n"
    if (! -e $parm_file);

die "Invalid movie lengt. Expected integer, got $movie_length\n"
    if (! IsDecimal($movie_length) && ! IsInteger($movie_length) );

$movie_size = "700,700";
$create_pdb = $create_pics = 1;
{
    print "--==== START ===--\n\n";
    print "Step 1: Getting Parameters...";
    GetParms();
    print "Done\nStep 2: Creating PDB Files...";
    CreatePdbFiles();
    print "Done\n";
    print "Step 3: Creating Pictures...";
    CreatePictures();
    print "Done\n";
    print "Step 4: Creating Movie...";
    CreateSanderMovie();
    print "Done\n\n--=== All Tasks Completed Sucessfully ===--\n";
}

sub GetParms() {
# This sub will open the paramater file and extract the relevant parameters
# It will also check the parameters for validity

    my ($in_data);

    $interval = 1;
    $keep_water = $keep_ions = 0;
    $has_Mg = $has_Na = 1;

    open INFILE, $parm_file or die "FAILURE: Cannot open $parm_file: $!\n";
    while (<INFILE>) {
	chomp;
	$in_data = $_;
	if ( $in_data =~ /^Trajectory: (\S+)/ ) {
	    $trajectory = $1;
	}elsif ( $in_data =~ /^Topology: (\S+)/ ) {
	    $topology = $1;
	}elsif ( $in_data =~ /^Start Frame: (\d+)/ )  {
	    $start_frame = $1;
	}elsif ( $in_data =~ /^End Frame: (\d+)/ ) {
	    $end_frame = $1;
	}elsif ( $in_data =~ /^Host: (\w+)/ ) {
	    $host = $1;
	}elsif ( $in_data =~ /^Has Mg: ([0|1])/ ) {
	    $has_Mg = $1;
	}elsif ( $in_data =~ /^Has Na: ([0|1])/ ) {
	    $has_Na = $1;
	}elsif ( $in_data =~ /^Interval: (\d+)/ ) {
	    $interval = $1;
	}elsif ( $in_data =~ /^Keep Water: ([1|0])/ ) {
	    $keep_water = $1;
	}elsif ( $in_data =~ /^Keep Ions: ([0|1])/ ) {
	    $keep_ions = $1;
	}elsif ( $in_data =~ /^Total Bases: (\d+)/ )  {
	    $total_bases = $1;
	}elsif ( $in_data =~ /^Info File: (\S+)/ ) {
	    $info_file = $1;
	}elsif ( $in_data =~ /^Create PDB Files: ([0|1])/ ) {
	    $create_pdb = $1;
	}elsif ( $in_data =~ /^Create Pictures: ([0|1])/ ) {
	    $create_pics = $1;
	}elsif ( $in_data =~ /^Trajectory Step: (\d+)/ ) {
	    $trj_step = $1;
	}elsif ( $in_data =~ /^Movie Size: (\S+)/ ) {
	    $movie_size = $1;
	}
    }

#    print "data: $trajectory $start_frame $end_frame $host $topology $total_bases $info_file\n";
    die "FAILURE: Invalid paramater file $parm_file\n"
	if (! $trajectory || ! $start_frame || ! $end_frame || ! $host || ! $topology || ! $total_bases || ! $info_file || ! $trj_step);
    die "FAILURE: Cannot locate trajectory file $trajectory: $!\n"
	if (! -e $trajectory);
    die "FAILURE: Cannot locate topology file $topology: $!\n"
	if (! -e $topology);
    die "FAILURE: Cannot locate information file $info_file: $!\n"
	if (! -e $info_file);
    ($start_frame, $end_frame) = ($end_frame, $start_frame)
	if ($start_frame > $end_frame);

    $Sander_Info = ParseSanderInfoFile($info_file);
}

sub CreatePdbFiles() {
# This sub will execute Ptraj and create the pdbfiles

    my ($out_string, $i, $curr_file, $counter, $n_step);

    system "mkdir -p pdbfiles";

    $out_string = "trajin $trajectory $start_frame $end_frame $interval\n";
    $out_string .= "center :1-" . $total_bases . " mass origin\n";
    $out_string .= "image origin center\n";
    $out_string .= "strip :WAT\n"
	if (! $keep_water);
    $out_string .= "strip :Na+\n"
	if (! $keep_ions && $has_Na);
    $out_string .= "strip :Mg+\n"
	if (! $keep_ions && $has_Mg);
    $out_string .= "trajout pdbfiles/$movie_name pdb\n";

    open OUTFILE, "> ptraj_in" or die "FAILURE: Error creating temporary file ptraj_in: $!\n";
    print OUTFILE $out_string;
    close OUTFILE;

    if ($create_pdb) {
	if (system ("/ul/maiti/ptraj-6.3/linux/ptraj $topology < ptraj_in >& junk")) {
	    die "FAILURE: A serious error occured when running ptraj. See the junk file for details\n";
	} else {
	    system "rm -f junk";
	}
    }

    $counter = 0;

    for $i (0 .. int(($end_frame - $start_frame)/$interval)) {
	$curr_file = "./pdbfiles/$movie_name" . "." . $i;
	if (-e $curr_file) {
	    $n_step = ($start_frame + ($interval * $i)) * $trj_step;
	    $PDB_Files{$counter} = (
				    {
					"FILE_LOC" => "$movie_name" . "." . $i,
					"TIME" => $Sander_Info->{$n_step}->{"TIME"},
				    }
				    );
	    $counter++;
	}
    }
	    
    $start_frame = 0;
    $end_frame = $i;

    die "FAILURE: Ptraj error: No PDB Files created"
	if ($counter == 0);
}

sub CreatePictures() {
# This sub will create the png files from the PDB Files
    my ($counter, $curr_file, @namot_cmd, $curr_pic, $time_string, $jpeg_file);

    system "mkdir -p pictures/png";
    system "mkdir -p pictures/jpeg";
    for $counter (sort Numerically keys %PDB_Files) {
	$curr_file = "./pdbfiles/" . $PDB_Files{$counter}{"FILE_LOC"};
	$curr_pic = "./pictures/png/" . $PDB_Files{$counter}{"FILE_LOC"} . ".png";
	push @namot_cmd, "load pdb rigid " . $curr_file;
	push @namot_cmd, "write png $curr_pic";
	push @namot_cmd, "close";
	if ($create_pics) {
	    DoNamotCmd(\@namot_cmd);
	}
	if (-e $curr_pic) {
	    $time_string = "Time " . $PDB_Files{$counter}{"TIME"} . " ps";
	    $jpeg_file = "./pictures/jpeg/" . $PDB_Files{$counter}{"FILE_LOC"} . ".jpg";
	    if ($create_pics) {
		AnnotatePic($curr_pic, $time_string, 20, 20, 0, 3, 2);
		ConvertPicture($curr_pic, $jpeg_file);
	    }
	    if (-e $jpeg_file) {
		push @Pictures, $jpeg_file;
	    }
	}
     }
}

sub CreateSanderMovie() {
    my $rec = (
               {
                   "start" => 1,
                   "end"   => ($#Pictures + 1) ,
                   "filebase" => "pictures/jpeg/" . $movie_name . ".",
                   "extension" => "jpg",
                   "host"      => $host,
                   "user"      => $ENV{USER},
                   "savename"  => $movie_name . "_movie.avi",
                   "framerate" => $movie_length/($#Pictures + 1),
                   "moviesize" => $movie_size,
               }
               );
    
    if (! CreateMovie($rec)) {
	die "Error while creating movie\n";
    }
 
}

sub Numerically {
    ($a <=> $b);
}

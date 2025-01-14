#!/usr/bin/perl -w
use strict;
use File::Basename;
# This program will create the p5namot directory
# and will compile all of the C source files stored
# in the file "flelist.txt" and stores them in the
# p5namot directory.
# Then it creates a Makefile for the module, then executes it

#--Start

die "Usage: $0 install_dir file_list\n"
    if (! @ARGV || $#ARGV < 1);

my ($Install_loc, $flename) = @ARGV;
my ($count, $outstring, $compilename);
my (@mylist, @outlist, $delete_files);

#die "Cannot find gnu libraries. LDFLAGS environment variable not set!\n"
#    if (! $ENV{LDFLAGS} );
#die "Cannot find gnu source. CPPFLAGS environment variable not set!\n"
#    if (! $ENV{CPPFLAGS} );

die "Invalid directory: $Install_loc\n"
    if (! -d $Install_loc);
die "Cannot locate file $flename: $!\n"
    if (! -e $flename);

# Check for the "flelist" file

if (!open(FLELIST, $flename)) {
   die "$!\nExecution Terminated.\n";
} else {
    while (<FLELIST>) {
	chomp;
	push @mylist, $Install_loc . "/" . $_;
    }
    close FLELIST;
}

$count =0;
print "\nCompiling " . ($#mylist + 1) . " files....\n";
for $count (0 .. $#mylist) {
    $compilename =  substr(basename($mylist[$count]),0,-2) . ".o";
    $outstring = 'gcc -DHAVE_CONFIG_H -I. -I. -I..  -fpic ' .
                 $ENV{CPPFLAGS} . " -I${Install_loc} " .
				 ' -DLIB_HOME="\"/home/tpascal/codes/namot/2.2.0-pre4//share/namot\"" -DHELP_FILE_DIR="\"/home/tpascal/codes/namot/2.2.0-pre4//share/namot\"" -I/usr/include/python2.7 -I/usr/lib/python2.7/config';
    $outstring .= " -c $mylist[$count] -o p5namot-" . $compilename;

    if (! system($outstring)) {
	push @outlist, "p5namot-" . $compilename;
    }
    print $outstring . "\n";
}

$count = 0;

#now create the Makefile.PL

print "Done...\n\nCreating MakeFile.PL....\n";
open MAKEFLE, "> Makefile.PL" or die "Error! Cannot Create MakeFile.PL";
print MAKEFLE "use ExtUtils::MakeMaker;\n" .
              "     WriteMakefile(\n" .
              "             'NAME'    => 'p5namot'," .
              "                 \# Name of package\n" .
			  "             'INC' => '-I${Install_loc} $ENV{CPPFLAGS} $ENV{CFLAGS}',  \n" .
              "             'LIBS'    => ['-L/usr/X11R6/lib $ENV{LDFLAGS}" . 
			  "-lreadline -lncurses -lutil -ldl -lpthread -lpng -lm -lXm -lXt -lXext -lgsl -lgslcblas -lX11  -lSM -lICE'],  " .
              "\# Name of custom libraries\n" . 
              "             'OBJECT'  => 'p5namot_wrap.o ";
for $count (0 .. $#outlist) {
    print MAKEFLE $outlist[$count] . " ";
}
print MAKEFLE "'  \# Object files\n" . 
              "     );\n";
close MAKEFLE;

print "Done...\n\nCompiling Makefile.PL.....";
system "perl Makefile.PL";
system "make";

print "Done...\n\nRemoving Temporary Files.....";
$delete_files = "rm -fr p5namot-*";
system $delete_files;
print "Done..\n";



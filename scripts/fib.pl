#!/usr/bin/perl -w
use strict;
use Errno qw(EAGAIN);
use POSIX ":sys_wait_h";

sub calcFib;

die "usage: $0 number_of_sequence\n"
    if (! @ARGV);

my ($numSequence) = $ARGV[0]; # get arguments
my ($pid, $kid); #declare variables
$numSequence =~ s/\s+//g; # remove white spaces

die "ERROR: Expected integer for number_of_sequence got $numSequence\n"
    if ($numSequence !~ /^\d+$/ || ! $numSequence); # check for validity

if ($numSequence == 1) {
    print "0\n";
} elsif ($numSequence == 2) {
    print "0 1\n";
} else {
    print "0 1 ";

  FORK: {
      if ($pid = fork) {
	  # this is the parent 
	  # the child process is available in $pid
	  # so wait for child process to complete
	  do {
	      $kid = waitpid(-1,&WNOHANG);
	  } until $kid == -1;
      } elsif (defined $pid) {
	  # this is the child process
	  # so do the fibonaci sequence
	  $numSequence -= 2;
	  calcFib(0, 1, $numSequence);
      } elsif ($! == EAGAIN) {
	  # recoverable fork error so wait for a while and try again
	  sleep 5;
	  redo FORK;
      } else {
	  # something weird happened
	  die "ERROR: Something weird happened while trying to fork!\n";
      }
  }
}

sub calcFib {
    my ($i, $j, $num) = @_;
    my ($newFibNum);

    while ($num > 0) {
	$newFibNum = $i + $j;
	print "$newFibNum ";
	$i = $j;
	$j = $newFibNum;
	$num--;
    }
    print "\n";
}

#!/usr/bin/perl -w
#
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Math::Macopt;

&main();
  
sub main
  {
        # Some settings
        my $N = 10;
        my $epsilon = 0.001;
  
        # Initialize the Macopt 
        my $macopt = new Math::Macopt::Base($N, 1);
  
        # Setup the function and its gradient
        my $func = sub {
                my $x = shift;
  
                my $size = $macopt->size();
                my $sum = 0;
                foreach my $i (0..$size-1) {
                        $sum += ($x->[$i]-$i)**2;
                }
                
                return $sum;
        };
        my $dfunc = sub {
                my $x = shift;
  
                my $size = $macopt->size();
                my $g = ();
                foreach my $i (0..$size-1) {
                        $g->[$i] = 2*($x->[$i]-$i); 
                }
   
                return $g;
        };
        $macopt->setFunc(\&$func);
        $macopt->setDfunc(\&$dfunc);
  
        # Optimizer using macopt 
        my $x = [(1)x($N)];
        $macopt->maccheckgrad($x, $N, $epsilon, 0) ;
        $macopt->macoptII($x, $N);

        # Display the result
        printf "[%s]\n", join(',', @$x);
}

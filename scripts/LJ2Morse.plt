#!/usr/bin/gnuplot --persist

if (ARGC < 2) {
	print "usage: ",ARG0," d0 r0"
	exit(1)
}

set fit quiet
d0 = ARG1
r0 = ARG2
print "Reading parameters... d0 = ",d0," r0 = ",r0
lj(x) = 4*d0*((r0/x)**12-(r0/x)**6)
morse(x) = D*(exp(-2*alpha*(x-r))-2*exp(-alpha*(x-r)))
xstart = 2
xstop = r0 + 6
set samples 1000
set table "_lj.dat"
pl [xstart:xstop] lj(x)
unset table
weight(x)=exp((x+d0)/8.314/.298) #Boltzman weighting 
D = d0
r = r0
alpha = 12
fit morse(x) '_lj.dat' u 1:2:(weight($2)) yerrors via D,r,alpha
print "Result: Morse parameter. D = ",D," r = ",r," alpha = ",alpha
print "Creating comparison plot lj-morse.compare.png"

set encoding iso_8859_1
set term pngcairo size 1400,1050 enhanced color font "Heveltica,40" rounded crop lw 3
set style line 80 lt rgb "#000000" lw 2
set style func linespoints 

# Line style for grid
set style line 81 lt 0 lc rgb "#000000"  lw 2
# axes
set style line 11 lc rgb '#808080' lt 1 lw 2 # grey
set tics nomirror out scale 0.75 tc rgb "black"
set border 3 front ls 11 back ls 80

set out 'lj-morse.compare.png'
set xlabel "dist [\305]"
set ylabel "energy"
set xtics 1
set mxtics 5
set xrange [xstart:xstop]
pl lj(x) w l lt 2 lw 1.5 t 'LJ', morse(x) w l lt 3 lw 1.5 t 'Morse', 0 w l lt 0 lw 1.5 not
system("rm -fr _lj.dat fit.log")

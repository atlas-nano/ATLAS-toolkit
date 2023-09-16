#!/usr/bin/gnuplot --persist

if (ARGC < 3) {
	print "usage: ",ARG0," D r alpha"
	exit(1)
}

set fit quiet
D = ARG1
r = ARG2
alpha = ARG3
print "Reading parameters... D = ",D," r = ",r," alpha = ",alpha
lj(x) = 4*d0*((r0/x)**12-(r0/x)**6)
morse(x) = D*(exp(-2*alpha*(x-r))-2*exp(-alpha*(x-r)))
xstart = r - 0.75
xstop = r + 4
set samples 1000
set table "_morse.dat"
pl [xstart:xstop] morse(x)
unset table
weight(x)=exp((x+D)/8.314/.298) #Boltzman weighting 
d0 = D
r0 = r
#fit lj(x) '_morse.dat' u 1:2:(weight(column(2))) yerrors via d0,r0
fit lj(x) '_morse.dat' via d0,r0
print "Result: Morse parameter. D = ",d0," r = ",r0
print "Creating comparison plot morse-lj.compare.png"

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

set out 'morse-lj.compare.png'
set xlabel "dist [\305]"
set ylabel "energy"
set xtics 1
set mxtics 5
set xrange [xstart:xstop]
pl morse(x) w l lt 2 lw 1.5 t 'Morse', lj(x) w l lt 3 lw 1.5 t 'LJ', 0 w l lt 0 lw 1.5 not
system("rm -fr _morse.dat fit.log")

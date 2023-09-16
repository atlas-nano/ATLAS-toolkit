#!/usr/bin/gnuplot --persist

if (ARGC < 2) {
	print "usage: ",ARG0," d0 r0"
	exit(1)
}

set fit quiet
lj(x) = 4*d0*((r0/x)**12-(r0/x)**6)
x6(x) = A*exp(-x/rho)-C/x**6
d0 = ARG1
r0 = ARG2
xstart = 3
set samples 1000
set table "_lj.dat"
pl [xstart:6] lj(x)
unset table
A = 100000
rho = 0.25
C = 100
fit x6(x) '_lj.dat' via A,rho,C
load '/home/tpascal/scripts/gnuplot_header_small.plt'
set out 'ljx6.compare.png'
set xlabel "dist [\305]"
set ylabel "energy"
set xtics 1
set mxtics 5
set xrange [xstart:6]
f(x) = 0
pl lj(x) w l lt 2 lw 1.5 t 'LJ', x6(x) w l lt 3 lw 1.5 t 'X6', f(x) w l lt 0 lw 1.5 not
set print "ljx6.compare.dat"
print "A = ",A," rho = ",rho," C = ",C

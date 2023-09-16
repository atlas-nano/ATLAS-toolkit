#!/bin/tcsh
#!/bin/csh

if ($#argv < 3) then
	echo "usage: $0 A rho C"
	exit(1)
endif

gnuplot <<DATA
set fit quiet
lj(x) = 4*d0*((r0/x)**12-(r0/x)**6)
x6(x) = A*exp(-x/rho)-C/x**6
A = $1
rho = $2
C = $3
xstart = 2.5
xstop = 6
set samples 1000
set table "_x6.dat"
pl [xstart:xstop] x6(x)
unset table
d0 = 0.6
r0 = 3
fit lj(x) '_x6.dat' via d0,r0
print "d0 = ",d0," r0 = ",r0
load '/home/tpascal/scripts/gnuplot_header_small.plt'
set out 'x6lj.compare.png'
set xlabel "dist [\305]"
set ylabel "energy"
set xtics 1
set mxtics 5
set xrange [xstart:xstop]
f(x) = 0
pl lj(x) w l lt 2 lw 1.5 t 'LJ', x6(x) w l lt 3 lw 1.5 t 'X6', f(x) w l lt 0 lw 1.5 not
DATA

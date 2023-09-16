#!/bin/tcsh
#!/bin/csh

if ($#argv < 3) then
	echo "usage: $0 A rho C"
	exit(1)
endif

gnuplot <<DATA
set fit quiet
morse(x) = d0*(exp(-2*alpha*(x-r0))-2*exp(-alpha*(x-r0)))
x6(x) = A*exp(-x/rho)-C/x**6
A = $1
rho = $2
C = $3
xstart = 1.5
xstop = 6
set samples 1000
set table "_x6.dat"
pl [xstart:xstop] x6(x)
unset table
d0 = GPVAL_DATA_Y_MIN
r0 = 2
alpha = 4
fit morse(x) '_x6.dat' via d0,r0,alpha
print "d0 = ",d0," alpha = ",alpha," r0 = ",r0
load '/home/tpascal/scripts/gnuplot_header_small.plt'
set out 'x6morse.compare.png'
set xlabel "dist [\305]"
set ylabel "energy"
set xtics 1
set mxtics 5
set xrange [xstart:xstop]
f(x) = 0
pl morse(x) w l lt 2 lw 1.5 t 'Morse', x6(x) w l lt 3 lw 1.5 t 'X6', f(x) w l lt 0 lw 1.5 not
DATA

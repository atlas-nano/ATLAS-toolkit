#!/usr/bin/gnuplot --persist

if (ARGC < 3) {
	print "usage: ",ARG0," A rho C"
	exit(1)
}

set fit quiet
A = ARG1
rho = ARG2
C = ARG3
print "Reading parameters... A = ",A," rho = ",rho," C = ",C
lj(x) = 4*d0*((r0/x)**12-(r0/x)**6)
x6(x) = A*exp(-x/rho)-C/x**6
dxX6(x)=-A*exp(-x/rho)/rho+6*C/x**7
#here we perform the 1st derivative test to obtain the X6 catastrophe point, an the minima
dx = 0.025 #x interval
i=0.1;
while (i<20) {
	j = i + dx
	k = i + 2*dx
	e_i=x6(i)
	e_j=x6(j)
	e_k=x6(k)
	d_ij = (e_j-e_i)/dx
	d_jk = (e_k-e_j)/dx
	if (d_ij>0.0 && d_jk<0.0) { #at the inflection point
		x6start=k
	} 
	if(exists("x6start") && !exists("xstart") && e_i<0.0 && d_ij<0.0 && d_jk < 0.0) {
		xstart=i-dx
		if(x6start<(i-dx*25)) {
			xstart=i-dx*25
		}
	} 
	if(d_ij<0.0 && d_jk>0.0) { #minima
		r0=j
		d0=-1*e_j	
		i=21	
	}
	i = i + dx
}
print "x6start ",x6start
if (i==20) {
	print "ERROR: Could not find the X6 catastrophe point..."
	exit(1)
}
if (!exists("xstart")) {
	print "ERROR: Count not figure out the fitting starting x value"
	exit(1)
}
print "Found inflection point at ",x6start
print "Starting fitting from ",xstart," energy minima at r0 = ",r0," with value d0 = ",d0
print "Starting fit from ",xstart
xstop = r0 + 6
set samples 1000
set table "_x6.dat"
pl [xstart:xstop] x6(x)
unset table
weight(x)=exp((x+d0)/8.314/.298) #Boltzman weighting 
fit lj(x) '_x6.dat' u 1:2:(weight($2)) yerrors via d0,r0
#fit lj(x) '_x6.dat' u 1:2 via d0,r0
print "Result: LJ parameter. d0 = ",d0," r0 = ",r0*2**(1/6)
print "Creating comparison plot x7lj.compare.png"

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

set out 'x6lj.compare.png'
set xlabel "dist [\305]"
set ylabel "energy"
set xtics 1
set mxtics 5
xstart = xstart + 0.5
set xrange [xstart:xstop]
f(x) = 0
pl lj(x) w l lt 2 lw 1.5 t 'LJ', x6(x) w l lt 3 lw 1.5 t 'X6', f(x) w l lt 0 lw 1.5 not
system("rm -fr _x6.dat fit.log")

#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 evsv_file [save_prefix]"
	exit(1)
endif

if !(-e $1) then
	echo "ERROR: Cannot locate $1"
	exit(1)
endif

set prefix = "${1}.EvsV"
if ($#argv > 1) set prefix = $2

#get some vals
set min = (`cat $1 | awk 'BEGIN{ymin=999999999999999999;}{if($2<ymin) {ymin=$2;xmin=$1}}END{print xmin,ymin}'`)
echo $min
set start = `grep -n " $min[2]" $1 | awk '{print $1-15}'`
set start = `head -${start} $1 | tail -1 | awk '{print $1}'`
set stop = `grep -n " $min[2]" $1 | awk '{print $1+15}'`
set stop = `head -${stop} $1 | tail -1 | awk '{print $1}'`

gnuplot<<DATA;
f(x)=E0+B0*V0*((1/dB0/(dB0-1))*(x/V0)**(1-dB0)+x/dB0/V0-1/(dB0-1))
#f(x) = E0 + (9*V0*B0/16)*((((V0/x)**(2/3))-1)**2*(6+dB0*(((V0/x)**(2/3))-1)-4**(2/3)))
E0=$min[2]
V0=$min[1]
dB0=4
B0=0.5
set fit quiet
fit [${start}:${stop}] f(x) '$1' u 1:2 via E0,V0,B0,dB0
set xlabel "Volume [\305^3]"
set ylabel "Energy [eV/atom]"
title_p = sprintf("VASP energies Fit to Murnaghan EOS\nB_0=%.3f GPa, E_0=%.3f B_0'=%.3f V_0=%.3f",B0*160.02176,E0,dB0,V0)
set key top left
load '/home/tpascal/scripts/gnuplot_header.plt'
set encoding iso_8859_1
set out '${prefix}.png'
set pointsize 2.5
set title title_p
set multiplot
set size 1.0,1.0
set origin 0.0,0.0
pl '${1}' u 1:2 w lp ls 1 t "data", '' u 1:2 w p pt 6 ps 2.0 lc rgb "white" not, f(x) w l t "fit"
set origin 0.4,0.3
set size 0.5,0.5
set xrange [${start}:${stop}]
set xtics 3
set ytics 0.3
unset key
unset title
unset xlabel
unset ylabel
rep
unset multiplot
DATA

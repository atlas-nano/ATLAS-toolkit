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

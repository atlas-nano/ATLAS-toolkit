set encoding iso_8859_1
set term pngcairo size 1400,1050 enhanced color font "Verdana,40" rounded crop lw 2
set style line 80 lt -1 lc rgb "#000000" lw 3
set style func linespoints 

# Line style for grid
#set style line 81 lt 0  # dashed
#set style line 81 lt rgb "#000000"  # grey
# axes
set style line 11 lc rgb '#808080' lt 1
set border 3 front ls 11
set tics nomirror out scale 0.75
# grid
set style line 12 lc rgb'#808080' lt 0 lw 1
#set grid back ls 12
#set grid xtics ytics mxtics

#set grid back linestyle 81

set border 31 back linestyle 80 # Remove border on top and right.  These
             # borders are useless and make it harder
             # to see plotted lines near the border.
    # Also, put it in grey; no need for so much emphasis on a border.
#set xtics nomirror
#set ytics nomirror

#set log x
set mxtics 2    # Makes logscale look good.
set mytics 2    # Makes logscale look good.

# Line styles: try to pick pleasing colors, rather
# than strictly primary colors or hard-to-see colors
# like gnuplot's default yellow.  Make the lines thick
# so they're easy to see in small plots in papers.
set style line 1 lt rgb "#000000" lw 1 pt 7 
set style line 2 lt 1 lc rgb '#FB9A99' # light red
set style line 3 lt 1 lc rgb '#A6CEE3' # light blue
#set style line 4 lt rgb "fuchsia" lw 2 pt 5 
set style line 5 lt rgb "#025214" lw 2 pt 9 
set style line 6 lt rgb "#031A49" lw 2 pt 11
set style line 7 lt rgb "#6D0D03" lw 2 pt 13
set style line 8 lt rgb "#E69F17" lw 2 pt 4
set style line 9 lt rgb "#8F234D" lw 2 pt 8 ps 0.2
set style line 10 lt rgb "#E62B17" lw 2 pt 10 
set style line 11 lt rgb "#8F643F" lw 2 pt 12 
set style line 12 lt rgb "#2F6C3D" lw 2 pt 5
set style line 13 lt rgb "#8F463" lw 2 pt 7


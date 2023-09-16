load '/home/tpascal/scripts/gnuplot_header_small.plt'
set tmargin at screen .7
set bmargin at screen .3
set lmargin at screen .3
set rmargin at screen .7
set term pngcairo size 3000,2400 transparent enhanced color font "Verdana,40" rounded crop lw 3
set palette rgb 21,22,23 negative
set view map
set contour base
#set dgrid3d 200 200
f_list = "`ls *.oh-oh.sym.dat | tr '\n' ' '`"
set xtics 0.5 in mirror
set ytics 0.5 in mirror
set mxtics 3
set mytics 3
#set xlabel "P[OH_1]"
#set ylabel "P[OH_2]"
unset key
set border 15
#set colorbox horizontal user origin 0.1,0.0 size 0.5,0.05
#set mcbtics 2
set grid cb
set isosamples 1000
set samples 1000
f(x) = x*1E3-5
do for [j=1:words(f_list)] {
	ncol=system(sprintf("tail -2 %s | head -1 | awk '{print NF-2}'",word(f_list,j)))
	t_list = "surface "
	if(ncol>1) {
		if(ncol>2) {
			t_list = "surface subsurface"
		}
		do for [i=3:(ncol-1)] {
			t_list = sprintf("%s layer%d",t_list,i)
		}
		t_list = sprintf("%s bulk",t_list)
	}
	print word(f_list,j)," ",ncol,"cols: ",t_list
	do for [i=1:words(t_list)] {
		print "	",word(t_list,i)
		set out system(sprintf("echo %s | sed 's/.dat/.%s.png/'",word(f_list,j),word(t_list,i)))
		if (ncol>2) { 
#set cbrange [-1:1]
			set palette defined (0 "red", 1 "white", 2 "blue")	
			spl word(f_list,j) u 1:2:(column(i+2)/column(ncol+2)) w pm3d lw 2 dt 2
		} else {
			set palette color
#set cbrange [0:]
			spl word(f_list,j) u 1:2:(f(column(i+2))) w pm3d lw 2 dt 2
		}
	}
}

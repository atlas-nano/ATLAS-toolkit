#!/bin/tcsh
#!/bin/csh

if ($#argv < 3) then
	echo "usage: $0 datafile(s) col(s)[use '+' to average, ',' to group, space for multiple, eg '1+5,2 4'] col_header(s) grp_header(s) [save prefix]"
	exit(1)
endif

set list = (`ls $1`)
if ($#list == 0) then
	echo "ERROR: No valid files found while searching '$1'"
	exit(1)
endif

set sstr = ($2)
set sprefix = `basename $list[1]`
set sprefix = $sprefix:r
if ($#argv > 4) set sprefix = $5

@ f = 1
set ctitle = ()
set ftitle = ()
set gcounter = ()
set colheaders = ($3)
#test grp headers
@ hasgrpheaders = 0
set grpheaders = ()
if ($#argv > 3) then
	set hasgrpheaders = 1
	set grpheaders = ($4)
	set grps = ($sstr)
	if !($#grpheaders == $#grps) then
		echo "ERROR: Number of grps ($#grps) does not equal number of grp headers ($#grpheaders)"
		exit(1)
	endif
	echo $grpheaders
endif

#test col headers
@ tcols = $#colheaders
foreach i ($sstr)
	set nc = (`echo $i | sed 's/,/ /g'`)
	if ($#nc <> $tcols && $tcols>1) then
		echo "ERROR: Number of column headers ($tcols) is not equal to '$i'($#nc)"
		exit(1)
	endif
end

foreach i ($list)
	if !(-e $i) continue
	set mol = `basename $i`
	#set tmp = (`echo $mol | sed 's/\./ /g'`)
	set tmp = (`echo $mol | sed 's/\.hbond.DAStats.dat//g'`)
	if ($#tmp > 2) then
		set mol = "$tmp[1].$tmp[2]"
	else
		set mol = $tmp[1]
	endif
	echo "$mol"
	echo "#index $colheaders" > ${mol}.stats.dat
	set ftitle = ($ftitle $mol)
	@ g = 1
	foreach j ($sstr)
		set grps = (`echo $j | sed 's/,/ /g'`)
		@ c = 0
		set clist = ()
		set xpos = ()
		set ypos = ()
		set stats = ()
		set std = ()
		set alist = ()
		echo "	$grpheaders[$g]"
		echo "#x y">labels.max.pos.dat
		foreach k ($grps)
			if !($hasgrpheaders) set grpheaders = ($grpheaders $g)
			set cols = (`echo $k | sed 's/\+/ /g' | sed 's/[^0-9\. ]//g'`)
			@ l = $c + 1
			echo "		$colheaders[$l] $cols"
			egrep '^\s*[0-9]' $i | awk -v var="$cols" 'BEGIN{split(var,cols," ")}{ c = 0; avg = 0; for(i in cols) {avg += $cols[i]; c++} avg /= c; print NR,avg}' > __tmp.dat
			cat __tmp.dat | awk '{v = sprintf("%.0f\n",$2/0.001); hist[v]++;}END{for (i in hist) printf "%s %s\n",i*0.001,hist[i]}' | sort -n > __f${f}.g${g}.c${c}.dat
			csh -f ~/scripts/fit_gaussian.csh __f${f}.g${g}.c${c}.dat __f${f}.g${g}.c${c}.fit.dat >& /dev/null
			set aval = `egrep '^a ' fit.log | grep '=' | tail -1 | awk '{print $3}'`
			set bval = `egrep '^b ' fit.log | grep '=' | tail -1 | awk '{printf "%.5f\n", $3}'`
			set cval = `egrep '^c ' fit.log | grep '=' | tail -1 | awk '{std = $3; if(std<0) std = -$3; printf "%.5f\n", std}'`
			#fix for nd
			if (`echo $bval | awk '{if($1==0)print 1;else print 0}'` && `echo $cval | awk '{if($1==0)print 1; else print 0}'`) then
				sed -i 's/^0 .*$//' __f${f}.g${g}.c${c}.dat
				csh -f ~/scripts/fit_gaussian.csh __f${f}.g${g}.c${c}.dat __f${f}.g${g}.c${c}.fit.dat '[-0.01:1]' >& /dev/null
				set aval = `egrep '^a ' fit.log | grep '=' | tail -1 | awk '{print $3}'`
				set bval = `egrep '^b ' fit.log | grep '=' | tail -1 | awk '{printf "%.5f\n", $3}'`
				set cval = `egrep '^c ' fit.log | grep '=' | tail -1 | awk '{std = $3; if(std<0) std = -$3; printf "%.5f\n", std}'`
			endif
			
			rm -fr __f${f}.g${g}.c${c}.dat __tmp.dat fit.log
			if !(-e __f${f}.g${g}.c${c}.fit.dat && -s __f${f}.g${g}.c${c}.fit.dat) then
				set stats = ($stats 0 0)
				set std = ($std 0)
				set alist = ($alist 0)
				set clist = ($clist __blank.dat)
				echo "0 0 0" > __blank.dat
				echo "0 0" | awk '{print $1,$2*$1}' >> labels.max.pos.dat
				rm -fr __f${f}.g${g}.c${c}.fit.dat
			else
				echo "$bval $aval" | awk '{print $1,$2*$1}' >> labels.max.pos.dat
				set stats = ($stats $bval $cval)
				set std = ($std `echo $cval | awk '{printf "%.2f\n",$1}'`)
				set clist = ($clist __f${f}.g${g}.c${c}.fit.dat)
				egrep '^\s*[0-9]'  __f${f}.g${g}.c${c}.fit.dat | awk '{print $1,$2*'$bval'}' > _tmp.dat
				mv _tmp.dat  __f${f}.g${g}.c${c}.fit.dat
				set alist = ($alist $aval)
			endif
			@ c += 1
		end
		@ n = $f * $g
		#nd fix
		if ($tcols == 3) set stats = (`echo $stats | awk '{nd=1-$3-$5;print nd; for(i=2;i<=NF;i++) print $i}'`)
		set grplabel = $grpheaders[$g]
		echo "$grplabel $stats" >> ${mol}.stats.dat
		set tlist = ($colheaders)
		set gcounter = ($gcounter $c)
		set ctitle = ($ctitle $tlist)
		gnuplot <<DATA;
load '/home/tpascal/scripts/gnuplot_header.plt'
load '/home/tpascal/scripts/gnuplot_line_styles_pub_new.plt'
set term pngcairo size 1400,1050 enhanced color font "Verdana,32" rounded crop lw 2 
tlist_str="$tlist"
clist_str="$clist"
std_str="$std"
clist(n)=sprintf("%s",word(clist_str,n))
title(n)=sprintf("%s",word(tlist_str,n))
slabel(index,n)=sprintf("%.2f (%s)%",index,word(std_str,n))
avg(n)=sprintf("%s",word(avg_str,n))
set ylabel "frequency [arb. units]"
set xrange [0:1]
set key top left horiz out
set ytics nomirror
set xtics nomirror
set border 3
set style data lines 
set style fill transparent pattern 4 bo
set out "${mol}.${grplabel}.png"
#pl for [ii=1:$c] clist(ii) u 1:(column(2)) w filledcu lt ii lw 3 t title(ii), 'labels.max.pos.dat' u 1:2:(slabel(column(1),ii)) w labels center offset 0,1 not
pl for [ii=1:$c] clist(ii) u 1:(column(2)) w filledcu lt ii lw 3 t title(ii)
DATA

		rm -fr $clist labels.max.pos.dat
		@ g += 1
	end
	@ f += 1
end
set n = 1
@ nbars = ($f - 1) * ($g - 1) * ($c - 1)
set xoffset = `echo $nbars | awk '{printf "%.3f\n",1/($1+1)}'`
@ xsize = ($f - 1) * 1400
foreach f (`seq 1 $#list`)
	if !(-e $list[$f]) continue
	set mol = `basename $list[$f]`
	set tmp = (`echo $mol | sed 's/\./ /g'`)
	set tmp = (`echo $mol | sed 's/\.hbond.DAStats.dat//g'`)
	if ($#tmp > 2) then
		set mol = "$tmp[1].$tmp[2]"
	else
		set mol = $tmp[1]
	endif
	if ($#list == 1) then
		set pstr = "pl '${mol}.stats.dat' "
	else if ($n == 1) then
		set pstr = "pl newhistogram '$mol', '${mol}.stats.dat' "
	else 
		set pstr = "$pstr, newhistogram '$mol', '${mol}.stats.dat' "
	endif
	foreach i (`seq 1 $tcols`)
		if ($i > 1) set pstr = "$pstr, ''" 
		@ j = $i * 2
		@ k = $j + 1
		if ($n == 1) then
			#set pstr = "$pstr u (column($j)):(column($k)):xtic(1) w histogram fs pattern $i lt $i t '$colheaders[$i]', '' u (column(0)):(column($j)):(plabel(column($j))) w labels center offset 0,1 not"
			set pstr = "$pstr u (column($j)):(column($k)):xtic(1) w histogram fs pattern $i lt $i t '$colheaders[$i]'"
		else
			#set pstr = "$pstr u (column($j)):(column($k)):xtic(1) w histogram fs pattern $i lt $i not, '' u (column(0)):(column($j)):(plabel(column($j))) w labels center offset 0,1 not"
			set pstr = "$pstr u (column($j)):(column($k)):xtic(1) w histogram fs pattern $i lt $i not"
		endif
	end
	@ n += 1
end

set xs = `echo $n | awk '{val=sprintf("%.0f",$1/2); print 1400*val}'`
cat >${sprefix}.plt <<DATA;
load '/home/tpascal/scripts/gnuplot_header.plt'
load '/home/tpascal/scripts/gnuplot_line_styles_pub_new.plt'
set term pngcairo size $xs,1050 enhanced color font "Verdana,40" rounded crop lw 2
set auto x
set style data histograms
set style histogram errorbars linewidth 1 gap 3 title offset character 0, 0, 0
set style fill transparent solid 0.9 bo
plabel(index)=sprintf("%.2f",index)
set xtics nomirror rotate by -45 font ",30"
set ytics nomirror 0.3
set yrange [0:]
set border 3
set bars 0.3 front
set key top center horiz out
set out "${sprefix}.png"
$pstr
DATA

gnuplot < ${sprefix}.plt

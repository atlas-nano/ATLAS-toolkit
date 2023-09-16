#!/bin/tcsh
#!/bin/csh

if ($#argv < 3) then
	echo "usage: $0 bgf_file lmp_trj_file 'fftype(s)' [sindex] [save_name]"
	exit(1)
endif

if !(-e $1) then
	echo "ERROR: Cannot locate $1"
	exit(1)
endif
if !(-e $2) then
	echo "ERROR: Cannot locate $2"
	exit(1)
endif

set fftypes = ($3)
set sindex = 2
if ($#argv > 3) set sindex = $4
set prefix = `basename $1`
set prefix = $prefix:r
if ($#argv > 4) set prefix = $5

rm -fr __tmp1.dat __tmp2.dat __header.dat
foreach i (`seq 1 $#fftypes`)
	set t = $fftypes[$i]
	echo 'set t'$i' [atomselect top "type '$t'"]' >> __tmp1.dat
	foreach j (`seq 1 $#fftypes`)
		set u = $fftypes[$j]
		echo >> __tmp2.dat
		echo 'set t'${i}'t'${j}' [measure gofr $t'${i}' $t'${j}' rmax 10 delta 0.1  usepbc 1 selupdate 0 first '$sindex' step 1]'  >> __tmp2.dat
		cat >> __tmp2.dat <<DATA
set outfile [open ${prefix}.gofr.${t}.${u}.dat w]
set r [lindex \$t${i}t${j} 0]
set gr2 [lindex \$t${i}t${j} 1]
set igr [lindex \$t${i}t${j} 2]
set i 0
foreach j \$r k \$gr2 l \$igr {
	puts \$outfile "\$j \$k \$l"
}
close \$outfile
DATA
	end
end

cat > __header.dat <<DATA
mol new $1 type bgf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile $2 type lammpstrj first 1 last -1 step 1 filebonds 1 autobonds 1 waitfor all
DATA

cat __header.dat __tmp1.dat __tmp2.dat > __cmd.vmd
/home/tpascal/codes/bin/vmd -dispdev none -nt < __cmd.vmd > /dev/null
#rm -fr __tmp1.dat __tmp2.dat __header.dat __cmd.vmd

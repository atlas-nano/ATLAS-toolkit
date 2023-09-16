#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
	echo "usage: $0 'thermo_file(s)' density_file [save_prefix]"
	exit(1)
endif

set thermo_files = ($1)
if ($#thermo_files == 0) then
	echo "ERROR: No valid file found while searching $1"
	exit(1)
endif

if !(-e $2) then
	echo "ERROR: Cannot locate $2"
	exit(1)
endif

set surface_vals = (`csh -f /home/tpascal/scripts/gnuplot_find_surface.csh $2 0.99`)

set fields = (Aq_ Eq_ ZPE Emd Cvq Cvc Sq_ fluid Diff)
set types = (d l t r v tot)

#get the number of fields in each file
set min_fields = `grep nmolecules $1 | awk 'BEGIN{ max_i = 0; } { c[NF]++; if(NF>max_i) max_i = NF;}END{max=0; i=1; for(i in c) { if(c[i]>max) { max=c[i]; max_i = i; }} print max_i;}'`
set thermo_files = (`grep nmolecules $1 | awk '{if(NF=='$min_fields') print $1}' | sed 's/:$//'`)
set fn_prefix = `set string1="$thermo_files[1]"; set string2="$thermo_files[$#thermo_files]"; printf "%s\n%s\n" "$string1" "$string2" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/' | sed 's/\.\w*$//'`
set fn_prefix = `basename $fn_prefix`
set ngrps = `grep property $thermo_files[1] | awk '{print $NF}' | sed 's/Tot\[G//' | sed 's/\]//' | awk '{print $1+0}'`
@ nfields = ($min_fields - 1) / $ngrps
@ start_field = 1
echo "${fn_prefix}: Found $#thermo_files files. Each file has $ngrps groups and $nfields fields per group"
echo "surface vals $surface_vals"
set prefix = $fn_prefix
if ($#argv > 2) set prefix = $3

#extract data and save to file
foreach f (nmol $fields)
	set c = 0
	set flist = ()
	foreach i (`seq $start_field $nfields`)
		@ c = $c + 1
		@ j = $i + 2
		set p = $types[$c]
		set fname = "${prefix}.${f}${p}.full.dat"
		grep $f $1 | sed "s/.*$fn_prefix.//" | sed 's/\..*://' | awk '{ if(NF=='$min_fields' && $1 >= 1000000) { j = 0; for(i=3;i<NF;i++) { if($i>100000000 || $i < -10000000) { j++; } } if (j==0) { printf "%s ",$1; i = '$j'; while(i<=(NF-'$nfields')) { printf "%s ", $i; i += '$nfields';} printf "\n"} } }' | sort -n > $fname
		if ($f == "nmol") then
			#get stats for nmol
			set stats_file = $fname
			set stats = (`awk -f ~/scripts/calcStats.csh $stats_file`)
			#save the stats in a gnuplot readable file
			awk -f ~/scripts/calcStats.csh $stats_file | awk '{i = 1; while(i<NF) { j = i + 1; printf "%s %s\n",$i,$j; i += 2; } }' > __${i}.dat
			echo $stats | awk '{ printf "%s ","#AVG"; i = 1; while(i<=NF) { printf "%s ",$i; i += 2; } printf "\n";}' >> $stats_file
			echo $stats | awk '{ printf "%s ","#DEV"; i = 2; while(i<=NF) { printf "%s ",$i; i += 2; } printf "\n";}' >> $stats_file
			set flist = ($flist "__${i}.dat")
		endif
	end
	if ($f == "nmol") then
		paste $flist > __tmp.dat
		echo "#group $types" > ${prefix}.${f}.stats.dat
		echo "#surface_vals $surface_vals $ngrps" >> ${prefix}.${f}.stats.dat
		cat __tmp.dat | awk '{printf "%s ", NR; for(i=1;i<=NF;i++) printf "%s ",$i; printf "\n"; }' >> ${prefix}.${f}.stats.dat
		rm -fr $flist __tmp.dat
	endif
end

#normalize data
foreach f ($fields)
	set c = 0
	set flist = ()
	foreach i (`seq $start_field $nfields`)
		@ c = $c + 1
		set p = $types[$c]
		set nmol_file = "${prefix}.nmol${p}.full.dat"
		if ($i < 3) set nmol_file = "${prefix}.nmolt.full.dat"
		set thermo_file = "${prefix}.${f}${p}.full.dat"
		set norm_file = "${prefix}.${f}${p}.norm.dat"
		join $thermo_file $nmol_file > __tmp.dat
		cat __tmp.dat | awk '{printf "%s ",$1; i = 2; j = 0; while(j<NF) { j = i + '$ngrps' - 1; printf "%s ",($i/$j); i++; } printf "\n";}' > $norm_file
		if ($f == "fluid" || $f == "Diff") cp $thermo_file $norm_file
		rm -fr __tmp.dat
		#calculate statistics for norm
		set stats_file = $norm_file
		set stats = (`awk -f ~/scripts/calcStats.csh $stats_file`)
		#save the stats in a gnuplot readable file
		awk -f ~/scripts/calcStats.csh $stats_file | awk '{i = 1; while(i<NF) { j = i + 1; printf "%s %s\n",$i,$j; i += 2; } }' > __${i}.dat
		echo $stats | awk '{ printf "%s ","#AVG"; i = 1; while(i<=NF) { printf "%s ",$i; i += 2; } printf "\n";}' >> $stats_file
		echo $stats | awk '{ printf "%s ","#DEV"; i = 2; while(i<=NF) { printf "%s ",$i; i += 2; } printf "\n";}' >> $stats_file
		set flist = ($flist "__${i}.dat")
	end
	paste $flist > __tmp.dat
	echo "#group $types" > __junk.dat
	echo "#surface_vals $surface_vals $ngrps" >> __junk.dat
	cat __tmp.dat | awk '{printf "%s ", NR; for(i=1;i<=NF;i++) printf "%s ",$i; printf "\n"; }' >> __junk.dat
	join __junk.dat ${prefix}.nmol.stats.dat > ${prefix}.${f}.stats.dat
	rm -fr $flist __tmp.dat __junk.dat
end

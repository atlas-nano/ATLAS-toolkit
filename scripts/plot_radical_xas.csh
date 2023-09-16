#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 prefix"
	exit(1)
endif

foreach i ("-UP" "-DN")
	set list = ("${1}.*.xas-${i}_1")
	if ($#list == 0) then
		echo "ERROR: Cannot locate any valid files while searching '${1}.*.${i}*")
		exit(1)
	endif
end

foreach i (1 2)
	foreach j (1 2)
		awk -f ~/scripts/calc_col_avg.awk ${1}.*.xas-UP_${i} ~/scripts/calc_col_avg.awk ${1}.*.xas-DN_${j} > ${1}.Spectrum-Avg_${i}_${j}
	end
end

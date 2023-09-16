#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
	echo "usage: $0 prefix num (has_spin=1) (save_prefix)"
	exit(1)
endif

set prefix = "$1"
@ num = $2 
set sprefix = "avg"
set has_spin = 1
if ($#argv > 2) set has_spin = $3
if ($#argv > 3) set sprefix = $4

#foreach i (`seq -f '%01.0f' 45 $num`)
foreach i (`seq -f '%03.0f' 1 $num`)
	foreach j (`seq 1 2`)
		foreach k (UP DN)
			set lfile = (`find . -name "${prefix}.*.S${i}-XCH.xas.5.xas-${k}_${j}"`)
			if !($has_spin) set lfile = (`find . -name "${prefix}.S${i}-XCH.xas.5.xas"`)
			if ($#lfile == 0) then
				echo "ERROR: Cannot locate file ${prefix}.*.S${i}-XCH.xas.5.xas-${k}_${j}"
				exit(1)
			endif
		end
	end
end

set term_atoms = (`echo 1 $num | awk '{printf "%03.0f %03.0f\n",$1,$2}'`)
@ num_m1 = $num - 1
set inter_atoms = (`seq -f '%03.0f' 2 $num_m1`) 

if ($has_spin) then
	set list = ""
	foreach i ($term_atoms)
		set list = "$list ${prefix}.*.S${i}-XCH.xas.5.xas-UP_1"
	end
	awk -f ~/scripts/calc_col_avg.awk $list > ${sprefix}.terminal.alpha.Spectrum-Avg-S
	set list = ""
	foreach i ($term_atoms)
		set list = "$list ${prefix}.*.S${i}-XCH.xas.5.xas-DN_2"
	end
	awk -f ~/scripts/calc_col_avg.awk $list > ${sprefix}.terminal.beta.Spectrum-Avg-S
	set list = ""
	foreach i ($term_atoms)
		set list = "$list ${prefix}.*.S${i}-XCH.xas.5.xas-DN_2 ${prefix}.*.S${i}-XCH.xas.5.xas-UP_1"
	end
else
	set list = ""
	foreach i ($term_atoms)
		set list = "$list ${prefix}.*.S${i}-XCH.xas.5.xas"
	end
endif
awk -f ~/scripts/calc_col_avg.awk $list > ${sprefix}.terminal.Spectrum-Avg-S

if ($has_spin) then
	set list = ""
	foreach i ($inter_atoms)
		set list = "$list ${prefix}.*.S${i}-XCH.xas.5.xas-UP_1"
	end
	awk -f ~/scripts/calc_col_avg.awk $list > ${sprefix}.central.alpha.Spectrum-Avg-S
	set list = ""
	foreach i ($inter_atoms)
		set list = "$list ${prefix}.*.S${i}-XCH.xas.5.xas-DN_2"
	end
	awk -f ~/scripts/calc_col_avg.awk $list > ${sprefix}.central.beta.Spectrum-Avg-S
	set list = ""
	foreach i ($inter_atoms)
		set list = "$list ${prefix}.*.S${i}-XCH.xas.5.xas-DN_2 ${prefix}.*.S${i}-XCH.xas.5.xas-UP_1"
	end
else
	set list = ""
	foreach i ($inter_atoms)
		set list = "$list ${prefix}.*.S${i}-XCH.xas.5.xas"
	end
endif
awk -f ~/scripts/calc_col_avg.awk $list > ${sprefix}.central.Spectrum-Avg-S

set list = "${prefix}.*.S*-XCH.xas.5.xas-UP_1 ${prefix}.*.S*-XCH.xas.5.xas-DN_2" 
if !($has_spin) set list = "${prefix}.*.S*-XCH.xas.5.xas ${prefix}.*.S*-XCH.xas.5.xas" 
awk -f ~/scripts/calc_col_avg.awk ${list} > ${sprefix}.Spectrum-Avg-S

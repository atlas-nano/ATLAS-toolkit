#!/bin/tcsh
#!/bin/csh

if ($#argv < 1) then
	echo "usage: $0 cif_file|bgf_file [mol_prefix] [cell_replication=2x2x2] [MD_temp=298K] [nodes=6] [ppn=8] [cp2k_template_file] [pbs_template_file]"
	exit(1)
endif

set cif_file = $1
if !(-e $cif_file) then
	echo "ERROR: Cannot locate CIF file $1"
	exit(1)
endif
set prefix = `basename $cif_file`
set prefix = $prefix:r
if ($#argv > 1) set prefix = $2

set replication = "2 2 2"
if ($#argv > 2) set replication = `echo $3 | sed 's/[a-z]*/ /gi'`
set rarray = ($replication)
if ($#rarray != 3) then
	echo $#rarray
	echo "ERROR: Expected 'i j k' for replication. Got $replication"
	exit(1)
endif
set repstr = `echo $replication | sed 's/ /x/g'`
if ($repstr == "1x1x1") then
	set repstr = ""
else
	set repstr = ".${repstr}"
endif

set temp = 298
if ($#argv > 3) set temp = `echo $4 | sed 's/[a-z]*//gi'`

set nodes = 6
if ($#argv > 4) set nodes = $5

set ppn = 8
if ($#argv > 5) set ppn = $6

@ nprocs = $nodes * $ppn

set cp2k_template_file = ~tpascal/scripts/cp2k.npt.dftd3.template.in
if ($#argv > 6 && -e $7 && ! -d $7) set cp2k_template_file = $7

set pbs_template_file = ~tpascal/scripts/template.solids.cp2k.shirley.slurm
if ($#argv > 7 && -e $8 && ! -d $8) set pbs_template_file = $8

echo
echo "=================================="
echo "Running with the following options"
echo "=================================="
echo "INPUT FILE: $1"
echo "SAVE PREFIX: $2"
echo "CELL REPLICATION: $3"
echo "MD TEMPERATURE: $temp"
echo "#NODES: $nodes"
echo "#PPN: $ppn"
echo "CP2K TEMPLATE FILE: $cp2k_template_file"
echo "QSCRIPT TEMPLATE FILE: $pbs_template_file"
echo

set ext = $cif_file:e
if ($ext == "cif") then 
	csh -f ~tpascal/scripts/cif2bgf.csh $cif_file ${prefix}.bgf > /dev/null || goto error
	if (`echo $rarray | awk '{if ($1 == 1 && $2 == 1 && $3 == 1) { print 0; }else { print 1;} }'`) then
		~tpascal/scripts/replicate.pl -b ${prefix}.bgf -d "$replication" -s ${prefix}${repstr}.bgf >& /dev/null || goto error
	else
		cp ${prefix}.bgf ${prefix}${repstr}.bgf
	endif
	set bgf_file = ${prefix}${repstr}.bgf
else
	set bgf_file = $cif_file
endif

babel -ibgf $bgf_file -oxyz ${prefix}${repstr}.xyz >& /dev/null || goto error
set cell = (`grep "^CRYSTX" $bgf_file | awk '{print $2,$3,$4,$5,$6,$7}'`)
set bg = `basename $bgf_file`
sed -i "s/^${bg}/${prefix}${repstr} CELL: $cell/" ${prefix}${repstr}.xyz
cat ${prefix}${repstr}.xyz | awk '{if(NR>2) { printf "%-4s %6.5f %6.5f %6.5f\n",$1,$2,$3,$4} }' > ${prefix}${repstr}.cp2k.xyz
cp $cp2k_template_file ${prefix}${repstr}.${temp}K.cp2k.in
sed -i "s/LA_HERE/$cell[1]/" ${prefix}${repstr}.${temp}K.cp2k.in
sed -i "s/LB_HERE/$cell[2]/" ${prefix}${repstr}.${temp}K.cp2k.in
sed -i "s/LC_HERE/$cell[3]/" ${prefix}${repstr}.${temp}K.cp2k.in
sed -i "s/ALPHA_HERE/$cell[4]/" ${prefix}${repstr}.${temp}K.cp2k.in
sed -i "s/BETA_HERE/$cell[5]/" ${prefix}${repstr}.${temp}K.cp2k.in
sed -i "s/GAMMA_HERE/$cell[6]/" ${prefix}${repstr}.${temp}K.cp2k.in
sed -i "s/TEMP_HERE/$temp/" ${prefix}${repstr}.${temp}K.cp2k.in
sed -i "s/PREFIX_HERE/${prefix}${repstr}/" ${prefix}${repstr}.${temp}K.cp2k.in

echo "" > ${prefix}.kind.dat
rm -fr __tmp.dat ${prefix}.kind.*.dat
foreach j (~tpascal/codes/cp2k/cp2k/data/GTH_BASIS_SETS ~tpascal/codes/cp2k/cp2k/data/BASIS_MOLOPT)
	set valid = 1;
	set basis_set_file = $j
	foreach i (`cat ${prefix}${repstr}.cp2k.xyz | awk '{print $1}'`)
		set list = (`egrep "^\s*${i} " $basis_set_file | grep GTH | grep DZ`)
		if ($#list == 0) then
			set valid = 0
			break
		endif
	end
	if ($valid) break
end
if !($valid) then
	echo "ERROR: Cannot find DZVP-GTH entry for $i in $basis_set_file"
	exit(1)
endif
foreach j (~tpascal/codes/cp2k/cp2k/data/GTH_POTENTIALS ~tpascal/codes/cp2k/cp2k/data/POTENTIAL)
	set valid = 1
	set pseudo_file = $j
	foreach i (`cat ${prefix}${repstr}.cp2k.xyz | awk '{print $1}'`)
		set list = (`egrep "^\s*${i} " $pseudo_file | grep GTH | grep PBE`)
		if ($#list == 0) then
			set valid = 0
			break
		endif
	end
	if ($valid) break
end
if !($valid) then
	echo "ERROR: Cannot find GTH-PBE entry for $i in $pseudo_file"
	exit(1)
endif
foreach i (`cat ${prefix}${repstr}.cp2k.xyz | awk '{print $1}'`)
	if (-e ${prefix}.kind.${i}.dat) continue
	set basis_set = (`egrep "^\s*${i} " $basis_set_file | grep GTH | grep DZ`)
	set pseudo = (`egrep "^\s*${i} " $pseudo_file | grep PBE | grep GTH`)
	cat > ${prefix}.kind.${i}.dat <<DATA;
&KIND $i
 BASIS_SET $basis_set[2]
 POTENTIAL $pseudo[2]
&END KIND
DATA

	cat ${prefix}.kind.dat ${prefix}.kind.${i}.dat > __tmp.dat
	mv __tmp.dat ${prefix}.kind.dat
end

set basis_set_file = `echo $basis_set_file | sed 's/.*\///'`
set pseudo_file = `echo $pseudo_file | sed 's/.*\///'`
sed -i "/KIND_HERE/r ${prefix}.kind.dat" ${prefix}${repstr}.${temp}K.cp2k.in
sed -i '/KIND_HERE/d' ${prefix}${repstr}.${temp}K.cp2k.in
sed -i "s/GTH_BASIS_SETS/$basis_set_file/" ${prefix}${repstr}.${temp}K.cp2k.in
sed -i "s/GTH_POTENTIALS/$pseudo_file/" ${prefix}${repstr}.${temp}K.cp2k.in
cp $pbs_template_file ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i 's/PBS_O_WORKDIR/SLURM_SUBMIT_DIR/' ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
set shirley_input_block = "~tpascal/scripts/Input_Block.in"
cp $shirley_input_block ${prefix}${repstr}.Input_Block.in
sed -i "s/nodes_here/$nodes/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i "s/nprocs_here/$nprocs/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i "s/ppn_here/$ppn/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i "s/rtemp_here/$temp/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i "s/fprefix_here/$prefix/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i "s/cell_p_here/$repstr/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i "s/c_a_here/$cell[1]/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i "s/c_b_here/$cell[2]/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i "s/c_c_here/$cell[3]/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i "s/c_alpha_here/$cell[4]/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i "s/c_beta_here/$cell[5]/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
sed -i "s/c_gamma_here/$cell[6]/" ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
if ($repstr == "") then
	sed -i 's/\.\${cell_p}//g' ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
	sed -i 's/\${cell_p}//g' ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
	sed -i '/^cell_p=/d' ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript
endif

if (`echo $pbs_template_file | grep -c '\.sp\.'`) then
	mv ${prefix}${repstr}.${temp}K.cp2k.shirley.qscript ${prefix}${repstr}.sp.cp2k.shirley.qscript
	mv ${prefix}${repstr}.${temp}K.cp2k.in ${prefix}${repstr}.sp.cp2k.in
endif

rm -fr __tmp.dat ${prefix}.kind*.dat 

exit:
exit(0)

error:
echo "ERROR occurred"
exit(1)

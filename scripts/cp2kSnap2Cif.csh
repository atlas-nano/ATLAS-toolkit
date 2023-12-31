#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
  echo "usage: $0 prefix snapshot [save_name]"
  exit(1)
endif

set prefix = $1
set snap = $2
set save_name = `basename $1`.snap${2}.cif
if ($#argv > 2) set save_name = $3

if !(-e ${prefix}.out) then
  echo "ERROR: Cannot locate output file ${prefix}.out"
  exit(1)
endif
if !(-e ${prefix}-pos-1.xyz) then
  echo "ERROR: Cannot locate XYZ file ${prefix}-pos-1.xyz"
  exit(1)
endif
grep "^ CELL LNTHS" ${prefix}.out | awk '{print NR,$4*.529,$5*.529,$6*.529}' > __cell_length.dat
set tot_cell = `wc -l __cell_length.dat | awk '{print $1}'`
grep -n " i = " ${prefix}-pos-1.xyz | sed 's/://' | awk '{print $1-2}' > __xyz_pos.dat
set tot_xyz = `wc -l __xyz_pos.dat | awk '{print $1}'`
if !(`echo $snap $tot_cell $tot_xyz | awk '{if($1>$2||$1>$3) print 0; else print 1}'`) then
	echo "ERROR: Requested snapshot $snap exceeds either number of xyz snap $tot_xyz or number of cell parameters $tot_cell"
	exit(1)
endif

set xyz_start = `head -${snap} __xyz_pos.dat | tail -1`
if(`echo $snap $tot_xyz | awk '{if($1==$2)print 1; else print 0}'`) then
	set xyz_stop = `wc -l ${prefix}-pos-1.xyz | awk '{print $1}'`
else
	@ i = $snap + 1
	set xyz_stop = `head -${i} __xyz_pos.dat | tail -1`
endif
@ xyz_offset = $xyz_stop - $xyz_start + 2

echo "start $xyz_start stop $xyz_stop offset $xyz_offset"

head -${xyz_stop} ${prefix}-pos-1.xyz | tail -${xyz_offset} > __tmp.xyz
set cell = (`cat __cell_length.dat | awk '{if($1=='$snap') print $2,$3,$4}'`)
grep 'CELL_TOP| Angle' ${prefix}.out | awk '{printf "%s ",$6; if(NR%3==0)printf "\n"}' > __cell_angle.dat
set tot_cell_angle = `wc -l __cell_angle.dat | awk '{print $1}'`
if ($tot_cell_angle == 1) then
	set cell_angle = (`cat __cell_angle.dat`)
else
	set cell_angle = (`cat __cell_angle.dat | awk '{if(NR=='$snap') print}'`)
endif

cat > $save_name <<DATA;
data_block_1
_audit_creation_data            '`date`'
_audit_creation_method          'generated by ~tpascal/scripts/cp2kSnap2cif.csh'

loop_
_cell_length_a $cell[1]
_cell_length_b $cell[2]
_cell_length_c $cell[3]
_cell_angle_alpha $cell_angle[1]
_cell_angle_beta $cell_angle[1]
_cell_angle_gamma $cell_angle[1]

loop_
_atom_site_label
_atom_site_cartn_x
_atom_site_cartn_y
_atom_site_cartn_z
DATA

cat __tmp.xyz | awk '{printf "%2s %10.5f %10.5f %10.5f\n",$1,$2,$3,$4}' >> ${save_name}
echo "copy ${save_name} ${save_name}" > __gdis_cmd
~tpascal/codes/bin/gdis_cmdline < __gdis_cmd > /dev/null
rm -fr __tmp.xyz __xyz_pos.dat __cell_length.dat __cell_angle.dat __gdis_cmd


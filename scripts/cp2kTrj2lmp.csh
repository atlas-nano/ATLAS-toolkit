#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
  echo "usage: $0 cp2k_output_fike cp2k_pos_file [save_name] [cp2k_vel_file] [cp2k_mulliken_file]"
  exit(1)
endif

set output_file = $1
if !(-e $output_file) then
  echo "ERROR: Cannot find $1"
  exit(1)
endif

set pos_file = $2
if !(-e $pos_file) then
  echo "ERROR: Cannot find $2"
  exit(1)
endif

set save_name = `basename $pos_file`
set save_name = $save_name:r
set save_name = "${save_name}.lammps"
if ($#argv > 2) set save_name = $3

set vel_file = ""
if ($#argv > 3 && -e $4) set vel_file = $4

set charge_file = ""
#if ($#argv > 4 && -e $5) set charge_file = $5

set num_entries_output = `grep -c "^ CELL LNTHS" $output_file`
set num_entries_trj = `grep -c "^ i =" $pos_file`
set num_entries_vel = 99999999999999
set num_entried_chg = 99999999999999
echo "nentries_output: $num_entries_output num_entries_trj: $num_entries_trj"
if (${?vel_file} > 0) set num_entries_vel = `grep -c "^ i =" $vel_file`
#if (${?charge_file} > 0) set num_entries_chg = `grep -c "^ MULLIKEN POPULATION ANALYSIS" $charge_file`
echo "HERE"

set nsnaps = 9999999999999999999
foreach i ($num_entries_output $num_entries_trj $num_entries_vel $num_entries_chg)
  set nsnaps = `echo $nsnaps $i | awk '{if($1>$2) print $2; else print $1;}'`
end

set nsnaps = 23614
set num_entries_output = 23614
echo "CP2K Options"
echo "============="
echo "OUTPUT $output_file $num_entries_output" | awk '{printf "%10sFILE: %s (%d snapshots)\n",$1,$2,$3}'
echo "COORDINATE $pos_file $num_entries_trj" | awk '{printf "%10sFILE: %s (%d snapshots)\n",$1,$2,$3}'
if (${?vel_file} > 0) echo "VELOCITY $vel_file $num_entries_vel" | awk '{printf "%10sFILE: %s (%d snapshots)\n",$1,$2,$3}'
if (${?charge_file} > 0) echo "CHARGE $charge_file $num_entries_vel" | awk '{printf "%10sFILE: %s (%d snapshots)\n",$1,$2,$3}'
echo "SAVE_TRJ $save_name" | awk '{printf "%10sFILE: %s\n",$1,$2}'

echo "using $nsnaps snapshots"
#set offset_out = (`grep -n "^ CELL LNTHS" $output_file | sed 's/:.*//' | awk '{print $1}'`)
set offset_pos = (`grep -n "^ i =" $pos_file | sed 's/:.*$//'`)
if (${?vel_file} > 0) set offset_vel = (`grep -n "^ i =" $vel_file | sed 's/:.*$//' | awk '{print $1}'`)
if (${?charge_file} > 0) set offset_chg = (`grep -n "^ MULLIKEN POPULATION ANALYSIS" $charge_file | sed 's/:.*//' | awk '{print $1}'`)

set lmp_atom_str = "id type xu yu zu"
if (${?vel_file} > 0) set lmp_atom_str = "$lmp_atom_str vx vy vz"
if (${?charge_file} > 0) set lmp_atom_str = "$lmp_atom_str q"
echo > $save_name

#set offset_out = ($offset_out `wc -l $output_file | awk '{print $1}'`)
#set offset_pos = ($offset_pos `wc -l $pos_file | awk '{print $1+2}'`)
#if (${?vel_file} > 0) set offset_vel = ($offset_vel `wc -l $vel_file |  awk '{print $1+2}'`)
#if (${?charge_file} > 0) set offset_chg = ($offset_chg `wc -l $charge_file |  awk '{print $1+5}'`)

foreach i (`seq 1 $nsnaps`)
  set per = `echo $i $nsnaps | awk '{printf "%.2f\n",$1*100/$2}'`
  /bin/echo -en "${i}/${nsnaps} (${per}%)\r"
  set tstep = `head -$offset_pos[$i] $pos_file | tail -1 | sed 's/,//' | awk '{print $3}'`
  @ o = $offset_pos[$i] - 1
  set natom = `head -${o} $pos_file | tail -1`
  #set cell = (`head -$offset_out[$i] $output_file | tail -1 | awk '{print $4*.529,$5*.529,$6*.529}'`)
  set cell = (11.17480 9.70490 48.28836)
  cat >> $save_name <<DATA
ITEM: TIMESTEP
$tstep
ITEM: NUMBER OF ATOMS
$natom
ITEM: BOX BOUNDS x y z
0 $cell[1]
0 $cell[2]
0 $cell[3]
ITEM: ATOMS $lmp_atom_str
DATA

  @ j = $i + 1
  @ o = $offset_pos[$j] - 2
  head -${o} $pos_file | tail -${natom} | awk '{print NR,1,$2,$3,$4}' > __pos.dat
  set paste_str = "__pos.dat"
  if (${?vel_file} > 0) then
    @ o = $offset_vel[$j] - 2
    head -${o} $vel_file | tail -${natom} | awk '{print $2,$3,$4}' > __vel.dat
    set paste_str = "${paste_str} __vel.dat"
  endif
  if (${?charge_file} > 0) then
    @ o = $offset_chg[$j] - 5
    head -${o} $charge_file | tail -${natom} | awk '{print $5}' > __chg.dat
    set paste_str = "${paste_str} __chg.dat"
  endif
  paste $paste_str > __atom.dat
  cat __atom.dat | awk '{for(i=1;i<=NF;i++) printf "%s ", $i; printf "\n"}' >> $save_name
  rm -fr $paste_str __atom.dat
end
echo "${nsnaps}/${nsnaps} (100%)...Done"

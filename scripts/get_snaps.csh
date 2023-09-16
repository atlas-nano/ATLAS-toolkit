#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
  echo "usage: $0 prefix [num_snaps=3] [save_prefix]"
  exit(1)
endif

set prefix = $1
set num_snaps = 3
if ($#argv > 1) set num_snaps = $2
set save_prefix = `basename $prefix`
set save_prefix = ${save_prefix}.equil
if ($#argv > 2) set save_prefix = $3

if !(-e ${prefix}.out) then
  echo "ERROR: Cannot locate output file ${prefix}.out"
  exit(1)
endif
if !(-e ${prefix}-pos-1.xyz) then
  echo "ERROR: Cannot locate XYZ file ${prefix}-pos-1.xyz"
  exit(1)
endif
grep "^ CELL LNTHS" ${prefix}.out | awk '{print NR,$4*.529,$5*.529,$6*.529}' > ${save_prefix}.cell_length.dat
set tot = `wc -l ${save_prefix}.cell_length.dat | awk '{print $1}'`
set half = `echo $tot | awk '{print int($1/2);}'`
set stats = (`tail -${half} ${save_prefix}.cell_length.dat | awk '{print $1,$2}' | awk -f ~/scripts/computeTotals.awk | grep "^Mean" | awk '{print $3,$7}'`)
set stats = ($stats `tail -${half} ${save_prefix}.cell_length.dat | awk '{print $1,$3}' | awk -f ~/scripts/computeTotals.awk | grep "^Mean" | awk '{print $3,$7}'`)
set stats = ($stats `tail -${half} ${save_prefix}.cell_length.dat | awk '{print $1,$4}' | awk -f ~/scripts/computeTotals.awk | grep "^Mean" | awk '{print $3,$7}'`)
grep "^ CELL LNTHS" ${prefix}.out | awk '{print ($4*.529-'$stats[1]'+$5*.529-'$stats[3]'+$6*.529-'$stats[5]')/3,NR,$4*.529,$5*.529,$6*.529}' | sed 's/-//' | sort -n | awk '{print $2,$3,$4,$5,$1}' > ${save_prefix}.cell_length.stats.dat
set found = 0
set file_n = 1
foreach i (`cat ${save_prefix}.cell_length.stats.dat | awk '{print $5}'`)
  set count = `egrep -c " $i" ${save_prefix}.cell_length.stats.dat`
  if ($count >= $num_snaps) then
    set found = 1
    set cell_str = `echo $i`
    break
  endif
  @ file_n = $file_n + 1
end

if !($found) then
  echo "ERROR: Could not find $num_snaps snapshots with same cell length"
  exit(1)
endif

set cell_length = (`grep " $cell_str" ${save_prefix}.cell_length.stats.dat | head -1 | awk '{print $2,$3,$4}'`)
set cell_angles = (`grep "^ CELL ANGLS" ${prefix}.out | awk '{ if (NR=='$file_n') printf "%f %f %f\n",$5,$6,$7; }'`)
echo "$count $cell_length $cell_angles" 
egrep " $cell_str" ${save_prefix}.cell_length.stats.dat | awk '{print $1}' > ${save_prefix}.snaps.dat
set tot = `wc -l ${prefix}-pos-1.xyz | awk '{print $1}'`
grep -n "i =" ${prefix}-pos-1.xyz | awk '{print $4,$1}' | sed 's/://' | sed 's/,//' | awk '{print $1,$2,'$tot'-$2}' > ${prefix}.xyz_offset.dat
set csize = `tail -1  ${prefix}.xyz_offset.dat | awk '{print $3+2}'`
foreach i (`head -${num_snaps} ${save_prefix}.snaps.dat`)
  set offset = `grep "^$i " ${prefix}.xyz_offset.dat | awk '{print $3+2}'`
  tail -${offset} ${prefix}-pos-1.xyz | head -${csize} > ${save_prefix}.snap${i}.xyz
  sed -i '2 s/$/ cell = '"${cell_length} ${cell_angles}"'/' ${save_prefix}.snap${i}.xyz
end

#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
  echo "usage: $0 'data_file1 data_file2 ...' ('column(s)') (save_prefix) (search_column=1)"
  exit(1)
endif

set data_files = ($1)
if ($#data_files < 2) then
  echo "ERROR: Need at least 2 data file"
  exit(1)
endif

set tot_cols = `tail -1 $data_files[1] | awk '{print NF}'`
set cols = (`seq 2 $tot_cols`)
if ($#argv > 1) set cols = ($2)

set save_prefix = `basename $data_files[1]`
set save_prefix = $save_prefix:r
if ($#argv > 2) set save_prefix = $3

set search_col = 1
if ($#argv > 3) set search_col = $4


set tot_lines = `wc -l $data_files[1] | awk '{print $1}'`
echo "columns: $cols, prefix: $save_prefix, search_column: $search_col tot_line: $tot_lines"
foreach i ($cols)
  if (`echo $i | grep -io '[a-z]' | wc -l`) then
    echo "ERROR: Expected integer for cols. Got $i"
    exit(1)
  endif
  rm -fr __${save_prefix}.tmp.joined.${i}.avg.dat __${save_prefix}.tmp.joined.${i}.std.dat
end
set avg_files = ()
set std_files = ()
rm -fr __${save_prefix}.tmp.index.dat
foreach i (`seq 1 $tot_lines`)
  if (`head -${i} $data_files[1] | tail -1 | grep -io '[a-z]' | wc -l`) continue
  if !(`head -${i} $data_files[1] | tail -1 | grep -io '[1-9]' | wc -l`) continue
  set index = `head -${i} $data_files[1] | tail -1 | awk '{print $'$search_col'}'`
  set per = `echo $i $tot_lines | awk '{printf "%.2f\n",$1*100/$2}'`
  /bin/echo -en "${i}/${tot_lines} (${per}%)\r"
  echo $index >> __${save_prefix}.tmp.index.dat
  foreach j ($cols)
    set vals = (`egrep "^[[:blank:]]*$index " $1 | awk '{print $1,$'$j'}' | awk -f ~/scripts/computeTotals.awk | grep "^Mean" | awk '{print $3,$7}'`)
    echo $vals[1] >> __${save_prefix}.tmp.joined.${j}.avg.dat
    echo $vals[2] >> __${save_prefix}.tmp.joined.${j}.std.dat
  end
end
/bin/echo "${i}/${tot_lines} (${per}%)"

paste __${save_prefix}.tmp.index.dat __${save_prefix}.tmp.joined.*.avg.dat > ${save_prefix}.avg.dat
paste __${save_prefix}.tmp.index.dat __${save_prefix}.tmp.joined.*.std.dat > ${save_prefix}.std.dat
rm -fr __${save_prefix}.tmp.index.dat 
rm -fr __${save_prefix}.tmp.joined.*.avg.dat
rm -fr __${save_prefix}.tmp.joined.*.std.dat

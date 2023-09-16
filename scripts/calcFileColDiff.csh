#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
  echo "usage: $0 'data_file1 data_file2' (column(s)=all) (save_name) (search_column=1) (scale_1 scale_2)"
  exit(1)
endif

set data_files = ($1)
if ($#data_files <> 2) then
  echo "ERROR: Need exactly 2 data file"
  exit(1)
endif

set tot_cols = `tail -1 $data_files[1] | wc -l`
set cols = (`seq 2 $tot_cols`)
if ($#argv > 1) set cols = ($2)

set save_name = `basename $data_files[1]`
set save_name = $save_name:r
set save_name = "${save_name}.diff.dat"
if ($#argv > 2) set save_name = $3

set search_col = 1
if ($#argv > 3) set search_col = $4

set scale = (1 1)
if ($#argv > 4) set scale = ($5)
if ($#scale < 2) then
  echo "Expected 2 scaling factors"
  exit(1)
endif

foreach i ($cols)
  if (`echo $i | grep -io '[a-z]' | wc -l`) then
    echo "ERROR: Expected integer for cols. Got $i"
    exit(1)
  endif
end
rm -fr __${save_name}.tmp.index.dat
set tot_lines = `wc -l $data_files[1] | awk '{print $1}'`
set c = 1
foreach i (`seq 1 $tot_lines`)
  if (`head -${i} $data_files[1] | tail -1 | grep -io '[a-z]' | wc -l`) continue
  if !(`head -${i} $data_files[1] | tail -1 | grep -io '[1-9]' | wc -l`) continue
  set index = `head -${i} $data_files[1] | tail -1 | awk '{print $'$search_col'}'`
  set l2 = `cat $data_files[2] | awk '{ if($'$search_col'=='$index') print NR;}'`
  if !(${?l2}) continue
  set per = `echo $i $tot_lines | awk '{printf "%.2f\n",$1*100/$2}'`
  /bin/echo -en "${i}/${tot_lines} (${per}%)\r"
  echo $index >> __${save_name}.tmp.index.dat
  foreach j ($cols)
    set val1 = `head -${i} $data_files[1]  | tail -1 | awk '{print $'$j'/'$scale[1]'}'`
    set val2 = `head -${l2} $data_files[2] | tail -1 | awk '{print $'$j'/'$scale[2]'}'`
    echo $val1 $val2 | awk '{ print $1-$2}' >> __${save_name}.tmp.${j}.diff.dat
  end
end
/bin/echo "${tot_lines}/${tot_lines} (100.00%)"

paste __${save_name}.tmp.index.dat __${save_name}.tmp.*.diff.dat > $save_name
rm -fr __${save_name}.tmp.index.dat __${save_name}.tmp.*.diff.dat

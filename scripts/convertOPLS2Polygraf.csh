#!/bin/tcsh
#!/bin/csh

if ($#argv < 3) then
  echo "usage: bgf_file(s) forcefield_dir rot_symmetry_file (machine)"
  exit(1)
endif

set bgfs = ($1)
if ($#bgfs == 0) then
    echo "ERROR: No valid bgf file found in $1"
    exit(1)
endif

set ff_dir = $2
if !(-e $ff_dir) then
    echo "ERROR: Cannot locate $ff_dir"
    exit(1)
endif

set rSymm = $3
if !(-e $rSymm) then
    echo "ERROR: Cannot locate $rSymm"
    exit(1)
endif

set machine = ""
if ($#argv > 3) then
    set machine = $4
endif

foreach b ($bgfs)
  set mol = `basename $b .bgf`
  set ff = ${ff_dir}/${mol}.ff
  if !(-e $ff ) continue 
  set name = $mol
  if !(`grep -c "^$name " $rSymm`) then
    set name = $name:r
    if !(`grep -c "^$name " $rSymm`) continue
  endif
  set rotSym = `grep "^$name " $rSymm | awk '{print $2}'`
  echo $mol
  ~/scripts/removeBGFCellInfo.pl -b $b -s ${mol}.bgf > /dev/null || continue
  ~/scripts/opls2polygrafBGF.pl -b ${mol}.bgf > /dev/null || continue
  ~/scripts/ceriusff2biograf.pl -f $ff -s ${mol}.polygraf.par -n 1 -b ${mol}.bgf -p 1> /dev/null || continue
  ~/scripts/createPolygrafThermoMacro.pl -b ${mol}.polygraf.bgf -r $rotSym -s ${mol}.polygraf.macro > /dev/null || continue
  if($#argv == 3) continue
  set c_host = $HOST
  set c_dir = $PWD
  tar -cof ss.tar ${mol}.polygraf.*
  scp ss.tar ${machine}:/ul/tpascal/tmp/
  ssh $machine 'mkdir -p ~/tmp/_polygraf_'${mol}';cd ~/tmp/_polygraf_'${mol}'; mv ~/tmp/ss.tar ./; tar -xof ss.tar; bgver 400 ; echo "/biograf/bio_msc batpoly '${mol}'.polygraf '${mol}'.polygraf.macro "> ~/'${mol}'.polygraf.out; rm -fr ss.tar; cd ../; rm -fr ~/tmp/_polygraf_'${mol}
  rm -fr ss.tar
end

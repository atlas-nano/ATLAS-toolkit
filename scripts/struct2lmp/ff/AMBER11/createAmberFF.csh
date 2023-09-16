#!/bin/tcsh
#!/bin/csh

foreach f ( /usr/local/amber11/dat/leap/parm/parm* /usr/local/amber11/dat/leap/parm/frcmod.* )
  set mol = `basename $f .dat`
  set prefix = $mol:r
  if ($prefix == frcmod) set mol = $mol:e
  echo "$f -> ${mol}.AMBER.ff"
  ~/scripts/amberff2cerius.pl -f $f -s ${mol}.AMBER > /dev/null || continue
  cat ../AMBERFF_header.ff ${mol}.AMBER > ${mol}.AMBER.ff
  rm -fr ${mol}.AMBER
end

#!/bin/tcsh
#!/bin/csh

if ($#argv < 3) then
  echo "usage: $0 bgf_file force_field required_density [save_name]"
  exit(1)
endif

set bgf = $1
set ff = $2
set density = $3
set save_name = `basename $bgf`
set save_name = $save_name:r
set save_name = "${save_name}.${density}_g_cc.bgf"
if ($#argv > 3) set save_name = $4

echo "Creating $save_name with box to make density $density"
/home/tpascal/scripts/removeBGFBox.pl -b $bgf -s __tmp.bgf > /dev/null
/home/tpascal/scripts/addBoxToBGF.pl __tmp.bgf __tmp.bgf > /dev/null || goto error
set cell = (`grep CRYSTX __tmp.bgf | tail -1 | awk '{print $2,$3,$4}'`)
set cell_vol = `echo $cell | awk '{ print $1*$2*$3}'`
set mass = `/home/tpascal/scripts/bgfmass.pl -b $bgf -f $ff | grep "Total Mass" | awk '{i = NF-1; print $i}'`
set vol = `echo $mass $density | awk '{print ($1/0.6023/$2)'}`
set scale_f = `echo $vol $cell_vol | awk '{print ($1/$2)**(1/3)}'`
set new_cell = (`echo $scale_f $cell | awk '{ for(i=2;i<=NF;i++) print $1*$i; }'`)
/home/tpascal/scripts/updateBGFBox.pl -b __tmp.bgf -s $save_name -c "$new_cell" > /dev/null || goto error

exit:
rm -fr __tmp.bgf
echo "Done"
exit(1)

error:
echo "Error occurred"
exit(1)

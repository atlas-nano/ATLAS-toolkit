#!/bin/tcsh
#!/bin/csh

if ($#argv < 3) then
    echo "usage: $0 jaguar_output_file jag_jag_resp [save_prefix] [charge=0] [multipliticy=1] [fftype=gaff|amber]"
    exit(1)
endif

set jag_out = $1
if (! -e $jag_out) then
    echo "Error: Cannot read jaguar output file $jag_out"
    exit(1)
endif

set jag_resp = $2
if (! -e $jag_resp) then
    echo "Error: Cannot read jaguar RESP file $jag_resp"
    exit(1)
endif

set sprefix = $3

set charge = 0
if ($#argv > 3) then
    set charge = $4
endif

set multiplicity=0
if ($#argv > 4) then
    set multiplicity = $5
endif

set resname = `echo $sprefix | cut -c1-3 | tr '[a-z]' '[A-Z]'`
set lc_res = `echo $resname | tr '[A-Z]' '[a-z]'`

set fftype = "gaff"
if ($#argv > 5) then
    set fftype = $6
    if !($fftype == "amber" || $fftype == "gaff") then
	echo "Error: Expected amber|gaff for forcefield. Got \"${fftype}\""
	exit(1)
    endif
    set sprefix = "${sprefix}.${fftype}"
endif

cat <<DATA;
jag out file:        $jag_out
jag resp file:       $jag_resp
molecular charge:    $charge
multiplicity (2s+1): $multiplicity
prepi resname:       $resname
prepi resid:         $lc_res 
save_prefix:         $sprefix
DATA

echo "Step 1. Creating .ac file"
${AMBERHOME}/bin/antechamber -i $jag_out -fi jout -o _tmp.ac -fo ac -nc $charge -m ${multiplicity} -rf ${resname} -rf ${lc_res} -at ${fftype} -pf y -dr no >& /dev/null || goto error
echo "Step 2. Creating Respgen files"
${AMBERHOME}/bin/respgen -i _tmp.ac -o _tmp.respin1 -f resp1 > /dev/null || goto error
${AMBERHOME}/bin/respgen -i _tmp.ac -o _tmp.respin2 -f resp2 > /dev/null || goto error
echo "Step 3. Resp stage 1"
cat $jag_resp | awk '{if(NR==1) printf "%5d %s\n",$1,$2; else print}' > __tmp.resp
${AMBERHOME}/bin/resp -O -i _tmp.respin1 -o _tmp.respout1 -e __tmp.resp -t _tmp_qout_stage1  > /dev/null || goto error
echo "Step 4. Resp stage 2"
${AMBERHOME}/bin/resp -O -i _tmp.respin2 -o _tmp.respout2 -e __tmp.resp -q _tmp_qout_stage1 -t _tmp_qout_stage2 > /dev/null || goto error
echo "Step 5. Antechamber with correct charges"
${AMBERHOME}/bin/antechamber -i _tmp.ac -fi ac -o ${sprefix}.prepi -fo prepi -c rc -cf _tmp_qout_stage2 -pf y -s 2 -at ${fftype} -rf ${resname} -rn ${lc_res} > /dev/null || goto error
${AMBERHOME}/bin/antechamber -i _tmp.ac -fi ac -o ${sprefix}.resp.mol2 -fo mol2 -c rc -cf _tmp_qout_stage2 -pf y -s 2 -at ${fftype} -rf ${resname} -rn ${lc_res} > /dev/null || goto error
/home/tpascal/scripts/mol2bgf.pl -m ${sprefix}.resp.mol2 > /dev/null || goto error
mv ${sprefix}.resp.bgf ${sprefix}.gaff.resp.bgf

echo "Stap 6. ESP/Mulliken charges with Dreiding Atomtypes"
/home/tpascal/scripts/opls2polygrafBGF.pl  -b ${sprefix}.gaff.resp.bgf -s ${sprefix}.resp.bgf > /dev/null
/home/tpascal/scripts/autoType.pl -i ${sprefix}.resp.bgf -s ${sprefix}.resp.bgf -f DREIDING > /dev/null
/home/tpascal/scripts/jag2bgf.pl -b ${sprefix}.resp.bgf -j ${jag_out} -f CHARGE -c mul -s ${sprefix}.mul.bgf > /dev/null 
/home/tpascal/scripts/jag2bgf.pl -b ${sprefix}.resp.bgf -j ${jag_out} -f CHARGE -c esp -s ${sprefix}.esp.bgf > /dev/null

rm -fr _tmp.ac _tmp.respin1 _tmp.respin2 _tmp.respout1 _tmp_qout_stage1 _tmp.respout2 _tmp_qout_stage2 punch esout _tmp.mol2 __tmp.resp
exit:
echo "All tasks completed"
exit(0)

error:
echo "Error occurred"
exit(1)

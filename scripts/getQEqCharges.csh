#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 bgf_file [savename]"
	exit(1)
endif

if !(-e $1) then
	echo "ERROR: Cannot locate $1"
	exit(1)
endif

set savename = `basename $1`
set savename = $savename:r
set savename = "${savename}.qeq.bgf"
if ($#argv > 1) set savename = $2

echo "Step 1: Calculating QEq charges"
/home/tpascal/scripts/modifyAtomData.pl -s $1 -w __test.bgf -a "index>0" -f "CHARGE:0" > /dev/null || goto error
sed -i '/^DISP/d' __test.bgf
sed -i '/^ORDER/d' __test.bgf
/home/tpascal/codes/bin/pqeq /home/tpascal/ff/pqeq_1_1.ff 100 __test.bgf __test.bgf >/dev/null || goto error
echo "Step 2: Saving charges to $savename"
sed -i '/^END/d' __test.bgf
sed "1,/^FORMAT CONECT/d" $1 > __bonds.dat
echo "FORMAT CONECT (a6,12i6)" >> __test.bgf
cat __bonds.dat >> __test.bgf
/home/tpascal/scripts/updateBGF.pl -b $1 -r __test.bgf -t "CHARGE" -s ${savename} > /dev/null || goto error
rm -fr __test.bgf __bonds.dat
echo "Done"
exit(0)

error:
echo "ERROR occurred"
exit(1)

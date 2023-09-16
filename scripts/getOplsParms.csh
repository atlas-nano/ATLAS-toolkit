#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
    echo "usage: $0 structure_file(s)|directory"
    exit(1)
endif

set files = `ls $1 | egrep '\.(bgf|mol2|xyz|pdb|cif)'`
if ($#files == 0) then
    echo "ERROR: No valid bgf found while searching $1"
    exit(1)
endif
echo "files $files"

set headers = (TotEng Bonds Angles Torsions Impropers VDW Elec)
set mm = (0 0 0 0 0 0 0 0)
set lmp = (0 0 0 0 0 0 0 0)
set diff = (0 0 0 0 0 0 0 0)
set diffp = (0 0 0 0 0 0 0 0)

set curr = $PWD

mkdir -p bgfs ffs lammps macromodel_data

foreach i ($files)
    cd $curr
    set suffix = `echo $i | awk -F . '{print $NF}'`
    set mol = `basename $i .${suffix}`
    set mol = "${mol}.opls"
    echo "Entered calculation for $mol"
    echo "	1. Assigning OPLS charges/atom types and creating FF...";
    if ( $suffix != "bgf") then
        $SCHRODINGER/utilities/babel -i${suffix} $i -obgf bgfs/__test.bgf >& /dev/null
    else
        cp $i bgfs/__test.bgf
    endif
    cd bgfs
    /ul/tpascal/scripts/removeBGFBox.pl -b __test.bgf -s __test.bgf > /dev/null
    set natoms = `egrep -c '^(HETATM|ATOM) ' __test.bgf`
    /ul/tpascal/scripts/addBoxToBGF.pl __test.bgf __test.bgf > /dev/null
    awk '{if($1 ~ /CRYSTX/) printf "CRYSTX %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n",$2+10,$3+10,$4+10,$5,$6,$7; else print}' __test.bgf > __tmp.bgf
    mv __tmp.bgf __test.bgf
    /ul/tpascal/scripts/replicate.pl -b __test.bgf -s __test.bgf -d "1 1 2" > /dev/null
    /ul/tpascal/scripts/removeBGFBox.pl -b __test.bgf -s __test.bgf > /dev/null
    cd ../
    cd macromodel_data/
    /ul/tpascal/scripts/opls2cerius.pl -i ../bgfs/__test.bgf -t bgf -s ${mol}.ff > ${mol}.opls2cerius.screen.out || goto error
    mv ${mol}.ff ../ffs/
    mv ${mol}.bgf ../bgfs/
    cd ../
    rm -fr bgfs/__test.bgf
    echo "	2. Calculating energy with LAMMPS..."
    cd lammps
	/ul/tpascal/scripts/createLammpsInput.pl -b ../bgfs/${mol}.bgf -f ../ffs/${mol}.ff -s ${mol} > ${mol}.createLammpsTrj.screen.out || continue
	sed -i 's/boundary        p p p/boundary        f f f/' in.${mol}_singlepoint
	sed -i 's/pair_style\s\+\(lj\S\+\).*/pair_style      lj\/cut\/coul\/cut 48.0 48.0/' in.${mol}_singlepoint
	sed -i 's/kspace_style.*/kspace_style    none/' in.${mol}_singlepoint
	sed -i 's/.*xhi/-140 140 xlo xhi/' data.${mol}
	sed -i 's/.*yhi/-140 140 ylo yhi/' data.${mol}
	sed -i 's/.*zhi/-140 140 zlo zhi/' data.${mol}
	~/codes/bin/lmp_serial -in in.${mol}_singlepoint -screen none -log ${mol}_singlepoint.log > /dev/null
    cd ../
    /ul/tpascal/scripts/getBGFAtoms.pl -b bgfs/${mol}.bgf -s  bgfs/${mol}.bgf -o "index<=$natoms" > /dev/null
    echo "	3. Comparing Energies..."
    #get macromodel energies in kcal/mol
    set mm[1] = `grep "    Total Energy =" macromodel_data/${mol}.log | awk '{printf "%14.5f", ($4/4.184)}'`
    set mm[2] = `grep "         Stretch =" macromodel_data/${mol}.log | awk '{printf "%14.5f", ($3/4.184)}'`
    set mm[3] = `grep "            Bend =" macromodel_data/${mol}.log | awk '{printf "%14.5f", ($3/4.184)}'`
    set mm[4] = `grep "         Torsion =" macromodel_data/${mol}.log | awk '{printf "%14.5f", ($3/4.184)}'`
    set mm[5] = `grep "Improper Torsion =" macromodel_data/${mol}.log | awk '{printf "%14.5f", ($4/4.184)}'`
    set mm[6] = `grep "             VDW =" macromodel_data/${mol}.log | awk '{printf "%14.5f", ($3/4.184)}'`
    set mm[7] = `grep "   Electrostatic =" macromodel_data/${mol}.log | awk '{printf "%14.5f", ($3/4.184)}'`
    #now get lammps energies
    set lmp[1] = `grep "^TotEng "  lammps/${mol}_singlepoint.log | awk '{print $3}'`
    set lmp[2] = `grep "^PotEng "  lammps/${mol}_singlepoint.log | awk '{print $6}'`
    set lmp[3] = `grep "^PotEng "  lammps/${mol}_singlepoint.log | awk '{print $9}'`
    set lmp[4] = `grep "^E_dihed " lammps/${mol}_singlepoint.log | awk '{print $3}'`
    set lmp[5] = `grep "^E_dihed " lammps/${mol}_singlepoint.log | awk '{print $6}'`
    set lmp[6] = `grep "^E_dihed " lammps/${mol}_singlepoint.log | awk '{print $9}'`
    set lmp[7] = `grep "^E_coul "  lammps/${mol}_singlepoint.log | awk '{print $3}'`
    #now calculate differences
    foreach i (`seq 1 7`)
	set diff[$i] = `echo $mm[$i] $lmp[$i] | awk '{printf "%14.5f", ($1-$2)}'`
	set diffp[$i] = `echo $mm[$i] $lmp[$i] | awk '{printf "%14.5f", 100*($1-$2)/($1+.00000000001)}'`
    end
    #print values
    echo "quantity MacroModel LAMMPS Diff Diff(%)" | awk '{printf "%10s%-10s %10s %10s %10s %10s\n","",$1,$2,$3,$4,$5}'
    foreach i (`seq 1 7`)
	set check_str = ""
	if($i == 1 && `echo $diffp[$i] | awk '{if($1> 0.01 || ($1<0 && $1<-0.01)) print 1; else print 0}'`) set check_str = "***"
	echo "$headers[$i] $mm[$i] $lmp[$i] $diff[$i] $diffp[$i] $check_str" | awk '{printf "%10s%-10s %10.5f %10.5f %10.5f %10.5f%3s\n","",$1,$2,$3,$4,$5,$6}'
    end
end

exit:
echo "All tasks completed"
exit(0)

error:
echo "Error occurred"
exit(1)

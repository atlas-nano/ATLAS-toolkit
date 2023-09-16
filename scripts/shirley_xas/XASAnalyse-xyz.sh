#!/bin/bash
# This script XASAnalyze.sh combines spectra from various atoms in various
# xyz-files with atomic alignment as average spectra by type of atom.
# Note that everything is controlled in the Input_Block.in - from the
# xyz-files considered to the atoms considered.
# You can run this script interactively without arguments.
# Since there is a parallel component for running xas_para.x, then you
# may need to request an interactive shell (qsub -I) and submit from there.
# All spectra are collected together in XAS/Spectra-{type}

# davegp, Dec 20, 2011, TMF

#SBATCH --job-name=test
#SBATCH --partition=vulcan 
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --time=01:00:00
#SBATCH --qos=normal
#SBATCH --account=vulcan


# change to working directory if PBS job
if [[ -n $PBS_O_WORKDIR ]]; then
  cd $PBS_O_WORKDIR
elif [[ -n $SLURM_JOBID ]] ; then
  export PBS_O_WORKDIR=$SLURM_SUBMIT_DIR
  export PBS_JOBID=$SLURM_JOBID
  cd $PBS_O_WORKDIR
else
  export PBS_O_WORKDIR=`pwd`
fi

. ./Input_Block.in

# spectral details
SIGMA=0.2   # Gaussian broadening
NENER=1000  # Number of energy points
ELOW=-10    # Lower limit for spectrum
EHIGH=25    # Upper limit for spectrum
NOPr=08    # Processors
TEMP=0.0    # Electron temp
STRETCH=1.0 #DFT Stretching Factor 

if [ -z $CRAY_SITE_LIST_DIR ]; then
  GLOBAL_PARA_PREFIX="$PARA_PREFIX -np $NOPr"
else
  GLOBAL_PARA_PREFIX="$PARA_PREFIX -n $NOPr"
fi
echo $GLOBAL_PARA_PREFIX

ExecDir="$SHIRLEY_ROOT/bin"
BibDir="$BIB_ROOT/Bibliothek"
xyz2inp="${BibDir}/../xyz2inp.sh"
MyDir=`pwd`

REF=atom
CALCGS=GS

ixyz=0
TYPES=''
for xyzfile in $XYZFILES
do

    if [[ -f $xyzfile ]]; then
	echo " working on xyz-file $xyzfile"
    else
	echo " unable to find xyz-file $xyzfile - stopping"
	exit
    fi

  # This defines the name of the directory where calculation results are stored
  # This should be linked to an xyz file
    CALC=`echo $xyzfile | sed 's/\.xyz$//'`

    if [[ ! -d $MyDir/XAS/$CALC ]]; then
	echo directory $MyDir/XAS/$CALC is missing - stopping
	exit
    fi
    
    ixyz=$((ixyz + 1))
    inpblk[$ixyz]=$MyDir/XAS/$CALC/Input_Block.in
    cp $MyDir/Input_Block.in ${inpblk[$ixyz]}
    $xyz2inp $xyzfile $XYZUNIT $XASELEMENTS >> ${inpblk[$ixyz]} &
    
done
wait
nxyz=$ixyz

echo

 SPIN_CALC=0 ;
# if [[ -n $SPIN  &&  -n $GS_MAG && $GS_MAG -gt 0 ]]  ; then
 if [[ -n $SPIN  &&  -n $GS_MAG  ]]  ; then
     SPIN_CALC=1;
 fi

TYPES=''
# now loop over excited atoms of this type
for ixyz in `seq 1 $nxyz`
do

  # source the right file
    . ${inpblk[$ixyz]}
    dir=`dirname ${inpblk[$ixyz]}`

    Atoms=""

  #find all excited atoms
    for i in `seq 1 $NAT` ; do
	IA=${ATOM_SYMBOL[$i]}
	if [ ${IND_EXCITATION[$i]} -eq 1 ]; then
	    Atoms="$Atoms $i"
	fi
    done
    
  #shift spectra
    for atom in $Atoms ; do 
	IA=${ATOM_SYMBOL[$atom]}
	elem=`echo $IA | sed 's/[0-9]*$//'`
	cd $MyDir
	
	echo "Reference atom calculations from directory:
  XAS/$REF/$elem"
	
	if [[ ! -f XAS/$REF/$elem/atom.$elem-GS.scf.out ]]; then
	    echo "missing atom GS file: XAS/$REF/$elem/atom.$elem-GS.scf.out"
	    exit
	fi
	ERefGS0=`grep '^!' XAS/$REF/$elem/atom.$elem-GS.scf.out | awk '{print $5}' | tail -1`
	
	if [[ ! -f XAS/$REF/$elem/atom.$elem-${CHAPPROX}.scf.out ]]; then
	    echo "missing atom ${CHAPPROX} file: XAS/$REF/$elem/atom.$elem-${CHAPPROX}.scf.out"
	    exit
	fi
	ERefES0=`grep '^!' XAS/$REF/$elem/atom.$elem-${CHAPPROX}.scf.out | awk '{print $5}' | tail -1`

	SpecDir=$MyDir/XAS/Spectrum-$elem
	if test ! -d $SpecDir ; then
	    mkdir $SpecDir
	fi
	cd $SpecDir
	
    # if this is a new type, create a new Spectra file and zero counter
	n=`echo "$TYPES" | grep -c "$elem"`
	if [[ $n == 0 ]]; then
	    :> Spectra
	    TYPES="$TYPES $elem"
	    echo TYPES="$TYPES"
	    
	    declare "NAVE_$elem"=0
	fi
	
	atomw=`seq -w $atom $NAT | head -1`

	if [ $SPIN_CALC -eq 1 ] ; then
	    for SP in UP DN ; do
		if [[ $SP == UP || ( $SP == DN && $GS_MAG -gt 0 ) ]] ; then
		XASPrefix=$dir/${ATOM_SYMBOL[$atom]}${atom}_${SP}/$MOLNAME.${ATOM_SYMBOL[$atom]}${atomw}-$CHAPPROX
	    
    #get the specific energies for shift
		if [[ ! -f ${XASPrefix}.scf.out ]]; then
		    echo "missing SCF output: ${XASPrefix}.scf.out"
		    exit
		fi
		EXCH=`grep '^!' ${XASPrefix}.scf.out | awk '{print $5}' | tail -1`

		if [[ ! -f $dir/$CALCGS/$MOLNAME.scf.out ]]; then
		    echo "missing GS output: $dir/$CALCGS/$MOLNAME.scf.out"
		    exit
		fi
		EGS=`grep '^!' $dir/$CALCGS/$MOLNAME.scf.out | awk '{print $5}' | tail -1`
		
    # find the Fermi level (and CBM)
#		EFERM=`grep "the Fermi energy is"  ${XASPrefix}.scf.out | awk '{print $5}' | tail -1`
		NELXCH=`echo "$NELEC+1" | bc -l` 
		TOMAG=`grep 'total magnetization       =' ${XASPrefix}.scf.out | tail -1 | awk '{print $4}'`
		echo TOMAG = $TOMAG
		$GLOBAL_PARA_PREFIX $ExecDir/efermi.x $PARA_POSTFIX $TEMP ${XASPrefix}.xas.$XAS_ARG $NELXCH $TOMAG > ${XASPrefix}.xas.${XAS_ARG}-fermi.out_XCH
		if [[  $SP  =~ "UP" ]]  ; then
		    EFERM=`grep 'VBM ='  ${XASPrefix}.xas.${XAS_ARG}-fermi.out_XCH | head -1 |  awk '{print $3}'`
		else
		    EFERM=`grep 'VBM ='  ${XASPrefix}.xas.${XAS_ARG}-fermi.out_XCH | tail -1 |  awk '{print $3}'`
		fi

		TOMAG=$GS_MAG
		$GLOBAL_PARA_PREFIX $ExecDir/efermi.x $PARA_POSTFIX $TEMP ${XASPrefix}.xas.$XAS_ARG $NELEC $TOMAG > ${XASPrefix}.xas.${XAS_ARG}-fermi.out
		
#		ELUMO=`tail ${XASPrefix}.xas.${XAS_ARG}-fermi.out | grep 'CBM =' | awk '{print $3}'`
		if [[  $SP  =~ "UP" ]]  ; then
		    ELUMO=`grep 'CBM ='  ${XASPrefix}.xas.${XAS_ARG}-fermi.out | head -1 |  awk '{print $3}'`
		else
		    ELUMO=`grep 'CBM ='  ${XASPrefix}.xas.${XAS_ARG}-fermi.out | tail -1 |  awk '{print $3}'`
		fi
		
		E=`echo | awk -v exch=$EXCH -v erefes0=$ERefES0 -v egs=$EGS -v erefgs0=$ERefGS0 -v elumo=$EFERM '{print (exch-erefes0-egs+erefgs0)*13.6056923 - elumo}'`
		
		echo "i : (XCH - XCH0) - (GS - GS0) - LUMO = Delta"
		echo $atom : $EXCH $ERefES0 - $EGS $ERefGS0 $EFERM = $E
		
    # check for previous .xas file and save
		for file in `ls $dir/${ATOM_SYMBOL[$atom]}${atom}_${SP}/*.xas 2> /dev/null`; do
		    mv $file ${file}-raw
		done
		
    # run xas_para.x
		$GLOBAL_PARA_PREFIX $ExecDir/xas_para.x $PARA_POSTFIX $ELOW $EHIGH $NENER $SIGMA $E $ELUMO $STRETCH ${XASPrefix}.xas.$XAS_ARG 2> /dev/null
	    
		XASFile="`echo $dir | sed 's/^.*\///'`-`basename ${XASPrefix}.xas.$XAS_ARG.xas`"
		XASFile="`echo ${XASFile}-${SP}`";
		echo "doing ${XASPrefix}.xas.$XAS_ARG.xas $XASFile "
		mv ${XASPrefix}.xas.$XAS_ARG.xas $XASFile
		
		paste Spectra $XASFile > $$
		mv $$ Spectra    
		fi
	    done
	else


	    XASPrefix=$dir/${ATOM_SYMBOL[$atom]}${atom}/$MOLNAME.${ATOM_SYMBOL[$atom]}${atomw}-$CHAPPROX
	    
    #get the specific energies for shift
	    if [[ ! -f ${XASPrefix}.scf.out ]]; then
		echo "missing SCF output: ${XASPrefix}.scf.out"
		exit
	    fi
	    EXCH=`grep '^!' ${XASPrefix}.scf.out | awk '{print $5}' | tail -1`
	    
	    if [[ ! -f $dir/$CALCGS/$MOLNAME.scf.out ]]; then
		echo "missing GS output: $dir/$CALCGS/$MOLNAME.scf.out"
		exit
	    fi
	    EGS=`grep '^!' $dir/$CALCGS/$MOLNAME.scf.out | awk '{print $5}' | tail -1`
	    
    # find the Fermi level (and CBM)
#	    EFERM=`grep "the Fermi energy is"  ${XASPrefix}.scf.out | awk '{print $5}' | tail -1`
	    echo "NELEC=" $NELEC
	    NELXCH=`echo "$NELEC+1" | bc -l` 
	    $GLOBAL_PARA_PREFIX $ExecDir/efermi.x $PARA_POSTFIX $TEMP ${XASPrefix}.xas.$XAS_ARG $NELXCH > ${XASPrefix}.xas.${XAS_ARG}-fermi.out_XCH
	    EFERM=`grep 'CBM =' ${XASPrefix}.xas.${XAS_ARG}-fermi.out_XCH | sort -gr -k3,3 | tail -1 | awk '{print $3}'`

	    $GLOBAL_PARA_PREFIX $ExecDir/efermi.x $PARA_POSTFIX $TEMP ${XASPrefix}.xas.$XAS_ARG $NELEC > ${XASPrefix}.xas.${XAS_ARG}-fermi.out
	    
	    ELUMO=`grep 'CBM =' ${XASPrefix}.xas.${XAS_ARG}-fermi.out | sort -gr -k3,3 | tail -1 | awk '{print $3}'`
	    
	    E=`echo | awk -v exch=$EXCH -v erefes0=$ERefES0 -v egs=$EGS -v erefgs0=$ERefGS0 -v elumo=$EFERM '{print (exch-erefes0-egs+erefgs0)*13.6056923 - elumo}'`
	    
	    echo "i : (XCH - XCH0) - (GS - GS0) - LUMO = Delta"
	    echo $atom : $EXCH $ERefES0 - $EGS $ERefGS0 $EFERM = $E
	    
    # check for previous .xas file and save
	    for file in `ls $dir/${ATOM_SYMBOL[$atom]}${atom}/*.xas 2> /dev/null`; do
		mv $file ${file}-raw
	    done
	    
    # run xas_para.x
	    echo "executing $GLOBAL_PARA_PREFIX $ExecDir/xas_para.x $PARA_POSTFIX $ELOW $EHIGH $NENER $SIGMA $E $ELUMO $STRETCH ${XASPrefix}.xas.$XAS_ARG 2> /dev/null"
	    $GLOBAL_PARA_PREFIX $ExecDir/xas_para.x $PARA_POSTFIX $ELOW $EHIGH $NENER $SIGMA $E $ELUMO $STRETCH ${XASPrefix}.xas.$XAS_ARG 2> /dev/null
	    
	    XASFile="`echo $dir | sed 's/^.*\///'`-`basename ${XASPrefix}.xas.$XAS_ARG.xas`"
	    mv ${XASPrefix}.xas.$XAS_ARG.xas $XASFile
	    
	    paste Spectra $XASFile > $$
	    mv $$ Spectra    
	fi

	declare "NAVE_$elem"=$(( NAVE_$elem + 1 ))
	
    done
    
done

# loop over types and generate average spectra
for t in $TYPES; do
    
    if [[ ! -d $MyDir/XAS/Spectrum-$t ]]; then
	echo "Unable to find Spectrum directory for type $t
$MyDir/XAS/Spectrum-$t"
	exit
    fi
    cd $MyDir/XAS/Spectrum-$t
    
    NAVE=$((NAVE_$t))
    
    AWKSTR=''
    AVESTR='(0'
    for a in `seq 1 $NAVE`
    do
	b=`echo $a | awk '{print $1*5-3}'`
	AWKSTR="${AWKSTR}, \$$b"
	AVESTR="${AVESTR}+\$$b"
    done
    AVESTR="${AVESTR})/$NAVE"
    #echo "{print \$1, ${AVESTR}${AWKSTR}}"
    awk "{print \$1, ${AVESTR}${AWKSTR}}" Spectra > Spectrum-Ave-$t
    perl -pi -e "s/#SPIN.*$/\n#SPIN/" Spectrum-Ave-$t
    echo output located in 
    echo   $SpecDir/Spectrum-Ave-$t
done
exit

#!/bin/bash
# Computes the reference GS calculations for each xyz-file in XYZFILES and
# computes the atomic GS and XCH reference energies.
# Assumption that cell dimesions and atomic species do not vary among xyz-files
#
# If atomic convergence is an issue, then alter TMP_ELEC_MIXING_BETA
#
# davegp

##PBS -q lr_batch
##PBS -l nodes=9:ppn=16:lr3
##PBS -l walltime=02:00:00

#SBATCH --job-name=test
#SBATCH --partition=vulcan 
#SBATCH --nodes=2
#SBATCH --tasks-per-node=8
#SBATCH --time=01:00:00
#SBATCH --qos=normal
#SBATCH --account=vulcan


#module unload intel
#module load intel/2011.11.339 mkl fftw openmpi

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

# Load user variables
. ./Input_Block.in

# number of processors
if [[ -n $PBS_NODEFILE ]]; then
  # Maybe we are at NERSC
  if [[ -n $CRAY_SITE_LIST_DIR ]]; then
    PPN=`qstat -f $PBS_JOBID | grep Resource_List.mppwidth | awk '{print $3}'`

    export PBS_NODEFILE="$PBS_O_WORKDIR/PBS_NODEFILE-$PBS_JOBID"
    :> $PBS_NODEFILE
    for i in `seq 1 $PPN`; do
      echo $i >> $PBS_NODEFILE
    done
  else
  # otherwise this is just a regular PBS job with all nodes listed in nodefile
    PPN=`wc -l $PBS_NODEFILE | awk '{print $1}'`
  fi
elif [[ -n $SLURM_NODELIST ]]; then
    export PBS_NODEFILE="$PBS_O_WORKDIR/PBS_NODEFILE-$PBS_JOBID" 
    :> $PBS_NODEFILE
    echo $SLURM_NODELIST > slurm.nodes ;
    perl -pi -e "s/,/ /g" slurm.nodes ;
    ncpu=`echo $SLURM_CPUS_ON_NODE`;
    PPN=0;
    for node in `cat slurm.nodes` ; do
	for((i=1;i<=$ncpu;i++));do
	    echo $node >> $PBS_NODEFILE ;
	    ((PPN++));
	done
    done
else
  # otherwise this is an interactive job without a scheduler (we must make a fake scheduler)
  if [ $# -lt 1 ]; then
    echo " interactive usage: ./XAS.sh N"
    echo "                    where N is the total number of processors you want to use"
    exit
  fi
  PPN=$1
  export PBS_JOBID=fake
  export PBS_NODEFILE="$PBS_O_WORKDIR/PBS_NODEFILE-$PBS_JOBID"
  :> $PBS_NODEFILE
  for i in `seq 1 $PPN`; do
    echo localhost >> $PBS_NODEFILE
  done
fi
echo " total number of processors = $PPN"

# divide by number of simultaneous atomic calculations
PPN=`echo $PPN / $NJOB | bc`
echo " number of independent calculations (NJOB) = $NJOB"
echo " number of processors given to each independent calculation = $PPN"

# Location of executables and script libraries
BibDir="$BIB_ROOT/Bibliothek"
export BibDir
ExecDir="$SHIRLEY_ROOT/bin"

# Calculation-specific file info
MyDir=`pwd`
export NodePool="${MyDir}/Nodes"
ScriptName=`echo "$0" | sed -e 's/\.\///g'`

# useful script shortcuts
ReVa="${BibDir}/ResetVariables.sh"
VP="${BibDir}/VarPen.sh" #Executable to change the file TMP_INPUT.in
Re="${BibDir}/Reverend.sh"
GeNo="${BibDir}/GetNode.sh"
xyz2inp="${BibDir}/../xyz2inp.sh"

#end of declarations

################################################################################
# Chop up the total number of processors into chunks defined by PPN
${BibDir}/Chop.sh $PBS_NODEFILE $PPN        #decompose parent node file
################################################################################

#XAS start

if [[ ! -d $MyDir/XAS ]]; then
mkdir $MyDir/XAS
fi

for xyzfile in $XYZFILES
do

  if [[ -f $xyzfile ]]; then
    echo " working on xyz-file $xyzfile"
  else
    echo " unable to find xz-file $xyzfile - skipping"
    continue
  fi

  # This defines the name of the directory where calculation results are stored
  # This should be linked to an xyz file
  CALC=`echo $xyzfile | sed 's/.xyz$//'`

  if [[ ! -d $MyDir/XAS/$CALC ]]; then
    mkdir $MyDir/XAS/$CALC
  fi

  # Do a GS calculation
  dir=$MyDir/XAS/$CALC/GS
  if test ! -d $dir ; then
    mkdir $dir
  fi
  SCFOUT=$dir/${MOLNAME}.scf.out
  if [ ! -f $SCFOUT ] || ! grep -q 'convergence has been achieved' $SCFOUT
  then

    # Make a copy of Input_Block.in in the working directory
    inpblk=$MyDir/XAS/$CALC/Input_Block.in${PBS_JOBID}
    cp $MyDir/Input_Block.in $inpblk
    echo " converting xyz to input..."
    $xyz2inp $xyzfile $XYZUNIT $XASELEMENTS >> $inpblk
    echo " ... done"
    . $inpblk

    cp $inpblk $dir/Input_Block.in

#    if [[ -n $SPIN  &&  -n $GS_MAG && $GS_MAG -gt 0 ]]  ; then
    if [[ -n $SPIN  &&  -n $GS_MAG  ]]  ; then
	SP2="tot_magnetization=$GS_MAG";
	SPIN=`echo -e $SPIN '\n' $SP2` ;	
	echo "#----redefine spin------" >> $dir/Input_Block.in
	echo SPIN="'$SPIN'" >> $dir/Input_Block.in ;
    fi

    #Create file: TMP_INPUT.in
    inp=$dir/TMP_INPUT.in${PBS_JOBID}

    $ReVa $dir
    ${VP} $inp TMP_BIN_DIR=$ExecDir TMP_wf_collect=.TRUE. TMP_disk_io='low' TMP_CHAPPROX=''

    ${GeNo} $dir $NodePool || continue
    # refresh file system
    ls $dir/Node-${PBS_JOBID}* > /dev/null

    # Run the calc
    echo " Ground State Calculation"
    echo " - working directory: $dir"
    ${BibDir}/SCF.sh $dir 1 &

    sleep 1
    ls $NodePool > /dev/null

  else
    echo "GS SCF completed for $dir"
  fi
done

# Now do atomic calculations - assuming that volume is the same for all xyz's
# atom GS and XAS
# subdirectory where calculations will appear - can be anything you want
CALC=atom

if test ! -d $MyDir/XAS/${CALC} ; then
mkdir $MyDir/XAS/${CALC}
fi
CalcDir=${MyDir}/XAS/${CALC}

# Make a copy of Input_Block.in in the working directory
# use the first xyzfile as template
#xyzfile=`echo $XYZFILES | awk '{print $1}'`
#inpblk=$MyDir/XAS/$CALC/Input_Block.in${PBS_JOBID}
#cp $MyDir/Input_Block.in $inpblk
#echo " converting xyz to input..."
#$xyz2inp $xyzfile $XYZUNIT $XASELEMENTS >> $inpblk
#echo " ... done"
#. $inpblk

TYPES=''
for xyzfile in $XYZFILES; do
  TYPES=`(echo "$TYPES" ; cat $xyzfile | awk '{if(NR>2){print $1}}') | sed 's/[0-9]*$//' | sort | uniq` 
done
NTYP=`echo $TYPES | wc -w | awk '{print $1}'`
echo TYPES = $TYPES
echo NTYP = $NTYP
echo XASEL = $XASELEMENTS
# Get periodic table info
. $BibDir/periodic.table $PSEUDO_FUNCTIONAL

# loop over types
for TYP_SYMBOL in $TYPES; do

  DO_TYP=0
  for E in `echo $XASELEMENTS` ; do
    E=`echo $E | sed 's/[0-9]*$//'`
    if [ $E == $TYP_SYMBOL ] ; then
      DO_TYP=1
    fi
  done

  # if we are interested in this atom type then do an atomic calc
  if [ $DO_TYP -eq 1 ] ; then
    echo "type=" $TYP_SYMBOL

    AN=`get_atomicnumber $TYP_SYMBOL`
    TYP_Z=${PSEUDO[$AN]}
    echo $PSEUDO_DIR/$TYP_Z
    if [ ! -f $PSEUDO_DIR/$TYP_Z ] ; then
      echo Pseudopotential not found: $PSEUDO_DIR/$TYP_Z
      exit
    fi
    NELEC=`grep 'Z valence' $PSEUDO_DIR/$TYP_Z | awk '{print $1}'`
    # Determine the number of bands
    NBND=`echo "($NELEC*0.5*1.2)/1" | bc`
    NBND1=`echo "($NELEC*0.5)/1+4" | bc`
    if [[ $NBND1 -gt $NBND ]] ; then
      NBND=$NBND1
    fi
    NBND2=`echo "$NBND*$PPP" | bc`
    if [[ $NBND2 -lt $PPN ]] ; then
      PPN=$PPP
################################################################################
# Chop up the total number of processors into chunks defined by PPN
${BibDir}/Chop.sh $PBS_NODEFILE $PPN        #decompose parent node file
################################################################################
    fi


    TMP_ATOMIC_SPECIES="ATOMIC_SPECIES
$TYP_SYMBOL ${MASS[$AN]} ${PSEUDO[$AN]}"
    TMP_ATOMIC_POSITIONS="ATOMIC_POSITIONS (angstrom)
$TYP_SYMBOL 0.0 0.0 0.0"

    if [[ -n $SPIN ]]; then
#      echo "Yo!" $SPIN
      TMP_SPIN=`echo $SPIN | sed -e "s/start.*//"`
      TMP_SPIN=`echo '"' $TMP_SPIN starting_magnetization\(1\)=1 '"'`
      #sedstr="/starting_magnetization($ityp)/s/.*starting_magnetization($ityp) *= *\([+-]*[0-9]*\.*[0-9]*\).*/starting_magnetization(1)=\1/p"
      #TMP_SPIN=`echo $SPIN | sed -n "$sedstr"`
      #echo "Yo1!" $TMP_SPIN
      #if [[ -n $TMP_SPIN ]]; then
        #TMP_SPIN=`echo $SPIN | sed -n '/nspin/s/.*\(nspin *= *[0-9]*\).*/\1/p'
#echo $TMP_SPIN`
      #fi
    fi
#    echo "Yo!" $TMP_SPIN
    if [[ -n $LDAU ]]; then
      sedstr="/Hubbard_U($ityp)/s/.*Hubbard_U($ityp) *= *\([+-]*[0-9]*\.*[0-9]*\).*/Hubbard_U(1)=\1/p"
      TMP_LDAU=`echo $LDAU | sed -n "$sedstr"`
      if [[ -n $TMP_LDAU ]]; then
        TMP_LDAU=`echo $LDAU | sed -n '/lda_plus_u/s/.*\(lda_plus_u *= *\.[a-zA-Z]*\.\).*/\1/p'
echo $TMP_LDAU`
      fi
    fi

    echo $TYP_SYMBOL
    echo $NELEC $NBND
    echo "$TMP_ATOMIC_SPECIES"
    echo "$TMP_ATOMIC_POSITIONS"
    echo "$TMP_SPIN"
    echo "$TMP_LDAU"

    dir=${CalcDir}/${TYP_SYMBOL}
    if test ! -d $dir ; then
      mkdir $dir
    fi

    SCFOUT=$dir/atom.${TYP_SYMBOL}-GS.scf.out
    if [ ! -f $SCFOUT ] || ! grep -q 'convergence has been achieved' $SCFOUT
    then

      cp $MyDir/Input_Block.in $dir/Input_Block.in

      #Create file: TMP_INPUT.in
      inp=$dir/TMP_INPUT.in${PBS_JOBID}

      $ReVa $dir  #Create file: TMP_INPUT.in
      TMP_ELEC_MIXING_BETA=0.1
      ${VP} $inp TMP_BIN_DIR=$ExecDir TMP_MOLNAME=atom.${TYP_SYMBOL}-GS TMP_wf_collect=.TRUE. TMP_disk_io='low' TMP_NAT=1 TMP_NELEC=$NELEC TMP_NBND=$NBND TMP_ATOMIC_SPECIES="$TMP_ATOMIC_SPECIES" TMP_ATOMIC_POSITIONS="$TMP_ATOMIC_POSITIONS" TMP_CHAPPROX='' TMP_ELEC_MIXING_BETA=$TMP_ELEC_MIXING_BETA TMP_SPIN=$TMP_SPIN TMP_LDAU=$TMP_LDAU TMP_tot_charge=0 

      ${GeNo} $dir $NodePool || continue
      ls $dir/Node-${PBS_JOBID}*

      ${BibDir}/SCF.sh $dir 1 &

      sleep 1
      ls $NodePool

    fi

    SCFOUT=$dir/atom.${TYP_SYMBOL}-${CHAPPROX}.scf.out
    if [ ! -f $SCFOUT ] || ! grep -q 'convergence has been achieved' $SCFOUT
    then

      cp $MyDir/Input_Block.in $dir/Input_Block.in

      #Create file: TMP_INPUT.in
      inp=$dir/TMP_INPUT.in${PBS_JOBID}

      $ReVa $dir  #Create file: TMP_INPUT.in
      TMP_ELEC_MIXING_BETA=0.1
      ${VP} $inp TMP_BIN_DIR=$ExecDir TMP_MOLNAME=atom.${TYP_SYMBOL}-${CHAPPROX} TMP_wf_collect=.TRUE. TMP_disk_io='low' TMP_NAT=1 TMP_NELEC=$NELEC TMP_NBND=$NBND TMP_ATOMIC_SPECIES="$TMP_ATOMIC_SPECIES" TMP_ATOMIC_POSITIONS="$TMP_ATOMIC_POSITIONS" TMP_CHAPPROX=${CHAPPROX} TMP_ELEC_MIXING_BETA=$TMP_ELEC_MIXING_BETA TMP_SPIN=$TMP_SPIN TMP_LDAU=$TMP_LDAU TMP_tot_charge=0
      ${BibDir}/Labeler_atom.sh 1 $inp

      ${GeNo}  $dir $NodePool || continue
      ls $dir/Node-${PBS_JOBID}*

      ${BibDir}/SCF.sh $dir 1 &

      sleep 1
      ls $NodePool

    fi

  fi
done

wait
exit


#!/bin/bash
# XAS-xyz.sh - this script manages XAS calculations based on a set of xyz-files
# provided as input via the Input_Block.in. This should be an improvement over
# previous scripts which required users to update Input_Block.in themselves
# using the xyz2inp.sh script (now incorporated directly into this script).
# Makes use of XYZFILES variable defined in Input_Block.in
# Note that this script should be code-version agnostic - version information 
# should be provided in Input_Block.in via the SHIRLEY_ROOT variable
#
# davegp, Dec 16, 2011, TMF

##PBS -q lr_batch
##PBS -l nodes=1:ppn=8:lr1
##PBS -l walltime=24:00:00
##PBS -N fe2uc10
##PBS -V

#SBATCH --job-name=test
#SBATCH --partition=vulcan 
#SBATCH --nodes=2
#SBATCH --tasks-per-node=8
#SBATCH --time=01:00:00
#SBATCH --qos=normal
#SBATCH --account=vulcan

#
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

# check that PPP divides PPN
if [ $((PPN % PPP)) -ne 0 ] ; then
  echo " procs for each job $PPN should be divisible by procs per pool $PPP for shirley_xas.x"
  echo " the remainder is $((PPN % PPP))"
  exit
fi
echo " number of processors per pool used by shirley_xas.x = $PPP"
export PPP

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
#xyz2inp="./xyz2inp.sh"
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
  CALC=`echo $xyzfile | sed 's/\.xyz$//'`

  if [[ ! -d $MyDir/XAS/$CALC ]]; then
    mkdir $MyDir/XAS/$CALC
  fi

  # Make a copy of Input_Block.in in the working directory
  inpblk=$MyDir/XAS/$CALC/Input_Block.in${PBS_JOBID}
  cp $MyDir/Input_Block.in $inpblk
  echo " converting xyz to input..."
  $xyz2inp $xyzfile $XYZUNIT $XASELEMENTS >> $inpblk
  echo " ... done"
  . $inpblk

  SPIN_CALC=0 ;
#  if [[ -n $SPIN  &&  -n $GS_MAG && $GS_MAG -gt 0 ]]  ; then
  if [[ -n $SPIN  &&  -n $GS_MAG  ]]  ; then
      echo "SPIN and GS_MAG greater than 0...will do spin calc" ;
      MAG_UP=`echo "$GS_MAG+1" | bc` ;
      MAG_DN=`echo "$GS_MAG-1" | bc` ;
      SPIN_CALC=1;
  fi

  # loop over atoms
  for i in `seq 1 $NAT` ; do

      if [ $SPIN_CALC -eq 1 ] ; then

	  if [ ${IND_EXCITATION[$i]} -eq 1 ] ; then
#------------------UP---------------------------------------------------------
	      dir=$MyDir/XAS/$CALC/${ATOM_SYMBOL[$i]}${i}_UP
	      if test ! -d $dir ; then
		  mkdir $dir
	      fi
	      cp $inpblk $dir/Input_Block.in

	      SPIN=`echo $SPIN | sed -e s/tot_magnetization=.*//g` ;
	      SP2="tot_magnetization=$MAG_UP";
#	      SPIN=`echo -e $SPIN '\n' $SP2` ;
#	      echo "#----redefine spin------" >> $dir/Input_Block.in
#	      echo SPIN="'$SPIN'" >> $dir/Input_Block.in ;
	      echo "#----define tmag------" >> $dir/Input_Block.in
	      echo TMAG="'$SP2'" >> $dir/Input_Block.in ;
	      
	      
	      CATOM=`seq -w $i $NAT | head -1`
	      
      #Create file: TMP_INPUT.in
	      inp=$dir/TMP_INPUT.in${PBS_JOBID}
	      $ReVa $dir
	      ${VP} $inp TMP_BIN_DIR=$ExecDir TMP_MOLNAME=${MOLNAME}.${ATOM_SYMBOL[$i]}${CATOM}-${CHAPPROX} TMP_wf_collect=.TRUE. TMP_disk_io='high'
	      ${BibDir}/Labeler.sh $i $inp
	      
      # Get nodes - or skip atom if error
	      ${GeNo} $dir $NodePool || continue
      # refresh file system
	      ls $dir/Node-${PBS_JOBID}* > /dev/null
	      
      # Run the calc
	      echo " treating atom ${ATOM_SYMBOL[$i]}${i}_UP"
	      echo " - working directory: $dir"
#	      exit
	      ${BibDir}/DoXAS.sh $dir $i ${ATOM_SYMBOL[$i]} $NodePool &
	      sleep 1
      # refresh file system
	      ls $NodePool > /dev/null

#------------------DN---------------------------------------------------------
	      if [ $GS_MAG -gt 0 ] ;  then
	      dir=$MyDir/XAS/$CALC/${ATOM_SYMBOL[$i]}${i}_DN
	      if test ! -d $dir ; then
		  mkdir $dir
	      fi
	      cp $inpblk $dir/Input_Block.in

	      SPIN=`echo $SPIN | sed -e s/tot_magnetization=.*//g` ;
	      SP2="tot_magnetization=$MAG_DN";
#	      SPIN=`echo -e $SPIN '\n' $SP2` ;
#	      echo "#----redefine spin------" >> $dir/Input_Block.in
#	      echo SPIN="'$SPIN'" >> $dir/Input_Block.in ;
	      echo "#----define tmag------" >> $dir/Input_Block.in
	      echo TMAG="'$SP2'" >> $dir/Input_Block.in ;
	      
	      CATOM=`seq -w $i $NAT | head -1`
	      
      #Create file: TMP_INPUT.in
	      inp=$dir/TMP_INPUT.in${PBS_JOBID}
	      $ReVa $dir
	      ${VP} $inp TMP_BIN_DIR=$ExecDir TMP_MOLNAME=${MOLNAME}.${ATOM_SYMBOL[$i]}${CATOM}-${CHAPPROX} TMP_wf_collect=.TRUE. TMP_disk_io='high'
	      ${BibDir}/Labeler.sh $i $inp

	      
      # Get nodes - or skip atom if error
	      echo here1
	      ${GeNo} $dir $NodePool || continue
	      echo here2
      # refresh file system
	      ls $dir/Node-${PBS_JOBID}* > /dev/null
	      
      # Run the calc
	      echo " treating atom ${ATOM_SYMBOL[$i]}${i}_DN"
	      echo " - working directory: $dir"
#	      exit
	      ${BibDir}/DoXAS.sh $dir $i ${ATOM_SYMBOL[$i]} $NodePool &
	      
	      sleep 1
      # refresh file system
	      ls $NodePool > /dev/null
	      fi
	  fi
#----------------------------------------------------------------------------------
      else

	  if [ ${IND_EXCITATION[$i]} -eq 1 ] ; then

	      dir=$MyDir/XAS/$CALC/${ATOM_SYMBOL[$i]}${i}
	      if test ! -d $dir ; then
		  mkdir $dir
	      fi
	      cp $inpblk $dir/Input_Block.in
	      
	      CATOM=`seq -w $i $NAT | head -1`
	      
      #Create file: TMP_INPUT.in
	      inp=$dir/TMP_INPUT.in${PBS_JOBID}
	      $ReVa $dir
	      ${VP} $inp TMP_BIN_DIR=$ExecDir TMP_MOLNAME=${MOLNAME}.${ATOM_SYMBOL[$i]}${CATOM}-${CHAPPROX} TMP_wf_collect=.TRUE. TMP_disk_io='high'
	      ${BibDir}/Labeler.sh $i $inp
	      
      # Get nodes - or skip atom if error
	      ${GeNo} $dir $NodePool || continue
      # refresh file system
	      ls $dir/Node-${PBS_JOBID}* > /dev/null
	      
      # Run the calc
	      echo " treating atom ${ATOM_SYMBOL[$i]}$i"
	      echo " - working directory: $dir"
	      ${BibDir}/DoXAS.sh $dir $i ${ATOM_SYMBOL[$i]} $NodePool &
	      
	      sleep 1
      # refresh file system
	      ls $NodePool > /dev/null
	      
	  fi
      fi

  done
done
wait
exit


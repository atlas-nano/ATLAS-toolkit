#!/bin/bash

job=$1 ;
user=pemmaraj ;

file="XAS-xyz_$job.sh" ;
echo $file ;

nmax=5 ;

i=0 ;

while [ $i -le $nmax ]

do
	name=`qstat -u $user | grep $job | awk '{print $4}'` ;

	if [ -z "$name" ] ; then

		echo "job " $name " not found...re-submitting" ;

		qsub $file ;

		(( i += 1 )) ;


	else 

		echo "job " $name " found" ;

	fi

	sleep 300s

done

echo "exiting!"

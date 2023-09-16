#!/bin/bash
echo "usage: mkbandsdat.sh nspin";
nsp=$1;
rm bands.dat 
i=0
for item in `cat hse.eig` ; do
let i++;
done
nbnd=`echo "$i/$nsp" | bc` ;
echo $nbnd | awk '{printf "%8i\n",$1}' > bands.dat ;
echo $nsp | awk '{printf "%8i\n",$1}' >> bands.dat ;

i=1
for item in `cat hse.eig` ; do 
echo $i $item | awk '{printf "%5i%10.5f%10.5f%10.5f%10.5f\n",$1,$2,$2,$2,$2}' >> bands.dat ; 
let i++ ;
if [[ $i -gt $nbnd ]] ; then
i=1 ;
fi 
done


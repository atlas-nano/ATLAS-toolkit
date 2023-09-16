#!/bin/gawk

BEGIN{
	split(type_str,type," ");
}
{
	if(NF>ml && $1 ~ /^[0-9]/){
		indx=int($q/dq); 
		val[$2][indx]++
		tot[$2]++
		if(!($2 in max) || max[$2] < indx)
			max[$2] = indx
		if(!($2 in min) || min[$2] > indx)
			min[$2] = indx
	}
}
END{
	for (i in val){
		fname=sprintf("%s.%s.qdistrib.dat",prefix,type[i]); 
		print "" > fname; 
		j=min[i]
		while(j<=max[i]) {
			v=0
			if(j in val[i])
				v=val[i][j]
			print j*dq,v,v/tot[i]>>fname
			j++
		}
	}
}

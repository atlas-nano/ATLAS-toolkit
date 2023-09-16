#!/bin/gawk

{
	if ($1 !~ /^[0-9]/) next
	count++	
	i=1
	col=1 
	while(i<NF) { 
		for(j=1;j<8;j++) { 
			k=j+i 
			val[col][j] += $k 
			val2[col][j]+=$k*$k 
		} 
		i += 7 
		col++
	}
}
END{
	printf "%-8s %8s %7s %8s %7s %8s %7s %8s %7s %8s %7s %8s %7s %8s %7s\n","shell","nmol","std","nHb/mol","std","%NoHb","std","%ND","std","%SD_","std","%SD|","std","%DD","std";
	for(i=1;i<col;i++) {
		j = 1
		val[i][j] /= count
		val2[i][j] /= count	
		std = sqrt(val2[i][j]-val[i][j]*val[i][j]);
		printf "%-8d %8d %7d ",i,val[i][j],std;
		j = 2
		val[i][j] /= count
		val2[i][j] /= count
		std = sqrt(val2[i][j]-val[i][j]*val[i][j]);
		printf "%8.3f %7.3f ",val[i][j]/val[i][1],std/val[i][1]
		for(j=3;j<8;j++) {
			val[i][j] /= count;
			val2[i][j] /= count;
			std = sqrt(val2[i][j]-val[i][j]*val[i][j]);
			printf "%8.3f %7.3f ",val[i][j],std
		}
		printf "\n";
	}
}

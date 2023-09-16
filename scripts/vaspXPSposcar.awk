#!/bin/gawk

{
	if(NR==1){
		for(i=1;i<=NF;i++){
			ele[i]=$i
		}
		tot=NF;
	} else if(NR==6){
		j=0;
		for(i=1;i<=tot;i++){
			natom[i]=$i
			atomStart[i]=j+1
			if(i>0)
				atomStop[i-1]=j
			k=1;
			while(k<=$i) {
				atomMapEle[++j] = ele[i];
				atomMapTypeID[j] = i;
				k++;
			}
		}
	}
}
END {
	iele = atomMapEle[idx];
	itypID = atomMapTypeID[idx];
	printf "" > "_xps_str";
	printf "" > "_xps_idx";
	for(j=1;j<=tot;j++) {
		printf "%s ",ele[j] >> "_xps_str";
		if(j==itypID) {
			if(natom[j] == 1) {
				printf "%4d",natom[j]>> "_xps_idx";
				id=itypID
			} else if(atomStart[j]==idx) { #case 1: atomID is 1st atom in list
				id=itypID
				printf "%s ",ele[j] >> "_xps_str";
				printf " %3d %3d",1,natom[j]-1 >> "_xps_idx";
			}else if (atomStop[j]==idx) { #case 2: atomID is last atom in list
				id=itypID+1
				printf "%s ",ele[j] >> "_xps_str";
				printf " %3d %3d",natom[j]-1,1 >> "_xps_idx";
			} else { #case 3: atomID is in middle of list
				id=itypID+1
				printf "%s %s ",ele[j],ele[j] >> "_xps_str";
				printf " %3d %3d %3d",idx-atomStart[j],1,natom[j]-(idx-atomStart[j])-1 >> "_xps_idx";
			}
		} else {
			printf "%4d",natom[j]>> "_xps_idx";
		}
	}
	printf "\n" > "_xps_str";
	printf "\n" > "_xps_idx";
	print iele,id
}

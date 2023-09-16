#!/bin/gawk
{
	if(NF==12 && $2 ~ /^core/ && $4 ~ /^core/){
		gsub("O2","O_2CA"); 
		gsub("O1","O_3CA"); 
		gsub("C ","C_CA"); 
		printf " %-11s %-11s %-11s %11.9E %11.5f %11.9E\n",$1,$3,"BUCKINGHAM",$5*23.06,1/$6,$7*23.06;
	}
}

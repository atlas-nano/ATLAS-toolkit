#!/bin/gawk
BEGIN{
	fle=0;
	iread=0;
}
{
	if(FNR==1){
		fle++;
	}
	if(fle==1) {
		if($2 ~ /^atoms$/) {
			iread=1
		}else if ($2 ~ /^bonds$/) {	
			iread=0
		}else if (iread==1 && $2 ~ /^opls/ && NF>7) {
			sub(/opls_/,"",$2);
			sub(/[0-9].*/,"",$5);
			atom[$1][0] = $2;
			atom[$1][1] = $5; #element
			atom[$1][2] = $7; #charge
		}
	} else if(fle==2) {
		if($1 ~ /ATOM|HETATM/ && $2 ~ /^[0-9]/ && NF>9) {
			$10=atom[$2][1] atom[$2][0];
			$13=atom[$2][2];
			printf "%-6s %5d %5s %3s %1s %5s%10.5f%10.5f%10.5f %5s%5d%2d %8.5f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13;
		} else {
			print
		}
	}
}

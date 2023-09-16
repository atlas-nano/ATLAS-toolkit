#!/bin/tcsh
#!/bin/csh

set list = (01 02 03 04 05 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 28 29 30 )
foreach l ($list)
 echo node7-${l}
 echo tpascal
 ssh node7-${l} ls /temp1/tpascal/
 echo
 echo ch121ah 
 ssh node7-${l} ls /temp1/ch121ah 
 echo
end

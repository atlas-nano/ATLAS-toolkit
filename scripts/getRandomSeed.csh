#!/bin/tcsh
#!/bin/csh

set seed = `awk 'BEGIN {srand();myRandomSeed=int(rand()*systime());print myRandomSeed}'`
set slen = `expr length $seed`
if($slen > 6) then
    set seed = `echo $seed | cut -c1-6`
endif
echo $seed

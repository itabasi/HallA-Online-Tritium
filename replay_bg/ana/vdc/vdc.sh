#! /bin/sh


init=111130
end=111500

i=`expr $init`
j=`expr $end`
    echo start run: $i
    echo end run: $j


    root="root -l"
    Gmp=".L GmpVDCt0Calib.cc+"

    eval $root
    echo $Gmp
    eval $Gmp
    
    while [ $i -le $j ]
do
	echo run number: $i


	
	VDC_R=".x vdct0_Right.C($i)"
	echo $VDC_R
	eval $VDC_R


	i=$((i+50));
	echo $i
    done
    	    

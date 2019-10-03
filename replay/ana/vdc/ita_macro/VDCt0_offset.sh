#! /bin/sh

	    
init=111160
end=111600
k=10 # num of branch 
i=`expr $init`
j=`expr $end`
    echo start run: $i
    echo end run: $j


    while [ $i -le $j ]
    do
	  echo run number: $i
	  #	  offset="./bin/VDCt0 -i $i"
	  offset="./bin/vdct0_off -i $i -n $k"	  
#	  echo $offset
	  eval $offset
	i=$((i+k-1));
    done


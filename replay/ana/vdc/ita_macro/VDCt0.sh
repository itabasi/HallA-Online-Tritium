#! /bin/sh

	    
init=111140
end=111840
k=10 # num of branch 
i=`expr $init`
j=`expr $end`
    echo start run: $i
    echo end run: $j


    while [ $i -le $j ]
    do
	  echo run number: $i
	  #	  offset="./bin/VDCt0 -i $i"
	  offset="./bin/VDCt0_raw -i $i -n $k"	  
#	  echo $offset
	  eval $offset
	i=$((i+k-1));
    done


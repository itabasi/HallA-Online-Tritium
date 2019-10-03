#! /bin/sh

	    
init=111160
end=111600

k=10 # Branch Num
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
	i=$((i+k));
    done


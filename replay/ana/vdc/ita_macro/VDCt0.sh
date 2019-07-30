#! /bin/sh

	    
init=111160
end=111479

i=`expr $init`
j=`expr $end`
    echo start run: $i
    echo end run: $j


    while [ $i -le $j ]
    do
	  echo run number: $i
	  #	  offset="./bin/VDCt0 -i $i"
	  offset="./bin/VDCt0_raw -i $i"	  
#	  echo $offset
	  eval $offset
	i=$((i+10));
    done


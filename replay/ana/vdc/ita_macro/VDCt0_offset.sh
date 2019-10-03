#! /bin/sh

	    
init=111160
end=111600

i=`expr $init`
j=`expr $end`
    echo start run: $i
    echo end run: $j


    while [ $i -le $j ]
    do
	  echo run number: $i
	  #	  offset="./bin/VDCt0 -i $i"
	  offset="./bin/vdct0_off -i $i"	  
#	  echo $offset
	  eval $offset
	i=$((i+2));
    done


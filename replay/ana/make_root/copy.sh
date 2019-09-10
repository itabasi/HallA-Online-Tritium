#! /bin/sh


###nnL_small_2 #######
### 111369 -111412 ###

init=111738
end=111830

i=`expr $init`
j=`expr $end`
    echo start run: $i
    echo end run: $j


    while [ $i -le $j ]
    do
	echo run number: $i
	copy="./bin/copy -f /data1/root/tritium\_$i.root -w /data1/nnL_smallroot/tritium\_$i.root"
	#copy="./bin/copy -f /data1/rootfiles/tritium\_$i.root -w /data1/nnL_smallroot/tritium\_$i.root"
	eval $copy
	for((k=1;k<7;k++))
	do
	    copy_sub="./bin/copy -f /data1/root/tritium\_$i\_$k.root -w /data1/nnL_smallroot/tritium\_$i\_$k.root"
	 #   copy_sub="./bin/copy -f /data1/rootfiles/tritium\_$i\_$k.root -w /data1/nnL_smallroot/tritium\_$i\_$k.root"
	    eval $copy_sub
	done
	i=$((i+1));
    done


delete="find /data1/nnL_small/*.root -size -1000 -delete"

echo $delete
eval $delete #delete empty files
    

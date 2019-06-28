#! /bin/sh

	    
init=111480
end=111840

i=`expr $init`
j=`expr $end`
    echo start run: $i
    echo end run: $j


    while [ $i -le $j ]
    do
	echo run number: $i
	copy="./bin/copy_VDC -f /data1/root/tritium\_$i.root -w /data/opt_small/VDC_tuning/tritium\_$i.root"
	eval $copy
	for((k=1;k<7;k++))
	do
	    copy_sub="./bin/copy_VDC -f /data1/root/tritium\_$i\_$k.root -w /data/opt_small/VDC_tuning/tritium\_$i\_$k.root"
	    eval $copy_sub
	done
	i=$((i+1));
    done


delete="find \/data\/opt_small/*.root \-size \-512 \-delete"

echo $delete
eval $delete #delete empty files
    

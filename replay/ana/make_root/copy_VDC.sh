#! /bin/sh

	    
init=111160
end=111840

i=`expr $init`
j=`expr $end`
    echo start run: $i
    echo end run: $j


    while [ $i -le $j ]
    do
	echo run number: $i
	copy="./bin/copy_VDC -f /data2/opt_small/VDC/initial/tritium\_$i.root -w /data2/opt_small/VDC/small/tritium\_$i.root"
#		copy="./bin/copy_VDC -f /data/VDC/root/20190707/tritium\_$i.root -w /data/VDC/small/20190707/tritium\_$i.root"
	eval $copy
#	for((k=1;k<7;k++))
#	do
	    #	    copy_sub="./bin/copy_VDC -f /data/VDC/root/20190707/tritium\_$i\_$k.root -w /data/VDC/small/20190707/tritium\_$i\_$k.root"
#	    	    copy_sub="./bin/copy_VDC -f /data2/opt_small/VDC/initial/tritium\_$i\_$k.root -w /data2/opt_small/VDC/small/tritium\_$i\_$k.root"
#	    eval $copy_sub
#	done
	i=$((i+1));
    done


delete="find \/data2\/opt_small\/VDC\/small\/*.root \-size \-1000 \-delete"

echo $delete
eval $delete #delete empty files
    

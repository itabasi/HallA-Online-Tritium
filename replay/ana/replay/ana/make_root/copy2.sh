#! /bin/sh

	    
init=111190
end=111199

i=`expr $init`
j=`expr $end`
    echo start run: $i
    echo end run: $j


    while [ $i -le $j ]
    do
	echo run number: $i
	copy="./bin/copy -f /data/VDC/root/LU1t0\_2ch/11119X_/tritium\_$i.root -w /data/VDC/small/LU1t0\_2ch/tritium\_$i.root"
	eval $copy
	for((k=1;k<7;k++))
	do
	copy_sub="./bin/copy -f /data/VDC/root/LU1t0\_2ch/11119X_/tritium\_$i\_$k.root -w /data/VDC/small/LU1t0\_2ch/tritium\_$i\_$k.root"	    
	    eval $copy_sub
	done
	i=$((i+1));
    done


delete="find \/data/VDC/small/LU1t0\_2ch/*.root \-size \-512 \-delete"

echo $delete
eval $delete #delete empty files
    

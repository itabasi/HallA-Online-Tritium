#! /bin/sh


###nnL_small_2 #######
### 111369 -111412 ###

init=111160
end=111840

i=`expr $init`
j=`expr $end`
    echo start run: $i
    echo end run: $j


    while [ $i -le $j ]
    do
	echo run number: $i
	# 	copy="./bin/copy -f /data1/root/tritium\_$i.root -w /data1/small/tritium\_$i.root"
#	copy="./bin/copy -f /data1/root/tritium\_$i.root -w /data2/test/small/tritium\_$i.root"
	copy=f"./bin/copy -f /data2/test/root/tritium\_$i.root -w /data2/small/small_Ole/tritium\_$i.root"
		#  	copy="./bin/copy -f /data1/root/tritium\_$i.root -w /data2/rootfiles/small/tritium\_$i.root"
	#	copy="./bin/copy -f /data1/root/tritium\_$i.root -w /data2/opt_small/tritium\_$i.root"
#	copy="./bin/copy -f /data2/opt_small/VDC/t0tuned/t0tuned_all/tritium\_$i.root -w /data2/small/tritium\_$i.root"
	#copy="./bin/copy -f /data2/opt_small/VDC/t0tuned/t0tuned_all/tritium\_$i.root -w /data2/test/tritium\_$i.root" #test

	
	eval $copy
	for((k=1;k<7;k++))
	do
	    #	    	    copy_sub="./bin/copy -f /data1/root/tritium\_$i\_$k.root -w /data1/small/tritium\_$i\_$k.root"
	    #	    	    	    copy_sub="./bin/copy -f /data1/root/tritium\_$i\_$k.root -w /data2/test/small/tritium\_$i\_$k.root"
	    copy_sub=f"./bin/copy -f /data2/test/root/tritium\_$i\_$k.root -w /data2/small/small_Ole/tritium\_$i\_$k.root"
	    #	    copy_sub="./bin/copy -f /data1/root/tritium\_$i\_$k.root -w /data2/rootfiles/small/tritium\_$i\_$k.root"
#	    copy_sub="./bin/copy -f /data1/root/tritium\_$i\_$k.root -w /data2/opt_small/tritium\_$i\_$k.root"	    
#	    copy_sub="./bin/copy -f /data2/opt_small/VDC/t0tuned/t0tuned_all/tritium\_$i\_$k.root -w /data2/small/tritium\_$i\_$k.root"
	  #copy_sub="./bin/copy -f /data2/opt_small/VDC/t0tuned/t0tuned_all/tritium\_$i\_$k.root -w /data2/test/tritium\_$i\_$k.root" #test
	    eval $copy_sub
	done
	i=$((i+1));
    done


#delete="find /data1/nnL_small/*.root -size -1000 -delete"
    #delete="find /data2/rootfiles/small/*.root -size -1000 -delete"
   delete="find /data2/small/small_Ole/*.root -size -1000 -delete"
#echo $delete
eval $delete #delete empty files




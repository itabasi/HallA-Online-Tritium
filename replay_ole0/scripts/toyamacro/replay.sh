#!/bin/sh -f

for seg in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
  echo "$seg"
  ./bin/mk_ -R -f runlist/test.txt -w rootfiles/RS2_${seg}.root -n 1000000
done

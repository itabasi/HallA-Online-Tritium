#!/usr/local/bin/python3

# Toshiyuki Gogami
# Mar 11, 2020

import sys
import time, os.path
import subprocess
from subprocess import call
from subprocess import check_output
from subprocess import run
from subprocess import Popen
from logging import StreamHandler, Formatter, INFO, getLogger
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures.process import ProcessPoolExecutor 


#target = "dummy"
#target = "T2"
#target = "h2"
#target = "h22"
target = "He3"

thisfile = "charge.py"


def main():
    runfile = "../" + target + "_replay.dat"
    outfile = "charge_" + target + ".dat"
    
    inputfile  = open(runfile,"r")
    outputfile = open(outfile,"w")

    lines = inputfile.readlines()
    for line in lines:
        data = line.split()
        com = "./charge " + data[0] + " " + outfile
        time.sleep(0.3)
        proc = subprocess.run(com ,shell=True)



stime = time.time()
main()


print("\n Jobs were done in %.0f sec \n" % float(time.time()-stime))

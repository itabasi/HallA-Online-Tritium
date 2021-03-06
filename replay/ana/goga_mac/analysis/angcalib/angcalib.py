#!/usr/bin/python

# Toshiyuki Gogami
# Nov 2, 2018

import sys
import time, os.path
from subprocess import call
#import concurrent.futures
from logging import StreamHandler, Formatter, INFO, getLogger
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures.process import ProcessPoolExecutor 
import numpy as np

nworkers=2
inpfile = "tunepar.dat"
thisfile = "angcalib.py"

def nude_start(command):
    time.sleep(1.0)
    call(command,shell=True)

def main():
    comlist = []
    inputfile = open(inpfile,"r")
    lines = inputfile.readlines()
    for line in lines:
        data = line.split()
        com = "./angcalib " + data[0]+ " " +data[1]
        comlist.append(com)
    with ProcessPoolExecutor(max_workers=nworkers) as executor:
        executor.map(nude_start,comlist)

stime = time.time()
main()
print("\n Jobs were done in %.0f sec \n" % float(time.time()-stime))


#!/usr/bin/python3
####################################
## Make Missing Mass ##############
###################################



import sys
import time, os.path
from subprocess import call
#import concurrent.futures
from logging import StreamHandler, Formatter, INFO, getLogger
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures.process import ProcessPoolExecutor 
import numpy as np
import shlex
import subprocess
import re
import sys
import os,os.path
import glob
import datetime
nworkers=7


#####################################

fname  ='param/Altuning.param'
#input_matrix ='../matrix/matrix_new.list'
input_matrix ='./matrix/momcalib_matrix_test.list'
input_root   ='../rootfiles/momcalib/momcalib_wAl_init/All.root'
root_dir     ='../rootfiles/momcalib/Altuning'
matrix_dir   ='./matrix/Altuning'
weight=[]
mean=[]
width=[]
#####################################

def ana(i,p0,p1,p2):            
    file_name =f'momcalib_5th_wAl_w{p0}_m{p1}_r{p2}'
    ROOT =file_name +'.root'
    Matrix = file_name
    cmd = f'./bin/momcalib -s {input_root} -m {input_matrix} -r {root_dir}/{ROOT} -w {matrix_dir}/{Matrix} -Al -n 5 -i 1 -0 {p0} -1 {p1} -2 {p2}'
    print(cmd)
    subprocess.run([cmd],shell=True)
    
def main():
    with ProcessPoolExecutor(max_workers=nworkers) as executor:
        weight,mean,width = np.loadtxt(fname,unpack=True)
        nparam = len(weight)
        for i in range(nparam):
            p0 = int(weight[i])
            p1 = int(mean[i])
            p2 = int(width[i])
            print(i,p0,p1,p2)
            mappings = {executor.submit(ana,i,p0,p1,p2)}
        
main()


print()
print(f"=================< COMMENT >=================")
print(f"Input Root directory : {input_root}")
print(f"Input matrix directory : {input_matrix}")
print(f"New Root directory : {root_dir}")
print(f"New matrix directory : {matrix_dir}")
print(f"==============================================")

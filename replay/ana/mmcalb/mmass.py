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
nworkers=3


#####################################
date=datetime.date.today()
matrix="../matrix/matrix_new.list"
root_dir="../rootfiles/mmass/ana_Lambda/"+str(date)
pdf_dir ="../pdf/mmass/ana_Lambda/"+str(date)
name=[]
param=[]
mode=[]
num=1
name.append("Lambda_small_optH1")
mode.append("H1")
name.append("nnL_small_opt1")
mode.append("T")
name.append("Lambda_small_optH2")
mode.append("H1")
name.append("Lambda_small_optT")
mode.append("H2")
name.append("nnL_small_opt2")
mode.append("T")
param.append("f1_Lambda_twc.param")
param.append("f1_Lambda_phase2_tuned.param")

#####################################

if os.path.exists(root_dir):
    while  num<100:
        if os.path.exists(root_dir+"_"+str(num)):
            num+=1
        else :
            root_dir=root_dir+"_"+str(num)
            pdf_dir=pdf_dir+"_"+str(num)
            break
        
print(f"Maked new root directory {root_dir}")
os.mkdir(root_dir)
print(f"Maked new pdf directory {pdf_dir}")
os.mkdir(pdf_dir)


def ana():
    if i<=2 :
        pfile=param[0]
    else :
        pfile=param[1]
            
    names=name[i]
    MODE=mode[i]
    cmd=f'./ana_Lambda -f ../run_list/nnlambda/{names}.list -p param/{pfile} -r {root_dir}/{names}.root -w {pdf_dir}/{names}.pdf -m {matrix} -{MODE}'
#    print(cmd)
    subprocess.run([cmd],shell=True)

def main():


    with ProcessPoolExecutor(max_workers=nworkers) as executor:
        for i in range(5):
            mappings = {executor.submit(ana,i)}

    
main()
#ana()

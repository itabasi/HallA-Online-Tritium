#!/usr/bin/python


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

nworkers=5

init=111160
end=111840

#####################################
root_dir="/data3/root_ole/root"
file_init="tritium_"
small_dir="/data3/root_ole/small"
#####################################

def sub_root(RUNNUM):
    i=0
    for f in glob.glob(f"{root_dir}/tritium_{RUNNUM}*"):
        i=i+1
    return i



def copy(RUNNUM):    
    print(f"go making small root #{RUNNUM}")
    copy =f' ./bin/copy -f {root_dir}/tritium_{RUNNUM}.root -w {small_dir}/tritium_{RUNNUM}.root'
    print(copy)
    subprocess.run([copy],shell=True)
    for sub in range(1,sub_root(RUNNUM)) :
        copy_sub=f' ./bin/copy -f {root_dir}/tritium_{RUNNUM}_{sub}.root -w {small_dir}/tritium_{RUNNUM}_{sub}.root'
        print(copy_sub)
        subprocess.run([copy_sub],shell=True)

def main():
    with ProcessPoolExecutor(max_workers=nworkers) as executor:
        for RUNNUM in range(init, end):
            mappings = {executor.submit(copy,RUNNUM)}
    
main()
cmd(111180)


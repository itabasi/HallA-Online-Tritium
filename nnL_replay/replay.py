#!/usr/bin/python3



# Toshiyuki Gogami
# Nov 2, 2018
# Modified by itabashi

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
nworkers=5
init=111200
end=111300


def cmd(RUNNUM):
    rep =f' analyzer -l \"replay_coinc_new.C({RUNNUM})\" '
    subprocess.run([rep],shell=True)

def main():
    for RUNNUM in range(init, end):
#    cmd(str(RUNNUM))
        with ProcessPoolExecutor(max_workers=nworkers) as executor:
            mappings = {executor.submit(cmd,RUNNUM)}
#        executor.map(cmd,RUNNUM)
    
main()


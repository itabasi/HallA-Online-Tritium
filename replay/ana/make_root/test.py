#!/usr/bin/python3


import sys
import time, os.path
from subprocess import call
import concurrent.futures
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

nworkers=1

init=111220
end=111320


#####################################
root_dir="/data3/root_ole/root2"
file_init="tritium_"
small_dir="/data3/root_ole/small2"
#####################################

def main():
    print('test')
main()


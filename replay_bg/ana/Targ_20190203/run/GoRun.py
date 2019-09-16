#!/usr/bin/python3

import time, threading, os
from subprocess import call

def jobreq(name,runnum):
    cmd = " "
    #cmd1 = "./TARG macro/run.mac input/param.in "
    cmd1 = "./TARG macro/run.mac input/param_l.in "
    #cmd1 = "./TARG macro/run.mac input/param_r.in "

    for i in range(runnum):
        cmd2 = cmd1 + str(i+1) + "; "
        cmd = cmd + cmd2
    print(cmd)
    call(cmd, shell=True)

stime = time.time()
jobreq("a",5)
print("Jobs completed in %.0f sec" % float(time.time()-stime))


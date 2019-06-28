#!/usr/bin/python3

# --- "CreateSeed.py"                  ---
# --- for seed file creation           ---
# --- Toshiyuki Gogami, April 30, 2018 ---

import time, threading, os
import os.path
import concurrent.futures
from subprocess import call


filedir = "../root"
#files = [ "tritium_wcell_noraster.root_1", "tritium_wcell_noraster.root_2"
#         ,"tritium_wcell_noraster.root_3", "tritium_wcell_noraster.root_4"
#         ,"tritium_wcell_noraster.root_5"]
#files = [  "carbon_2.1GeV.root_1", "carbon_2.1GeV.root_2"
#          ,"carbon_2.1GeV.root_3", "carbon_2.1GeV.root_4"
#          ,"carbon_2.1GeV.root_5"]
files = [  "helium4_3Htar.root_1", "helium4_3Htar.root_2"
          ,"helium4_3Htar.root_3", "helium4_3Htar.root_4"
          ,"helium4_3Htar.root_5"]

ofile    = "helium4_3Htar"
seedflag = 400
nevents  = 61000
#nevents  = 10000
#pcent    = 2220.0
pcent    = 600.0
pbite    = 5.5
thcent   = 0.235
thbite   = 0.1
phicent  = 0.0
phibite  = 0.4

com = "./KINEMA_EEK"

nfiles = len(files)


def FileCreation(k):
    inputfile = open("input.dat","w")
    ofile_this = ofile + "/" + ofile + "_" + str(k)
    inputfile.write(filedir + "/" + files[i] + "\n")
    inputfile.write(ofile_this + "\n")
    inputfile.write(str(seedflag) + "\n")
    inputfile.write(str(nevents) + "\n")
    inputfile.write(str(pcent) + " " + str(pbite) + "\n")
    inputfile.write(str(thcent) + " " + str(thbite) + "\n")
    inputfile.write(str(phicent) + " " + str(phibite)) 


ok=os.path.isdir(ofile)
if ok==False:
    com1 = "mkdir " + ofile
    call(com1,shell=True)

stime = time.time()
for i in range(len(files)):
    FileCreation(i+1)
    com2 = "cp input.dat " + "input_" + ofile + "_" + str(i+1) + ".dat"
    #print(com2)
    call(com,shell=True)
    #print(files[i])


print("\n Jobs were done in %.0f sec \n" % float(time.time()-stime))


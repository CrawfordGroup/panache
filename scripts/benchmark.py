#!/usr/bin/env python3

import argparse
import subprocess
import os
import re


argparser = argparse.ArgumentParser()


testdir = "/home/ben/programming/psi4/libpanache/testfiles"
molbasis = "h2o-sto3g"
runtestpath = "runtest/runtest"

ntotalproc = 1
ntotalthread = 1

filepath = os.path.join(testdir, molbasis)


def MPIRun(nthreads, nproc, cyclops, disk, blocksize):
    cmd = ["mpirun", "-n", str(nproc), runtestpath, "-v"]

    if(cyclops):
        cmd.append("-c")
    elif(disk):
        cmd.append("-d")
    
    if(blocksize > 0):
        cmd.extend(["-b",str(blocksize)])

    cmd.append(filepath)

    outpipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    output = outpipe.stdout.read().splitlines()

    timingline = 0
    for i in range(len(output)-1, 0, -1):
        if(re.search(r"LibPANACHE DF Tensor Timings", output[i])):
            timingline = i
            break

    if(timingline == 0):
        raise RuntimeError("Error - can't find timing information")

    timingoutput = output[timingline+6:timingline+11]

    info = {}
    for line in timingoutput:
        m = re.match(r'(\w+) *(\d+) *\( *(\d+) *\) *(\d+) *\( *(\d+) *\) *(\d+) *\( *(\d+) *\)', line)
        if not m:
            raise RuntimeError("Error parsing output!")

        dictinfo = { 'Gen':m.group(2), 'GetBatch':m.group(4), 'GetQBatch':m.group(6) }
        info[m.group(1)] = dictinfo
    
    return info





results = {}

for nproc in range(1,ntotalproc+1):
    nthreadres = {}
    for nthread in range(1,ntotalthread+1):
        resmem  = MPIRun(nthread, nproc, False, False, 0)
        rescyc  = MPIRun(nthread, nproc, True, False, 0)
        resdisk = MPIRun(nthread, nproc, False, True, 0)
        nthreadres[nthread] = { 'Mem':resmem, 'Cyclops':rescyc, 'Disk':resdisk }

    results[nproc] = nthreadres

#print(results[1][1]['Cyclops']['QSO']['Gen'])


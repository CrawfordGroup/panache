#!/usr/bin/env python3

import argparse
import subprocess
import os
import re


argparser = argparse.ArgumentParser()
argparser.add_argument("-m", help="Path to molecule/basis directory to use", required=True)
argparser.add_argument("-r", help="Path to runtest executable", required=True)
argparser.add_argument("-o", help="Output directory", required=True)
argparser.add_argument("-t", help="Max number of threads per process", type=int,  required=True)
argparser.add_argument("-c", help="Test cyclops", action='store_true')
argparser.add_argument("-M", help="Enable testing with MPI", action='store_true')
argparser.add_argument("-n", help="Max number of processes to use (if using mpi)", type=int)
args = argparser.parse_args()

nproc = 1
if (args.M and args.n) and args.c:
    nproc = args.n




def RunTest(mpi, nthreads, nproc, cyclops, disk, blocksize):
    print("NTHREADS: {} NPROC: {}".format(nthreads, nproc))

    os.environ['OMP_NUM_THREADS'] = str(nthreads)

    cmd = []

    if mpi:
        cmd.extend(["mpirun", "-n", str(nproc)])

    cmd.extend([args.r, "-v"])


    if(cyclops):
        cmd.append("-c")
    elif(disk):
        cmd.append("-d")
    
    if(blocksize > 0):
        cmd.extend(["-b",str(blocksize)])

    cmd.append(args.m)

    print("Command: {}".format(cmd))

    outpipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    output = outpipe.stdout.read().splitlines()

    timingline = 0
    for i in range(len(output)-1, 0, -1):
        if(re.search(r'LibPANACHE DF Tensor Timings', output[i])):
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



def PrintRes(tensor, timer):
    column1 = '{:<12} '
    datacolumn = ' {:>13} '
    print("")
    print("Tensor: {}  Timing: {}".format(tensor, timer))

    headstr = column1.format('Type')
    for nthread in range(1,args.t+1):
        headstr += datacolumn.format('{}_Threads'.format(nthread))

    for np in range(1,nproc+1):
        if args.c:
            cycstr = column1.format('Cyclops_{}'.format(np))
    
        if np == 1:
            memstr = column1.format('Memory'.format(np))
            diskstr = column1.format('Disk'.format(np))
    
        for nthread in range(1,args.t+1):
            if args.c:
                cycstr += datacolumn.format(results[np][nthread]['Cyclops'][tensor][timer])
    
            if np == 1:
                memstr += datacolumn.format(results[np][nthread]['Mem'][tensor][timer])
                diskstr += datacolumn.format(results[np][nthread]['Disk'][tensor][timer])
    
    
        if np == 1:
            print(headstr)
            print(memstr)
            print(diskstr)
    
        if args.c:
            print(cycstr)




results = {}



for np in range(1,nproc+1):
    nthreadres = {}
    for nthread in range(1,args.t+1):
        nthreadres[nthread] = {}
        if np == 1:
            nthreadres[nthread]['Mem'] = RunTest(args.M, nthread, np, False, False, 0) 
            nthreadres[nthread]['Disk'] = RunTest(args.M, nthread, np, False, True, 0)
        if args.c:
            nthreadres[nthread]['Cyclops'] = RunTest(args.M, nthread, np, True, False, 0)

    results[np] = nthreadres

PrintRes('QSO','Gen')
PrintRes('QMO','Gen')




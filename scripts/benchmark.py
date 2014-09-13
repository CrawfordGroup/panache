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





results = {}

nproc = 1
if (args.M and args.n) and args.c:
    nproc = args.n


for nproc in range(1,nproc+1):
    nthreadres = {}
    for nthread in range(1,args.t+1):
        nthreadres[nthread] = {}
        if nproc == 1:
            resmem  = RunTest(args.M, nthread, nproc, False, False, 0)
            resdisk = RunTest(args.M, nthread, nproc, False, True, 0)
            nthreadres[nthread]['Mem'] = resmem
            nthreadres[nthread]['Disk'] = resdisk
        if args.c:
            rescyc  = RunTest(args.M, nthread, nproc, True, False, 0)
            nthreadres[nthread]['Cyclops'] = rescyc

    results[nproc] = nthreadres



# Generation of Qso
column1 = '{:<12} '
datacolumn = ' {:>13} '
for nproc in range(1,nproc+1):
    headstr = column1.format('Type')

    if args.c:
        cycstr = column1.format('Cyclops_{}'.format(nproc))

    if nproc == 1:
        memstr = column1.format('Memory'.format(nproc))
        diskstr = column1.format('Disk'.format(nproc))

    for nthread in range(1,args.t+1):
        headstr += datacolumn.format('{}_Threads'.format(nthread))

        if args.c:
            cycstr += datacolumn.format(results[nproc][nthread]['Cyclops']['QSO']['Gen'])

        if nproc == 1:
            memstr += datacolumn.format(results[nproc][nthread]['Mem']['QSO']['Gen'])
            diskstr += datacolumn.format(results[nproc][nthread]['Disk']['QSO']['Gen'])

    print(headstr)
    print(memstr)
    print(diskstr)
    if args.c:
        print(cycstr)

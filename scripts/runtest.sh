#!/bin/bash

RUNTEST=$1
TESTDIR=/home/ben/programming/psi4/libpanache/testfiles

echo "===================================================================="
echo "Testing ${RUNTEST}"
echo "===================================================================="

for T in ${TESTDIR}/*
do
  for N in `seq 1 4`; do
    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:      :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST}          ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:  B   :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST}    -b    ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:    D :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST}       -d ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:  B D :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST}    -b -d ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:T     :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -t       ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:T B   :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -t -b    ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:T   D :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -t    -d ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:T B D :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -t -b -d ${T} | grep OVERALL | awk '{print $3}'`"
  done
  echo 
done

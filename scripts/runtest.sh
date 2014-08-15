#!/bin/bash

RUNTEST=$1
TESTDIR=/home/ben/programming/psi4/libpanache/testfiles

echo "===================================================================="
echo "Testing ${RUNTEST}"
echo "===================================================================="

for T in ${TESTDIR}/*
do
  for N in `seq 1 4`; do
  for B in `seq 0 5`; do
    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:${B}:      :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -b ${B}          ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:${B}:  D   :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -b ${B}    -d    ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:${B}:T     :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -b ${B} -t       ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:${B}:T D   :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -b ${B} -t -d    ${T} | grep OVERALL | awk '{print $3}'`"


    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:${B}:    Q :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -b ${B}       -q ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:${B}:  D Q :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -b ${B}    -d -q ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:${B}:T   Q :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -b ${B} -t    -q ${T} | grep OVERALL | awk '{print $3}'`"

    PREFIX=`printf "%-20s %s" "$(basename $T)" ":${N}:${B}:T D Q :"`
    echo "${PREFIX} `OMP_NUM_THREADS=${N} ${RUNTEST} -b ${B} -t -d -q ${T} | grep OVERALL | awk '{print $3}'`"
  done
  done
  echo 
done

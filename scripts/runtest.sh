#!/bin/bash

RUNTEST=$1
TESTDIR=/home/ben/programming/psi4/libpanache/testfiles

echo "===================================================================="
echo "Testing ${RUNTEST}"
echo "===================================================================="

for T in ${TESTDIR}/*
do
  PREFIX=`printf "%-20s %s" "$(basename $T)" ":     :"`
  echo "${PREFIX} `${RUNTEST}       ${T} | grep OVERALL | awk '{print $3}'`"

  PREFIX=`printf "%-20s %s" "$(basename $T)" ": B   :"`
  echo "${PREFIX} `${RUNTEST} -b    ${T} | grep OVERALL | awk '{print $3}'`"

  PREFIX=`printf "%-20s %s" "$(basename $T)" ":   D :"`
  echo "${PREFIX} `${RUNTEST}    -d ${T} | grep OVERALL | awk '{print $3}'`"

  PREFIX=`printf "%-20s %s" "$(basename $T)" ": B D :"`
  echo "${PREFIX} `${RUNTEST} -b -d ${T} | grep OVERALL | awk '{print $3}'`"

  echo 
done

#!/bin/bash

RUNTESTSH=$1

for R in *intel* *gcc* *clang*
do
  ${RUNTESTSH} ${R}/bin/runtest
done

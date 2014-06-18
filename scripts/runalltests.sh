#!/bin/bash

RUNTESTSH=$1

for R in *intel* *gcc*
do
  ${RUNTESTSH} ${R}/bin/runtest
done

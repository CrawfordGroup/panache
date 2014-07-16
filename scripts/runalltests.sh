#!/bin/bash

RUNTESTSH=$1

(
for I in *gcc49* *clang-lapack*
do
  ${RUNTESTSH} ${I}/bin/runtest
done
)


. /opt/tag/gcc-4.8.3.tag

(

for I in *gcc48*
do
  ${RUNTESTSH} ${I}/bin/runtest
done
)



(
. /opt/tag/intel-2013-sp1-update3.tag

for I in *intel14*
do
  ${RUNTESTSH} ${I}/bin/runtest
done
)


(
. /opt/tag/intel-2015-beta.tag

for I in *intel15*
do
  ${RUNTESTSH} ${I}/bin/runtest
done
)


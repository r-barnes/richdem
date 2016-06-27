#!/bin/bash

#set -x

rm -f output*asc
for p in `seq 1 2`;
do
  for i in `seq 1 10`;
  do
    echo "Performing Test #$i with $((p + 1)) processors"
    mpirun -n $((p + 1)) ./parallel_d8_accum.exe one @evict "tests/dev/testdem${i}.d8" /z/out-%n.tif -w 5 -h 5 &> "test_output$i.log"
    ./assemble_ascii.exe /z/out-%n.tif 1 1 > /z/out-combined
    if [ -a tests/dev/output${i} ]
    then
      diff -w tests/dev/output${i} /z/out-combined
    else
      echo "No good values against which to check Test #$i"
    fi
  done
done

#!/bin/bash

rm -f output*asc
for p in `seq 1 2`;
do
  for i in `seq 1 10`;
  do
    echo "Performing Test #$i with $((p + 1)) processors"
    mpirun -n $((p + 1)) ./a.out "tests/testdem${i}.d8" &> "test_output$i.log"
    if [ -a tests/output${i}_0 ]
    then
      cat output*asc | diff -w tests/output${i}_0 -
    else
      echo "No good values against which to check Test #$i"
    fi
  done
done
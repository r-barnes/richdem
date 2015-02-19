#!/bin/bash

for i in `seq 1 10`;
do
  echo "Performing Test #$i"
  mpirun -n 3 ./a.out "tests/testdem${i}.d8" &> "test_output$i.log"
  if [ -a tests/output${i}_0 ]
  then
    diff output0 tests/output${i}_0
    diff output1 tests/output${i}_1
  else
    echo "No good values against which to check Test #$i"
  fi
done
#!/bin/bash

#set -x

if [[ $1 = "test" ]]; then
  echo "Running output test on '$2'"
  if [ -a tests/dev/$2 ]
  then
    mpirun -n 2 ./parallel_d8_accum.exe one @evict "tests/dev/$2" /z/out-%n.tif -w 5 -h 5 &> "/z/test_output$2.log"
    ./assemble_ascii.exe /z/out-%n.tif 2 2 > /z/out-combined
    cat /z/out-combined
    exit 0
  else
    echo "Test file not found!"
    exit 1
  fi
fi


for p in `seq 1 2`;
do
  for i in `ls tests/dev/*d8 | sed 's/.*\///'`;
  do
    echo "Performing Test '$i' with $((p + 1)) processors"
    mpirun -n $((p + 1)) ./parallel_d8_accum.exe one @evict "tests/dev/$i" /z/out-%n.tif -w 5 -h 5 &> "/z/test_output$i.log"
    ./assemble_ascii.exe /z/out-%n.tif 1 1 > /z/out-combined
    name=$(echo ${i} | cut -f 1 -d '.')
    if [ -a tests/dev/${name}.out ]
    then
      diff -w tests/dev/${name}.out /z/out-combined
    else
      echo "No good values against which to check Test '$i'"
    fi
  done
done

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
    rm -f /z/tout-*.tif
    mpirun -n $((p + 1)) ./parallel_d8_accum.exe one @evict "tests/dev/$i" /z/tout-%n.tif -w 5 -h 5 &> "/z/test_output$i.log"
    width=`ls /z/tout-*tif | sed 's/.*-//' | sed 's/_.*//' | sort -k1 -n | tail -1`
    height=`ls /z/tout-*tif | sed 's/.*_//' | sed 's/\..*//' | sort -k1 -n | tail -1`
    ./assemble_ascii.exe /z/tout-%n.tif "$width" "$height" > /z/out-combined
    name=$(echo ${i} | cut -f 1 -d '.')
    if [ -a tests/dev/${name}.out ]
    then
      diff -w tests/dev/${name}.out /z/out-combined
    else
      echo "No good values against which to check Test '$i'"
    fi
  done
done

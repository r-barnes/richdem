#!/bin/bash

#set -x

TESTLOC=../../tests/flow_accum
TEMPDIR=temp

function DoTest {
  rm -f $TEMPDIR/tout-*tif
  mpirun -n $1 ./parallel_d8_accum.exe one @evict $TESTLOC/$2.d8 $TEMPDIR/tout-%n.tif -w 5 -h 5 &> "$TEMPDIR/test_output_$i.log"
  gdal_merge.py -o $TEMPDIR/tout-combined.tif $TEMPDIR/tout-*tif
  if [ -a $TESTLOC/$2.out ]
  then
    ../../apps/rd_compare.exe $TESTLOC/$2.out $TEMPDIR/tout-combined.tif 1>/dev/null 2>/dev/null
    if [ $? -ne 0 ]; then
      echo "Test '$2' failed!"
    fi
  else
    echo "No good values against which to check Test '$i'"
  fi
}

if [[ $1 = "test" ]]; then
  echo "Running output test on '$2'"
  if [ -a $TESTLOC/$2.out ]
  then
    name=$(echo ${i} | cut -f 1 -d '.')
    DoTest $((p + 1)) $name
    exit 0
  else
    echo "Test file not found!"
    exit 1
  fi
fi

for p in `seq 1 2`;
do
  for i in `ls $TESTLOC/*d8 | sed 's/.*\///'`;
  do
    name=$(echo ${i} | cut -f 1 -d '.')
    echo "Performing Test '$name' with $((p + 1)) processors"
    DoTest $((p + 1)) $name
  done
done

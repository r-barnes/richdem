#!/bin/bash

#This loop reads in a file specified on the command line, one line at a time
while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "Text read from file: $line"

  exec 5>&1
  CMD_OUT=$(eval "$1 install -v $line" 2>&1 |tee >(cat - >&5); exit ${PIPESTATUS[0]})
  exit_code=$?
  exec 5>&-

  if [ $exit_code -ne 0 ]; then #Bitter failure... but wait! We can repipinate this sucker!
    best_version=`echo "$CMD_OUT" | grep Found | sed 's/\s*Found link //' | sed -r 's/\([^)]+\)//' |  sed 's/, version: /|/' | awk -F '|' '{print $2" "$1}' | grep -v rc | grep cp | sort -t. -k 1,1n -k 2,2n -k 3,3n -k 4,4n | tail -n 1 | sed 's/ .*//'`
    vline_stripped=`echo "$line" | sed -r 's/^([A-Za-z]+).*$/\1/'`
    eval "$1 install $vline_stripped==$best_version"
  fi

done < "$2" #The file to be read in.

#!/usr/bin/env python3

import sys

if len(sys.argv)!=2:
  print("{0} <LAYOUT FILE>".format(sys.argv[0]))
  sys.exit(-1)

data = open(sys.argv[1],'r').readlines()

data  = [x.strip().split(',') for x in data]
data  = [ ['#' if len(x.strip())>0 else ' ' for x in line] for line in data]

for line in data:
  print(''.join(line))
#!/usr/bin/env python3

import sys
import numpy as np

def LargestSquareOfOnes(mat):
  mat     = mat.copy()
  bestval = -1
  bestloc = None
  for y in range(data.shape[0]-2,-1,-1):
    for x in range(data.shape[1]-2,-1,-1):
      if mat[y,x]>0:
        mat[y,x] = min(mat[y+1,x],mat[y,x+1],mat[y+1,x+1])+1
        if mat[y,x]>bestval:
          bestval = mat[y,x]
          bestloc = (y,x)
  return bestval,bestloc

if len(sys.argv)!=3:
  print("{0} <prefix/onlydim> <LAYOUT FILE>".format(sys.argv[0]))
  print("\tprefix  = Prefix for layout files")
  print("\tonlydim = Only output the largest square dimension")
  sys.exit(-1)


data = open(sys.argv[2],'r').readlines()

#a = "ned/ned_tiles.layout"
#data  = open(a,'r').readlines()

#data = np.random.rand(10,10)
#data = (data>0.3).astype(np.int)

data  = [x.strip().split(',') for x in data]
fgrid = data.copy()
data  = [ [1 if len(x.strip())>0 else 0 for x in line] for line in data]
data  = np.array(data)

bestval, bestloc = LargestSquareOfOnes(data)

if sys.argv[1]=="onlydim":
  print("{0:43} {1}".format(sys.argv[2],bestval))
else:
  for size in range(1,bestval+1):
    with open(sys.argv[1]+str(size)+".layout",'w') as fout:
      for y in range(bestloc[0],bestloc[0]+size):
        fout.write(','.join(fgrid[y][bestloc[1]:bestloc[1]+size])+"\n")
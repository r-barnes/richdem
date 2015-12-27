#!/usr/bin/env python3

import sys
import glob
import os
import subprocess
import time

def doRaw(cmd):
  print(cmd)
  p          = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  p.wait()
  #while p.poll() is None:
  #  l = p.stdout.readline() # This blocks until it receives a newline.
  #  print(l.decode("utf-8"),end="")
  # When the subprocess terminates there might be unconsumed output
  # that still needs to be processed.
  #print(p.stdout.read().decode("utf-8"),end="")
  return p

def GetFilenames(the_glob):
  #print('Globbing {0}'.format(the_glob))
  return ' '.join(glob.glob(the_glob))

def FillAndTest(authoritative,inpfile,width,height,intermed,tempglob,n=4):
  strat = intermed
  if strat=='@saveall':
    intermed = 'temp/s{width}_{height}'.format(width=width,height=height)
  elif strat=='@retainall' and n<100:
    n = 100
  a = doRaw('time mpirun -n {n} ./parallel_pit_fill.exe one {intermed} {file} temp/{strat}_{width}_{height} -w {width} -h {height}'.format(file=inpfile,width=width,height=height,n=n,strat=strat,intermed=intermed))
  if a.returncode!=0:
    return "BadTileSize"
  a = doRaw('gdal_merge.py -o out/{file}_{strat}_{width}_{height}.tif -of GTiff -n -9999 -a_nodata -9999 {outfiles}'.format(file=inpfile,width=width,height=height,outfiles=GetFilenames('temp/'+strat+'_'+tempglob),strat=strat))
  a = doRaw('gdal_calc.py --NoDataValue=-9999 -A {authoritative} -B out/{file}_{strat}_{width}_{height}.tif --calc="abs(A-B)" --outfile=out/{file}_{strat}_{width}_{height}-diff.tif'.format(file=inpfile,width=width,height=height,authoritative=authoritative,strat=strat))
  a = doRaw('gdalinfo -mm out/{file}_{strat}_{width}_{height}-diff.tif | grep Max'.format(file=inpfile,width=width,height=height,strat=strat))
  comp = a.stdout.readlines()[0]
  print(comp)
  return comp==b'    Computed Min/Max=0.000,0.000\n'

if len(sys.argv)!=2:
  print('Syntax: {0} <TEST GLOB>'.format(sys.argv[0]))
  sys.exit(-1)

files = glob.glob(sys.argv[1])

if not files:
  print('Glob did not match any file!')
  sys.exit(-1)

print('Making directories "out" and "temp"')
if not os.path.exists('out'):
  os.makedirs('out')
if not os.path.exists('temp'):
  os.makedirs('temp')

print('Generating authoritative copy')

sizes = [500,600,700]
name  = os.path.splitext(os.path.basename(sys.argv[1]))[0]

doRaw('mpirun -n 4 ./parallel_pit_fill.exe one @offloadall {file} out/singlecore-'.format(file=sys.argv[1]))

for width in sizes:
  for height in sizes:
    for intermed in ['@offloadall','@retainall','@saveall']:
      outglob = "{width}_{height}{file}-*-fill.tif".format(width=width,height=height,file=name)

      ret = FillAndTest(
        authoritative = 'out/singlecore-{file}-0-fill.tif'.format(file=name),
        inpfile       = sys.argv[1],
        width         = width,
        height        = height,
        intermed      = intermed,
        tempglob      = outglob
      )

      print("\033[96m{intermed}: {width}x{height} - ".format(intermed=intermed,width=width,height=height),end="")
      if ret=="BadTileSize":
        print("\033[91mBad tile size!",end="")
      elif ret==True:
        print("Good",end="")
      else:
        print("\033[91mBad result comparison!",end="")
      print("\033[39m")
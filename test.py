#!/usr/bin/env python3

#In the case of a layout file being used, the `test.py` script will merge all of
#the tiles together. This merged file, or, in the case of a single input file
#being used, that file, will be depression filled using the algorithm in a
#single-core mode. This generates an authoritative answer against which
#correctness is checked. The program then iterates over many tile sizes to
#ensure that they all compare correctly against this authoritative answer.
#See README.md for further details.

import sys
import glob
import os
import subprocess
import time
from osgeo import gdal

def doRaw(cmd):
  print(cmd+"\n")
  p          = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  p.wait()

  #TODO: Mark what to un-comment in order to debug issues

  #while p.poll() is None:
  #  l = p.stdout.readline() # This blocks until it receives a newline.
  #  print(l.decode("utf-8"),end="")
  # When the subprocess terminates there might be unconsumed output
  # that still needs to be processed.
  #print(p.stdout.read().decode("utf-8"),end="")
  return p

def FileInfo(filename):
  """Returns the NoData value and data type of the specified file"""
  src_ds  = gdal.Open( filename )
  srcband = src_ds.GetRasterBand(1)
  nodata  = srcband.GetNoDataValue()
  dtype   = gdal.GetDataTypeName(srcband.DataType)

  return nodata,dtype

def FillAndTest(
  authoritative,
  n,
  many_or_one,
  strat,
  inpfile,
  width  = -1,
  height = -1
):
  for file in glob.glob('temp/manycore-*'):
    os.remove(file)
  for file in glob.glob('temp/test-*dat'):
    os.remove(file)
  if os.path.exists('temp/diff.tif'):
    os.remove('temp/diff.tif')


  ######################
  #Generate filled tiles
  ######################
  a = doRaw("""mpirun -n {n} ./parallel_pf.exe {manyone} {strat} {inpfile} temp/manycore-%n.tif -w {width} -h {height}""".format(
      n        = n,
      manyone  = many_or_one,
      strat    = strat,
      inpfile  = inpfile,
      width    = width,
      height   = height
    ))

  if a.returncode!=0:
    print("\033[96m{strat} {manyone}: {width}x{height} - \033[91mError in fill!\033[39m".format(strat=strat,manyone=many_or_one,width=width,height=height))
    return


  ######################
  #Merge filled tiles
  ######################
  outfiles = glob.glob('temp/manycore-*.tif')

  nodata, dtype = FileInfo(outfiles[0])

  a = doRaw("""gdal_merge.py -o temp/manycore_merged.tif -of GTiff -ot {dtype} -n {nodata} -a_nodata {nodata} {outfiles}""".format(
              dtype    = dtype,
              nodata   = nodata,
              outfiles = ' '.join(outfiles)
      ))
  if a.returncode!=0:
    print("\033[96m{strat} {manyone}: {width}x{height} - \033[91mError in merge!\033[39m".format(strat=strat,manyone=many_or_one,width=width,height=height))
    return

  ######################
  #Difference merged filled tiles with authoritative answer
  ######################
  a = doRaw("""gdal_calc.py --NoDataValue=-9999 -A {authoritative} -B temp/manycore_merged.tif --calc="abs(A-B)" --outfile=temp/diff.tif""".format(
                  authoritative = authoritative
      ))
  if a.returncode!=0:
    print("\033[96m{strat} {manyone}: {width}x{height} - \033[91mError in diff!\033[39m".format(strat=strat,manyone=many_or_one,width=width,height=height))
    return

  ######################
  #Determine if there were differences
  ######################
  a = doRaw('gdalinfo -mm temp/diff.tif | grep Max')
  if a.returncode!=0:
    print("\033[96m{strat} {manyone}: {width}x{height} - \033[91mError in gdalinfo!\033[39m".format(strat=strat,manyone=many_or_one,width=width,height=height))
    return

  comp = a.stdout.readlines()[0]
  print(comp)
  print("\033[96m{strat} {manyone}: {width}x{height} - ".format(strat=strat,manyone=many_or_one,width=width,height=height),end="")
  if comp==b'    Computed Min/Max=0.000,0.000\n':
    print("Good",end="")
  else:
    print("\033[91mBad result comparison!",end="")
  print("\033[39m")


#TODO: Add full debug info switch
def main():
  if len(sys.argv)!=3:
    print('Syntax: {0} <many/one> <Test DEM or Layout File>'.format(sys.argv[0]))
    sys.exit(-1)

  many_or_one = sys.argv[1]
  inputfile   = sys.argv[2]
  basename    = os.path.basename(inputfile)   #Filename
  filename    = os.path.splitext(basename)[0] #Filename without extension

  if not many_or_one in ['many','one']:
    print("Must specify 'many' or 'one'")
    sys.exit(-1)

  if not os.path.exists(inputfile):
    print("Input file '{0}' does not exist!".format(inputfile))
    sys.exit(-1)

  print('Ensuring directory "temp" exists')
  if not os.path.exists('temp'):
    os.makedirs('temp')


  if many_or_one=='many' and not os.path.exists('temp/merged.tif'):
    print('Merging tile files to make dataset for authoritative answer')
    #Merge the many files into one to create an authoritative copy
    with open(inputfile) as fin:
      #Turn the layout file into a list of filenames
      filetiles = fin.read()
      filetiles = filetiles.replace(',',' ').replace('\n',' ').replace('\r',' ')
      filetiles = filetiles.split(' ')

      #Files are assumed to be stored at a location relative to the layout file,
      #so append the layout file's path to each filename in the list
      filetiles = [os.path.join(os.path.dirname(inputfile),x) for x in filetiles]

      #Get NoData value and data type of the input files
      nodata, dtype = FileInfo(filetiles[0])

      #Merge all of the file tiles together using GDAL
      a = doRaw("""gdal_merge.py -o temp/merged.tif -of GTiff \\
                             -ot {dtype} -n {nodata} -a_nodata {nodata} \\
                             {filetiles}""".format(
                              dtype     = dtype,
                              nodata    = nodata,
                              filetiles = ' '.join(filetiles)))
      if a.returncode!=0:
        print('Error merging!')
        sys.exit(-1)


  #Determine what the authoritative input should be
  auth_input = "temp/merged.tif" if many_or_one=='many' else inputfile


  print('Generating authoritative answer')
  if not os.path.exists('temp/singlecore-0.tif'):
    a = doRaw('mpirun -n 4 ./parallel_pf.exe one @offloadall {file} temp/singlecore-%n.tif'.format(file=auth_input))
    if a.returncode!=0:
      print('Error generating authoritative answer!')
      sys.exit(-1)
  else:
    print('Using pre-existing authoritative answer')

  print('Authoritative answer is at: temp/singlecore-0.tif')


  if many_or_one=='many':
    sizes = [-1]
  else:
    sizes = [500,600,700]


  for width in sizes:
    for height in sizes:
      for strat in ['@offloadall','@retainall','temp/test-%n-']:
        ret = FillAndTest(
          authoritative = 'temp/singlecore-0.tif',
          n             = 100 if strat == '@retainall' else 3,
          many_or_one   = many_or_one,
          strat         = strat,
          inpfile       = inputfile,
          width         = width,
          height        = height
        )

main()
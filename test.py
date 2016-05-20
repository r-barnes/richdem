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
import argparse
import multiprocessing
from osgeo import gdal

VERBOSE = False

def doRaw(cmd):
  if VERBOSE:
    print('      '+cmd)
  p      = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  p.wait()
  output = p.stdout.readlines()
  err    = p.returncode
  if VERBOSE:
    print(''.join(map(lambda x: '      '+x.decode("utf-8"),output)))
    print('')
  return output, err

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

  print("\033[96m###{strat} {manyone}: {width}x{height}###\033[39m".format(strat=strat if strat[0]=='@' else '@cache', manyone=many_or_one, width=width, height=height))

  ######################
  #Generate filled tiles
  ######################
  output, err = doRaw("""mpirun -n {n} ./parallel_pf.exe {manyone} {strat} {inpfile} temp/manycore-%n.tif -w {width} -h {height}""".format(
      n        = n,
      manyone  = many_or_one,
      strat    = strat,
      inpfile  = inpfile,
      width    = width,
      height   = height
    ))

  if err!=0:
    print("\033[91mError in fill! Tiles probably too small; does not indicate a serious issue.\033[39m".format(strat=strat,manyone=many_or_one,width=width,height=height))
    return


  ######################
  #Merge filled tiles
  ######################
  outfiles = glob.glob('temp/manycore-*.tif')

  nodata, dtype = FileInfo(outfiles[0])

  output, err = doRaw("""gdal_merge.py -o temp/manycore_merged.tif -of GTiff -ot {dtype} -n {nodata} -a_nodata {nodata} {outfiles}""".format(
              dtype    = dtype,
              nodata   = nodata,
              outfiles = ' '.join(outfiles)
      ))
  if err!=0:
    print("\033[91mError in merge!\033[39m".format(strat=strat,manyone=many_or_one,width=width,height=height))
    return

  ######################
  #Difference merged filled tiles with authoritative answer
  ######################
  output, err = doRaw("""gdal_calc.py --NoDataValue=-9999 -A {authoritative} -B temp/manycore_merged.tif --calc="abs(A-B)" --outfile=temp/diff.tif""".format(
                  authoritative = authoritative
      ))
  if err!=0:
    print("\033[91mError in diff!\033[39m".format(strat=strat,manyone=many_or_one,width=width,height=height))
    return

  ######################
  #Determine if there were differences
  ######################
  output, err = doRaw('gdalinfo -mm temp/diff.tif | grep Max')
  if err!=0:
    print("\033[91mError in gdalinfo!\033[39m".format(strat=strat,manyone=many_or_one,width=width,height=height))
    return

  if output[0].strip()==b'Computed Min/Max=0.000,0.000':
      print("\033[96m!!!Good.\033[39m")
  else:
    print("\033[91m!!!Bad result comparison!\033[39m")
    print("Got: '{0}'".format(output))
    print("Expected: '{0}'".format('Computed Min/Max=0.000,0.000'))


def is_valid_file(parser, arg):
  if not os.path.exists(arg):
    parser.error("Input file '{0}' does not exist!".format(arg))
  else:
    return arg

#TODO: Add full debug info switch
def main():
  global VERBOSE
  parser = argparse.ArgumentParser()
  action = parser.add_mutually_exclusive_group(required=True)
  action.add_argument('--many', action='store_true', help='Test a multiple file DEM. A layout file is expected.')
  action.add_argument('--one',  action='store_true', help='Test a single file DEM. A GDAL readable DEM is expected.')
  parser.add_argument('-v', '--verbose', action='store_true', help='Print all output information: use if something goes wrong.')
  parser.add_argument('--evict',  action='store_true', help='Test the EVICT strategy. (Default is to test all.)')
  parser.add_argument('--cache',  action='store_true', help='Test the CACHE strategy. (Default is to test all.)')
  parser.add_argument('--retain', action='store_true', help='Test the RETAIN strategy. (Default is to test all.)')
  parser.add_argument('--cores',   action='store', type=int, help='Number of cores to use', default=multiprocessing.cpu_count())
  parser.add_argument('inputfile', help='Layout or DEM file', type=lambda x: is_valid_file(parser, x) )
  args = parser.parse_args()

  VERBOSE = args.verbose
  print(VERBOSE)

  print('Ensuring directory "temp" exists')
  if not os.path.exists('temp'):
    os.makedirs('temp')


  if args.many and not os.path.exists('temp/merged.tif'):
    print('Merging tile files to make dataset for authoritative answer')
    #Merge the many files into one to create an authoritative copy
    with open(args.inputfile) as fin:
      #Turn the layout file into a list of filenames
      filetiles = fin.read()
      filetiles = filetiles.replace(',',' ').replace('\n',' ').replace('\r',' ')
      filetiles = filetiles.split(' ')

      #Files are assumed to be stored at a location relative to the layout file,
      #so append the layout file's path to each filename in the list
      filetiles = [os.path.join(os.path.dirname(args.inputfile),x) for x in filetiles]

      #Get NoData value and data type of the input files
      nodata, dtype = FileInfo(filetiles[0])

      #Merge all of the file tiles together using GDAL
      output,err = doRaw("""gdal_merge.py -o temp/merged.tif -of GTiff \\
                             -ot {dtype} -n {nodata} -a_nodata {nodata} \\
                             {filetiles}""".format(
                              dtype     = dtype,
                              nodata    = nodata,
                              filetiles = ' '.join(filetiles)))
      if err!=0:
        print('Error merging!')
        sys.exit(-1)


  #Determine what the authoritative input should be
  auth_input = "temp/merged.tif" if args.many else args.inputfile


  print('Generating authoritative answer')
  if not os.path.exists('temp/singlecore-0.tif'):
    output,err = doRaw('mpirun -n 4 ./parallel_pf.exe one @evict {file} temp/singlecore-%n.tif'.format(file=auth_input))
    if err!=0:
      print('Error generating authoritative answer!')
      sys.exit(-1)
  else:
    print('Using pre-existing authoritative answer')

  print('Authoritative answer is at: temp/singlecore-0.tif')


  if args.many:
    sizes = [-1]
  elif args.one:
    sizes = [500,600,700]

  strats_to_test = []
  if args.evict:
    strats_to_test.append('@evict')
  if args.retain:
    strats_to_test.append('@retain')
  if args.cache:
    strats_to_test.append('temp/test-%n-')
  if len(strats_to_test)==0:
    strats_to_test = ['@evict','@retain','temp/test-%n-']

  for width in sizes:
    for height in sizes:
      for strat in strats_to_test:
        ret = FillAndTest(
          authoritative = 'temp/singlecore-0.tif',
          n             = args.cores,
          many_or_one   = 'many' if args.many else 'one',
          strat         = strat,
          inpfile       = args.inputfile,
          width         = width,
          height        = height
        )

main()
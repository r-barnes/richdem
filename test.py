#!/usr/bin/env python3

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
  #while p.poll() is None:
  #  l = p.stdout.readline() # This blocks until it receives a newline.
  #  print(l.decode("utf-8"),end="")
  # When the subprocess terminates there might be unconsumed output
  # that still needs to be processed.
  #print(p.stdout.read().decode("utf-8"),end="")
  return p

def FileInfo(filename):
  src_ds  = gdal.Open( filename )
  srcband = src_ds.GetRasterBand(1)
  nodata  = srcband.GetNoDataValue()
  dtype   = gdal.GetDataTypeName(srcband.DataType)

  return nodata,dtype

def FillAndTest(
  authoritative,
  n,
  many,
  strat,
  inpfile,
  width  = -1,
  height = -1,
  flipV  = False,
  flipH  = False
):
  intermed=strat
  if strat=='@saveall':
    intermed = 'temp/s{width}_{height}'.format(width=width,height=height)

  filename = os.path.splitext(os.path.basename(inpfile))[0]

  many  = 'many'    if many   else 'one'
  flipV = '--flipV' if flipV  else ''
  flipH = '--flipH' if flipH  else ''

  #Do the stuff
  a = doRaw("""time mpirun -n {n} ./parallel_pf.exe {manyone} {intermed} \\
    {inpfile} temp/{strat}_{width}_{height}_{filename}_ -w {width} -h {height} {flipV}     \\
    {flipH}""".format(
      n        = n,
      manyone  = many,
      intermed = intermed,
      inpfile  = inpfile,
      strat    = strat,
      width    = width,
      height   = height,
      filename = filename,
      flipV    = flipV,
      flipH    = flipH
    ))

  if a.returncode!=0:
    print("\033[96m{strat} {manyone}: {width}x{height} - \033[91mError in fill!\033[39m".format(strat=strat,manyone=many,width=width,height=height))
    return

  outglob = 'temp/{strat}_{width}_{height}_{filename}_*-fill.tif'.format(
      strat    = strat,
      width    = width,
      height   = height,
      filename = filename
    )

  outfiles = glob.glob(outglob)

  nodata, dtype = FileInfo(outfiles[0])

  a = doRaw("""gdal_merge.py -o out/{file}_{strat}_{width}_{height}.tif \
            -of GTiff -ot {dtype} -n {nodata} -a_nodata                 \
            {nodata} {outfiles}""".format(
              file     = filename,
              strat    = strat,
              width    = width,
              height   = height,
              dtype    = dtype,
              nodata   = nodata,
              outfiles = ' '.join(outfiles)
      ))
  if a.returncode!=0:
    print("\033[96m{strat} {manyone}: {width}x{height} - \033[91mError in merge!\033[39m".format(strat=strat,manyone=many,width=width,height=height))
    return

  a = doRaw("""gdal_calc.py --NoDataValue=-9999 -A {authoritative}           \
                -B out/{file}_{strat}_{width}_{height}.tif --calc="abs(A-B)" \
                --outfile=out/{file}_{strat}_{width}_{height}-diff.tif""".format(
                  authoritative = authoritative,
                  file          = filename,
                  strat         = strat,
                  width         = width,
                  height        = height
      ))
  if a.returncode!=0:
    print("\033[96m{strat} {manyone}: {width}x{height} - \033[91mError in calc!\033[39m".format(strat=strat,manyone=many,width=width,height=height))
    return

  a = doRaw('gdalinfo -mm out/{file}_{strat}_{width}_{height}-diff.tif | grep Max'.format(
    file   = filename,
    strat  = strat,
    width  = width,
    height = height
  ))
  if a.returncode!=0:
    print("\033[96m{strat} {manyone}: {width}x{height} - \033[91mError in gdalinfo!\033[39m".format(strat=strat,manyone=many,width=width,height=height))
    return

  comp = a.stdout.readlines()[0]
  print(comp)
  print("\033[96m{strat} {manyone}: {width}x{height} - ".format(strat=strat,manyone=many,width=width,height=height),end="")
  if comp==b'    Computed Min/Max=0.000,0.000\n':
    print("Good",end="")
  else:
    print("\033[91mBad result comparison!",end="")
  print("\033[39m")



def main():
  if len(sys.argv)!=3:
    print('Syntax: {0} <many/one> <Test DEM or Layout File>'.format(sys.argv[0]))
    sys.exit(-1)

  if not sys.argv[1] in ['many','one']:
    print("Must specify 'many' or 'one'")
    sys.exit(-1)

  if not os.path.exists(sys.argv[2]):
    print("Input file '{0}' does not exist!".format(sys.argv[2]))
    sys.exit(-1)

  print('Ensuring directories "out" and "temp" exist')
  if not os.path.exists('out'):
    os.makedirs('out')
  if not os.path.exists('temp'):
    os.makedirs('temp')


  basename = os.path.basename(sys.argv[2])
  filename = os.path.splitext(os.path.basename(sys.argv[2]))[0]


  if sys.argv[1]=='many' and not os.path.exists('out/{filename}_merged.tif'.format(filename=filename)):
    print('Merging tile files to make dataset for authoritative answer')
    #Merge the many files into one to create an authoritative copy
    with open(sys.argv[2]) as fin:
      filetiles = fin.read()
      filetiles = filetiles.replace(',',' ').replace('\n',' ').replace('\r',' ')
      filetiles = filetiles.split(' ')
      filetiles = [os.path.dirname(sys.argv[2])+'/'+x for x in filetiles]

      nodata,dtype = FileInfo(filetiles[0])

      #Merge all of the file tiles together
      a = doRaw("""gdal_merge.py -o out/{filename}_merged.tif -of GTiff \\
                             -ot {dtype} -n {nodata} -a_nodata {nodata} \\
                             {filetiles}""".format(
                              filename  = filename,
                              dtype     = dtype,
                              nodata    = nodata,
                              filetiles = ' '.join(filetiles)))
      if a.returncode!=0:
        print('Error merging!')
        sys.exit(-1)

    auth_input = "out/{filename}_merged.tif".format(filename=filename)
  else:
    auth_input = sys.argv[2]


  if sys.argv[1]=='many':
    filename  += '_merged'
    sizes = [-1]
  else:
    sizes = [500,600,700]


  print('Generating authoritative answer')
  doRaw('mpirun -n 4 ./parallel_pf.exe one @offloadall {file} out/singlecore-'.format(file=auth_input))

  print('Authoritative: '+'out/singlecore-{filename}-0-fill.tif'.format(filename=filename))

  for width in sizes:
    for height in sizes:
      for strat in ['@offloadall','@retainall','@saveall']:
        ret = FillAndTest(
          authoritative = 'out/singlecore-{filename}-0-fill.tif'.format(filename=filename),
          n             = 100 if strat=='@retainall' else 3,
          many          = (sys.argv[1]=='many'),
          strat         = strat,
          inpfile       = sys.argv[2],
          width         = width,
          height        = height,
          flipV         = False,
          flipH         = False
        )

main()
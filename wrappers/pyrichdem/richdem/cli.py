#!/usr/bin/env python

import richdem as rd
import sys
import argparse
from argparse import RawTextHelpFormatter
import numpy as np

def DepressionFilling():
  parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='RichDEM Depression Filling')

  parser.add_argument('dem',     type=str,                     help='Elevation model')
  parser.add_argument('outname', type=str,                     help='Name of output file')
  parser.add_argument('-g', '--gradient', action='store_true', help='Ensure that all cells are at least an epsilon above their downstream cell. This ensures that each cell has a defined flow direction.')
  parser.add_argument('-v', '--version',  action='version', version=rd._RichDEMVersion())
  args = parser.parse_args()

  dem = rd.LoadGDAL(args.dem)
  rd._AddAnalysis(dem, ' '.join(sys.argv))
  rd.FillDepressions(dem, epsilon=args.gradient, in_place=True)
  rd.SaveGDAL(args.outname, dem)



def BreachDepressions():
  parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""RichDEM Depression Breaching""")

  parser.add_argument('dem',                    type=str,                help='Elevation model')
  parser.add_argument('outname',                type=str,                help='Name of output file')
  parser.add_argument('-v', '--version',  action='version', version=rd._RichDEMVersion())
  args = parser.parse_args()

  dem = rd.LoadGDAL(args.dem)
  rd._AddAnalysis(dem, ' '.join(sys.argv))
  rd.BreachDepressions(dem)
  rd.SaveGDAL(args.outname, dem)



def FlowAccumulation():
  parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""RichDEM Flow Accumulation

A variety of methods are available.

Method            Note                           Reference
Tarboton          Alias for Dinf.
Dinf              Alias for Tarboton.
Quinn             Holmgren with exponent=1.
Holmgren(E)       Generalization of Quinn.
Freeman(E)        TODO
FairfieldLeymarie Alias for Rho8.
Rho8              Alias for FairfieldLeymarie.
OCallaghan        Alias for D8.                  10.1016/S0734-189X(84)80011-0
D8                Alias for OCallaghan.          10.1016/S0734-189X(84)80011-0

Methods marked (E) require the exponent argument.
""")
  parser.add_argument('dem',              type=str,                help='Elevation model')
  parser.add_argument('outname',          type=str,                help='Name of output file')
  parser.add_argument('-m', '--method',   type=str, required=True, help='Flow accumulation method to use')
  parser.add_argument('-e', '--exponent', type=float,              help='Some methods require an exponent')
  parser.add_argument('-v', '--version',  action='version', version=rd._RichDEMVersion())
  args = parser.parse_args()

  dem = rd.LoadGDAL(args.dem)
  rd._AddAnalysis(dem, ' '.join(sys.argv))
  accum = rd.FlowAccumulation(dem, method=args.method, exponent=args.exponent)
  rd.SaveGDAL(args.outname, accum)



def TerrainAttribute():
  parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""RichDEM Terrain Attribute

A variety of methods are available.

Parameters:
dem      -- An elevation model
attrib   -- Terrain attribute to calculate. (See below.)
zscale   -- How much to scale the z-axis by prior to calculation

Method:
slope_riserun
slope_percentage
slope_degrees
slope_radians
aspect
curvature
planform_curvature
profile_curvature
""")
  parser.add_argument('dem',              type=str,                help='Elevation model')
  parser.add_argument('outname',          type=str,                help='Name of output file')
  parser.add_argument('-a', '--attrib',   type=str, required=True, help='Terrain attribute to calculate')
  parser.add_argument('-z', '--zscale',   type=float, default=1.0, help='Scale elevations by this factor prior to calculation')
  parser.add_argument('-v', '--version',  action='version', version=rd._RichDEMVersion())
  args = parser.parse_args()

  dem = rd.LoadGDAL(args.dem)
  rd._AddAnalysis(dem, ' '.join(sys.argv))
  tattrib = rd.TerrainAttribute(dem, attrib=args.attrib, zscale=args.zscale)
  rd.SaveGDAL(args.outname, tattrib)



def RdInfo():
  parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""RichDEM Dataset Information

Parameters:
rda      -- A dataset
""")
  parser.add_argument('rda',              type=str,                help='Elevation model')
  parser.add_argument('-s','--show',      action='store_true',     help='Show the model')
  parser.add_argument('--cmap',           type=str, default='jet', help='Colormap (Default: jet)')
  args = parser.parse_args()

  rda = rd.LoadGDAL(args.rda)

  print('File      = {0}    '.format(args.rda                 ))
  print('Data type = {0}    '.format(rda.dtype                ))
  print('Width     = {0}    '.format(rda.shape[1]             ))
  print('Height    = {0}    '.format(rda.shape[0]             ))
  print('Shape     = {0}x{1}'.format(rda.shape[1],rda.shape[0]))
  print('Metadata:')

  for k,v in rda.metadata.items():
    if k=='PROCESSING_HISTORY':
      continue
    print('\t{0} = {1}'.format(k,v))

  print('Processing History:')
  print('-------------------')
  if 'PROCESSING_HISTORY' in rda.metadata:
    for ph in rda.metadata['PROCESSING_HISTORY'].split('\n'):
      ph = ph.strip()
      if len(ph)==0:
        continue
      print(ph)
  print('-------------------')

  if args.show:
    rd.rdShow(rda, cmap=args.cmap)

def RdCompare():
  parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""RichDEM Dataset Comparison

Parameters:
rda      -- A dataset
""")
  parser.add_argument('rda1',              type=str,            help='Elevation model')
  parser.add_argument('rda2',              type=str,            help='Elevation model')
  args = parser.parse_args()

  ds1 = rd.LoadGDAL(args.rda1)
  ds2 = rd.LoadGDAL(args.rda2)

  if ds1.no_data!=ds2.no_data:
    print("NoData differs")

  if ds1.geotransform!=ds2.geotransform:
    print("Geotransform differs")

  if ds1.projection!=ds2.projection:
    print("Projection differs")    

  if np.any(np.isnan(ds1)):
    print("NaN in '{0}'".format(filename))

  if np.any(np.isnan(ds2)):
    print("NaN in '{0}'".format(filename))    

  diff = np.array(ds1-ds2)

  print("Absolute Max difference: {0:10.6}".format(np.nanmax(np.abs(diff))))
  print("Absolute Min difference: {0:10.6}".format(np.nanmin(np.abs(diff))))
  print("Absolute Avg difference: {0:10.6}".format(np.mean  (np.abs(diff))))
  print("RMS difference:          {0:10.6}".format(np.sqrt(np.mean(np.square(diff)))))

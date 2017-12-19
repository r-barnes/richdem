#!/usr/bin/env python

import richdem as rd
import sys
import argparse

def DepressionFilling():
  parser = argparse.ArgumentParser(description='RichDEM Depression Filling')

  parser.add_argument('dem',    type=str,                      help='Elevation model to depression-fill')
  parser.add_argument('outdem', type=str,                      help='Name of output file')
  parser.add_argument('-g', '--gradient', action='store_true', help='Ensure that all cells are at least an epsilon above their downstream cell. This ensures that each cell has a defined flow direction.')
  parser.add_argument('-v', '--version',  action='version', version=rd._RichDEMVersion())
  args = parser.parse_args()

  dem = rd.LoadGDAL(parser.dem)
  rd._AddAnalysis(dem, ' '.join(sys.argv))
  rd.FillDepressions(dem, epsilon=parser.gradient, in_place=True)
  rd.SaveGDAL(parser.outdem, dem)

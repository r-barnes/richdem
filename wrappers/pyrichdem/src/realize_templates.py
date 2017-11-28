#!/usr/bin/env python3

dtypes = [
  'float',
  'double',
  'int8_t',
  'int16_t',
  'int32_t',
  'uint8_t',
  'uint16_t',
  'uint32_t'
]

fin  = open("Array2D_wrapper.hpp.template","r").readlines()
fout = open("Array2D_wrapper.hpp", "w")

for Ttype in dtypes:
  for line in fin:
    lineout = line.replace("@T@", Ttype)
    if "@U@" in lineout:
      for Utype in dtypes:
        ulineout = lineout.replace("@U@", Utype)
        fout.write(ulineout)
    else:
      fout.write(lineout)

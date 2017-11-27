#!/bin/bash
cat Array2D_wrapper.hpp.template | sed 's/@T@/float/g'    >  Array2D_wrapper.hpp
cat Array2D_wrapper.hpp.template | sed 's/@T@/double/g'   >> Array2D_wrapper.hpp
cat Array2D_wrapper.hpp.template | sed 's/@T@/int8_t/g'   >> Array2D_wrapper.hpp
cat Array2D_wrapper.hpp.template | sed 's/@T@/int16_t/g'  >> Array2D_wrapper.hpp
cat Array2D_wrapper.hpp.template | sed 's/@T@/int32_t/g'  >> Array2D_wrapper.hpp
cat Array2D_wrapper.hpp.template | sed 's/@T@/uint8_t/g'  >> Array2D_wrapper.hpp
cat Array2D_wrapper.hpp.template | sed 's/@T@/uint16_t/g' >> Array2D_wrapper.hpp
cat Array2D_wrapper.hpp.template | sed 's/@T@/uint32_t/g' >> Array2D_wrapper.hpp
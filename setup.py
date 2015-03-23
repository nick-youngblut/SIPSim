#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

#setup(name="SIPSim_cpp",
#      ext_modules=[
#          Extension('SIPSimCpp', ['./src/SIPSimCpp.cpp'],
#                    libraries = ["boost_python"]),
#      ]
#)

setup(
    ext_modules = cythonize('./lib/SIPSimCython.pyx'),# './lib/helloworld.pyx'),
#    ext_modules = cythonize('./lib/helloworld.pyx'),
    include_dirs = [numpy.get_include()]

)





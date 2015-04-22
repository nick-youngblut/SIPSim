#!/usr/bin/env python

# SETUP:
## python setup.py build_ext --inplace


from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

setup(name="SIPSim_cpp",
      ext_modules=[
          Extension('SIPSimCpp', ['./src/SIPSimCpp.cpp'],
                    libraries = ["boost_python"]),
      ]
)

setup(
    name = 'SIPSim',
    version = '0.1',
    description = 'Simulate Hi Res Stable Isotope Datasets',
    author = 'Nick Youngblut',
    author_email = 'nyoungb2@gmail.com',
    ext_modules = cythonize('./lib/SIPSimCython.pyx'),
    include_dirs = [numpy.get_include()]
)





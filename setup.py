#!/usr/bin/env python

# SETUP:
## python setup.py build_ext --inplace

from setuptools import setup, find_packages
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


install_reqs = [
    'docopt',
    'cython',
    'dill',
    'numpy',
    'pandas',
    'scipy',
    'pathos',
    'configobj',
    'biopython'
]

#scripts = [
#    'scripts/SIPSim',
#    'scripts/SIPSimR'
#]

# note: include R scripts

## install cpp 
setup(name="SIPSim_cpp",
      ext_modules=[
          Extension('SIPSimCpp', ['./src/SIPSimCpp.cpp'],
                    libraries = ["boost_python"]),
      ]
)

## install main application
desc = 'Simulate High Resolution Stable Isotope Probing Datasets'
setup(
    name = 'SIPSim',
    version = '0.2',
    description = desc,
    long_description = desc + '\n See README for more information.',
    author = 'Nick Youngblut',
    author_email = 'nyoungb2@gmail.com',
    entry_points={
        'console_scripts': [
            'SIPSim = SIPSim.__main__:main'
        ]
    },
    ext_modules = cythonize('./SIPSim/*.pyx'),
    install_requires = install_reqs,
    include_dirs = [numpy.get_include()],
    license = "MIT license",
    packages = find_packages(),
    package_dir={'SIPSim':
                 'SIPSim'},
#    scripts = scripts,
    url = 'https://github.com/nick-youngblut/SIPSim'
)





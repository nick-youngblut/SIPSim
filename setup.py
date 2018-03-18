#!/usr/bin/env python

from setuptools import setup, find_packages, Extension
#from distutils.core import setup
from distutils.extension import Extension
#from Cython.Build import cythonize
import numpy

try:
    from Cython.Build import cythonize
except ImportError:
    cythonize = lambda x: None

# dependencies
install_reqs = [
    'docopt>=0.6.2',
    'cython>=0.25',
    'dill>=0.2',
    'numpy>=1.10',
    'pandas>=0.18',
    'scipy>=0.17',
    'sympy>=1.0',
    'multiprocess>=0.70',
    'pox>=0.2',
    'pathos>=0.2.0',
    'configobj>=5.0.6',
    'biopython>=1.68',
    'dendropy>=4.2.0',
    'intervaltree>=2.1',
    'pyfasta>=0.5',
    'matplotlib>=2.1',
    'SIPSim_pymix>=0.1'
]

## install main application
desc = 'Simulate High Resolution Stable Isotope Probing Datasets'
setup(
    name = 'SIPSim',
    version = '0.3.1',
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
    setup_requires=[
        'cython>=0.27',
    ],
    url = 'https://github.com/nick-youngblut/SIPSim'
)





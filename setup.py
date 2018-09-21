#!/usr/bin/env python
from setuptools import setup, find_packages
from distutils.extension import Extension
import os
import glob
import numpy

# cython
## don't cythonize if no cython dependency
cmdclass = {}
try:
    from Cython.Build import cythonize
    ext = '.pyx'
    cmdclass = {'build_ext' : Cython.Distutils.build_ext}
except ImportError:
    ext = '.c'

## cython extensions
ext_modules = [ ]
cython_files = glob.glob('./SIPSim/*Cython.pyx')
for f in cython_files:
    x = os.path.splitext(os.path.split(f)[1])
    x = 'SIPSim.{}'.format(x[0])
    y = os.path.splitext(f)[0]
    ext_modules.append(Extension(x, [y + ext]))
if ext == '.pyx':
    ext_modules = cythonize(ext_modules)
    
    
# dependencies
install_reqs = [
    'docopt>=0.6.2',
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
    'matplotlib>=2.1,<3',
    'SIPSim_pymix>=0.1.2',
    'MFEprimer_linux>=2.1.1'
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
    cmdclass = cmdclass,
    ext_modules = ext_modules,
    install_requires = install_reqs,
    include_dirs = [numpy.get_include()],
    license = "MIT license",
    packages = find_packages(),
    package_dir={'SIPSim':
                 'SIPSim'},
    url = 'https://github.com/nick-youngblut/SIPSim'
)





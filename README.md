[![Build Status](https://travis-ci.org/nick-youngblut/SIPSim.svg?branch=master)](https://travis-ci.org/nick-youngblut/SIPSim)

SIPSim
======

SIPSim is a toolset for simulating data from high resolution 
stable isotope probing (HR-SIP) experiments.

>Note: currently SIPSim is only Python 2.7 compatable,
mainly because MFEprimer is used for simulating sequences from genomes.
We recommend using an anaconda environment,
which will also help with installing the dependencies.


#### Sections

- [REFERENCE](#reference)
- [INSTALLATION](#installation)
- [TUTORIALS](#tutorials)
- [SIMULATION WORKFLOW](#simulation_workflow)
- [CHANGE LOG](#changelog)
- [LICENSE](#license)


# REFERENCE

[[top](#sections)]

If you use SIPSim, please cite:

> Youngblut, Nicholas D., Samuel E. Barnett, and Daniel H. Buckley. 2018. "SIPSim: A Modeling Toolkit to Predict Accuracy and Aid Design of DNA-SIP Experiments." Frontiers in Microbiology 9: 570. doi: 10.3389/fmicb.2018.00570


# INSTALLATION

[[top](#sections)]

## DEPENENCIES

* [SIPSim_pymix](https://github.com/nick-youngblut/SIPSim_pymix)
  * This is a modified version of [pymix](http://www.pymix.org/pymix/)
     * It has been modified to fix bugs resulting from the integration with SIPSim

* [MFEprimer_linux](https://github.com/nick-youngblut/MFEprimer_linux)
  * This is a modified version of [MFEprimer-2.0](https://github.com/quwubin/MFEprimer)
    * It has been modified for installation into a linux environment via `python setup.py install`

### Dependency install issues (using Anaconda)

* boost-python
  * install via conda or see [boost.org](http://www.boost.org/doc/libs/1_64_0/libs/python/doc/html/index.html) on other methods to install
    * Note: there's currently no official conda channel for boost-python
* scipy libgrfortran issues
  * See https://github.com/ilastik/ilastik-build-conda/issues/17
* scipy MKL issues
  * See https://github.com/BVLC/caffe/issues/3884
  * MKL can be shut down. See [this blog post](https://www.continuum.io/blog/developer-blog/anaconda-25-release-now-mkl-optimizations)
    * This can be done by: `conda install nomkl`
    * NOTE: OpenBlas will try to use all threads. To limit threads, use `export OMP_NUM_THREADS=N`, where `N` is the number of threads to use.

## Installation of SIPSim

### Clone the repo

See the "before_install:" and "install:" sections in `.travis.yml` file for installation instructions.
If the travisCI tests are passing, then these instructions should work on a linux machine.

## Installation of SIPSimR

[SIPSimR](https://github.com/nick-youngblut/SIPSimR) contains R scripts for data
analysis and plotting of data produced by SIPSim. See the SIPSimR README for more information.


# TUTORIALS

[[top](#sections)]

* [An example with 3 genomes](./ipynb/example/1_dataset.ipynb)
* [Recreating Fig 1 from Lueders et al., 2004](./ipynb/example/Lueders2004.ipynb)


# SIMULATION_WORKFLOW

[[top](#sections)]

![simulation pipeline](img/simulation_pipeline.png)


# CHANGELOG

[[top](#sections)]

## v0.3

* Added many unit tests and created SIPSim-specific dependency packages.

## v0.2

* Restructered SIPSim to be fully installable via `setup.py`. This involved spliting the
software into 3 repositories (SIPSim, SIPSimR, and MFEprimer-linux).


# LICENSE

[[top](#sections)]

* Free software: MIT license
SIPSim
======
SIPSim is a toolset for simulating data from high resolution 
stable isotope probing (HR-SIP) experiments.


#### Sections

- [REFERENCE](#reference)
- [INSTALLATION](#installation)
- [TUTORIALS](#tutorials)
- [SIULATION WORKFLOW](#simulation workflow)
- [CHANGE LOG](#changelog)
- [LICENSE](#license)


# REFERENCE

[[top](#sections)]

If you use SIPSim, please cite:

> Youngblut, ND, Buckley DH. Evaluating the accuracy of DNA stable isotope probing. doi: https://doi.org/10.1101/138719


# INSTALLATION

[[top](#sections)]

## DEPENENCIES

* [pathos](https://github.com/uqfoundation/pathos)
  * [pathos external dependencies](https://github.com/uqfoundation/pathos/tree/master/external)
* [pyfasta](https://pypi.python.org/pypi/pyfasta/)
* [intervaltree](https://github.com/chaimleib/intervaltree)

### Dependency install issues

* scipy libgrfortran issues
  * See https://github.com/ilastik/ilastik-build-conda/issues/17
  

## Installation of SIPSim

### Clone the repo

~~~
git clone https:github.com/nyoungb2/SIPSim.git
cd SIPSim
~~~

### Compile C code; set up paths; add bash completion

~~~
python setup.py build
python setup.py install --prefix=~
echo 'source '`pwd`'/sourceMe' >> ~/.bashrc
~~~


# TUTORIALS

[[top](#sections)]

* [An example with 3 genomes](./ipynb/example/1_dataset.ipynb)
* [Recreating Fig 1 from Lueders et al., 2004](./ipynb/example/Lueders2004.ipynb)


# SIULATION WORKFLOW

[[top](#sections)]

![simulation pipeline](img/simulation_pipeline.png)


# CHANGELOG

[[top](#sections)]


# LICENSE

[[top](#sections)]

* Free software: MIT license
SIPSim
======
SIPSim is a pipeline for simulating data from high resolution 
stable isotope probing (HR-SIP).


INSTALLATION
============

## DEPENENCIES

* [pathos](https://github.com/uqfoundation/pathos)
  * [pathos external dependencies](https://github.com/uqfoundation/pathos/tree/master/external)
* [pyfasta](https://pypi.python.org/pypi/pyfasta/)
* [intervaltree](https://github.com/chaimleib/intervaltree)

## Install of SIPSim

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
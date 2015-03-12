SIPSim
======
SIPSim is a pipeline for simulating data from high resolution 
stable isotope probing (HR-SIP).


INSTALLATION
============

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
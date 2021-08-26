Installation
=========================

Prerequisites
----------------------------


  * Python 3.5 to 3.9
  * numpy, scipy, matplotlib, pyqt5, h5py

Python 3.5-3.7 on Ubuntu
-----------------------------

Prerequisites

 * libgl1-mesa-glx
 * libpython3.5, libpython3.6, or libpython3.7 (to match your Python version)
 * libgfortran3
 * DAE Tools version 1.9.0

Python 3.8/3.9 on Ubuntu
-----------------------------

Prerequisites

 * libgl1-mesa-glx
 * libpython3.8 or libpython3.9 (to match your Python version)
 * libgfortran5
 * libsuperlu5
 * DAE tools version 1.9.1 from github.com/v1kko/daetools

Install via PyPi
-----------------------------

MPET is available on PyPi, the Python Packaging Index and can be installed with :

``pip install mpet``

Install from source
----------------------------

You can also download the source code and install the latest version

 * clone the repository : ``git clone https://bitbucket.org/bazantgroup/mpet``
 * Enter the mpet directory : ``cd mpet``
 * install MPET using pip ``pip install -e .``

Test your installation
---------------------------
 To test your installation make sure to install MPET following :

 ``pip install -e .[test]``

 Then go into the test directory

 ``cd tests``

 and run the following command

 ``pytest --baseDir=ref_outputs --modDir=../bin/workdir/modified compare_tests.py``



Installation
=========================

Prerequisites
----------------------------

  * Python 3.5 to 3.7
  * DAE Tools version 1.9.0

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



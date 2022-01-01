Installation
=========================

Prerequisites
----------------------------


  * Python 3.6 to 3.10
  * numpy, scipy, matplotlib, pyqt5, h5py
  * daetools

MPET on Windows
-----------------------------

MPET on Windows needs to use python 3.6 or 3.7 because daetools on
windows is only available for those versions.

First make sure that you have a correct python version with (ana)conda for
example:


.. code-block:: bash

  conda create -n mpet python=3.7 pip
  conda activate mpet

Then install daetools via PyPi


.. code-block:: bash

  conda activate mpet
  pip install daetools


Then either clone the MPET repository if you want to work on the source code, or
install the mpet package through PyPi (see steps below)


MPET on Linux
-----------------------------

The easiest way to install MPET on Linux is via conda

Install the daetools dependency via conda (in a new environment called mpet)


.. code-block:: bash

  conda create -n mpet -c conda-forge daetools
  conda activate mpet

Then either clone the MPET repository if you want to work on the source code, or
install the mpet package through PyPi (see steps below)


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

 
 Then you can run the tests with

 ``./bin/mpettest.py``

Common Installation Bugs
---------------------------

One of the most common bugs is with QT plugins (it is not acutally a problem with MPET, but with one of the packages that MPET uses). The bug will usually cause plots to not be able to initialize and have the following error message:


.. code-block:: RST

    qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
    This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.
    Available platform plugins are: eglfs, linuxfb, minimal, minimalegl, offscreen, vnc, wayland-egl, wayland, wayland-xcomposite-egl, wayland-xcomposite-glx, webgl, xcb.
    Aborted (core dumped)``

If you get this bug, first, check to see that your X11 server is in use!
If not, try turning on debugging for QT plugins with 

``export QT_DEBUG_PLUGINS=1``. 

Often, the library issue that appears is 


.. code-block:: RST

    Cannot load library /.../lib/python3.7/site-packages/PyQt5/Qt5/plugins/platforms/libqxcb.so:
    (libxcb-xinerama.so.0: cannot open shared object file: No such file or directory)``.

If this is the issue, outside of an virtual environment, install the libxcb-xinerama0 package with

``sudo apt-get install libxcb-xinerama0``.

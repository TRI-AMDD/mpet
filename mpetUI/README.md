cbar# mpetUI -- User Interface for MPET (Multiphase Porous Electrode Theory)

mpetUI is a user interface (UI) developed to run mpet in either GUI mode (using the --gui option),or using the command-line interface by passing --parameters. MPET (Multiphase Porous Electrode Theory) is designed to run simulations of batteries with porous electrodes using porus electrode theory. MPET was developed by Raymond Smith [1]. The mpetUI also uses the VKML library, which was developed by Alex Bartol [2].

[1] Smith, R. B., and Bazant M. Z., Multiphase Porous Electrode Theory, [Journal of the Electrochemical Society] (https://doi.org/10.1149/2.0171711jes), 2017, 164 (11) E3291-E3310, [arXiv pre- print](https://arxiv.org/abs/1702.08432)

[2] Bartol, A., Ely, D. R., Guyer, J., García, R. E. The Virtual Kinetics of Materials Laboratory. 2015.


# Installation

Installation is similar to other Python libraries

The executive summary of steps is:
```
tar -xzvf mpetUI.tar.gz
cd mpetUI
python setup.py install
```

## Prerequisites

The current version of mpetUI has been tested to work on machines running Linux  and Mac OS X. The following 
- [Python (3.6 )](http://www.python.org) .
- Mac OS X users will need to install an [*X11 server*](http://xquartz.macosforge.org/).
- Download and install the [latest release  of MPET](https://bitbucket.org/bazantgroup/mpet/downloads/?tab=tags), (SEE MPET installation instructions.)
- `tkinter` module for python3 should be installed. 

## Procedure

Download `mpetUI.tar.gz` to a suitable folder. 
### 1.Unpack

  Unpack the .tar.gz file.  The usual way is to run `tar -xzvf` on the
  file you want to unpack.  This will create a subdirectory named
  "mpetUI-1.0" in the directory where you run tar.

### 2. Install

Switch to the newly-created directory after unpacking 'mpetUI-1.0/', 

```
cd mpetUI

```
To install mpetUI, run the command 

```
python setup.py install
```
This will install mpetUI in the standard location for Python
extensions on your system. We can also install mpetUI in a different location, like this:

```
python setup.py install --prefix=<prefix>
```
<prefix> is the directory where you want to install. If \<prefix>/bin is not in your Unix command path, you'll need to add it to the PATH environment variable. 

## Running mpetUI

1. After installation, navigate to any folder of your choice which will be the working directory. 
2. It is recommended to run `mpetUI` command once in the working directory before running any mpetUI commands. 
    This creates three files `params_a.cfg`, `params_c.cfg` and, `params_system.cfg` and runs a default mpet simulation. The simulation output is saved in a time-stamped subdirectory called `history`. The data contents of the most recent output will be copied to the directory called `sim_output`. The output directory should contain:
    - the output data (`.mat` file)
    - copies of the input parameters files defining the simulation
    - a copy of the daetools config parameters (e.g. solver tolerances)
    - information about the script used to run the simulation
    - information about the simulation (e.g. run time)
    - processed, dimensional and nondimensional parameters as
      Python-pickled dictionary objects
      
3. In this working directory, you can run a simulation through graphical user intreface mode of mpetUI using the command `mpetUI --gui`, or run in command line interface mode using the command `mpetUI --parameter1=value1 --parameter2=value2 ...` . The summary of all parameters and possible commands can be found by typing `mpetUI --help` or `man mpetUI`.  
4. The simulation output after every simulation is saved to the the `sim_output` folder in the working directory. 
5. The simulation results can be visualized using the `mpetviz --gui`. A summary of commands available to visualization and analysis of data can be obtaing using `mpetviz --help`. 

## Uninstalling

To uninstall mpetUI in the terminal use `pip uninstall mpetUI` command.

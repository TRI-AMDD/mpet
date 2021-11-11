# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.7] - 2021-10-01
### Added
- Example configs added for reproducing the results from classic Newman papers.
- Constant power (CP) charging/discharging implemented.

### Changed
- Refactored code for input file processing with better error checking and a schema to define the structure of config files.
- Regression test suite made more robust by comparing mean error.

### Fixed
- Resolved occasional plot crashes with mpetplot.py.
- Big performance and stability improvements when using noise.


## [0.1.6] - 2021-04-23
### Added
- Implemented hdf5 file output using new 'dataReporter' option, and added an 'hdf5Fast' option which saves fewer variables for smaller file size.
- Coupled ion electron transfer (CIET) added as an intercalation reaction rate.
- Input parameters 'specified_psd_a' and 'specified_psd_c' added for specifying particle radii.
- New '1C_current_density' parameter added for defining a nominal 1C current.
- Example config file added for the LiCoO2/graphite cell from the LIONSIMBA benchmark problem.

### Changed
- Fixed temperature now fully implemented (without heat generation). Activation energy added as a reaction input parameter.
- Code style now conforms to pep8 standards.

### Fixed
- Finite volumes discretization now properly accounts for variable cell size when computing face values. Harmonic mean only used for transport properties at faces.
- Improved stability during initializtion with larger MaxNumItersIC.
- Updated lambda parameter for graphite & LFP with more physical estimates.


## [0.1.5] - 2021-01-19
### Added
- Improved regression test suite based on pytest.
- Docker config file for running MPET in a container.
- Continuous integration implemented for Github and Bitbucket source code repositories.

### Fixed
- Port variables and redundant time variables are no longer saved as output, reducing file size.


## [0.1.4] - 2020-07-27
### Added
- Adds a separate package for mpet.plot.
- Ending condition message indicates why simulation stopped.

### Fixed
- Improved initialization of simulations.
- The prevDir restart option initializes correctly.
- Y-axis limits on electrolyte potential plots chosen automatically.

### Deprecated
- capFrac is now optional and defaults to 1.


## [0.1.3] - 2020-04-30
### Added
- Support for daetools 1.9.0 and Python 3.7.
- Missing plot types restored.

### Fixed
- Initialization of simulations is much more robust.
- More realistic parameter choices for the default params_system.cfg.
- setup.py now installs prerequisite packages.

### Deprecated


## [0.1.2] - 2019-12-19
### Added
- Adds a help flag to mpetrun and mpetplot.
- Implementation of CCsegments and CVsegments properly as discontinuous equations.
- The ramp variable "tramp" can now be set to zero.

### Fixed
- Performance improvements, including using daetools Compute Stack execution model when possible.
- Improved formatting of plots in mpetplot.

### Deprecated
- The variable "tramp" is no longer required and will be removed in future releases.


## [0.1.1] - 2019-07-18
### Added
- Package installer (setup.py) and simplified install instructions.
- Missing reference solutions for regression testing that were no longer available for download.

### Fixed
- Compatibility with DAE Tools 1.8
- mpettest.py now prints all failing variables

### Deprecated
- Some variables (Omega_a,Omega_b,Omega_c, and EvdW) in the [Material] section of electrode config files are no longer required to be defined if they are not being used.


## [0.1.0] - 2017-08-16
### Added
- Last commit from Ray Smith's PhD thesis

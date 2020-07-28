# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


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

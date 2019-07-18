# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.1] - 2019-07-18
### Added
- Package installer (setup.py) and simplified install instructions.
- Missing reference solutions for regression testing that were no longer available for download.

### Fixed
-Compatibility with DAE Tools 1.8
-mpettest.py now prints all failing variables

### Deprecated
-Some variables (Omega_a,Omega_b,Omega_c, and EvdW) in the [Material] section of electrode config files are no longer required to be defined if they are not being used.

## [0.1.0] - 2017-08-16
### Added
- Last commit from Ray Smith's PhD thesis

# pyArtaios

## Introduction

pyArtaios is a simple python interface to control ARTAIOS calculations. Without any knowledge of how to start quantum chemical calculations and to perform the post-processing using Artaios, it yields properties as the transmissions function for a given structure of a junction and a handfull of settings.

```
# import module
import pyartaios as pa

# define settings
settings = {'basis set': 'lanl2dz', 'functional': 'b3lyp'}

# read structure
natoms, atoms, coords = pa.read_xyz('file.xyz')

# initialize class object
transport = artaios(natoms, atoms, coords, settings)

# perform calculations
transport.calculate()

# get results
transmission = transport.get('transmission')
```

## Adding the library to the python path

There are two ways to make the pyArtaios library findable for python on your system. The first way is to add it to the `$PYTHONPATH` variable in the `.bashrc` by adding the line `export PYTHONPATH=/home/<username>/<path to your local git>`. But, and this is crucial: **We have to be careful with the file and functionnames!** If somebody creates a file 'numpy', this would overwrite all `numpy` functionalities, since now two python modules exist with that name.

The second way is to add the link at the beginning of each script by using
```
import sys
sys.path.append(<path to pyartaios library>)
```
## Settings

The default settings are set in the `_default_settings` function of the `pyartaios` class and are overwritten by the settings supplied by the user. The full list of settings can be found there, too. The keywords which cover most of the calculations are given below.

|keyword | default | description|
|--------|---------|-----------|
|basis set | 'lanl2dz' | basis set for the quantum chemical calculation|
|charge | 0 | charge of the molecule/system|
|elow | -8 | lower energy limit for the transmission function in eV|
|estep | 400 | number of steps of the transmissions function|
|eupp | -8 | lower energy limit for the transmission function in eV|
|fermi | -5 | Fermi energy used for transport calculation in eV|
|functional | 'bp86' | exchange-correlation functional|
|method |'dft' | method for quantum chemical calculation|
|multi | 1 | multiplicity|
|n cores | 1 | number of cores for calculation|
|nosymmetry | False | Use of "nosymmetry" keyword in Gaussian|
|qcprog | 'gaussian' | 'gaussian' or 'turbomole'|
|ri | False | switch for density fitting|

The settings are organized as a dictionary, as shown in the example above.

## Functions

With the initialization of the pyartaios class object, all necessary checks will be performed and the consistency of the chosen settings will be checked. If problems occur, the class will raise an error, so the user has to adjust the settings accordingly.

If the settings are not consistent, the `pyartaios.calculate` function will not execute. This function performs all necessary calculations in a subfolder called `artaios`. Existing folders with the same name will be moved to `artaios_backup_<i>`.

The available properties and results of the calculations can bet obtained with `pyartaios.get(<keyword>)`. The keywords can be shown with `transport.results.keys()`, which yields the keys of the internal "results" dictionary.

### List of Functions

The function which are needed for all basic calculations are already given in the example above. This is a full list of all function of the pyArtaios class. Under normal circumstances, they should not be used to perform calculations.

| function name | description |
| ---------| -------- |
| `calculate` | performs all necessary calcualtions |
| `get` | getter function for results |
| `plot_local` | plots the local transmissions |
| `_artaios_transport` | performs the artaios transport calculation |
| `_check` | performs consistency check |
| `_default_settings` | sets default settings |
| `_get_subsys_energies` | routine to perform the subsystem analysis with artaios |
| `_get_local_transmission` | get local transmissions from artaios calculation |
| `_get_mo_energies` | get mo energies from single point calculation |
| `_singlepoint` | performs single point calculation |
| `__init__` | initializes pyArtaios class |

### List of Helperfunctions

pyArtaios contains a bunch of scripts to facilitate the usage of python, see below for a list.

| function name | description |
| ---------| -------- |
| `read_xyz` | reads an xyz file|

### List of accessible variables

Even though for the results a getter function exists (`pyartaios.get('keyword')`), the necessity of adjusting a variable of the class can arise

|variable name| description|
|-------------|------------|
|natoms | number of atoms|
|atoms | list of atom symbols|
|coords | numpy array with cartesian coordinates|
|left | list of atoms for left electrode|
|central | list of atoms for central region|
|right | list of atoms for right electrode|
|results | dictionary with results|

## TODOs

- [x] add turbomole support
- [x] add calculation of subsystem MOs (closed shell)
- [x] add read in of subsystem MO energies (closed shell)
- [ ] add calculation of subsystem MOs (open shell)
- [ ] add read in of subsystem MO energies (open shell)
- [x] add calculation of local transmissions (closed shell)
- [x] add visualization of local transmission (closed shell)
- [ ] add calculation of local transmissions (open shell)
- [ ] add visualization of local transmission (open shell)
- [x] update readme
- [ ] add visualization tools
- [ ] add jgreen
- [ ] add switch for existing QC calculation, so existing calculations can be reused
- [ ] read MO energies for turbomole, cclib a problem with two-character elements, see https://github.com/cclib/cclib/pull/780
- [ ] find solution for printing information/status message - is that in order with python conventions?
- [ ] add config/path info for QC programs (gaussian, turbomole, artaios), so that this no only runs on our systems

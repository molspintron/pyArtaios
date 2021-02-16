import sys
# adjust path
sys.path.append("/home/deffner/code/pyartaios")
import pickle
import numpy as np
import pyArtaios as pa

# read in junction structure
natoms, atoms, coords = pa.read_xyz('junction.xyz')

# define settings
settings = {'basis set': 'lanl2dz', 'qcprog': 'gaussian', 'method' : 'dft', 'functional' : 'b3lyp', 'n cores': 10}

# create object
transport = pa.artaios(natoms, atoms, coords, settings)

# perform calculations
transport.calculate()

# save object
pickle.dump(transport, open("results.p", "wb"))

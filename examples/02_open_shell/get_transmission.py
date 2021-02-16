import sys
# adjust path
sys.path.append("/home/deffner/code/pyartaios")
import pickle
import numpy as np
import pyArtaios as pa

# example for simple open shell transport calculation with pyartaios

# load structure
natoms, atoms, coords = pa.read_xyz('junction.xyz')

# define settings
settings = {'charge': 1, 'multi': 2, 'spin': 2, 'basis set': 'lanl2dz', 'qcprog': 'gaussian', 'method' : 'dft', 'functional' : 'b3lyp', 'n cores': 10}

# create objext
transport = pa.artaios(natoms, atoms, coords, settings)

# perform calculations
transport.calculate()

# save class
pickle.dump(transport, open("results.p", "wb"))

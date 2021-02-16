import sys
# adjust path
sys.path.append("/home/deffner/code/pyartaios")
import pickle
import numpy as np
import pyArtaios as pa

# read in junction structure
natoms, atoms, coords = pa.read_xyz('junction.xyz')

# define settings
settings = {'basis set': 'def2-SVP', 'qcprog': 'turbomole', 'method' : 'dft', 'functional' : 'b3-lyp'}

# create object
transport = pa.artaios(natoms, atoms, coords, settings)

# perform calculations
transport.calculate()

print(transport.results.keys())

# print transmission
print(transport.get('transmission'))

# save object
#pickle.dump(transport, open("results.p", "wb"))

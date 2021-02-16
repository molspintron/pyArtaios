import sys
# adjust path
sys.path.append("/home/deffner/code/pyartaios")
import numpy as np
import pyArtaios as pa
import matplotlib.pyplot as plt

# read in junction structure
natoms, atoms, coords = pa.read_xyz('junction.xyz')

# define settings
settings = {'basis set': 'def2-SVP', 'qcprog': 'turbomole', 'method' : 'dft', 'functional' : 'bp86', 'ri': True, 'charge': -1, 'multi' : 2}

# create object
transport = pa.artaios(natoms, atoms, coords, settings)

# perform calculations
transport.calculate()

# plot transmission
plt.figure()
plt.semilogy(transport.get('transmission')[:, 0], transport.get('transmission')[:, 1])
plt.savefig('transmission.pdf')


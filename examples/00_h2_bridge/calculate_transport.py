import sys
# adjust path
sys.path.append("/home/deffner/code/pyArtaios")
import pyArtaios as pa
import matplotlib.pyplot as plt
import numpy as np

transmissions = []
mo_energies = []
subsys_energies = []

# read in junction structure
natoms, atoms, coords = pa.read_xyz('input.xyz')
# pa.write_coord(natoms, atoms, coords)

# define settings
settings = {'charge' : -1, 'multi' : 2, 'spin' : 2, 'basis set': 'def2svp', 'qcprog': 'gaussian', 'method' : 'dft', 'functional' : 'bp86', 'n cores': 1, 'elow': -12, 'eupp': 2, 'esteps' : 2}

# create object
transport = pa.artaios(natoms, atoms, coords, settings)

transport.left    = [2]
transport.right   = [3]
transport.central = [0, 1]

# perform calculations
transport.calculate()

transmissions.append(transport.get('transmission'))
mo_energies.append(transport.get('mo energies'))
subsys_energies.append(transport.get('subsys energies'))

transport.plot_local(-5)

'''

# define settings for turbomole calculation
transport.settings['basis set'] = 'def2-SVP'
transport.settings['qcprog']    = 'turbomole'

# perform calculations
transport.calculate()

transmissions.append(transport.get('transmission'))
#mo_energies.append(transport.get('mo energies'))
subsys_energies.append(transport.get('subsys energies'))

# plotting
plt.figure(1)

plt.plot(transmissions[0][:, 0], transmissions[0][:, 1], label = 'gaussian')
plt.plot(transmissions[1][:, 0], transmissions[1][:, 1], label = 'turbomole')

plt.scatter(mo_energies[0][0][:, 0], [1 for x in range(len(mo_energies[0][0][:, 0]))], label = 'mo energies, gaussian')
#plt.scatter(mo_energies[1][0][:, 0], [0.9 for x in range(len(mo_energies[1][0][:, 0]))], label = 'mo energies, turbomole')

plt.scatter(subsys_energies[0], [0.8 for x in range(len(subsys_energies[0]))], label = 'subsys energies, gaussian')
plt.scatter(subsys_energies[1], [0.7 for x in range(len(subsys_energies[1]))], label = 'subsys energies, turbomole')

plt.xlabel('Energy [eV]')
plt.xlim([-10, 0])
plt.legend()
plt.grid()
#plt.show()
plt.savefig('results.pdf')
plt.show()

#print(mo_energies[0][0][:, 0])
#print(mo_energies[1][0][:, 0])

#print(subsys_energies[0])
#print(subsys_energies[1])
'''


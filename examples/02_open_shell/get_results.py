import sys
# adjust path
sys.path.append("/home/deffner/code/pyartaios")
import pickle
import numpy as np
import pyArtaios as pa
import matplotlib.pyplot as plt

# example to load an old calculation and plot the data

# load object
transport = pickle.load(open('results.p', 'rb'))

# print available keys/results
print('available results:')
print(transport.results.keys())

# transport.get('mo energies')

# do some plotting
fig, axes = plt.subplots(2, 1, sharex = True, figsize = (6,4))

axes[0].semilogy(transport.results['transmission'][:, 0],
                 transport.results['transmission'][:, 1],
                 label = 'alpha')

axes[1].semilogy(transport.results['transmission'][:, 0],
                 transport.results['transmission'][:, 2],
                 label = 'beta')

axes[0].grid()
axes[0].set_xlabel('Energy [eV]')
axes[0].set_ylabel('Transmission [a.u.]')

axes[1].grid()
axes[1].set_xlabel('Energy [eV]')
axes[1].set_ylabel('Transmission [a.u.]')


plt.show()

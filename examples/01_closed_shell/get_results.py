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

transport.plot_local()

#print(transport.get('mo energies'))

'''
# do some plotting
fig, axes = plt.subplots(1, 1, figsize = (6,4))

axes.semilogy(transport.results['transmission'][:, 0],
                 transport.results['transmission'][:, 1],
                 label = 'transmission')

axes.grid()
axes.set_xlabel('Energy [eV]')
axes.set_ylabel('Transmission [a.u.]')

plt.show()
'''

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


# Import relevant modules
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants, optimize
from pyproj import Geod
import pickle
import pandas as pd
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize

# Load the node_config dictionary that contains information about the nodes
# (receivers) position and name
node_config = np.load("node_config_fixed.npy", allow_pickle=True).item()

# Load in waveforms
waveforms = np.load("exported_waveforms.npy", allow_pickle=True)

# For each waveform the following variables are available
# envelope      : The magnitude of the waveform envelope with time
# FT            : The Fast fourier transform of samples
# FT_recon      : The Fast fourire transform of the reconstructed time series
# node_unique_id: A 17 character unique identifier for a node
# peak_index    : The index in the raw data stream that the peak in envelope resides
# prod_time     : The time the waveform was extracted
# raw_data      : The raw time series the waveform was extracted from that has been high pass filtered at 1kHz
# recon_data    : The reconstructed complex time series the waveform was extracted from same as reconstructed samples
# samples       : The real part of extracted waveform
# site_name     : Name of site waveform was recorded at
# t1            : The sum of samples ^ 2
# tdv_penalty   : A time difference variance penalty 
# timestamp     : np.datetime64 object which contains the timestamp at index 256 of samples, envelope, raw_data and recon_data

# The sampling resolution of LEELA is 9.142 us


# Create a dictionary node_code with the key being the site and the value being the node code
node_code = {}
for key in node_config.keys():
    node_code[node_config[key]['Site']] = key


def TATD(waveforms, guess):
    '''
    Finds the Theoretical Arrival Time Difference for each node by calculating distance from a reference node 
    to each node using pyproj's Geod function. Adds these TATDs to one list.
    '''
    
    # Define the Geod function
    geod = Geod(ellps='WGS84')
    
    # Calculate the distance from the reference node to the event
    _, _, dist_GA = geod.inv(waveforms[0]['lon'], waveforms[0]['lat'], guess[0], guess[1])

    # Open list of distances from non-reference nodes to event
    dist_GBs = []
    # Append the distances from non-reference nodes to event to dist_GBs
    for i in range(1, 7):
        _, _, dist_GB = geod.inv(waveforms[i]['lon'], waveforms[i]['lat'], guess[0], guess[1])
        dist_GBs.append(dist_GB)
    
    # Create list of TATDs
    TATDs = []

    # Define the phase velocity
    vp = 2.99792458e8 * 1.00425
    
    # Calculate the TATD for each node and append to TATDs
    for i in range(0, 6):
        TATD = abs(dist_GBs[i] / vp) - abs(dist_GA / vp)
        TATDs.append(TATD)
    return TATDs


def RES(guess, waveforms, OATDs):
    '''
    Function to be minimised. Calculates the residual of how close the guess is 
    to the actual lightning strike location.
    '''
    # Use TATD function to find TATDs
    TATDs = TATD(waveforms, guess)
    # Find variance
    vari = np.var(OATDs)
    # Determine N from amount of values in OATDs
    N = len(OATDs)

    # Calculate the value of RES
    sumATD = 0
    for i in range(0, 5):
        sumATD = sumATD + (TATDs[i] - OATDs[i]) ** 2 / vari
    RES = np.sqrt((1 / (N - 2) * sumATD))
    
    # Print RES value
    print(RES)
    
    # Return RES value for each minimisation
    return RES


#Define the OATDs list
OATDs = []

# Add the latitude and longitude to the waveforms dictionary from node_config.
for i in range(0,7):
    waveforms[i]['lat'] = node_config[node_code[waveforms[i]['site_name']]]['Position']['lat']
    waveforms[i]['lon'] = node_config[node_code[waveforms[i]['site_name']]]['Position']['lon']

# Find the observed time difference by taking the value of the difference between the two timestamps.
for i in range(1, 7):
    OATDs.append((waveforms[i]['timestamp'] - waveforms[0]['timestamp']).astype(float)/1e9)
    
# Print the OATDs
for i in range(0, 6):
    print(OATDs[i])


# Initalise guess in London (Long=-0.118092 , Lat=51.509865)
guess = np.array([-0.118092, 51.509865])

# Call the minimise function from scipy, minimising the RES function until it is sufficiently small.
# (Can change the method of minimisation)
result = minimize(RES, guess, args=(waveforms, OATDs), method='Nelder-Mead')

# Print the result
print(result)

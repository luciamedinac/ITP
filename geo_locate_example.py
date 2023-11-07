import numpy as np
import matplotlib.pyplot as plt
from scipy import constants, optimize
from pyproj import Geod
import pickle
import pandas as pd
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize
# Load the node config dictionary this contains information about the nodes
# (receivers) position and name
node_config = np.load("node_config_fixed.npy", allow_pickle=True).item()

# Now load in waveforms
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



node_code = {}
for key in node_config.keys():
    node_code[node_config[key]['Site']] = key



'''
This function finds the Theoretical Arrival Time Difference. This is done through calculating
the distances between the guess and each node. distGA is the distance between the guess and 
the reference node. distGBs is a list of all the distances between the guess and all other nodes.
This is done using the geod function in pyproj which calculates distances between 2 points taking
into account the geodisic lines of the earth (WGS84 being the best representation of the earth).
Theoretical Arrival Time Difference is calculated using both distGA and distGB as well as the phase
velocity (vp) which will be a little over the speed of light. The list TATD is outputted to be utilised
in the RES function.
'''

def TATD(waveforms, guess):
    dist_GBs = []
    geod = Geod(ellps='WGS84')
    _, _, dist_GA = geod.inv(waveforms[0]['lon'], waveforms[0]['lat'], guess[0], guess[1])

    for i in range(1, 7):
        _, _, dist_GB = geod.inv(waveforms[i]['lon'], waveforms[i]['lat'], guess[0], guess[1])
        dist_GBs.append(dist_GB)
        # print(waveforms[x]['site_name'], waveforms[x]['lat'], waveforms[x]['lon'], dist_GB)
    TATDs = []

    vp = 2.99792458e8 * 1.00425

    for i in range(0, 6):
        TATD = abs(dist_GBs[i] / vp) - abs(dist_GA / vp)
        TATDs.append(TATD)
    return TATDs



'''
RES function takes in an intial guess of where the lightening point will be,
the waveforms dictionary containing data necessary to find the TATD (see TATD
function) and OATDs is the list of OATDs (see above). The guess will change after
each iteration of the minimisation function is complete. The point of this
fuction is to calculate the residual (res) of how close the guess is to the actual
lightening strike location. The smaller res gets, the closer the guess is to the
location of the lightening strike.
'''

def RES(guess, waveforms, OATDs):
    # Finds TATDs
    TATDs = TATD(waveforms, guess)
    # Finds variance
    vari = np.var(OATDs)
    # determines N from amount of values in OATDs
    N = len(OATDs)

    # Calculate the value of RES
    sumATD = 0
    for i in range(0, 5):
        sumATD = sumATD + (TATDs[i] - OATDs[i]) ** 2 / vari
    RES = np.sqrt((1 / (N - 2) * sumATD))
    print(RES)
    return RES


'''
This adds longitude and and latitude to the waveforms dictionary for ease to only need to pass 
through one dictionary to all functions. A list is also defined to hold all observed arrival
time differences. This is the time difference between our reference node (waveforms[0]) and 
another arbitrary node. This is obviously doesn't vary so it makes sense to not attach this
to a function. This reference node is consistent when calculating distances between 
nodes as well as the theoretical arrival time difference.
'''
OATDs = []
for i in range(0,7):
    waveforms[i]['lat'] = node_config[node_code[waveforms[i]['site_name']]]['Position']['lat']
    waveforms[i]['lon'] = node_config[node_code[waveforms[i]['site_name']]]['Position']['lon']
for i in range(1, 7):
    OATDs.append((waveforms[i]['timestamp'] - waveforms[0]['timestamp']).astype(float)/1e9)
for i in range(0, 6):
    print(OATDs[i])


#Initalise guess in london
'''
guess = [0, 0]
guess[0] = -0.118092 #longitude
guess[1] = 51.509865 #latitude
'''
guess = np.array([-0.118092, 51.509865])
'''
This calls the minimise function in scipy which will minimise RES function until RES is sufficenitly 
small. Minimise allows the use of a variety of diiferent methods to minimise but we are using the Nelder-Mead
method.
'''
result = minimize(RES, guess, args=(waveforms, OATDs), method='Nelder-Mead')

#blametest

print(result)


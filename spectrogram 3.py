# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 21:04:58 2024

@author: archi
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants, optimize, signal
from pyproj import Geod
from scipy.optimize import minimize
import os
import glob
import math
from matplotlib.colors import LogNorm
from scipy.signal import savgol_filter

def filter_t_series(t, samples, sos, thresh=np.timedelta64(20000, "ns")):
    """
    Filters time series accounting for discontinuities

    Parameters
    ----------
    t : np.array of np.datetime64
        time array
    samples : np.array of float
        samples
    sos : filter object
        sso filter object
    thresh : np.timedelta64, optional
        The minimal time gap needed for a discontinuity.
        The default is np.timedelta64(20000, "ns").

    Returns
    -------
    np.array float
        filtered array.

    """
    if len(t) > 0:
        dt = np.diff(t)
        new_samples = []
        current_idx = 0
        l1 = np.where(dt > thresh)[0]
        l1 = np.hstack((l1, len(samples) - 1))
        l1 += 1
    
        for end_idx in l1:
            new_sample_tmp = signal.sosfilt(sos,
                                            samples[current_idx:end_idx])
    
            current_idx = end_idx
    
            new_samples += list(new_sample_tmp)
    
        return new_samples
    else:
        return []


def atddict(long, lat, ref):
    '''
   Finds the Arrival Time Difference between a reference node and each other node
    based on an inputted longitude and latitude of a meteorite. Returns a dictionary
    with Site as key and ATD as value.
    '''
    
    #start atd dictionary
    atd = {}
    
    #define the speed of light in m/s (can change if we find it is different)
    c = 299792458
    
    # Loop between each node in node_config
    for node in node_config:
        geod = Geod(ellps='WGS84')
        
        # Use the geod function to find the dist to the for the node
        _, _, dist = geod.inv(node_config[node]['Position']['lon'], node_config[node]['Position']['lat'], long, lat)

        # Use the geod function to find the dist to the reference node
        _, _, ref_dist = geod.inv(node_config[ref]['Position']['lon'], node_config[ref]['Position']['lat'], long, lat)

        # Calculate the time difference
        timediff = (dist / c) - (ref_dist / c)

        # Add the ATD to the ATD dictionary in ns
        atd[node_config[node]['Site']] = np.timedelta64(int(timediff * 1e9),'ns')
        
    return atd

# Load the node config dictionary this contains information about the nodes
# (receivers) position and name
prog_start = np.datetime64("now")
node_config = np.load("node_config_fixed.npy", allow_pickle=True).item()
node_code = {}
datadict = {key: {"time": [],
                  "waveform": []} for key in node_config}
# LEELA sampling dt
dt = np.timedelta64(int(1e9 / 109375), "ns")
t = np.arange(0, 1024, 1) * dt

# Set parameter atds to be a dictionary with Watnall as the reference node
atds = atddict(-5.83,51.09,'569218Q0B001D0029')
window_dt = np.timedelta64(1000, "ms")
start_time = np.datetime64("2022-07-09T01:42:51.5")
end_time = np.datetime64("2022-07-09T01:42:52.5")

# Create filter here so only processed once:
lower_f = 2000
upper_f = 18000
sos = signal.butter(6, [lower_f,upper_f], 'bp', fs=1/9142e-9, output='sos')

for key in node_config.keys():

    # File path to data files (if have hard drive then: /Volumes/WD_ITP/Bath_VLF_2023_decoded)
    folder_path = 'D:/bath_students_2023/bath_students_2023/13/'
    # Create list of file paths for each file in the folder
    file_list = glob.glob(folder_path + "*" + key + '*20220709*.npy')

    # Code to create a list of all the data for each of the files
    for file_path in file_list:
        #print("Reading: " + file_path)
        chunks = np.load(file_path, allow_pickle=True).item()
        for chunk in chunks:
            datadict[key]["time"].append(chunks[chunk]["starttime"])
            datadict[key]["waveform"].append(chunks[chunk]["wvfmdata"])

print("data loaded")
#569218Q0B0025002A for camborne
#569218Q0B001D0029 for watnall
node_id = '569218Q0B0025002A'

datadict[node_id]["time"] = np.array(datadict[node_id]["time"],
                                              dtype=np.datetime64)
l1 = np.where((datadict[node_id]["time"] >= start_time) &
                          (datadict[node_id]["time"] <= end_time))[0]
samples = np.array([])
time = np.array([], dtype=np.datetime64)

#lower_f = 3000
#upper_f = 10000
sos = signal.butter(10, [lower_f,upper_f], 'bp', fs=1/9142e-9, output='sos')

for i in l1:
    samples = np.hstack((samples, datadict[node_id]["waveform"][i]))
    time = np.hstack((time, datadict[node_id]["time"][i] + t))

    filtered_samples = np.array(filter_t_series(time,samples,sos))
    


#smoothed_samples = savgol_filter(filtered_samples,window_length=15,polyorder=2)
 
# Matplotlib.pyplot.specgram() function to
# generate spectrogram
fig,ax = plt.subplots(figsize=(20,10))
spec = ax.specgram(samples, Fs=1/9.142e-6, cmap ='magma',vmin=-50, vmax=0)


ax.set_ylim(0, 18000)  # Set the y-axis limit from 0 to 18 kHz
yticks = np.arange(0, 18001, 2000)  # Define y-axis tick positions in Hz
yticks_khz = yticks / 1000  # Convert Hz to kHz for tick labels
ax.set_yticks(yticks)  # Set y-axis tick positions
ax.set_yticklabels(['{:.0f}'.format(y) for y in yticks_khz], fontsize=16)
 
# Set the title of the plot, xlabel and ylabel
# and display using show() function
plt.title(str(node_config[node_id]['Site']), fontsize=18)
cbar = fig.colorbar(spec[3], ax=ax)

# Set colorbar label
cbar.set_label('Intensity', fontsize=18)
ax.set_xlabel('Time', fontsize=18)
ax.set_ylabel("Frequency / kHz", fontsize=18)
ax.set_ylim(2000,18000)
plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 09:56:30 2024

@author: tjr2099
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import constants, optimize, signal
from pyproj import Geod
from scipy.optimize import minimize
import os
import glob

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
            new_sample_tmp = signal.sosfiltfilt(sos,
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

node_config = np.load("node_config_fixed.npy", allow_pickle=True).item()

# LEELA sampling dt
dt = np.timedelta64(int(1e9 / 109375), "ns")
t = np.arange(0, 1024, 1) * dt

pf = pd.read_excel('/Users/tjr2099/Library/CloudStorage/OneDrive-UniversityofBath/Physics w Astro Year 4/Industry Team Project/fireball_data/potential_fireballs.xlsx')
pf = pf.set_index('num')

for n in pf['num']:
    
    save_to = '/Users/tjr2099/Library/CloudStorage/OneDrive-UniversityofBath/Physics w Astro Year 4/Industry Team Project/fireballs/'
    os.mkdir(save_to + str(n)) 
    
    print(n)

    # Load the node config dictionary this contains information about the nodes
    # (receivers) position and name
    prog_start = np.datetime64("now")
    
    datadict = {key: {"time": [],
                      "waveform": []} for key in node_config}
   
    
    # Set parameter atds to be a dictionary with Watnall as the reference node
    atds = atddict(pf['long'][n],pf['lat'][n],'569218Q0B001D0029')
    window_dt = np.timedelta64(1000, "ms")
    start_time = np.datetime64(str(pf['date'][n].date())+'T'+str(pf['start_time'][n]))
    end_time = np.datetime64(str(pf['date'][n].date())+'T'+str(pf['end_time'][n]))
    
    # Create filter here so only processed once:
    lower_f = 5000
    upper_f = 12000
    sos = signal.butter(10, [lower_f,upper_f], 'bp', fs=1/9142e-9, output='sos')
    
    for key in node_config.keys():
    
        # File path to data files (if have hard drive then: /Volumes/WD_ITP/Bath_VLF_2023_decoded)
        folder_path = '/Users/tjr2099/Library/CloudStorage/OneDrive-UniversityofBath/Physics w Astro Year 4/Industry Team Project/fireball_data/'+str(n)+'/'
        # Create list of file paths for each file in the folder
        file_list = glob.glob(folder_path + "*" + key + '*.npy')
    
        # Code to create a list of all the data for each of the files
        for file_path in file_list:
            #print("Reading: " + file_path)
            chunks = np.load(file_path, allow_pickle=True).item()
            for chunk in chunks:
                datadict[key]["time"].append(chunks[chunk]["starttime"])
                datadict[key]["waveform"].append(chunks[chunk]["wvfmdata"])
    
    print("data loaded")
    
    study_times = np.arange(start_time, end_time, window_dt)
    # Create figure
    for study_time in study_times:
        study_time_end = study_time + window_dt
        fig, axs = plt.subplots(len(datadict)- 4,
                                figsize=(15, 20),
                                sharex=True,
                                constrained_layout=True)
        plt.suptitle("Fireball "+str(n), fontsize=14)
        
        ii = 0
        #print("Processing: " + study_time.astype(str))
        for site in datadict:
            if (node_config[site]["Site"] != 'Gibraltar') and \
                    (node_config[site]["Site"] != 'Valentia') and \
                    (node_config[site]["Site"] != 'Keflavik') and \
                    (node_config[site]["Site"] != 'Herstmonceux'):
    
                datadict[site]["time"] = np.array(datadict[site]["time"],
                                                  dtype=np.datetime64)
                l1 = np.where((datadict[site]["time"] >= study_time) &
                              (datadict[site]["time"] <= study_time_end))[0]
    
                samples = np.array([])
                time = np.array([], dtype=np.datetime64)
    
                for i in l1:
                    samples = np.hstack((samples, datadict[site]["waveform"][i]))
                    time = np.hstack((time, datadict[site]["time"][i] + t))
    
                filtered_samples = np.array(filter_t_series(time,
                                                            samples,
                                                            sos))
    
                # Plot the time_inc_atd against cleaned_waveform
                # You'll need to double check if this next line should be
                # time + atds or time - atds
                axs[ii].plot(time - atds[node_config[site]["Site"]],
                             filtered_samples)
                axs[ii].set_title(node_config[site]["Site"])
    
                # Set x and y limits to zoom in
                axs[ii].set_xlim(study_time, study_time_end)
                axs[ii].grid()
                if len(filtered_samples) > 0:
                    axs[ii].set_ylim(np.min(filtered_samples) - 25,
                            np.max(filtered_samples) + 25)
                #axs[ii].set_yscale('log')
                ii += 1
    
        # Save the figure
        # Convert study_time to a windows happy filename (i.e no :)
        filename_stub = study_time.astype(str).replace(":", "")
        filename_stub = filename_stub.replace(".", "_")
        
        save_as = save_to + str(n) + '/' + filename_stub + ".png"
        
        fig.savefig(save_as)
        
        plt.show()
        
        plt.close(fig)
    
    prog_end = np.datetime64("now")
    print("Seconds 2 run: " + (prog_end - prog_start).astype(str))

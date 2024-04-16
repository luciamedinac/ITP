#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:06:01 2024

@author: tjr2099
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import constants, optimize, signal
from pyproj import Geod
from scipy.optimize import minimize
import os
import datetime
import glob

# Import data into two lists
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

def atddict(long, lat):
    '''
   Finds the Arrival Time Difference between a reference node and each other node
    based on an inputted longitude and latitude of a meteorite. Returns a dictionary
    with Site as key and ATD as value.
    '''
    
    #start atd dictionary
    atd = {}
    
    #define the speed of light in m/s (can change if we find it is different)
    c = 299792458
    
    dists = {}
    
    # Loop between each node in node_config
    for node in node_config:
        geod = Geod(ellps='WGS84')
        
        # Use the geod function to find the dist to the for the node
        _, _, dist = geod.inv(node_config[node]['Position']['lon'], node_config[node]['Position']['lat'], long, lat)
        
        
        dists[node] = dist
        
        # Add the ATD to the ATD dictionary in ns
        atd[node] = np.timedelta64(int((dist/c) * 1e9),'ns')
    
    
    
    return atd, [key for key in dists if dists[key] == min(dists.values())][0]

node_config = np.load("node_config_fixed.npy", allow_pickle=True).item()

pf = pd.read_excel('/Users/tjr2099/Library/CloudStorage/OneDrive-UniversityofBath/Physics w Astro Year 4/Industry Team Project/fireball_data/potential_fireballs.xlsx')
pf = pf.set_index('num')
pf.drop(32, inplace=True)

# Create filter here so only processed once:
lower_f = 2000
upper_f = 18000
sos = signal.butter(10, [lower_f,upper_f], 'bp', fs=1/9142e-9, output='sos')

dt = np.timedelta64(int(1e9 / 109375), "ns")
t = np.arange(0, 1024, 1) * dt

for n in pf.index:
    
    print(f"Starting on fireball {n}")
    
    atds = atddict(pf['long'][n],pf['lat'][n])[0]
    
    close_node = atddict(pf['long'][n],pf['lat'][n])[1]

    # Find the indexes within t_m + t_ref +/- 6 seconds
    datadict = {key: {"time": [],
                      "waveform": []} for key in node_config}

    for node in node_config.keys():
            
        # File path to data files (if have hard drive then: /Volumes/WD_ITP/Bath_VLF_2023_decoded)
        folder_path = '/Volumes/WD_ITP/bath_students_2023/bath_students_2023/'+str(n)+'/'
        # Create list of file paths for each file in the folder
        file_list = glob.glob(folder_path + "*" + node + '*.npy')    
        
        
        
        #Code to create a list of all the data for each of the files
        for file_path in file_list:
            #print("Reading: " + file_path)
            chunks = np.load(file_path, allow_pickle=True).item()
            
            for chunk in chunks:
                
                datadict[node]["time"].append(chunks[chunk]["starttime"]+t)
                datadict[node]["waveform"].append(chunks[chunk]["wvfmdata"])
                
        datadict[node]["time"] = np.array(datadict[node]["time"], dtype=np.datetime64)
        datadict[node]['waveform'] = np.array(datadict[node]["waveform"])   
    
        datadict[node]['time'] = datadict[node]['time'].flatten()
        datadict[node]['waveform'] = datadict[node]['waveform'].flatten()
                
        node_config[node]['start_time'] = (np.datetime64(str(pf['date'][n].date())+'T'+str(pf['ref_time'][n])) 
                  + atds[node] 
                  - np.timedelta64(6, "s"))
        
        node_config[node]['end_time'] = (np.datetime64(str(pf['date'][n].date())+'T'+str(pf['ref_time'][n])) 
                  + atds[node] 
                  + np.timedelta64(6, "s"))
        
        l1 = np.where((datadict[node]["time"] >= node_config[node]['start_time']) &
                      (datadict[node]["time"] <= node_config[node]['end_time']))
    
        datadict[node]['time'] = datadict[node]['time'][l1]
        datadict[node]['waveform'] = datadict[node]['waveform'][l1]
        
        sorted_indices = np.argsort(datadict[node]['time'])
        
        # Arrange the arrays based on sorted_indices
        datadict[node]['time'] = datadict[node]['time'][sorted_indices]
        datadict[node]['waveform'] = datadict[node]['waveform'][sorted_indices]
        
        if len(datadict[node]['waveform']) > 0:
            datadict[node]['waveform'] = signal.sosfiltfilt(sos, datadict[node]['waveform'])
            
        #print(node_config[node]['Site'],node)
        #print(len(datadict[node]['waveform']))
    
    # Check they are the same lengths
    
    # Find m by dividing 25ms by dt
    
    m = int(np.timedelta64(25, "ms")/dt)
    
    
    
    
    for node in node_config:
        
        print(str(node_config[node]['Site']+' starting'))
        
        cc = []
        mdt = []
        
        compare_site = ''
        
        if (len(datadict[node]['waveform']) <= 1312501) and \
                (len(datadict[node]['waveform']) >= 1312499):
            
        
            if (len(datadict[close_node]['waveform']) <= 1312501) and \
                (len(datadict[close_node]['waveform']) >= 1312499) and \
                    node != close_node:
                
                    
                compare_site = node_config[close_node]['Site']
                
                num_zeros = np.abs(datadict[node]['waveform'].shape[0] - datadict[close_node]['waveform'].shape[0])
                
                if len(datadict[node]['waveform']) < len(datadict[close_node]['waveform']):
                    datadict[node]['waveform'] = np.concatenate((datadict[node]['waveform'],np.zeros(num_zeros)))

                    
                if len(datadict[node]['waveform']) > len(datadict[close_node]['waveform']):
                    datadict[close_node]['waveform'] = np.concatenate((datadict[close_node]['waveform'],np.zeros(num_zeros)))

                
                iii=0
                N = len(datadict[node]['waveform'])

                for i in range(-m,m):
                    sum_n = 0
    
                    temp = np.roll(datadict[close_node]['waveform'], i)
        
                    sum_n = np.sum(temp[0:N]*datadict[node]['waveform'][0:N])
                
                    cc.append(sum_n)
                    mdt.append(i*dt)
                    
                peak = mdt[cc.index(max(cc))]/np.timedelta64(1,'ms')

                if np.abs(peak) < 1:
                    print(f'Hold on a second... The peak in the correlation is at {peak} microseconds for fireball {n} at site{node_config[node]["Site"]}')         

                fig,ax = plt.subplots(1,1,figsize=(8, 6))
                 
                ax.plot(mdt/np.timedelta64(1, 'ms'),cc)
                ax.set_xlabel('Shift [ms]')
                ax.set_ylabel('Signal Strength [V/m]')
                ax.set_title(f'Fireball {n} Cross Correlation: {str(compare_site)} against {str(node_config[node]["Site"])}')
                ax.grid()
                plt.show()
                
                save_to = '/Users/tjr2099/Library/CloudStorage/OneDrive-UniversityofBath/Physics w Astro Year 4/Industry Team Project/fireballs/'+str(n)
                
                if not os.path.exists(save_to + '/cross_correlations'):
                    os.makedirs(save_to + '/cross_correlations')
                
                filename_stub = str(n)+'_cc_'+str(compare_site)+'v'+str(node_config[node]['Site'])
                
                save_as = filename_stub + ".png"
                
                fig.savefig(save_to + '/cross_correlations/' + save_as)
                
                plt.show()
                
                plt.close()
            
            elif (len(datadict['569218Q0B0025002A']['waveform']) <= 1312501) and \
                (len(datadict['569218Q0B0025002A']['waveform']) >= 1312499) and \
                    node != '569218Q0B0025002A' and node != close_node:
    
                compare_site = 'Camborne'
    
                num_zeros = np.abs(datadict[node]['waveform'].shape[0] - datadict['569218Q0B0025002A']['waveform'].shape[0])
            
                if len(datadict[node]['waveform']) < len(datadict['569218Q0B0025002A']['waveform']):
                    datadict[node]['waveform'] = np.concatenate((datadict[node]['waveform'],np.zeros(num_zeros)))
                    
                    
                if len(datadict[node]['waveform']) > len(datadict['569218Q0B0025002A']['waveform']):
                    datadict['569218Q0B0025002A']['waveform'] = np.concatenate((datadict['569218Q0B0025002A']['waveform'],np.zeros(num_zeros)))

                
                iii=0
                N = len(datadict[node]['waveform'])

                for i in range(-m,m):
                    sum_n = 0
    
                    temp = np.roll(datadict['569218Q0B0025002A']['waveform'], i)
        
                    sum_n = np.sum(temp[0:N]*datadict[node]['waveform'][0:N])
                
                    cc.append(sum_n)
                    mdt.append(i*dt)


                peak = mdt[cc.index(max(cc))]/np.timedelta64(1,'ms')
    
                if np.abs(peak) < 1:
                    print(f'Hold on a second... The peak in the correlation is at {peak} microseconds for fireball {n} at site{node_config[node]["Site"]}')         
                
                
                
                fig,ax = plt.subplots(1,1,figsize=(8, 6))
                     
                ax.plot(mdt/np.timedelta64(1, 'ms'),cc)
                ax.set_xlabel('Shift [ms]')
                ax.set_ylabel('Signal Strength [V/m]')
                ax.set_title(f'Fireball {n} Cross Correlation: {str(compare_site)} against {str(node_config[node]["Site"])}')
                ax.grid()
                plt.show()
                
                save_to = '/Users/tjr2099/Library/CloudStorage/OneDrive-UniversityofBath/Physics w Astro Year 4/Industry Team Project/fireballs/'+str(n)
                
                if not os.path.exists(save_to + '/cross_correlations'):
                    os.makedirs(save_to + '/cross_correlations')
                
                filename_stub = str(n)+'_cc_'+str(compare_site)+'v'+str(node_config[node]['Site'])
                
                save_as = filename_stub + ".png"
                
                fig.savefig(save_to + '/cross_correlations/' + save_as)
                
                plt.show()
                
                plt.close()
                
                
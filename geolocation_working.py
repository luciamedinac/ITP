#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:19:09 2024

@author: tjr2099
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants, optimize, signal
from pyproj import Geod
import glob
import pandas as pd
import math

from scipy.signal import find_peaks
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize

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
    
    #print(long,lat)
    
    # Loop between each node in node_config
    for node in node_config:
        
        if node_config[node]['Site'] not in avoid_list:
        
            #print(long,lat)
            
            geod = Geod(ellps='WGS84')
            
            # Use the geod function to find the dist to the for the node
            _, _, dist = geod.inv(node_config[node]['Position']['lon'], node_config[node]['Position']['lat'], long, lat)
    
            # Use the geod function to find the dist to the reference node
            _, _, ref_dist = geod.inv(node_config[ref]['Position']['lon'], node_config[ref]['Position']['lat'], long, lat)
            
            
            
            # Calculate the time difference
            timediff = (dist / c) - (ref_dist / c)
                
            
            atd[node_config[node]['Site']] = np.timedelta64(int(timediff * 1e9),'ns')
        
    return atd

#n = 24
#peak_time = '09:59.998'

pf = pd.read_excel('/Users/tjr2099/Library/CloudStorage/OneDrive-UniversityofBath/Physics w Astro Year 4/Industry Team Project/fireball_data/potential_fireballs.xlsx')
pf = pf.set_index('num')

fb_spike = pd.read_excel('/Users/tjr2099/Library/CloudStorage/OneDrive-UniversityofBath/Physics w Astro Year 4/Industry Team Project/fireball_data/fireball_spikes.xlsx')

#pot = []
d_lat = []
d_long= []
d_lat13 = []
d_long13 = []

for i,j in fb_spike.iterrows():
    
    #n=13
    #peak_time='42:51.972000'
    
    n = j['fireball']
    peak_time =  f"{fb_spike['spike_time'][i].minute:02d}:{fb_spike['spike_time'][i].second:02d}.{fb_spike['spike_time'][i].microsecond:06d}"
    
    print(n,peak_time)
    
    try:
        
        prog_start = np.datetime64("now")
        node_config = np.load("node_config_fixed.npy", allow_pickle=True).item()
        datadict = {key: {"time": [],
                          "waveform": []} for key in node_config}
        # LEELA sampling dt
        dt = np.timedelta64(int(1e9 / 109375), "ns")
        t = np.arange(0, 1024, 1) * dt
        
        # Set parameter atds to be a dictionary with Watnall as the reference node
        window_dt = np.timedelta64(10, "ms")
        
        peak_t = np.datetime64(str(pf['date'][n].date())+'T'+f"{pf['ref_time'][n].hour:02d}"+':'+str(peak_time))
        start_time = peak_t - np.timedelta64(5, 'ms')
        end_time = peak_t + np.timedelta64(5, 'ms')
        
        # Create filter here so only processed once:
        lower_f = 2000
        upper_f = 18000
        sos = signal.butter(10, [lower_f,upper_f], 'bp', fs=1/9142e-9, output='sos')
        '''
        for key in node_config.keys():
        
            # File path to data files (if have hard drive then: /Volumes/WD_ITP/Bath_VLF_2023_decoded)
            folder_path = '/Volumes/WD_ITP/bath_students_2023/bath_students_2023/'+str(n)+'/'
            # Create list of file paths for each file in the folder
            file_list = glob.glob(folder_path + "*" + key + '*.npy')
        
            # Code to create a list of all the data for each of the files
            for file_path in file_list:
                #print("Reading: " + file_path)
                chunks = np.load(file_path, allow_pickle=True).item()
                for chunk in chunks:
                    datadict[key]["time"].append(chunks[chunk]["starttime"])
                    datadict[key]["waveform"].append(chunks[chunk]["wvfmdata"])
        '''
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
                    
                    
                
        #print("data loaded")
        
        def geolocate(ref,avoid_list):
            
            '''
            Produce a set of OATDs using the find_peaks function to then be used
            with the TATDs to minimise the residual between the two lists
            '''
            
            
            
            
            # FOR THE REFERENCE:
            # Create a list of indices that are between the start and end times
           
            l1 = np.where((datadict[ref]["time"] >= start_time) &
                      (datadict[ref]["time"] <= end_time))
        
            datadict[ref]['time'] = datadict[ref]['time'][l1]
            datadict[ref]['waveform'] = datadict[ref]['waveform'][l1]
            
            sorted_indices = np.argsort(datadict[ref]['time'])
            
            # Arrange the arrays based on sorted_indices
            ref_time = datadict[ref]['time'][sorted_indices]
            datadict[ref]['waveform'] = datadict[ref]['waveform'][sorted_indices]
            
        
            ref_filtered_samples = signal.sosfiltfilt(sos, datadict[ref]['waveform'])
        
            # Assign time peak to ref_time_of_peak
            ref_time_of_peak = ref_time[100:][find_peaks(ref_filtered_samples[100:],max(ref_filtered_samples[100:])-0.1)[0][0]]
        
            time_of_peak={}
            oatds = []
        
            # FOR OTHER SITES:
            for site in datadict:
                    # Miss out sites with bad samples
                    if node_config[site]["Site"] not in avoid_list:
        
                        # Rest is the same as reference above...            
        
                        l1 = np.where((datadict[site]["time"] >= start_time) &
                                      (datadict[site]["time"] <= end_time))
                
                        datadict[site]['time'] = datadict[site]['time'][l1]
                        datadict[site]['waveform'] = datadict[site]['waveform'][l1]
                        
                        sorted_indices = np.argsort(datadict[site]['time'])
                        
                        # Arrange the arrays based on sorted_indices
                        time = datadict[site]['time'][sorted_indices]
                        datadict[site]['waveform'] = datadict[site]['waveform'][sorted_indices]
                        
        
                        filtered_samples = signal.sosfiltfilt(sos, datadict[site]['waveform'])
                        
                        #if len(filtered_samples[100:]) != 0:
                        time_of_peak[node_config[site]['Site']] = time[100:][find_peaks(filtered_samples[100:],max(filtered_samples[100:])-0.1)[0][0]]
                        
                        # Calculate OATD by taking away reference peak time
                        oatds.append((time_of_peak[node_config[site]['Site']] - ref_time_of_peak)/ np.timedelta64(1, 's'))
            
            return oatds

        avoid_list = [node_config[node]['Site'] for node in node_config if not len(datadict[node]['waveform'])>13000000]
        
        if 11-len(avoid_list)<4:
            print('Not enough nodes')
            continue
        
        def RES(guess, OATDs):
            '''
            Function to be minimised. Calculates the residual of how close the guess is 
            to the actual lightning strike location.
            '''
        
            # Use atddict function to find TATDs
            atds = atddict(guess[0],guess[1],the_chosen_one)
        
            TATDs = []
            
            # Make list of TATDs by using atddict function and avoiding troubled sites
            for site in atds:
                if site not in avoid_list:
                    TATDs.append(atds[site]/ np.timedelta64(1, 's'))
                        
            
            # Calculate standard deviation of OATDs
            std = np.std(OATDs)
            # Determine N from amount of values in OATDs
            N = len(OATDs)
            
            print(OATDs,TATDs)
            
            # Calculate the value of RES
            sumATD = 0
            
            for i in range(0, 11-len(avoid_list)):
                sumATD = sumATD + (TATDs[i] - OATDs[i]) ** 2 / std
                
            res_value = np.sqrt((1 / (N - 2) * sumATD))
            
            # Return RES value for each minimisation
            return res_value
            
        
        guess = np.array([pf['long'][n],pf['lat'][n]])
        
        if 'Camborne' not in avoid_list:
            the_chosen_one = '569218Q0B0025002A'
        elif 'Watnall' not in avoid_list:
            the_chosen_one = '569218Q0B001D0029'
        elif 'Cabauw' not in avoid_list:
            the_chosen_one = '569218Q0B003D003E'
        elif 'Lerwick' not in avoid_list:
            the_chosen_one = '569218Q0B0019002A'
        else:
            the_chosen_one = avoid_list[0]
        
        
        
        # Call the minimise function from scipy, minimising the RES function until it is sufficiently small.
        # (Can change the method of minimisation)
        result = minimize(RES, guess, args=(geolocate(the_chosen_one,avoid_list)), method='Nelder-Mead') 
        '''
        if np.abs(pf['long'][n] - result['x'][0]) < 2 and np.abs(pf['lat'][n] - result['x'][1]) < 2:
            print('**************************')
            pot.append((n,peak_time))
            '''
        print(f"Result: {result['x'][0]:.2f},{result['x'][1]:.2f} with residual value {result['fun']:.10f}")
        print(f"Actual: {pf['long'][n]},{pf['lat'][n]}\n")
        
        
        if n != 13:
            d_lat.append(result['x'][1] - pf['lat'][n])
            d_long.append(result['x'][0] - pf['long'][n])
        else: 
            d_lat13.append(result['x'][1] - pf['lat'][n])
            d_long13.append(result['x'][0] - pf['long'][n])
    except:
        print("ERROR\n")
'''       
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1, 1, 1)
ax.plot(d_long,d_lat,'b.',ms=10,label='Other')
ax.plot(d_long13,d_lat13,'r.',ms=10,label='Spikes for 09/07/2022 Event')
ax.set_xlim(-30,30)
ax.set_ylim(-30,30)
ax.tick_params(axis='x', labelsize=15)
ax.tick_params(axis='y',labelsize=15)
#ax.set_xticklabels(['',-20,-10,0,10,20,''])
#ax.set_yticklabels(['',-20,-10,0,10,20,''])
ax.set_xlabel('Difference in longitude [degrees]',fontsize=15)
ax.set_ylabel('Difference in latitude [degrees]',fontsize=15)

ax.patch.set_edgecolor('black')  

ax.patch.set_linewidth(1)
plt.grid()
plt.legend(fontsize=14)
'''
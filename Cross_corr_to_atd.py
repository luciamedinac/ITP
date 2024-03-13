import numpy as np
import matplotlib.pyplot as plt
from scipy import constants, optimize, signal
from pyproj import Geod
from scipy.optimize import minimize
import os
import glob
geod = Geod(ellps='WGS84')

def middle_element(arr):
    # Calculate the middle index
    middle_index = len(arr) // 2

    # Find the index at which the middle of the array is
    closest_index = np.abs(np.arange(len(arr)) - middle_index).argmin()
    return closest_index

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
datadict = {key: {"time": [], "waveform": []} for key in node_config}
# LEELA sampling dt
dt = np.timedelta64(int(1e9 / 109375), "ns")
# Data points in each chunk in datadict
t = np.arange(0, 1024, 1) * dt

# Set parameter atds to be a dictionary with Watnall as the reference node variables we change
meteor_lat = -1.37
meteor_lon = 51.49
dt_e = np.timedelta64(10, "s")
meteor_time = np.datetime64("2023-01-31T00:01:18")

atds = atddict(meteor_lat,meteor_lon,'569218Q0B001D0029')
start_time = meteor_time - dt_e
end_time = meteor_time + dt_e
window_dt = end_time - start_time

# Create filter here so only processed once:
lower_f = 2000
upper_f = 18000
sos = signal.butter(6, [lower_f,upper_f], 'bp', fs=1/9142e-9, output='sos')
arrival_td = {}

for key in node_config.keys():

    # File path to data files (if have hard drive then: /Volumes/WD_ITP/Bath_VLF_2023_decoded)
    folder_path = '/Users/mylesjauncey/INDUSTRY_PROJECT/Bath_VLF_2023_decoded/2023-01-31/'
    # Create list of file paths for each file in the folder
    file_list = glob.glob(folder_path + "*" + key + '*20230131*.npy')

    # Code to create a list of all the data for each of the files
    for file_path in file_list:
        print("Reading: " + file_path)
        chunks = np.load(file_path, allow_pickle=True).item()
        for chunk in chunks:
            datadict[key]["time"].append(chunks[chunk]["starttime"])
            datadict[key]["waveform"].append(chunks[chunk]["wvfmdata"])

print("data loaded")

#tao = mdt, tao is cc domain
m = np.arange(0, 10000, 1)
tao = m*dt

#initiliase array
cc = np.zeros(len(m))

study_times = np.arange(start_time, end_time, window_dt)
# Create figure
for study_time in study_times:
    study_time_end = study_time + window_dt
    fig, axs = plt.subplots(len(datadict)-6,
                            figsize=(15, 20),
                            sharex=True,
                            constrained_layout=True)
    ii = 0
    print("Processing: " + study_time.astype(str))
    # Complete operations for ref node first ('569218Q0B001D0029')
    ref_node_key = '569218Q0B001D0029'

    #finds times in datadict within study times
    datadict[ref_node_key]["time"] = np.array(datadict[ref_node_key]["time"], dtype=np.datetime64)
    l1 = np.where((datadict[ref_node_key]["time"] >= study_time) & (datadict[ref_node_key]["time"] <= study_time_end))[0]
    
    #initialise arrays
    samples_ref = np.array([])
    time_ref = np.array([], dtype='datetime64[ms]')

    #filters reference node and creates time array
    for i in l1:
        samples_ref = np.hstack((samples_ref, datadict[ref_node_key]["waveform"][i]))
        time_ref = np.hstack((time_ref, datadict[ref_node_key]["time"][i] + t))

    filtered_samples_ref = np.array(filter_t_series(time_ref, samples_ref, sos))


    for site in datadict:
        if (node_config[site]["Site"] != 'Akrotiri') and \
                (node_config[site]["Site"] != 'Gibraltar') and \
                (node_config[site]["Site"] != 'Keflavik') and \
                (node_config[site]["Site"] != 'Lerwick') and \
                (node_config[site]["Site"] != 'Herstmonceux') and \
                (node_config[site]["Site"] != 'Valentia'):

            #finds times in datadict within study times
            datadict[site]["time"] = np.array(datadict[site]["time"], dtype=np.datetime64)
            l1 = np.where((datadict[site]["time"] >= study_time) & (datadict[site]["time"] <= study_time_end))[0]

            #initialise arrays
            samples = np.array([])
            time = np.array([], dtype='datetime64[ms]')

            #filters each node and creates time array
            for i in l1:
                samples = np.hstack((samples, datadict[site]["waveform"][i]))
                time = np.hstack((time, datadict[site]["time"][i] + t))

            filtered_samples = np.array(filter_t_series(time, samples, sos))
            
            
            # If the length of the filtered samples is different, then we need to
            # pad the begining and end with 0's
            if len(filtered_samples_ref) > len(filtered_samples):
                diff = len(filtered_samples_ref) - len(filtered_samples)
                temp_filtered_samples = np.pad(filtered_samples, (int(diff/2), diff-int(diff/2)), mode='constant')
                temp_filtered_samples_ref = filtered_samples_ref
            if len(filtered_samples_ref) < len(filtered_samples):
                diff = len(filtered_samples) - len(filtered_samples_ref)
                temp_filtered_samples_ref = np.pad(filtered_samples_ref, (int(diff/2), diff-int(diff/2)), mode='constant')
                temp_filtered_samples = filtered_samples
            if len(filtered_samples_ref) == len(filtered_samples):
                temp_filtered_samples = filtered_samples
                temp_filtered_samples_ref = filtered_samples_ref

            
            n=0
            iii=0
            N = len(temp_filtered_samples)
            #calculates cc for each m point (np.roll used instead of in sum due to shape error)
            for i in m:
                cc_temp_filtered_samples_ref = np.roll(temp_filtered_samples_ref, -i)
                cc[iii] = np.sum(cc_temp_filtered_samples_ref[n:N]*temp_filtered_samples[n:N])
                if i % 100 == 0: #troubleshooting
                    print(i)
                iii+=1
            cc_normalized = cc / np.max(np.abs(cc))
            
            
            axs[ii].set_title(f'Cross-Correlation between Watnall and {node_config[site]["Site"]}')
            plt.plot(tao,cc_normalized)
            # Set x and y limits to zoom in
            #axs[ii].set_xlim(study_time , study_time_end)
            axs[ii].grid()
            #if len(filtered_samples) > 0:
            #    axs[ii].set_ylim(25 + np.max(filtered_samples),
            #                     np.min(filtered_samples) - 25)
            plt.grid(True)
            plt.show()
            raise SystemExit
            ii += 1

    # Save the figure
    # Convert study_time to a string and then to a Windows-compatible filename (i.e., no ":" or ".")
    filename_stub = study_time.astype(str).replace(":", "").replace(".", "_")

    save_as = filename_stub + ".png"
    print("Saving figure: " + save_as)
    fig.savefig(save_as)

    plt.close(fig)

prog_end = np.datetime64("now")
print("Seconds 2 run: " + (prog_end - prog_start).astype(str))


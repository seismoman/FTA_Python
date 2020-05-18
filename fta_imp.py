'''
File of helper functions for the Cross Correlation files
cwh2020
'''
import numpy as np
from math import sin, cos, sqrt, atan2, radians
from numpy.fft import fft, ifft, fftfreq
import matplotlib.pyplot as plt
import scipy as sp
from scipy import signal
from scipy.signal import argrelextrema
from scipy.interpolate import interp1d
from matplotlib import gridspec
import matplotlib.ticker as ticker
import shutil
import os

# function to save output file from fta
def save_fta_output(output_list, filename="fta.csv"):
    master_string = [str(i) for i in output_list]
    line = ','.join(master_string)
    with open(filename,'a+') as f:
        f.write(line + '\n')
    f.close()

# function to compute distance (in km) between two coordinates (in degrees)
def compute_distance(lon1,lat1,lon2,lat2):
    from math import sin, cos, sqrt, atan2, radians
    # approximate radius of earth in km
    R = 6371.0
    # convert geographic coordinates to radians
    lon1 = radians(lon1)
    lat1 = radians(lat1)
    lon2 = radians(lon2)
    lat2 = radians(lat2)
    # get distance between points
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    # solve for arc between points
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    # return distance in kilometers
    distance = R * c
    return round(distance,3)

# function to compute snr and reject xcorrs that are too noisy
# reference: Stehly et al. 2009, "Tomography of the Alpine region from observations of seismic ambient noise"
def compute_snr(xcor, distance, fs, periods, order=4):
    # corner frequencies for filtering
    low = 1/max(periods)
    high = 1/min(periods)
    # velocities of interest
    minv = 1
    maxv = 5
    # define temporal window of interest to isolate surface waves
    start = int(distance/maxv*fs)
    stop = int(distance/minv*fs)
    # filter the traces
    b, a = signal.butter(order, [low, high], btype='bandpass', fs=fs)
    filtered = signal.filtfilt(b,a,xcor)
    t0 = int((len(xcor) - 1) / 2)
    # separate into components
    causal = filtered[t0 + 1:]
    acausal = filtered[:t0][::-1]
    symmetric = causal + acausal
    # separate into surface waves and noise
    sig_c = causal[start:stop]
    sig_a = acausal[start:stop]
    sig_s = symmetric[start:stop]
    noise_c = causal[stop:]
    noise_a = acausal[stop:]
    noise_s = symmetric[stop:]
    # get amplitudes of signal
    amp_c = max(abs(sig_c))
    amp_a = max(abs(sig_a))
    amp_s = max(abs(sig_s))
    # get standard deviations of noise
    std_noise_c = np.std(noise_c)
    std_noise_a = np.std(noise_a)
    std_noise_s = np.std(noise_s)
    # compute snr for each channel
    snr_c = amp_c / std_noise_c
    snr_a = amp_a / std_noise_a
    snr_s = amp_s / std_noise_s
    # return snr information
    return ((snr_c, snr_a, snr_s))

# function to compute snr for a trace that has been convolved with Gaussian window
def compute_snr_period(trace,distance,fs):
    # velocities of interest
    minv = 1
    maxv = 5
    # define temporal window of interest to isolate surface waves
    start = int(distance/maxv*fs)
    stop = int(distance/minv*fs)
    t0 = int((len(trace)-1)/2) # starting time for each component
    # separate into causal, acausal, and symmetric components
    causal = trace[t0+1:]
    acausal = trace[:t0][::-1]
    symmetric = causal + acausal
    # separate into surface waves and noise
    sig_c = causal[start:stop]
    sig_a = acausal[start:stop]
    sig_s = symmetric[start:stop]
    noise_c = causal[stop:]
    noise_a = acausal[stop:]
    noise_s = symmetric[stop:]
    # get amplitude of signal
    amp_c = max(abs(sig_c))
    amp_a = max(abs(sig_a))
    amp_s = max(abs(sig_s))
    # get standard deviation of the noise
    std_noise_c = np.std(noise_c)
    std_noise_a = np.std(noise_a)
    std_noise_s = np.std(noise_s)
    # compute snr for each channel
    snr_c = amp_c / std_noise_c
    snr_a = amp_a / std_noise_a
    snr_s = amp_s / std_noise_s
    # return snr information
    return((snr_c, snr_a, snr_s))

# compute alpha for narrow bandpas Gaussian filter for Frequency Time Analysis
# reference: personal communication with L. Stehlym 2019
# note I kept alpha in the numerator
def compute_alpha(distance):
    if distance >= 1000: alpha = 0.3
    if distance <  1000: alpha = 0.4
    if distance <   500: alpha = 0.6
    if distance <   200: alpha = 1.0
    return alpha


# Frequency Time Analysis for Dispersion measurement and Group Velocity calculation
# reference: Hermann, 1973
# reference: Bensen et al.,
# reference: Yang et al.,
def compute_fta(xcor, time, fs, distance, alpha, periods, vmin=1, vmax=5,
    snr_threshold=4.0):
    # functions used within Frequency Time Analysis routine

    def compute_instantaneous(xcor,fs):
    	analytic = signal.hilbert(xcor)
    	env = abs(analytic)
    	instantaneous_phase = np.unwrap(np.angle(analytic))
    	instantaneous_frequency = (np.diff(instantaneous_phase) /(2.0*np.pi) * fs)
    	instantaneous_frequency = np.append(instantaneous_frequency,0);  # to have the same number of points
    	return env,  instantaneous_frequency

    def split_trace(trace,time):
        n=int(np.floor((len(time))/2)) # index of 0 value on time axis
        matrix = np.zeros(shape=(4,n))
        matrix[0,:] = trace[n:-1] # casual
        matrix[1,:] = np.flipud(trace[1:n+1]) # acausal
        matrix[2,:] = matrix[0,:]+matrix[1,:] # symmetric
        matrix[3,:] = time[n:-1] # time
        return matrix

    def compute_max(trace,start,stop):
        max_index = np.argmax(trace)
        if max_index == (start-1) or max_index == (stop-1):
            local_index = argrelextrema(trace,np.greater)[0]
            local_value = trace[local_index]
            if len(local_index) > 1:
                sorted = np.argsort(local_value)[::-1]
                max_index = sorted[1]
        return max_index

    # define temporal window of interest to isolate surface waves
    start = int(distance/vmax*fs)
    stop = int(distance/vmin*fs)
    window_index = np.arange(start,stop)
    window_length = stop-start

    # matrices to store information on the group velocity and instantaneous period for plotting/saving
    amplitude = np.zeros(shape=(3,window_length,len(periods)))
    phase = np.zeros(shape=(3,window_length,len(periods)))
    trace_label = ["P", "N", "S"]

    # matrices to store group velocities and periods for plotting
    group_velocities = np.zeros(shape=(3,len(periods)))
    group_periods = np.zeros(shape=(3,len(periods)))

    # to be filled with acceptable group velocities
    output_matrix = np.zeros(shape=(3,len(periods)))*np.nan

    # compute fft for the entire xcor just once
    xfft = fft(xcor)
    freqs = fftfreq(len(xcor), d=fs)
    # freqs = np.linspace(0,fs,len(time))

    # loop over periods
    for k, T in enumerate(periods):  # enumerate to get index and value of interst
        # construct Gaussian filter
        filter_gaussian = np.exp(-((freqs-1./T)/(alpha*1./T))**2) # Gaussian window
        # multiply fft of signal by filter in the frequency domain (colvolution)
        xfilter = xfft*filter_gaussian
        # inverse Fourier transform back into time domain & keep only real components
        xfreal = ifft(xfilter).real
        # now check SNR for given period using Gaussian windowed trace
        snr_period = compute_snr_period(xfreal,distance,fs)
        # skip to next period if snr too low for any component
        if min(snr_period)<snr_threshold:
            continue
        # get instantaneous amplitudes (envelope) and frequencies
        x_env, x_phase = compute_instantaneous(xfreal,fs)
        env_3 = split_trace(x_env,time) # matrix of causal, acausal, symmetric parts
        per_3 = split_trace(x_phase,time) # matrix of phase info for 3 traces

        for j in range(3): #loop over causal, acausal, symmetric components
            # extract tempoeral window of interest for plotting
            amplitude[j,:,k] = env_3[j][window_index]
            # extract tempoeral window of interest for disperion curve
            phase[j,:,k] = per_3[j][window_index]
            # get group velocity information
            env_3[j][0:start] = 0.
            env_3[j][stop:] = 0.
            max_index = compute_max(env_3[j],start,stop)
            travel_time = env_3[3][max_index]
            if travel_time == 0:
                gvel = np.nan # group velocity measurement
                gper = np.nan # corresponding phase
            else:
                total_time = travel_time
                gvel = distance/total_time
            # save group velocity and group period
            group_velocities[j][k] = gvel
            group_periods[j][k] = T
    # prune based on wavelength
    for i in range(3):
        for period_index, T in enumerate(periods):
            velocity = group_velocities[i,period_index]
            if velocity == 0:
                continue
            # check min number of wavelengths
            wavelength = velocity*T
            if not wavelength > 0:
                continue
            n_waves = int(distance/wavelength)
            if n_waves < 2:
                continue
            output_matrix[i,period_index] = round(velocity,4)
    # prune based on mismatch between causal and acausal measurements
    for i in range(len(periods)):
        difference = abs(output_matrix[0,i]-output_matrix[1,i])
        if difference < 0.2:
            continue
        else:
            output_matrix[:,i] = np.nan
    # prune based on maximum velocity
    for i in range(len(periods)):
        if max(output_matrix[:,i]) >= vmax:
            output_matrix[:,i] = np.nan
    # return valid group velocity measurements
    return output_matrix



# end

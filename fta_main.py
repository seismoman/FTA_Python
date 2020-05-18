'''
main file to run the software
cwh2020
'''


import fta_classes as classes # class for hdf5 files
import fta_imp as imp # functions for reading and plotting hdf5 data
import fta_plotters as fplt # do stuff
import glob # task for managing files a-la the terminal
import os # like using the terminal to move around the computer
import matplotlib.pyplot as plt
import numpy as np


# directory where all the hdf5 files reside
target_dir = '/media/cwharris/Transcend/CC/Xcorr'
# go to the target directory
os.chdir(target_dir)
#check to see if FTA output file already exists and if so delete it
if os.path.exists('fta.csv'):
    print('fta.csv already exists')
    exit()
    os.remove('fta.csv')
print('computing new group velocity file')

# set a target string for the files of interest
target_string = '*h5'
# make a list of all the files that meet the search criteria
files = glob.glob(target_string)

# NOTE: ordering of time components: 1.Causal, 2.Acausal, 3.Symmetric
    # this NEVER CHANGES
components = ['C','A','S']
# loop over all the files
for i in files:
    # initate class object
    CC = classes.CC(i)
    # get timing information
    time = CC.time
    # sampling frequency
    fs = CC.fs
    # define periods of interest
    periods = [8,15,25,45,65,85,105,125]
    # periods = np.arange(5,100)
    # get location info
    lons, lats, zeds = CC.get_xyz()
    # get station pair names to navigate tree structure
    pairs = CC.ids
    # list of all possible channel xcorrs
    channel = ['ZZ','ZR','RZ','RR']
    # construct name of xcorr record in tree structure
    for index in range(len(pairs)):
        #get station names, need to convert from binary to ascii
        sta1 = pairs[index][0].decode('ascii')
        sta2 = pairs[index][1].decode('ascii')
        # compute distance
        # if distance is too small, reject
        lon1,lon2 = lons[index][0],lons[index][1]
        lat1,lat2 = lats[index][0],lats[index][1]
        distance = imp.compute_distance(lon1,lat1,lon2,lat2)
        min_dist = 1
        if distance <= min_dist:
            continue # skip over autocorrelations

        # loop over channels
        for chan in channel:
            # read in xcorr trace
            xcor = CC.get_xcorr(sta1,sta2,chan)

            # make sure xcor isn't all zeros
            if max(xcor) == 0:
                continue

            # perform SNR check, if too low, reject
            snr = imp.compute_snr(xcor,distance,fs,periods)
            min_snr = 7.0
            if min(snr) < min_snr:
                continue

            # get alpha for FTA
            alpha = imp.compute_alpha(distance)

            # perform FTAN on the current xcorr array
                # default output: vector of group velocity and error for periods of interest
                # optional output: fta plot of amplitudes and interpolated group velocities
            outfile = sta1+"."+sta2+"."+chan #  optional plot file name
            fta_output = imp.compute_fta(xcor,time,fs,distance,alpha,periods)
            # disregard erroneous output (forwhich fta_output == False)
            if fta_output is False:
                continue
            # count number of measurements kept, skipping intances of 0 kept measurements
            n_kept = np.count_nonzero(~np.isnan(fta_output))
            if n_kept == 0:
                continue

            for j in range(3): # get each time component on separate line
                output_list = list() # create new line
                output_list.append(components[j]) # time component
                for item in (chan,distance,lon1,lat1,lon2,lat2):
                    output_list.append(item)
                for gv in fta_output[j]:
                    output_list.append(gv)
                # save vector of FTAN results:comp,chan,dist,lon1,lat1,lon2,lat2,gv_i
                imp.save_fta_output(output_list)


# end

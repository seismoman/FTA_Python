'''
File of helper functions for the Cross Correlation files
cwh2019
'''

import shutil
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection
from matplotlib import cm
from matplotlib import colors
from mpl_toolkits.basemap import Basemap



# function to save geographic info for plotting
# if file exists then append
# if file doesn't exist then create
# save x,y,z values from current hdf5 file
def save_geofile(coords, file='latlonzed.csv'):
    # default output filename is latlonzed.csv
    # you can change this when you call the function
    #    #save plotting file in the right directory
    #    plotting_dir = '/home/cwharris/Bureau/Sola/Bruit/Laurent/doc/cwh'
    #    os.chdir(plotting_dir)
    with open(file, 'a+') as f:
        # write it one row at a time
        for i in coords:
            # strip the brackets from the tuples but keep the commas
            value = ','.join(map(str, i))
            f.write(str(value) + '\n')
        f.close()


# function to save distance info for plotting
# if file exists then append
# if file doesn't exist then create
# save x1,y1,x2,y2,d values from current hdf5 file
def save_distancefile(dist, file='distances.csv'):
    with open(file, 'a+') as f:
        # write 1 row at a time
        for i in dist:
            # strip the brackets from the tuples but keep the commas
            value = ','.join(map(str, i))
            f.write(str(value) + '\n')
        f.close()


# function to prune duplicates from plotting file
# run this before plotting
# this could be incorporated into the save_geofile function
# also tracks hitcount for each station
def prune_duplicates_coords(file='latlonzed.csv'):
    # make sure file exists
    if not os.path.exists(file):
        print("couldn't find: " + file)
        print("file containing station plotting info doesn't exist")
        print("run 'save_geofile' first")
        quit()
    # read in plotting info from geofile
    data = np.genfromtxt(file, delimiter=',')
    # new array for format purposes
    new_array = [tuple(row) for row in data]
    # prune
    newer = set(new_array)
    # replace with new file
    with open(file, 'w') as f:
        # write it one row at a time
        for i in newer:
            # add in hitcount to track how many traces exist for each station
            t = tuple((i[0], i[1], i[2], new_array.count(i)))
            # strip the brackets from the tuples but keep the commas
            value = ','.join(map(str, t))
            f.write(str(value) + '\n')
        f.close()


# function to prune duplicates from plotting file
# run this before plotting
def prune_duplicates_distance(file='distances.csv'):
    # make sure file exists
    if not os.path.exists(file):
        print("couldn't find: " + file)
        print("file containing station pairs doesn't exist")
        print("run 'save_distancefile' first")
        quit()
    # read in plotting info from distancefile
    data = np.genfromtxt(file, delimiter=',')
    # new array for format purposes
    new_array = [tuple(row) for row in data]
    # prune
    newer = set(new_array)
    # replace with new file
    with open(file, 'w') as f:
        # write it one row at a time
        for i in newer:
            # strip the brackets from the tuples but keep the commas
            value = ','.join(map(str, i))
            f.write(str(value) + '\n')
        f.close()


# function to plot the seismic network
# requires geographic info as input
def plot_network(geofile='latlonzed.csv', save=0, cutoff=1):
    # default filename can be changed if needed
    # default save=0 means display but not save, change to 1 to save
    # make sure geofile exists
    if not os.path.exists(geofile):
        print("couldn't find: " + geofile)
        print("file containing station plotting info doesn't exist")
        print("run 'save_geofile' first")
        quit()
    # if geofile exists, create the basemap
    map = Basemap(llcrnrlon=-8.5, urcrnrlon=48, llcrnrlat=29.5, urcrnrlat=63,
                  projection='lcc', lat_1=15., lat_2=65., lat_0=40, lon_0=10.,
                  resolution='l')
    # draw the map
    map.drawmapboundary(fill_color='aqua')
    # color the land/water/lakes
    map.fillcontinents(color='coral', lake_color='aqua', zorder=1)
    # draw coastlines & countries
    map.drawcoastlines()
    map.drawcountries()
    # add shaded relief for topography
    map.shadedrelief(zorder=1)
    # read in plotting info from geofile using pandas
    data = pd.read_csv(geofile, header=None)
    # prune if wanted based on hit count
    p = data[(data[3] >= cutoff)]
    x, y = map(p[0].values, p[1].values)
    # plot the stations as inverted triangles
    map.scatter(x, y, marker='v', s=15, c=p[3], cmap='hot', zorder=2)
    # add colorbar
    map.colorbar()
    # either display or save the figure
    if not save:
        plt.show()
    if save:
        outfile = 'array_20_map.pdf'
        with PdfPages(outfile) as pdf:
            pdf.savefig()


# function to plot station pairs and color the line between them by distance (km)
# requires the geofile so stations are plotted only once each
def plot_pairs(distancefile='distances.csv', geofile='latlonzed.csv', save=0):
    # default filename can be changed if needed
    # default save=0 means display but not save, change to 1 to save
    # make sure distancefile exists
    if not os.path.exists(distancefile):
        print("couldn't find: " + distancefile)
        print("file containing station pairs doesn't exist")
        print("run 'save_distancefile' first")
        quit()
    # make sure geofile exists
    if not os.path.exists(geofile):
        print("couldn't find: " + geofile)
        print("file containing station plotting info doesn't exist")
        print("run 'save_geofile' first")
        quit()
    # if both files exist, create the basemap
    map = Basemap(llcrnrlon=-8.5, urcrnrlon=48, llcrnrlat=29.5, urcrnrlat=63,
                  projection='lcc', lat_1=15., lat_2=65., lat_0=40, lon_0=10.,
                  resolution='l')
    # draw the map
    map.drawmapboundary(fill_color='aqua')
    # color the land/water/lakes
    map.fillcontinents(color='coral', lake_color='aqua', zorder=1)
    # draw coastlines & countries
    map.drawcoastlines()
    #    map.drawcountries()
    # add shaded relief for topography
    map.shadedrelief(zorder=1)
    # read in plotting info from geofile using pandas
    data = pd.read_csv(geofile, header=None)
    # prune if wanted based on hit count
    x, y = map(data[0].values, data[1].values)
    # plot the stations as inverted triangles
    map.scatter(x, y, marker='v', s=15, c='r', zorder=2)
    # plot lines linking stations that share cross-correlations
    data = pd.read_csv(distancefile, header=None)
    x1, y1 = map(data[0].values, data[1].values)
    x2, y2 = map(data[2].values, data[3].values)
    map.plot((x1, y1), (x2, y2), c='k')
    # display
    plt.show()

# function to plot group velocities as lines between stations colored by gv
def plot_gv_lines(inputfile,periods):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    # read in the fta file:
    '''
    Format: lon1, lat1, lon2, lat2, gv8, gv15, gv25, gv40, gv125
    This will be changed to include error estimates
    '''
    # cd to directory with dispersion plots and fta file
    target_dir = '/home/cwharris/Bureau/Sola/Bruit/Cwh/New/Disp_Curves/'
    os.chdir(target_dir)
    # read in the file with pandas
    table = pd.read_csv(inputfile,header=None)
    # make figure
    plt.figure()
    # plot stations
    for i, T in enumerate(periods):
        if i>3:
            plt.show()
            exit()
        line_segments = list()
        color_vals = list()
        sort = table.sort_values(4+i)
        data = sort.reset_index(drop=True) # reindex so we can access values in sorted order
        gv = data[4+i] # shift index because of 4 coordinate values
        for j in range(len(gv)):
            if np.isnan(gv[j]): # don't plot empty values
                continue
            x0,y0 = data[0][j], data[1][j]
            x1,y1 = data[2][j], data[3][j]
            line_segments.append(([(x0,y0),(x1,y1)]))
            color_vals.append((gv[j]))
        # create colorscale
        cnorm = colors.Normalize(min(color_vals),max(color_vals))
        colormapper = cm.ScalarMappable(norm=cnorm)
        # store all lines as collection to plot all at once. MUCH faster
        collection = LineCollection(line_segments,colors=colormapper.to_rgba(color_vals),linewidths=0.25)
        # plot geographic map
        ax = plt.subplot(2,2,i+1,projection=ccrs.PlateCarree())
        ax.set_extent([-10, 30, 25, 65], crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE)
        # overlay lines
        ax.add_collection(collection)
        ax.autoscale()
        # add title
        ax.set_title(str(T)+'s')
        # add colorbar
        num_ticks = 5
        ticks = np.linspace(0,1,num_ticks)
        bottom = int(min(color_vals))
        top = int(max(color_vals))
        labels = np.linspace(bottom,top,num_ticks)
        axcb = plt.colorbar(collection)
        axcb.set_ticks(ticks)
        axcb.set_ticklabels(labels)



        #
        #     outfile = str(T)+"s_gv_lineplot.png"
        #     if os.path.exists(outfile):
        #         os.remove(outfile)
        #         plt.savefig(outfile)
    # #     # set up plot
    #     trace_label = ["P","N","S"]
    #     fig, axs = plt.subplots(1,3)
    #     extent = (min(periods),max(periods),vmin,vmax)
    #     velocity_axis = [distance/(i+start) for i in range(stop-start)]
    # # plot amplitudes
    # if outfile:
    #     for i in range(3):
    #         function = interp1d(periods,group_velocities[i])
    #         vinterp = function(group_periods[i])
    #         norm_amplitude = amplitude[i].copy()
    #         # norm_amplitude = norm_amplitude/np.linalg.norm(norm_amplitude)
    #         for j in range(len(periods)):
    #            norm_amplitude[:,j] = norm_amplitude[:,j]/max(norm_amplitude[:,j])
    #         axs[i].pcolormesh(periods,velocity_axis,norm_amplitude,cmap='jet',alpha=0.8)
    #         if i == 0 or i == 1:
    #             axs[i].plot(periods,vinterp,c='k')
    #         if i == 2:
    #             axs[i].plot(periods,output_matrix[i],c='w')
    #         axs[i].set_xlim(min(periods),max(periods))
    #         axs[i].set_title(trace_label[i])
    #         axs
    #         if i == 0:
    #             axs[i].set_xlabel('Period s')
    #             axs[i].set_ylabel('Velocity km/s')
    #         if i == 1:
    #             axs[i].set_xlabel(outfile)
    #         if i == 2:
    #             axs[i].set_xlabel(str(round(distance,0))+' km')
    #         plt.show()

# end

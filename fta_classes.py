'''
Class to read in Cross Correlation files in hdf5 format
    This class helps extract info from the hdf5 files:
        -geographic info: lon, lat, zed
        -cross-correlation info
    This class is compatable with plotting scripts
        -to be developed: plot network
        -to be developed: plot cross-correlations
cwh2019
'''

#import libraries
import h5py

class CC(object):
    def __init__(self,filename):
        #get filename when you call the class
        self.filename = filename
        #read in the hdf5 file
        self.f = h5py.File(self.filename,'r')
        self.fs = self.f['md_c']['tau'][()]
        self.time = self.f['md_c']['t'][()]
        self.ids = self.f['md']['id'][()]

    def get_xyz(self):
        #extract the geographic data from the file
        x = self.f['md']['lon'][()]
        y = self.f['md']['lat'][()]
        z = self.f['md']['elev'][()]
        #these arrays have many duplicates so we need to prune
        #to do so I used zip, list comprehension, and set
#        l = [(i[0][0],i[1][0],i[2][0]) for i in zip(x,y,z)]
#        coords = [i for i in set(l)]
        return(x,y,z)


    def get_xcorr(self,sta1,sta2,channel):
        xcor = self.f['ref'][str(sta1)][str(sta2)][str(channel)][()]
        return xcor





#    def get_pairs(self):
#        #extract the geographic data from the file
#        x = self.f['md']['lon'][()]
#        y = self.f['md']['lat'][()]
#        #these arrays have many duplicates so we need to prune
#        #to do so I used zip, list comprehension, and set
#        l = [(i[0],i[1]) for i in zip(x,y)]
#        #make list to be filled with distance info for plotting
#        d = []
#        #get lon1,lat1,lon2,lat2 to find distance between points
#        #file is organized: x1,x2,y1,y2
#        for i in l:
#            lon1 = i[0][0]
#            lon2 = i[0][1]
#            lat1 = i[1][0]
#            lat2 = i[1][1]
#            distance = cc_helpers.compute_distance(lon1,lat1,lon2,lat2)
#            #don't save 0km distance for station compared with itself
#            if distance>0:
#                d.append(tuple((lon1,lat1,lon2,lat2,distance)))
#        return d




#end

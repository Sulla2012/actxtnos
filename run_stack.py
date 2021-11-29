
from pixell import enmap,utils, reproject, enplot
import numpy as np
import matplotlib.pyplot as plt
import os,sys
from scipy.interpolate import interp1d
import math
import pandas as pd
import pickle as pk
import h5py
import time
from astropy import wcs
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy.io import fits
from astropy.table import QTable
import re

def tnoStamp(ra, dec, kmap, frhs, width = 0.5):
    #Takes an individual stamp of the kmap and frhs, and takes
    #the ratio of the two to generate a stamp of the S/N per pixel
    #frhs is a matched filter map and kmap is the inverse variance per pixel
    #both maps are at a ~3 day cadence

    #Inputs: ra, dec in degrees, j2000
    #kmap, and frhs, described above. They must have the same wcs
    #but if you're using this code and sigurd's maps they will
    #width, the desired width of the stamp in degrees

    #Output: a stamp, centered on ra/dec, with width width, where
    #each pixel in the map records the S/N. This stamp is for one
    #object, one 3 day map. These must then be stacked for each object
    #and then the objects can be stacked together

    #Find the pixel 
    coords = np.deg2rad(np.array((dec,ra)))
    ypix,xpix = enmap.sky2pix(kmap.shape,kmap.wcs,coords)

    #Take the ratio of the maps        
    tile = frhs/np.sqrt(kmap)
    
    #nans are formed when try to form S/n for pixels with no hits
    #I just set the S/N to 0 which I think is safe enough
    for row in tile:
        for i in range(len(row)):
            if math.isnan(row[i]):
                row[i] = 0

    #Reproject will attempt to take a stamp at the desired location: if it can't
    #for whatever reason just return None. We don't want it throwing errors
    #while on the wall, however the skips should probably be checked after
    try:
        stamp = reproject.postage_stamp(tile, ra, dec, width*60, 0.5)
    except:
        return None
    
    return stamp

class OrbitInterpolator:
    def __init__(self, table):
        self.table = table
        self.targets = np.unique(table['targetname'])

        self._construct_dictionary()

    def _interpolate_radec(self, target):
        table = self.table[self.table['targetname'] == target]
        zero = np.min(table['datetime_jd'])
        ra_interp = interp1d(table['datetime_jd'] - zero, table['RA'])
        dec_interp = interp1d(table['datetime_jd'] - zero, table['DEC'])

        return zero, ra_interp, dec_interp 

    def _construct_dictionary(self):
        self.obj_dic = {}
    
        for i in self.targets:
            z, r, d = self._interpolate_radec(i)
            self.obj_dic[i] = {}
            self.obj_dic[i]['zero'] = z
            self.obj_dic[i]['RA'] = r
            self.obj_dic[i]['DEC'] = d


    def get_radec(self, target, time):
        time = time + 2400000.5
        
        t_intep = time - self.obj_dic[target]['zero']

        ra = self.obj_dic[target]['RA'](t_intep)
        dec = self.obj_dic[target]['DEC'](t_intep)

        return ra, dec

def tnoStacker(oribits, obj):
    #Returns a stack over ~~~1~~~ objects orbit
    
    #Inputs, orbits, and OrbitInterpolator instance, which contains all the orits of interest, i.e. for multiple objects
    #obj, the name/id of the object of interest
    
    #Path to the maps. Probably shouldn't be hard coded
    path = '/home/r/rbond/sigurdkn/scratch/actpol/planet9/20200801/maps/combined/'
    
    #Initialize stack/divisor
    stack = 0
    divisor = 0
    
    #we're now going to check each directory in the path, each of which corresponds to one
    #~3 day coadd
    for dirname in os.listdir(path=path):
        with h5py.File(path + dirname +"/info.hdf", "r") as hfile:
            #Find the (rough) mjd center of the map
            mjd_cent = hfile["mjd"][()]

        #Get the ra/dec of the object 
        ra, dec = orbits.get_radec(obj, mjd_cent)

        #Get wcs info from the kmap for this 3day coadd: for sigurd's maps the kmap and 
        #frhs map wcs info will be the same
        hdu = fits.open(path + dirname + '/kmap.fits')
        w = wcs.WCS(hdu[0].header)
        
        #Find pixel corresponding to our ra/dec. I actually can no longer recall why
        #I did this and it doesn't seem to do anything so it could probably be
        #removed
        c = SkyCoord(ra, dec, unit="deg")
        x, y = w.world_to_pixel(c)           
        
        #Read in the maps and take the stamp using tnoStamp
        kmap = enmap.read_map(path + dirname + '/kmap.fits')
        frhs = enmap.read_map(path + dirname + '/frhs.fits')

        stamp = tnoStamp(ra, dec, kmap, frhs)

        #See tnoStamp but if stamp is None it means something went wrong.
        #We also check to see if the maps are uniformly zero or if they contain
        #any infinities
        if stamp is None:
            continue
        if not np.any(stamp[0]):
            continue
        if np.any(np.isinf(stamp[0])):
            continue

        #Stack it up, divide and return
        stack += stamp[0]
        #print(stack)
        divisor += 1

    stack /= divisor

    return stack, divisor
tic = time.perf_counter()

#Set path to project dir
path = '/project/r/rbond/jorlo/act_tnos/indv_stamps/'

#Command line arg that tells us which tno we're looking at by index
#This makes it reasonably convenient to 'parallelize' this in a submission script
#but isn't how you want to do this if you're not running on a cluster
tno_index = int(sys.argv[1])

#We make the stacks along an individual objects orbit, and then the resulting
#stacked maps can be combined with whatever waiting we desire

#Load orbits
orbits = pk.load(open('orbits.p', 'rb'))
hdu = fits.open('/project/r/rbond/jorlo/tno_obs.fits')

#Get a list of names, removing duplicates in a way that has consistant ordering
names = hdu[1].data['targetname']
names = np.unique(names)
names.sort()

name = names[tno_index]

#Run the stacking code
stack, divisor = tnoStacker(orbits,name)

tno_dict = {'Name':name, 'stack':stack, 'weight':divisor}

#Make the individual tno dir if it doesn't exist
dir_name = name.replace(' ', '_')
dir_name = dir_name[1:-1]

print(dir_name)

full_path = os.path.join(path, dir_name)
if not os.path.exists(full_path):
    os.makedirs(full_path)
    print('Full path made')

pk.dump(tno_dict, open(path+dir_name+'/{}.p'.format(dir_name), 'wb'))

plt.imshow(stack, vmax = 5)
plt.title(name)
plt.savefig(path+dir_name+'/{}.pdf'.format(dir_name))
plt.close()

toc = time.perf_counter()
time_hours = round((toc-tic)/3600,2)
print('Job {} ran in {:0.2f} seconds'.format(tno_index, time_hours))



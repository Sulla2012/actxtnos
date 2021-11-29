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
    
    coords = np.deg2rad(np.array((dec,ra)))
    ypix,xpix = enmap.sky2pix(kmap.shape,kmap.wcs,coords)
        
    tile = frhs/np.sqrt(kmap)
    for row in tile:
        for i in range(len(row)):
            if math.isnan(row[i]):
                row[i] = 0
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
    #Returns a stack over ~~~1~~~ orbit
    
    path = '/home/r/rbond/sigurdkn/scratch/actpol/planet9/20200801/maps/combined/'
    
    stack = 0
    divisor = 0
    
    for dirname in os.listdir(path=path):
        with h5py.File(path + dirname +"/info.hdf", "r") as hfile:
            #Find the (rough) mjd center of the map
            mjd_cent = hfile["mjd"][()]
        #print(dirname)
        ra, dec = orbits.get_radec(obj, mjd_cent)
        hdu = fits.open(path + dirname + '/kmap.fits')
        w = wcs.WCS(hdu[0].header)
        
        c = SkyCoord(ra, dec, unit="deg")
        x, y = w.world_to_pixel(c)           
        
        kmap = enmap.read_map(path + dirname + '/kmap.fits')

        frhs = enmap.read_map(path + dirname + '/frhs.fits')

        stamp = tnoStamp(ra, dec, kmap, frhs)
        if stamp is None:
            continue
        if not np.any(stamp[0]):
            continue
        if np.any(np.isinf(stamp[0])):
            continue


        stack += stamp[0]
        #print(stack)
        divisor += 1

    stack /= divisor

    return stack, divisor
tic = time.perf_counter()

#Set path to project dir
path = '/project/r/rbond/jorlo/act_tnos/ps_stamps/'

#Command line arg that tells us which tno we're looking at by index
tno_index = int(sys.argv[1])

#Load orbits
orbits = pk.load(open('ps_orbits.p', 'rb'))
#hdu = fits.open('/project/r/rbond/jorlo/tno_obs.fits')

#Get a list of names, removing duplicates in a way that has consistant ordering

name = orbits.targets[tno_index]

#Run the stacking code
stack, divisor = tnoStacker(orbits,name)

tno_dict = {'Name':name, 'stack':stack, 'weight':divisor}

#Make the individual tno dir if it doesn't exist
dir_name = str(tno_index)


full_path = os.path.join(path, str(dir_name))
if not os.path.exists(full_path):
    os.makedirs(full_path)
    


plt.imshow(stack, vmax = 5)
plt.title(name)
plt.savefig(path+dir_name+'/{}.pdf'.format(dir_name))
plt.close()
pk.dump(tno_dict, open(str(path+dir_name+'/{}.p').format(dir_name), 'wb'))





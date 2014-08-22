# 
# august 2014
# Initial file from Matt.
# Now works with astropy (Loic's changes)
#



#
# program to write out the sample in a form suitable for downloading the sdss
# jpegs 
#
#import asciidata
from astropy.io import ascii
import numpy as np
#from astLib import astCoords
import funcfindingchart
import sys, shutil, glob, os

# Loic dowloaded telarchive from http://www.mpe.mpg.de/~erwin/code/telarchive/
# It seems great for downloading fits files from sdss and other sources
from telarchive import fetchsdss

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits

import aplpy
import matplotlib.pyplot as plt



#
# read in the data
#
highsSFRsample = ascii.read('NRTsamplewidezall')
MJD = highsSFRsample.field('MJD')
PLATEID = highsSFRsample.field('plateid')
FIBERID = highsSFRsample.field('fiberid')
PHOTID = highsSFRsample.field('photoid')
RAhrs= highsSFRsample.field('rahh')
RAmin = highsSFRsample.field('ramm')
RAsec = highsSFRsample.field('rasec')
Decdeg = highsSFRsample.field('decdd')
Decmin = highsSFRsample.field('decmm')
Decsec = highsSFRsample.field('decss')
#
# read in the ra and decs in decimal degress
#
origcat = 'NRTwideselectionSDSSdata.dat'
origdata = ascii.read(origcat)
SDSS = origdata.field('Survey')
PLATEIDorig = origdata.field('plateID')
MJDorig = origdata.field('MJD')
FIBERIDorig = origdata.field('Fiber')
RAdeg =  origdata.field('RA')
Decdeg =  origdata.field('Dec')
z = origdata.field('redshift')
zwarning = origdata.field('zwarning')
StoN = origdata.field('StoN')
#RAdeg =  origdata.field('RA').tonumpy()
#Decdeg =  origdata.field('Dec').tonumpy()
#z = origdata.field('redshift').tonumpy()
#zwarning = origdata.field('zwarning').tonumpy()
#StoN = origdata.field('StoN').tonumpy()
objclass= origdata.field('Class')
#
#
#
for ind in np.arange(len(MJD)):
    subname = str(MJD[ind])+str(PLATEID[ind])+str(FIBERID[ind])
    name = "images/" + subname + ".jpg"

    #funcfindingchart.getSDSSfindingchart(str(RAdeg[ind]),str(Decdeg[ind]),name)

    coords = SkyCoord(ra=RAdeg[ind]*u.degree, dec=Decdeg[ind]*u.degree, frame='icrs')

    #print RAdeg[ind]*u.degree
    #print str(coords.ra.hms)
    #print Decdeg[ind]*u.degree
    #print str(coords.dec.dms)

    strcoords = "--coords="+coords.to_string('hmsdms')
    filters = "gr"

    """
    fetchsdss.main(["", strcoords, filters, "--output="+subname+"_sdss_"])

    files = glob.iglob("*"+subname+"*fit*")
    for file in files:
        print file
        os.remove('images/'+file)
        shutil.move(file, 'images/')

    """

    # What filter to use ?

    filter = 'r'

    file = glob.iglob("images/*"+subname+"*_"+filter+"_*.fits.gz")
    for f in file:

        radius = 0.003

        fig = aplpy.FITSFigure(f)
        fig.recenter(RAdeg[ind], Decdeg[ind], radius=radius/3)
        fig.show_grayscale()
        y_min_pix, x_min_pix = fig.world2pixel(RAdeg[ind]-radius, Decdeg[ind]-radius)
        y_max_pix, x_max_pix = fig.world2pixel(RAdeg[ind]+radius, Decdeg[ind]+radius)
        plt.show()
        plt.close()

        y_center_pix, x_center_pix = fig.world2pixel(RAdeg[ind], Decdeg[ind])
        dist_pix = np.mean([(y_max_pix - y_min_pix)/2., (x_max_pix - x_min_pix)/2.])
        dist_pix = int(dist_pix)
        y_center_pix, x_center_pix = int(y_center_pix), int(x_center_pix)
        y_min_pix, x_min_pix = y_center_pix - dist_pix, x_center_pix - dist_pix
        y_max_pix, x_max_pix = y_center_pix + dist_pix, x_center_pix + dist_pix


        # version without AplPy
        hdulist = fits.open(f)
        #print hdulist.info()
        image = hdulist[0].data


        #subimage = image[170:352, 1610:1792]
        subimage = image[x_min_pix:x_max_pix, y_min_pix:y_max_pix]
        rotsubimage = np.rot90(np.rot90(subimage))


        print x_min_pix-x_max_pix, y_min_pix-y_max_pix

        print subimage

        fig = plt.figure()
        plt.imshow(abs(subimage-rotsubimage)/subimage, interpolation='none')
        plt.show()
        plt.close()

        fig = plt.figure()
        plt.imshow(subimage, interpolation='none')
        plt.show()
        plt.close()




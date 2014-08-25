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
from astropy.table import Table


import aplpy
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

from PyPDF2 import PdfFileMerger, PdfFileReader

from scipy import ndimage

pdfmerger = PdfFileMerger()


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

# Just do the first galaxies or all of them :
#n_galaxies_limit = len(MJD)
n_galaxies_limit = 30


# Makes a table to store results in.
galaxyIDs = np.zeros(n_galaxies_limit)
asymetries_corrected = np.zeros(n_galaxies_limit)
asymetries_object = np.zeros(n_galaxies_limit)
asymetries_void = np.zeros(n_galaxies_limit)

table_results = Table([galaxyIDs, asymetries_object, asymetries_void, asymetries_corrected], names=("ID", "asymetry object", "asymetry void", "asymetry corrected"), meta={'name': 'Results'})


# Make a table of the voids (just once, then I'll read the table directly)
voidfile = "void_positions.dat"
if not os.path.exists(voidfile):
    galaxyIDs = np.array([])
    void_xmin = np.array([])
    void_xmax = np.array([])
    void_ymin = np.array([])
    void_ymax = np.array([])
    table_voids = Table([galaxyIDs, void_xmin, void_xmax, void_ymin, void_ymax], names=("ID", "Xmin", "Xmax", "Ymin", "Ymax"), meta={'name': 'Table of void positions'})
    table_voids.write(voidfile, format='ascii.commented_header')
else:
    table_voids = Table.read(voidfile, format='ascii.commented_header')
    #print table_voids


# Loops on galaxies
for ind in np.arange(n_galaxies_limit):
    subname = str(MJD[ind])+str(PLATEID[ind])+str(FIBERID[ind])
    name = "images/" + subname + ".jpg"

    #funcfindingchart.getSDSSfindingchart(str(RAdeg[ind]),str(Decdeg[ind]),name)

    coords = SkyCoord(ra=RAdeg[ind]*u.degree, dec=Decdeg[ind]*u.degree, frame='icrs')

    #print RAdeg[ind]*u.degree
    #print str(coords.ra.hms)
    #print Decdeg[ind]*u.degree
    #print str(coords.dec.dms)

    strcoords = "--coords="+coords.to_string('hmsdms')
    filterstofetch = "r"

    """
    fetchsdss.main(["", strcoords, filterstofetch, "--output="+subname+"_sdss_"])

    files = glob.iglob("*"+subname+"*fit*")
    for file in files:
        print file
        if os.path.exists('images/'+file): os.remove('images/'+file)
        shutil.move(file, 'images/')

    """

    # What filter to use ?

    filter = 'r'

    file = glob.iglob("images/"+subname+"*_"+filter+"_*.fits.gz")
    for f in file:

        print "ici"
        print f

        radius = 0.002

        fig = aplpy.FITSFigure(f)
        fig.recenter(RAdeg[ind], Decdeg[ind], radius=radius/3)
        fig.show_grayscale()
        y_min_pix, x_min_pix = fig.world2pixel(RAdeg[ind]-radius, Decdeg[ind]-radius)
        y_max_pix, x_max_pix = fig.world2pixel(RAdeg[ind]+radius, Decdeg[ind]+radius)
        #plt.show()
        plt.close()


        # version without AplPy
        hdulist = fits.open(f)
        #print hdulist.info()
        image = hdulist[0].data

        y_center_pix, x_center_pix = fig.world2pixel(RAdeg[ind], Decdeg[ind])
        dist_pix = np.mean([(y_max_pix - y_min_pix)/2., (x_max_pix - x_min_pix)/2.])
        dist_pix = int(dist_pix)
        y_center_pix, x_center_pix = int(y_center_pix), int(x_center_pix)
        y_min_pix, x_min_pix = y_center_pix - dist_pix, x_center_pix - dist_pix
        y_max_pix, x_max_pix = y_center_pix + dist_pix, x_center_pix + dist_pix




        ########################
        # Petrosian radius     #
        ########################


        subimage = image[x_min_pix:x_max_pix, y_min_pix:y_max_pix]

        fig = plt.figure()
        plt.imshow(subimage, interpolation='none')
        plt.show()
        plt.close()

        sub_center_x = x_center_pix - x_min_pix
        sub_center_y = y_center_pix - y_min_pix


        y,x = np.ogrid[-sub_center_x:2*dist_pix-sub_center_x, -sub_center_y:2*dist_pix-sub_center_y]

        radial_luminosity = np.array([])
        radius = np.arange(dist_pix)
        for r in radius:
            circles = x*x + y*y <= r*r
            """
            fig = plt.figure()
            plt.imshow(circles, interpolation='none')
            plt.show()
            plt.close()
            print circles
            """
            radial_luminosity = np.append(radial_luminosity, np.sum(circles*subimage)/(np.pi*r*r))

        fig = plt.figure()
        plt.plot(radius, radial_luminosity)
        plt.show()
        plt.close()




        ########################
        # Asymetry measurement #
        ########################

        #TODO: delete sky from subimage for asymetry.
        #TODO: put blank sky closer to galaxy
        #TODO: what radius to take for asymetry ? The asymmetry (x 4.2) and clumpiness (x 4.3) parameters are also measured within the 1:5?? 0 radius (see CBJ00 for a discussion of the benefits of using this particular value). IN THE RELATIONSHIP BETWEEN STELLAR LIGHT DISTRIBUTIONS OF GALAXIES AND THEIR FORMATION HISTORIES Conselice 2003






        # Will test different centers and take the one
        # where the asymetry is minimal
        n_offset_x = 20
        n_offset_y = n_offset_x

        # Puts the first center at the lower left corner of the possible centers
        x_min_pix -= n_offset_x/2
        x_max_pix -= n_offset_x/2
        y_min_pix -= n_offset_y/2
        y_max_pix -= n_offset_y/2

        asymetry_array = np.zeros((n_offset_x, n_offset_y))


        # runs over all possible centers
        for offset_x in np.arange(n_offset_x):
            for offset_y in np.arange(n_offset_y):

                subimage = image[x_min_pix+offset_x:x_max_pix+offset_x, y_min_pix+offset_y:y_max_pix+offset_y]
                rotsubimage = np.rot90(np.rot90(subimage))

                diffrot = subimage-rotsubimage

                asymetry_array[offset_x, offset_y] = np.sum(np.abs(diffrot)) / np.sum(np.abs(subimage))

        asymetry = np.min(asymetry_array)
        center = np.where(asymetry_array == asymetry)
        print center

        subimage = image[x_min_pix+center[0]:x_max_pix+center[0], y_min_pix+center[1]:y_max_pix+center[1]]
        rotsubimage = np.rot90(np.rot90(subimage))
        diffrot = subimage-rotsubimage

        pdfpagename = 'results/results_'+subname+'.pdf'
        pdf = PdfPages(pdfpagename)

        fig = plt.figure()
        fig.text(.1, .7, subname)

        plt.subplot(231)
        plt.title("selected region")
        plt.imshow(subimage, interpolation='none')
        #plt.show()
        #pdf.savefig()

        plt.subplot(232)
        plt.title("difference")
        plt.imshow(diffrot, interpolation='none')
        #plt.show()
        #pdf.savefig()

        plt.subplot(233)
        plt.title("center search")
        plt.imshow(asymetry_array, interpolation='none')
        #plt.show()



        # Find a void region if it has not been already computed:
        voidfile_row = (np.where(table_voids['ID'] == int(subname)))[0]
        print voidfile_row
        if not len(voidfile_row):
            # in this case it computes the void and saves its parameters.

            n_regions_x = image.shape[0] - dist_pix * 2
            n_regions_y = image.shape[1] - dist_pix * 2

            regions_moment = np.zeros((n_regions_x, n_regions_y))

            for region_x in np.arange(n_regions_x):
                for region_y in np.arange(n_regions_y):
                    region = image[region_x:region_x+dist_pix*2, region_y:region_y+dist_pix*2]
                    regions_moment[region_x, region_y] = ndimage.standard_deviation(region)

            voidpos = np.where(regions_moment == np.min(regions_moment))

            void = image[voidpos[0]:voidpos[0]+dist_pix*2, voidpos[1]:voidpos[1]+dist_pix*2]

            table_voids.add_row([subname, voidpos[0], voidpos[0]+dist_pix*2, voidpos[1], voidpos[1]+dist_pix*2])
            table_voids.write(voidfile, format='ascii.commented_header')

        else:
            # in this case it just gets the void from the file.
            void = image[table_voids["Xmin"][voidfile_row]:table_voids["Xmax"][voidfile_row], table_voids["Ymin"][voidfile_row]:table_voids["Ymax"][voidfile_row]]



        rotvoid = np.rot90(np.rot90(void))
        diffrotvoid = void-rotvoid
        asymetry_void = np.sum(np.abs(diffrotvoid)) / np.sum(np.abs(subimage))

        asymetry_corrected = asymetry - asymetry_void
        print "asymetry_corrected = ", asymetry_corrected
        print "asymetry = ", asymetry
        print "asymetry_void = ", asymetry_void

        table_results["ID"][ind] = subname
        table_results["asymetry object"][ind] = asymetry
        table_results["asymetry void"][ind] = asymetry_void
        table_results["asymetry corrected"][ind] = asymetry_corrected




        fig.text(.1, .04, "asymetry_corrected = "+str(asymetry_corrected))
        fig.text(.1, .02, "asymetry = "+str(asymetry))
        fig.text(.1, .0, "asymetry_void = "+str(asymetry_void))

        plt.subplot(234)
        plt.title("void")
        plt.imshow(void, interpolation='none')

        pdf.savefig()
        plt.close()

        pdf.close()




        pdfmerger.append(PdfFileReader(pdfpagename, 'rb'))


table_results.write('results.html', format='html')


pdfmerger.write("results.pdf")


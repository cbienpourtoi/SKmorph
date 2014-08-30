#
# august 2014
# This file will initialize pymorph:
# - creer des dossiers aux noms des objets.
# - mettre des images decompressees dedans
# - CROPPER CES IMAGES ? on gagnerait du temps, et la psf serait forcement assez proche...
# - creer des listes contenant les identifiants, coordonneess et redshift de l obj
# - creer un fichier de config
# - passer le parametre 1 en 0 ?????????????????


from astropy.io import ascii
import numpy as np
import funcfindingchart
import sys, shutil, glob, os


# telarchive from http://www.mpe.mpg.de/~erwin/code/telarchive/
# It seems great for downloading fits files from sdss and other sources
from telarchive import fetchsdss

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, hstack, Column

import cutout

import aplpy
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

from PyPDF2 import PdfFileMerger, PdfFileReader

from scipy import ndimage


# Set to True if you wanna download sdss jpeg images
getjpg = True

# Set to True if you wanna download sdss fits images (slow but necessary the first time)
getfits = True

main_directory = "/home/loic/Projects/pymorph/"

launcher_name = main_directory + "launch_ID.sh"
launcher_file = open(launcher_name, "w")
launcher_file.write("#!/bin/bash\n")
launcher_file.close()


# These two files contain different (or sometimes redondant) informations about the same objects, in the same order:
NRTsamplewidezall = Table.read('NRTsamplewidezall', format='ascii.commented_header')
NRTwideselectionSDSSdata = Table.read('NRTwideselectionSDSSdata.dat', format='ascii.commented_header')

# Here I check that the files really contain the same stuff:

print "The 4 next lines must be [], else something is wrong with the files."
print np.where(NRTsamplewidezall['plateid'] != NRTwideselectionSDSSdata['plateID'])
print np.where(NRTsamplewidezall['MJD'] != NRTwideselectionSDSSdata['MJD'])
print np.where(NRTsamplewidezall['fiberid'] != NRTwideselectionSDSSdata['Fiber'])
print np.where(np.abs(NRTsamplewidezall['redshift'] - NRTwideselectionSDSSdata['redshift']) > 0.00006)

galaxy_table = hstack([NRTsamplewidezall, NRTwideselectionSDSSdata])

names = []
i = 0
IDs = []
for galaxy in galaxy_table:
    IDs.append("ID_"+str("%04d"%i))
    sdssname = str(galaxy['MJD_1']) + '_' + str(galaxy['plateid']) + '_' + str(galaxy['fiberid'])
    coordname = str(galaxy['rahh']) + '_' + str(galaxy['ramm']) + '_' + str(galaxy['rasec']) + '_' + str(galaxy['decdd']) + '_' + str(galaxy['decmm']) + '_' + str(galaxy['decss'])
    names.append(sdssname + '__' + coordname)
    i+=1
galaxy_table.add_column(Column(name='name', data=names))
galaxy_table.add_column(Column(name='ID', data=IDs), index=0)

print galaxy_table

"""
fig = plt.figure()
plt.plot(galaxy_table['RA'], galaxy_table['Dec'], 'x')
plt.show()
plt.close()
"""


print galaxy_table.columns

galaxy_table.write("galaxy_table.html")
galaxy_table.write("galaxy_table.txt", format='ascii.commented_header')


error_filename = "errors_galaxies.txt"
error_file = open(error_filename, "w")
error_file.close()

# Loops on galaxies
for galaxy in galaxy_table:

    error = False

    #if galaxy["ID"] == "ID_0120":
    if True:

        directory = main_directory + galaxy['ID'] + "/"
        if not os.path.exists(directory):
            os.makedirs(directory)

        # gets jpeg images from sdss
        if getjpg:
            jpgname = directory + galaxy['ID'] + ".jpg"
            if not os.path.exists(jpgname):
                funcfindingchart.getSDSSfindingchart(str(galaxy['RA']), str(galaxy['Dec']), jpgname)

        coords = SkyCoord(ra=galaxy['RA']*u.degree, dec=galaxy['Dec']*u.degree, frame='icrs')


        strcoords = coords.to_string('hmsdms')
        filter = "g"

        stamp_name = directory+galaxy['ID']+'_'+filter+'_small.fits'
        large_image_file = directory+galaxy['ID']+'_'+filter+'.fits'

        if getfits:
            if not os.path.exists(directory+galaxy['ID']+'_'+filter+'.fits') or not os.path.exists(stamp_name):
                if not os.path.exists(directory+galaxy['ID']+'_'+filter+'.fits'):
                    if len(glob.glob("*"+galaxy['ID']+"*"+filter+"*"+"*fit*gz*")) < 1:
                        fetchsdss.main(["", "--coords="+strcoords, filter, "--output="+galaxy['ID']+'_'])

                    files = glob.glob("*"+galaxy['ID']+"*"+filter+"*"+"*fit*gz*")

                    for file in files:
                        print file
                        hdulist = fits.open(file)
                        print hdulist[0].header['filter']
                        hdulist.writeto(large_image_file, output_verify='ignore', clobber = True)

                if not os.path.exists(stamp_name):
                    try:
                        cutout.cutout(large_image_file, galaxy['RA'], galaxy['Dec'], 0.035, 0.035, units='wcs', outfile=stamp_name, coordsys='celestial')
                    except:
                        error_file = open(error_filename, "a")
                        error_file.write(stamp_name+"\n")
                        error_file.close()
                        error = True

        if error == False:

            # Will construct the catalog which had the position and resdhift of the galaxy to analyze (called clus_cata in pymorph)
            catalog_file = directory+"clus_cata.cat"
            f = open(catalog_file, 'w')
            f.write("gal_id ra1 ra2 ra3 dec1 dec2 dec3 z\n")
            f.write(str(galaxy['ID'])+" "+coords.to_string('hmsdms', sep=" ")+" "+str(galaxy['redshift_1']))
            f.close()

            config_file_example = "config.py"
            f_example = open(config_file_example, 'r')
            f_config = open(directory+config_file_example, 'w')
            for line in f_example:
                line = line.replace('replaceheretherightfile.fit', os.path.basename(stamp_name))
                f_config.write(line)
            f_example.close()
            f_config.close()

            launcher_file = open(launcher_name, "a")
            launcher_file.write("cd "+directory+"\n")
            launcher_file.write("~/Projects/pymorph/pymorph.py\n")
            launcher_file.write("mkdir "+galaxy['ID']+"_phase1\n")
            launcher_file.write("cp -f /home/loic/Projects/pymorph/"+galaxy['ID']+"/*  /home/loic/Projects/pymorph/"+galaxy['ID']+"/"+galaxy['ID']+"_phase1/\n\n")


            launcher_file.close()





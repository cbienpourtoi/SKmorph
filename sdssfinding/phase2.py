#
# august 2014
# This file will change pymorpg config.py file from indicator 1 (psf check) to 0 (run)
# only if psf has been detected


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


main_directory = "/home/loic/Projects/pymorph/"

launcher_name = main_directory + "launch_ID_phase2.sh"
launcher_file = open(launcher_name, "w")
launcher_file.write("#!/bin/bash\n")
launcher_file.close()


galaxy_table = Table.read("galaxy_table.txt", format='ascii.commented_header')

print galaxy_table

for galaxy in galaxy_table:

    directory = main_directory + galaxy['ID'] + "/"
    if os.path.exists(directory+"psflist.list"):
        launcher_file = open(launcher_name, "a")
        launcher_file.write("cd "+directory+"\n")
        launcher_file.write("~/Projects/pymorph/pymorph.py\n")
        launcher_file.close()

        config_file_example = "config.py"
        config_file_phase2 = "config2.py"
        f_example = open(directory+config_file_example, 'r')
        f_config = open(directory+config_file_phase2, 'w')
        for line in f_example:
            line = line.replace('psfselect = 1', 'psfselect = 0')
            print line
            f_config.write(line)
        f_example.close()
        f_config.close()
        os.rename(directory+config_file_example, directory+config_file_example+".sauvegarde")
        os.rename(directory+config_file_phase2, directory+config_file_example)




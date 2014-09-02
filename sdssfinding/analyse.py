#
# sept 2014
# Analysis of the results output by pymorph


import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
import sys, shutil, glob, os

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, hstack, Column, vstack, join

#import funcfindingchart

# telarchive from http://www.mpe.mpg.de/~erwin/code/telarchive/
# It seems great for downloading fits files from sdss and other sources
#from telarchive import fetchsdss

#import cutout
#import aplpy

#from PyPDF2 import PdfFileMerger, PdfFileReader
#from scipy import ndimage


# Where the data output from pymorph is:
main_directory = "/home/loic/Projects/pymorph/"

# created table name:
name_largetable = "all_results.txt"


use_existing_table = True

if use_existing_table:
    large_table = Table.read(name_largetable, format='ascii.commented_header')
else:

    # My inital data about galaxies from M and W:
    galaxy_table = Table.read("galaxy_table.txt", format='ascii.commented_header')

    print "Total number of objects: " + str(len(galaxy_table))

    # make a selection on this data:
    limit_Mstar = 9.5
    selected_positions = np.where(galaxy_table['logMstar'] > limit_Mstar)[0]
    print "Number of selected objects: " + str(len(selected_positions))
    selection = np.zeros(len(galaxy_table), np.bool)
    selection[selected_positions] = True


    """
    # In case I want to select by sfr:
    limit_sSFR = -8.5
    selection = np.where(galaxy_table['logsSFR'] > limit_sSFR)[0]
    print "Number of selected objects: " + str(len(selection))
    not_selection = np.ones(len(galaxy_table), np.bool)
    not_selection[selection] = False
    """

    fig = plt.figure()
    plt.plot(galaxy_table['logMstar'][selection.__invert__()], galaxy_table['logsSFR'][selection.__invert__()], 'xb')
    plt.plot(galaxy_table['logMstar'][selection], galaxy_table['logsSFR'][selection], 'xr')
    plt.xlabel("log(Mstar)")
    plt.ylabel("log(sSFR)")
    plt.title("Selected objects")
    plt.savefig("Mstar_sSFR.png")
    #plt.show()
    plt.close()

    fig = plt.figure()
    plt.plot(galaxy_table['logMstar'][selection.__invert__()], galaxy_table['logSFR'][selection.__invert__()], 'xb')
    plt.plot(galaxy_table['logMstar'][selection], galaxy_table['logSFR'][selection], 'xr')
    plt.xlabel("log(Mstar)")
    plt.ylabel("log(SFR)")
    plt.title("Selected objects")
    plt.savefig("Mstar_SFR.png")
    #plt.show()
    plt.close()

    subsample = galaxy_table[selection]

    pymorph_result = Table.read(main_directory+subsample['ID'][0]+'/'+'result.csv', format='ascii.csv')
    pymorph_result.remove_row(0)

    n_errors = 0
    for galaxy in subsample:
        try:
            pymorph_result = vstack([pymorph_result,Table.read(main_directory+galaxy['ID']+'/'+'result.csv', format='ascii.csv')])
        except:
            print "No result.csv for "+str(galaxy['ID'])
            n_errors+=1
    print "Number of errors : " + str(n_errors)

    pymorph_result.rename_column('Name_1','ID')
    print pymorph_result

    for galaxy in pymorph_result:
        galaxy['ID'] = (galaxy['ID'])[1:]

    large_table = join(subsample, pymorph_result, keys='ID')
    large_table.write(name_largetable, format='ascii.commented_header')

    """
    print len(large_table)
    print len(subsample)
    print len(pymorph_result)

    print large_table
    """

#selection = np.where(np.abs(large_table['C_10']) < 50.)
#large_table = large_table[selection]

#selection = np.where(np.abs(large_table['A_12']) < 50.)
#large_table = large_table[selection]


colormap = np.array(['b', 'g', 'r'])

s_color = np.int8(large_table['S_14'].data*0+1)
s_color[np.where(large_table['S_14']<0.1)] = 2
s_color[np.where(large_table['S_14']>0.35)] = 0

c_color = np.int8(large_table['C_10'].data*0+1)
c_color[np.where(large_table['C_10']>4)] = 2
c_color[np.where(large_table['C_10']<3)] = 0

a_color = np.int8(large_table['A_12'].data*0+1)
a_color[np.where(large_table['A_12']<0.1)] = 2
a_color[np.where(large_table['A_12']>.35)] = 0


fig = plt.figure(figsize=[12,12])
plt.subplot(221)
plt.scatter(large_table['C_10'], large_table['A_12'], c=colormap[s_color], s=10, lw=0)
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.xlabel('C (Concentration)')
plt.ylabel('A (Asymetry)')
plt.xlim([5.1,1.8])
plt.ylim([1.,-0.1])
#plt.show()
#plt.close()

#fig = plt.figure()
plt.subplot(222)
plt.scatter(large_table['S_14'], large_table['A_12'], c=colormap[c_color], s=10, lw=0)
plt.gca().invert_yaxis()
plt.xlabel('S (Clumpiness)')
plt.ylabel('A (Asymetry)')
plt.xlim([-.1,1.])
plt.ylim([1.,-0.1])
#plt.show()
#plt.close()


#fig = plt.figure()
plt.subplot(223)
plt.scatter(large_table['C_10'], large_table['S_14'], c=colormap[a_color], s=10, lw=0)
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.xlabel('C (Concentration)')
plt.ylabel('S (Clumpiness)')
plt.xlim([5.1,1.8])
plt.ylim([1.,-0.1])

plt.savefig("Conselice.png")
plt.show()
plt.close()


#
# sept 2014
# Analysis of the results output by pymorph


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
#from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys, shutil, glob, os
from matplotlib import gridspec

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, hstack, Column, vstack, join

import Image

from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4, cm,landscape
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, TableStyle, Spacer
from reportlab.platypus import Table as reportTable
from reportlab.platypus import Image as reportImage
from reportlab.lib.enums import TA_LEFT, TA_CENTER
from reportlab.lib import colors
from reportlab.lib.units import inch

#import funcfindingchart

# telarchive from http://www.mpe.mpg.de/~erwin/code/telarchive/
# It seems great for downloading fits files from sdss and other sources
#from telarchive import fetchsdss

#import cutout
#import aplpy

#from PyPDF2 import PdfFileMerger, PdfFileReader
#from scipy import ndimage

def to_matrix(l, n):
    return [l[i:i+n] for i in xrange(0, len(l), n)]

# For matplotlib histograms in percents
"""
def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if plt.rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'
"""



# Where the data output from pymorph is:
main_directory = "/home/loic/Projects/pymorph/"

# Where the figures will be put:
figdir = "figures/"

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

print large_table.colnames


large_table = large_table[np.where(large_table['C_10'] < 10)]
large_table = large_table[np.where(large_table['A_12'] < 10)]




#################################
#### List of objects plots ######
####  selected by parameters ####
#################################


table_high_C = large_table[np.where(large_table['C_10'] > 3.2)]
table_high_C.meta['name'] = "High Concentration"
print table_high_C.meta['name']

table_low_C = large_table[np.where(large_table['C_10'] < 2.3)]
table_low_C.meta['name'] = "Low Concentration"
print table_low_C.meta['name']

print large_table['S_14']
large_table.sort('S_14')
print large_table['S_14']
print len(large_table)
large_table.meta['name'] = "All"








sys.exit()

list_of_tables = [large_table]

print len(table_low_C)
for table_sampled in list_of_tables:
    images = []
    for object in table_sampled:
        image = reportImage(main_directory+object['ID']+'/'+object['ID']+'.jpg')
        image.drawHeight = .4*inch
        image.drawWidth = .4*inch
        images.append(image)
        #print images
        #image.show()
    data=to_matrix(images, 25)
    doc = canvas.Canvas(figdir+"Images"+table_sampled.meta['name']+".pdf", pagesize=landscape(A4))
    table = reportTable(data, colWidths=.41*inch, rowHeights=.41*inch)
    table.setStyle(TableStyle([
                               ('INNERGRID', (0,0), (-1,-1), 0., colors.black),
                               ('BOX', (0,0), (-1,-1), 0., colors.black),
                               ('BACKGROUND',(0,0),(-1,2),colors.white)
                               ]))
    table.wrapOn(doc, 200, 400)
    table.drawOn(doc, 2, 5)
    #doc.append(table)
    doc.save()




list_of_tables = [table_high_C, table_low_C]

print len(table_low_C)
for table_sampled in list_of_tables:
    images = []
    for object in table_sampled:
        image = reportImage(main_directory+object['ID']+'/'+object['ID']+'.jpg')
        image.drawHeight = 1*inch
        image.drawWidth = 1*inch
        images.append(image)
        #print images
        #image.show()
    data=to_matrix(images, 10)
    doc = canvas.Canvas(figdir+"Images"+table_sampled.meta['name']+".pdf", pagesize=landscape(A4))
    table = reportTable(data, colWidths=1.1*inch, rowHeights=1.1*inch)
    table.setStyle(TableStyle([
                               ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),
                               ('BOX', (0,0), (-1,-1), 0.25, colors.black),
                               ('BACKGROUND',(0,0),(-1,2),colors.white)
                               ]))
    table.wrapOn(doc, 200, 400)
    table.drawOn(doc, 2, 5)
    #doc.append(table)
    doc.save()


#############################
####  Lee13 like plots   ####
#############################

# C vs A, with C and A on the sides:

fig = plt.figure(figsize=(6, 7))
gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[3, 2])
fig.subplots_adjust(hspace=0,wspace=0)

ax0 = plt.subplot(gs[0])
ax0.scatter(large_table['C_10'], large_table['A_12'], s=10, lw=0)
plt.ylabel('A (Asymetry)')
plt.xlim([1.9, 3.7])
plt.ylim([-.2, 0.7])
plt.setp(ax0.get_xticklabels(), visible=False)


ax1 = plt.subplot(gs[1])
ax1.hist(large_table['A_12'], bins=9, range=[-.2, 0.7], orientation='horizontal')
#plt.xlabel('A')
#plt.ylabel('#')
plt.ylim([-.2, 0.7])
plt.setp(ax1.get_yticklabels(), visible=False)

ax2 = plt.subplot(gs[2])
ax2.hist(large_table['C_10'], bins=18, range=[1.9, 3.7])
plt.xlabel('C (Concentration)')
#plt.ylabel('#')
plt.xlim([1.9, 3.7])

#plt.show()

plt.savefig(figdir+"Lee_AandC_Fig9.png")

plt.close()


#large_table = large_table[np.where(large_table['G_16'] < 10)]

# hist: G
fig = plt.figure(figsize=[6,4])
plt.hist(large_table['G_16'], bins=6, range=[0.1, 0.7])
plt.xlabel('G (Gini)')
plt.ylabel('#')
plt.xlim([0.1, 0.7])
#plt.show()
plt.savefig(figdir+"Lee_G_Fig5.png")
plt.close()


# hist: M20
fig = plt.figure(figsize=[6,4])
plt.hist(large_table['M_17'], bins=8, range=[-2.2, -.6])
plt.xlabel('M20')
plt.ylabel('#')
plt.xlim([-2.2, -.6])
#plt.show()
plt.savefig(figdir+"Lee_M_Fig5.png")
plt.close()



# hist: A
fig = plt.figure(figsize=[6,4])
plt.hist(large_table['A_12'], bins=9, range=[-.2, 0.7])
plt.xlabel('A')
plt.ylabel('#')
#plt.show()
plt.savefig(figdir+"Lee_A.png")
plt.close()


# C vs A
fig = plt.figure(figsize=[6,8])
plt.scatter(large_table['C_10'], large_table['A_12'], s=10, lw=0)
plt.xlabel('C (Concentration)')
plt.ylabel('A (Asymetry)')
plt.xlim([1.9, 3.7])
plt.ylim([-.2, 0.7])
plt.savefig(figdir+"Lee_CvsA.png")
#plt.show()
plt.close()


# G vs M20
fig = plt.figure(figsize=[6,6])
plt.scatter(large_table['M_17'], large_table['G_16'], s=10, lw=0)
plt.ylabel('G (Gini)')
plt.xlabel('M20')
plt.xlim([-2.3, -.6])
plt.ylim([.2, .7])
#plt.show()
plt.savefig(figdir+"Lee_GvsM20_Fig6.png")
plt.close()



# G vs Re
fig = plt.figure(figsize=[6,5])
plt.scatter(np.log10(large_table['re_kpc_29']), large_table['G_16'], s=10, lw=0)
plt.ylabel('G (Gini)')
plt.xlabel('Re [kpc]')
plt.xlim([-1.6, 1.])
plt.ylim([.2, .8])
plt.savefig(figdir+"Lee_RevsG.png")
#plt.show()
plt.close()


# M20 vs Re
fig = plt.figure(figsize=[6,8])
plt.scatter(np.log10(large_table['re_kpc_29']), large_table['M_17'], s=10, lw=0)
plt.ylabel('M20')
plt.xlabel('Re [kpc]')
plt.xlim([-.6, 1.])
plt.ylim([-2.1, -.9])
plt.savefig(figdir+"Lee_RevsM.png")
#plt.show()
plt.close()



#############################
#### Conselice like plot ####
#############################


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

plt.savefig(figdir+"Conselice.png")
#plt.show()
plt.close()


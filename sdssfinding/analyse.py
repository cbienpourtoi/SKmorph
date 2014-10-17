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

from scipy import interpolate
from scipy import interp

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


sort_index = 'C_10'

print large_table[sort_index]
large_table.sort(sort_index)
print large_table['ID'], large_table[sort_index]
print len(large_table)
large_table.meta['name'] = "All"

filter = 'g'


#Test concentration for 1 galaxy :
test_name = 'ID_1105'
test_name = 'ID_1532'
test_name = 'ID_2211'
test_name = 'ID_0641'
test_name = 'ID_0348'
test_name = 'ID_1891'
test_name = 'ID_1481'
test_name = 'ID_2196'
test_name = 'ID_1072'
test_name = 'ID_2131'
test_name = 'ID_2100'
test_name = 'ID_1414'
test_name = 'ID_0002'
test_name = 'ID_1299'
test_name = 'ID_2399'
test_name = 'ID_2689'
test_name = 'ID_2510'
test_name = 'ID_0818'
test_name = 'ID_2130'
test_name = 'ID_2572'
test_name = 'ID_2233'
test_name = 'ID_2576'
test_name = 'ID_0222'
test_name = 'ID_0452'
test_name = 'ID_0304'
test_name = 'ID_2427'
test_name = 'ID_2048'
test_name = 'ID_0112'
test_name = 'ID_2692'
test_name = 'ID_1258'
test_name = 'ID_0329'
test_name = 'ID_1732'
test_name = 'ID_0724'
test_name = 'ID_1805'
test_name = 'ID_1974'
test_name = 'ID_1805'
test_name = 'ID_1633'
test_name = 'ID_1852'
test_name = 'ID_2520'
test_name = 'ID_0381'
test_name = 'ID_1832'
test_name = 'ID_1533'
test_name = 'ID_0225'
test_name = 'ID_2312'
test_name = 'ID_1243'
test_name = 'ID_0428'
test_name = 'ID_0335'

hdulist = fits.open(main_directory+test_name+'/'+test_name+'_'+filter+'_mini.fits')
mini = hdulist[0].data
hdulist.close()

print large_table[np.where(large_table['ID'] == test_name)]['RA']
print large_table[np.where(large_table['ID'] == test_name)]['Dec']


table_test_objects = Table([[''], [2.], ['']], names=('ID', 'skyvalue', 'comment'), meta={'name': 'table of test objects'})
table_test_objects.add_row(['ID_1105', 1076.5, 'x'])

print table_test_objects

sys.exit()
skyvalue = 1076.5 #ID_1105
skyvalue = 1068 #ID_1532
skyvalue = 1057 #ID_2211
skyvalue = 1080 #ID_0641
skyvalue = 1081 #ID_0348 LOOK VERY FALSE, but also he is an "extremist"
skyvalue = 1060 #ID_1891
skyvalue = 1065 #ID_1481
skyvalue = 1061 #ID_2196
skyvalue = 1112 #ID_1072
skyvalue = 1070 #ID_2131 Is excentered. So I change the center positions in my code, and I have 2.5, for 2.3. So ok, actually.
skyvalue = 1060 #ID_2100
skyvalue = 1087 #ID_1414 Seems to be maybe a star, so I check the diff with and without. Pymorph gives 2.38. I have 2.61 deleting the star, and 2.31 keeping it: no big deal.
skyvalue = 1098 #ID_0002
skyvalue = 1055 #ID_1299 As a bright star compared to the galaxy. Pymorph: 2.97, Me with star: 3.8 (huge difference). Me without the big star: 3.2. If i also delete a smaller (but still important) star: 3.09, so close to pymorph. COnclusion : PyMorph does a good job (for this one at least)!
skyvalue = 1060 #ID_2399
skyvalue = 1086 #ID_2689 PyMorph : 4.45, me: 2.85.  I had to delete stars and recenter. Obviously here pymorph is wrong. But that's an extrem again.
skyvalue = 1083 #ID_2510. Extreme in pymorph = 4.34. I have 3.18, so a large difference. I don't understand why, seeing the object. At least, all values are high. But still...
skyvalue = 1074.6 #ID_0818 ANother extree value in PyMorph. I measure 2.1 instead of 4.3. Is that because I choose a wrong center ? Would it vary a lot with another center ?
skyvalue = 1061 #ID_2130 The object is blended with another galaxy (z=0.0367 vs z=0.0341). This companion (?) galaxy is very strong. If I don't zero it, I have 2.3. If I delete it, I have 3.49. Pymorph gives 4.24.
skyvalue = 1072 #ID_2572
skyvalue = 1062 #ID_2233. I kill just a small star, improves the result just a little, but with or without it is very good...
skyvalue = 1070 #ID_2576
skyvalue = 1088 #ID_0222
skyvalue = 1082 #ID_0452. I kill a small star, and values are perfect.
skyvalue = 1096 #ID_0304
skyvalue = 1064 #ID_2427
skyvalue = 1092 #ID_2048
skyvalue = 1064 #ID_0112
skyvalue = 1064 #ID_2692
skyvalue = 1053 #ID_1258
skyvalue = 1082 #ID_0329
skyvalue = 1064 #ID_1732 OK when star killed
skyvalue = 1075 #ID_0724 Wrong and I have no idea why...
skyvalue = 1069 #ID_1805
skyvalue = 1058 #ID_1974
skyvalue = 1069 #ID_1805
skyvalue = 1074 #ID_1633
skyvalue = 1095 #ID_1852
skyvalue = 1107 #ID_2520
skyvalue = 1111 #ID_0381
# from here random
skyvalue = 1086 #ID_1832
skyvalue = 1061 #ID_1533
skyvalue = 1096.5 #ID_0225
skyvalue = 1059 #ID_2312 Had to change center
skyvalue = 1053 #ID_1243 2 objects, a bit complex. i also deleted a star.
skyvalue = 1078 #ID_0428
skyvalue = 1104 #ID_0335







print np.mean(mini[0:20, 0:20])

#skyvalue = 1077

mini = mini - skyvalue

dist = mini *0.

if test_name == 'ID_1414':
    mini[9:29, 46:67] = mini[9:29, 46+20:67+20] #for ID_1414

if test_name == 'ID_1299':
    mini[72:87, 13:28] = mini[72+15:87+15, 13:28] #for ID_1299
    mini[92:105, 43:60] = mini[92+15:105+15, 43:60] #for ID_1299 (second star)


if test_name == 'ID_2689':
    mini[8:27, 75:95] = mini[8+20:27+20, 75:95] #for ID_2689
    mini[80:95, 75:88] = mini[80:95, 75+15:88+15] #for ID_2689 (second star)

if test_name == 'ID_0818':
    mini = mini[10:81, :] #for ID_0818 , band to kill star.

if test_name == 'ID_2130':
    mini[0:40, 50:-1] = 0. #ID_2130, kill companion galaxy.

if test_name == 'ID_2233':
    mini[40:55, 80:] = 0. #ID_2233, kills just a small star, improves the result just a little, but with or without it is very good...

if test_name == 'ID_0452':
    mini[40:55, :15] = 0. ##ID_0452 small star

if test_name == 'ID_2692':
    mini[60:80, :25] = 0. ###ID_2692 star

if test_name == 'ID_1732':
    mini[:20, 60:] = 0. ###ID_1732 star

if test_name == 'ID_1633':
    mini[:, 0:20] = 0. ###ID_1633 object at the border

if test_name == 'ID_1243':
    mini[40:55, 0:15] = 0. ####ID_1243 star at the border

if test_name == 'ID_0381':
    mini[34:53, 88:103] = 0. ###ID_0381 star
    mini[102:, 0:40] = 0. ###ID_0381 2nd order small object far




centerx = len(mini[:,0])/2
centery = len(mini[0,:])/2

if test_name == 'ID_2131':
    centery = 14.7 #for ID_2131
    centerx = 49.25 #for ID_2131

if test_name == 'ID_2689':
    centery = 66 #for ID_2689
    centerx = 102 #for ID_2689

if test_name == 'ID_1633':
    centery = 53 #for ID_1633
    centerx = 35 #for ID_1633

if test_name == 'ID_2312':
    centery = 26 #for #ID_2312
    centerx = 48 #for #ID_2312

for i in np.arange(len(mini[:, 0])):
    for j in np.arange(len(mini[0, :])):
        dist[i,j] = np.sqrt((i-centerx)*(i-centerx)+(j-centery)*(j-centery))

radii = np.linspace(1, min(len(mini[:, 0])/1.5, len(mini[0, :])/1.5))
integrated_light = np.array([])
for r in radii:
    integrated_light = np.append(integrated_light, np.sum(mini[np.where(dist <= r)]))

integrated_light = integrated_light/max(integrated_light)

inner_radius = 0.2
outer_radius = 0.8
Rin = interp(inner_radius, integrated_light, radii)
Rout = interp(outer_radius, integrated_light, radii)

plt.figure()
plt.plot(radii, integrated_light)
plt.plot([Rin, Rin], [0, 1])
plt.plot([Rout, Rout], [0, 1])
plt.show()
plt.close()


plt.figure()
plt.imshow(mini)
plt.contour(dist)
plt.show()
plt.close()


Concentration = 5. * np.log10(Rout/Rin)
print test_name
print "My method :", Concentration
print "PyMorph :\n", large_table[np.where(large_table['ID'] == test_name)][sort_index]


sys.exit()





Concentrations = []
IDs = []
pymorph_concentrations = []

for line in large_table:

    #print line['ID']
    ellfit = Table.read(main_directory+line['ID']+"/E__"+line['ID']+".txt", format="ascii")
    #print ellfit.columns

    integral = np.array([])
    for r in np.arange(len(ellfit)):
        integral = np.append(integral, np.sum(ellfit['inte'][0:r+1]))
    integral = integral/max(integral)
    f_inverse = interpolate.interp1d(integral, ellfit['sma'])

    #print "integral = ", integral

    inner_radius = 0.2
    outer_radius = 0.8

    Rin = interp(inner_radius, integral, ellfit['sma'])
    Rout = interp(outer_radius, integral, ellfit['sma'])

    '''
    r = np.linspace(1, max(ellfit['sma']), 0.1)
    plt.figure()
    #plt.plot(ellfit['sma'], ellfit['inte'])
    plt.plot(ellfit['sma'], integral, 'x-')
    plt.plot([Rin, Rin], [0, 1])
    plt.plot([Rout, Rout], [0, 1])
    plt.show()
    plt.close()
    '''

    Concentration = 5. * np.log10(Rout/Rin)
    #print Concentration
    #print line['ID']
    #print line[sort_index]

    Concentrations.append(Concentration)
    IDs.append(line['ID'])
    pymorph_concentrations.append(line[sort_index])

concentration_table = Table([IDs, pymorph_concentrations, Concentrations], names=('ID', 'PymorphC', 'MyC'), meta={'name': 'comparison of concentrations'})

print concentration_table

plt.figure()
#plt.plot(ellfit['sma'], ellfit['inte'])
plt.plot(concentration_table['PymorphC'], concentration_table['MyC'], 'x')
plt.plot([0,5], [0,5])
plt.show()
plt.close()

sys.exit()


number_of_objects = len(large_table)
leny = 2
lenx = int(np.floor(number_of_objects/leny+1)) #objects per row

print lenx, leny
print len(large_table)


plt.figure(figsize=[leny, lenx])
plt.suptitle(sort_index+'\nfrom '+str(large_table[0][sort_index])+'\nto '+str(large_table[-1][sort_index]), fontsize=14)
gs1 = gridspec.GridSpec(lenx, leny)
gs1.update(wspace=0., hspace=0.) # set the spacing between axes.
plt.set_cmap('Set1')

#large_table = large_table[10]

i=0
for line in large_table:

    hdulist = fits.open(main_directory+line['ID']+'/'+line['ID']+'_'+filter+'_mini.fits')
    mini = hdulist[0].data
    hdulist.close()


    if i < number_of_objects/2+1:
        ax1 = plt.subplot(gs1[2*i])
    else:
        ax1 = plt.subplot(gs1[number_of_objects-2*i+1])
    ax1.imshow(np.log10(mini))
    #ax1.text(0,-0.2, line[sort_index])
    ax1.annotate("{0:.2f}".format(line[sort_index]), xy=(1, 0), xycoords='axes fraction', fontsize=10,horizontalalignment='right', verticalalignment='bottom')
    ax1.annotate(line['ID'], xy=(0.8, 0.8), xycoords='axes fraction', fontsize=10, horizontalalignment='right', verticalalignment='bottom')
    #axarr[i, j].axes.get_xaxis().set_visible(False)
    #axarr[i, j].axes.get_yaxis().set_visible(False)

    #axarr[i, j].set(gca, 'XTickLabel', [],'XTick',[])

    plt.axis('off')
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_aspect('equal')

    i+=1

# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
#plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)

plt.savefig(figdir+sort_index+"_all_face2face.png")
plt.show()
plt.close()




sys.exit()


number_of_objects = len(large_table)
lenx = 20 #objects per line
leny = int(np.floor(number_of_objects/lenx+1))

plt.figure(figsize=[lenx*3, leny*3])
plt.suptitle(sort_index+' from '+str(large_table[0][sort_index])+' to '+str(large_table[-1][sort_index]), fontsize=140)
gs1 = gridspec.GridSpec(lenx, leny)
gs1.update(wspace=0., hspace=0.) # set the spacing between axes.
plt.set_cmap('Set1')

#large_table = large_table[10]

i=0
for line in large_table:
    hdulist = fits.open(main_directory+line['ID']+'/'+line['ID']+'_'+filter+'_mini.fits')
    mini = hdulist[0].data
    hdulist.close()

    ax1 = plt.subplot(gs1[i])
    ax1.imshow(np.log10(mini))
    #axarr[i, j].axes.get_xaxis().set_visible(False)
    #axarr[i, j].axes.get_yaxis().set_visible(False)

    #axarr[i, j].set(gca, 'XTickLabel', [],'XTick',[])

    plt.axis('off')
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_aspect('equal')

    i+=1

# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
#plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)

plt.savefig(figdir+sort_index+"_all.png")
plt.show()
plt.close()





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


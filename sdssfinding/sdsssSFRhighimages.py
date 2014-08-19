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
     name = "images/" + str(MJD[ind])+str(PLATEID[ind])+str(FIBERID[ind]) + ".jpg"
#     spRAsec = str(RAsec[ind]).split('.')
#     spDecsec = str(Decsec[ind]).split('.')
#     if (len(spRAsec[0]) < 2):
#          strsize = len(spRAsec[0])+len(spRAsec[1])+2
#     else:
#	  strsize = len(spRAsec[0])+len(spRAsec[1])+1
##
#     if (len(spDecsec[0]) < 2):
#          strsize2 = len(spDecsec[0])+len(spDecsec[1])+2
#     else:
#	  strsize2 = len(spDecsec[0])+len(spDecsec[1])+1
#   
##
#     RAstr = str(RAhrs[ind]).zfill(2)+':'+str(RAmin[ind]).zfill(2)+':'+ str(RAsec[ind]).zfill(strsize)
#     if (Decdeg[ind] < 0):
#          Decstr = str(Decdeg[ind]).zfill(3)+':'+str(Decmin[ind]).zfill(2)+':'+str(Decsec[ind]).zfill(strsize2)
#     else:
#          Decstr = str(Decdeg[ind]).zfill(2)+':'+str(Decmin[ind]).zfill(2)+':'+str(Decsec[ind]).zfill(strsize2)
     funcfindingchart.getSDSSfindingchart(str(RAdeg[ind]),str(Decdeg[ind]),name)

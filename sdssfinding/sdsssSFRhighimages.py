#
# program to write out the sample in a form suitable for downloading the sdss
# jpegs 
#
import asciidata
import numpy as np
from astLib import astCoords
import funcfindingchart
#
# read in the data
#
highsSFRsample = asciidata.open('../NRTsamplewidezall')
MJD = highsSFRsample[0]
PLATEID = highsSFRsample[1]
FIBERID = highsSFRsample[2]
PHOTID = highsSFRsample[3]
RAhrs= highsSFRsample[4]
RAmin = highsSFRsample[5]
RAsec = highsSFRsample[6]
Decdeg = highsSFRsample[7]
Decmin = highsSFRsample[8]
Decsec = highsSFRsample[9]
#
# read in the ra and decs in decimal degress
#
origcat = '../spectra/NRTwideselectionSDSSdata.dat'
origdata = asciidata.open(origcat)
SDSS = origdata[0]
PLATEIDorig = origdata[1]
MJDorig = origdata[2]
FIBERIDorig = origdata[3]
RAdeg =  origdata[4].tonumpy()
Decdeg =  origdata[5].tonumpy()
z = origdata[6].tonumpy()
zwarning = origdata[7].tonumpy()
StoN = origdata[8].tonumpy()
objclass= origdata[9]
#
#
#
for ind in np.arange(len(MJD)):
     name = str(MJD[ind])+str(PLATEID[ind])+str(FIBERID[ind])+'.jpg'
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

def getSDSSfindingchart(RA,Dec,name,Scale='0.39612',Width='512'):
#
# subroutine to get finding charts from the SDSS
# 
     import urllib2, urllib
#
#
#
     FindingChartURL = 'http://skyserver.sdss3.org/dr8/en/tools/chart/chart.asp'
     FindingChartNaviURL = 'http://skyserver.sdss3.org/dr8/en/tools/chart/navi.asp'
     FindingChartImURL = 'http://casjobs.sdss.org/ImgCutoutDR7/getjpeg.aspx'
#
#
#
     ImURL   = FindingChartImURL+'?ra='+RA+'&dec='+Dec+'&opt=&scale='+Scale+'&width='+Width+'&height='+Width
     urllib.urlretrieve(ImURL, filename=name)

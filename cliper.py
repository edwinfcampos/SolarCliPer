#!/opt/local/bin/python 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import time
import re, string, sys, os
import time as timing
import numpy as np
import cPickle as pickle
#import pprint 
from datetime import datetime,timedelta
import urllib2, urllib
import ConfigParser

#Config = ConfigParser.RawConfigParser()   
Config = ConfigParser.ConfigParser()
Config.read('cliper.ini');

AllSitesList=Config.sections()

print 'Number of sites to read: '+str(len(AllSitesList)) 
print '----------------------------------'
site_no=1
for iSite in AllSitesList:
    fcn_str=Config.get(iSite,'Function')
    print 'Site ('+str(site_no)+'/'+str(len(AllSitesList))+'): '+iSite+'; function: '+fcn_str
    execfile(fcn_str)
    print '----------------------------------'
    site_no+=1




#!/opt/local/bin/python 

import re, string, sys, os
import time as timing
import math
from numpy import linalg as LA 

import numpy as np
import cPickle as pickle
import pprint 
from datetime import datetime,timedelta
import urllib2, urllib
import csv


def DateToJulian(dt):
    dtuple=dt.timetuple()
    julday=dtuple.tm_yday

    ordinal=dt.toordinal()
    Dhr=dt.hour
    Dmin=dt.minute
    juldate=1721422.5 + ordinal + 2 + Dhr/24. + Dmin/(24*60.)

    return juldate


def JulianToDate(juldate):
    juldate_shift=juldate-0.5+(0.6/(24*3600))
    ordinal= int(np.floor(juldate_shift - 1721422.0 - 2))
    dt=datetime.fromordinal(ordinal)
    Dhour=int(np.floor((juldate_shift-np.floor(juldate_shift))*24)) 
    Dminute=int(np.floor((((juldate_shift-np.floor(juldate_shift))-(Dhour/24.))*(24*60))+0.5))

    if(Dminute==60):
        Dhour+=1
        Dminute=0
        

    dtx=datetime(dt.year,dt.month,dt.day,Dhour,Dminute)
    return dtx


def JulianToDateVec(juldate_in):
    juldate_in=np.array(juldate_in)
    vlen=len(juldate_in)
    dtx_out=list()
    for i in range(vlen):
        dtx_out.append(JulianToDate(juldate_in[i]))
    return dtx_out


def all_indices(value, qlist):
    indices = []
    idx = -1
    while True:
        try:
            idx = qlist.index(value, idx+1)
            indices.append(idx)
        except ValueError:
            break
    return indices

def all_indices3(value1,qlist1,value2,qlist2,value3,qlist3):
    indices = []
    idx = -1
    while True:
        try:
            idx = qlist1.index(value1, idx+1)
            idx = qlist2.index(value2, idx)
            idx = qlist3.index(value3, idx)
            if(qlist1[idx]==value1 and qlist2[idx]==value2 and qlist3[idx]==value3):
                indices.append(idx)
        except ValueError:
            break
    return indices

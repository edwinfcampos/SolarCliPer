#!/opt/local/bin/python 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.clf()

import math
import numpy as np
import os
import csv
import cPickle as pickle
from datetime import datetime,timedelta
from utility import JulianToDate,DateToJulian
import urllib
import time

class Climate_record(object):
    DateJulian=[]
    Dyear=[]
    Dmonth=[]
    Dday=[]
    Dhour=[]
    Dminute=[]
    ANLIrradiance=[] 
    ANLPressure=[]
    ANLTemperature=[]
    pass

data_file_txt=Config.get(iSite,'Data_file_txt')
data_file_bin=Config.get(iSite,'Data_file_bin')

months=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
url_anl_data='http://www.atmos.anl.gov/ANLMET/numeric/15Minute/'

start_year=Config.getint(iSite,'Start_year')
start_month=Config.getint(iSite,'Start_month')

end_year=Config.getint(iSite,'End_year')
end_month=Config.getint(iSite,'End_month')

print '----------------------------------'
print 'Argonne site: '+Config.get(iSite,'Function')
print 'Climatology record: '+str(start_month)+'/'+str(start_year)+' - '+str(end_month)+'/'+str(end_year)
print 'Forecast horizon: '+str(Config.getint(iSite,'Forecast_Hours'))+' hours'
#print '----------------------------------'

years=end_year-start_year

start_mo=start_month

if (years==0):
    end_mo=end_month
else:
    end_mo=12


def MeanThreshold(WW,Thr):
    n=len(WW)
    meanT=0.0
    meancnt=0
    for k in range(0,n):
        if(WW[k]>Thr):
            meanT=meanT+WW[k]
            meancnt+=1
    if(meancnt==0):
        rv=0
    else:
        rv=meanT/meancnt
    return rv

def GetIrradiaceClimatology(rightnow,start_year,end_year,CR):
    from utility import all_indices3
    
    arrdate=np.array([],dtype=float)
    arrirr=np.array([],dtype=float)
    
    for iyears in range(start_year,end_year+1):
        idx=all_indices3(iyears,CR.Dyear,rightnow.month,CR.Dmonth,rightnow.day,CR.Dday)
        if(len(idx)>0):
            if(arrdate.size==0):
                if(CR.Dday[idx[0]]==1 and CR.Dmonth[idx[0]]==1):
                    arrdate=np.array([DateToJulian(datetime(CR.Dyear[idx[0]],CR.Dmonth[idx[0]],CR.Dday[idx[0]],0,0))])
                    arrdate=np.concatenate((arrdate,[CR.DateJulian[i] for i in idx]))               
                else:
                    arrdate=np.concatenate((arrdate,[CR.DateJulian[i] for i in idx]))
            else:
                stack_v=np.array([CR.DateJulian[i] for i in idx])
                arrdate=np.vstack((arrdate,stack_v))
            if(arrirr.size==0):
                if(CR.Dday[idx[0]]==1 and CR.Dmonth[idx[0]]==1):
                    arrirr=np.array([0])
                    arrirr=np.concatenate((arrirr,[CR.ANLIrradiance[i] for i in idx]))
                else:
                    arrirr=np.concatenate((arrirr,[CR.ANLIrradiance[i] for i in idx]))
            else:
                arrirr=np.vstack((arrirr,[CR.ANLIrradiance[i] for i in idx]))
    
    arrirrm=np.zeros(arrirr.shape[1])
    arrirrmax=np.zeros(arrirr.shape[1])
    arrirrcnt=np.zeros(arrirr.shape[1])
    
    for icomp in range(arrirr.shape[1]):
        #cnt=0
        for ivec in range(arrirr.shape[0]):
            if(arrirr[ivec,icomp]<>99999):
                if(arrirr[ivec,icomp]<1):
                    arrirr[ivec,icomp]=0
                arrirrcnt[icomp]+=1
                arrirrm[icomp]+=arrirr[ivec,icomp]
                arrirrmax[icomp]=max(arrirrmax[icomp],arrirr[ivec,icomp])
        arrirrm[icomp]= arrirrm[icomp]/arrirrcnt[icomp]
    
    #replace bad values with the expected value
    for icomp in range(arrirr.shape[1]):
        for ivec in range(arrirr.shape[0]):
            if(arrirr[ivec,icomp]==99999):
                arrirr[ivec,icomp]=arrirrm[icomp]
                
    return (arrdate,arrirr)


def GetIrradiaceANLobservations(TodayDOY):
    url_anl_data='http://www.atmos.anl.gov/ANLMET'
    data48file='anltower.48'
    string_wget=url_anl_data+'/'+data48file
    #os.system('wget  -q ' + string_wget)
    urllib.urlretrieve(string_wget, filename=data48file)
    #if(not os.path.isfile(data48file)):
    #   raise NameError('HiThere')
    
    
    ANLTime=[]
    ANLIrradiance=[]
    
    fid=open(data48file,'r')
    p=csv.reader(fid,delimiter=" ")
    t=p.next()
    t=p.next()
    while 1:
        t=p.next()
        if (t[0]=='JDA'):
            break
        if (int(t[0])<>TodayDOY):
            continue
        g=''    
        for i in range(len(t)):
            g=g+t[i]+'^'
        g=g.replace('^^^^^^',' ')
        g=g.replace('^^^^^',' ')
        g=g.replace('^^^^',' ')
        g=g.replace('^^^',' ')
        g=g.replace('^^',' ')
        g=g.replace('^',' ')
        
        t=g.split(' ')
        
        
        JDA=int(t[0])
        sr=t[1]
        julhr=int(float(sr[sr.find('+')+1:sr.find('+')+3]))
        julmin=int(float(sr[sr.find(':')+1:]))
        sr=t[16]
        ANLTime.append(julhr+julmin/60.)
        irr=max(float(sr),0)
        if irr<1:
            irr=0.0
        ANLIrradiance.append(irr)
    
    fid.close()
    #This works for UNIX
    #os.system('rm -fr anltower.48*')
    #This works for MS-Windows
    #os.system('del anltower.48*')
    os.remove('anltower.48')

    return (ANLTime,ANLIrradiance)



if(os.path.isfile(data_file_bin)):
    print 'The climatology data is already available.'
else:
    print 'The climatology data is not available; downloading data:'
    str_file_out=data_file_txt

    fout=open(str_file_out,'w')
    pid=csv.writer(fout,delimiter=" ")

    CR=Climate_record

    for iyears in range(start_year,end_year+1):
        if (iyears==end_year):
            end_mo=end_month
        for imonths in range(start_mo-1,end_mo):
            string_fname=months[imonths]+str(iyears)[2:4]+'met.dat'
            string_wget=url_anl_data+str(iyears)+'/'+string_fname
            #os.system('wget -q ' + string_wget)
            urllib.urlretrieve(string_wget, filename=string_fname)
            sys.stdout.write('+')
            sys.stdout.flush()
            fid=open(string_fname,'a')
            fid.write('EOF')
            fid.write('\n')       
            fid.close()

            fid=open(string_fname,'r')
            p=csv.reader(fid,delimiter=" ")

            while 1:
                try:
                    t=p.next()
                except StopIteration:
                    print 'found exception'
                    break

                if (t[0]=='EOF'):
                    break
                if (t[1]!='M'):
                    continue
                g=''    
                for i in range(len(t)):
                    g=g+t[i]+'^'
                g=g.replace('^^^^^',' ')
                g=g.replace('^^^^',' ')
                g=g.replace('^^^',' ')
                g=g.replace('^^',' ')
                g=g.replace('^',' ')

                t=g.split(' ')

                rightnow=datetime.now()
                current_year_shift=datetime(iyears,1,1)

                srday=t[2]  
                Dday=int(float(srday[srday.find('+')+1:]))
                srmo=t[3] 
                Dmonth=int(float(srmo[srmo.find('+')+1:]))
                sr=t[5] 
                Dhr=int(float(sr[sr.find('+')+1:sr.find('+')+3]))
                Dmin=int(float(sr[sr.find('+')+1+2:sr.find('+')+3+2]))           

                if (Dhr==24):
                    Dhr=0
                    ddelta = timedelta(days=+1)
                else:
                    ddelta = timedelta(days=+0)

                Ddate=datetime(iyears,Dmonth,Dday,Dhr,Dmin)
                Ddate=Ddate+ddelta

                DJuldate=DateToJulian(Ddate)

                t1=p.next()
                g=''
                for i in range(len(t1)):
                    g=g+t1[i]+'^'
                g=g.replace('^^^^^',' ')
                g=g.replace('^^^^',' ')
                g=g.replace('^^^',' ')
                g=g.replace('^^',' ')
                g=g.replace('^',' ')

                t1=g.split(' ')
                g=''
                t2=p.next()
                for i in range(len(t2)):
                    g=g+t2[i]+'^'

                g=g.replace('^^^^^',' ')
                g=g.replace('^^^^',' ')
                g=g.replace('^^^',' ')
                g=g.replace('^^',' ')
                g=g.replace('^',' ')

                t2=g.split(' ')

                allfields=t+t1+t2

                CR.DateJulian.append(DJuldate)
                CR.Dyear.append(Ddate.year)
                CR.Dmonth.append(Ddate.month)
                CR.Dday.append(Ddate.day)
                CR.Dhour.append(Ddate.hour)
                CR.Dminute.append(Ddate.minute)
                CR.ANLIrradiance.append(float(allfields[23]))
                srt=allfields[28]  
                temp=float(srt[srday.find('+')+1:])
                CR.ANLPressure.append(temp)

                srt=allfields[18]  
                temp=float(srt[srday.find('+')+1:])

                CR.ANLTemperature.append(temp)
                pid.writerow([DJuldate,Ddate.year,Ddate.month,Ddate.day,Ddate.hour,Ddate.minute,float(allfields[23]),float(allfields[23+5]),(allfields[18])])
                dd=JulianToDate(DJuldate)
            fid.close()
            #os.system('rm *met.dat')
            os.remove(string_fname)
        sys.stdout.write('$'+str(iyears)+'$\n')
        sys.stdout.flush()
        start_mo=1
        
    fout.close()


    outpcl = open(data_file_bin, 'wb')
    pickle.dump(CR.DateJulian,outpcl)#,protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(CR.Dyear,outpcl)
    pickle.dump(CR.Dmonth,outpcl)
    pickle.dump(CR.Dday,outpcl)
    pickle.dump(CR.Dhour,outpcl)
    pickle.dump(CR.Dminute,outpcl)
    pickle.dump(CR.ANLIrradiance,outpcl)
    pickle.dump(CR.ANLPressure,outpcl)
    pickle.dump(CR.ANLTemperature,outpcl)
    outpcl.close()

    


print 'Reading data from file'
CR=Climate_record
indpcl = open('data_anl_001.pcl', 'rb')
CR.DateJulian=pickle.load(indpcl)
CR.Dyear=pickle.load(indpcl)
CR.Dmonth=pickle.load(indpcl)
CR.Dday=pickle.load(indpcl)
CR.Dhour=pickle.load(indpcl)
CR.Dminute=pickle.load(indpcl)
CR.ANLIrradiance=pickle.load(indpcl)
indpcl.close()
print 'Done reading the data'

rightnow=datetime.now()
if time.localtime().tm_isdst==1:
    rightnow=rightnow+timedelta(hours=-1)

dt=rightnow.timetuple()
TodayDOY=dt.tm_yday
OneDay=timedelta(days=1)

Yesterday=rightnow-OneDay
YesterdayT=Yesterday.timetuple()
YesterdayDOY=YesterdayT.tm_yday

Tomorrow=rightnow+OneDay
TomorrowT=Tomorrow.timetuple()
TomorrowDOY=TomorrowT.tm_yday

(arrdate,arrirr) = GetIrradiaceClimatology(rightnow,start_year,end_year,CR)
(arrdateY,arrirrY) = GetIrradiaceClimatology(Yesterday,start_year,end_year,CR)
(arrdateT,arrirrT) = GetIrradiaceClimatology(Tomorrow,start_year,end_year,CR)

arrirrm=arrirr.mean(0)
arrirrmax=arrirr.max(0)

arrirrYm=arrirrY.mean(0)
arrirrYmax=arrirrY.max(0)

arrirrTm=arrirrT.mean(0)
arrirrTmax=arrirrT.max(0)

ANLClimIrr=np.concatenate([arrirrYm.ravel(),arrirrm.ravel(),arrirrTm.ravel()])
ANLClimIrrMax=np.concatenate([arrirrYmax.ravel(),arrirrmax.ravel(),arrirrTmax.ravel()])
ANLClimTime=np.concatenate([arrdateY[0][:].ravel(),arrdate[0][:].ravel(),arrdateT[0][:].ravel()])

(ANLTime,ANLIrradiance) = GetIrradiaceANLobservations(TodayDOY)
(ANLTimeY,ANLIrradianceY) = GetIrradiaceANLobservations(YesterdayDOY)

ANLObsTime=np.concatenate([np.array(ANLTimeY).ravel()-24.0,np.array(ANLTime).ravel()])
ANLObsIrr=np.concatenate([np.array(ANLIrradianceY).ravel(),np.array(ANLIrradiance).ravel()])


LenANLObs=len(ANLObsTime)

ANLObsClimIrr=ANLObsIrr-ANLClimIrr[0:LenANLObs]

idxObsS=([])
idxObsE=([])
DistObs=([])
IrrAutocorr=np.zeros((11,16))
IrrObsClimAutocorr=np.zeros((11,16))

ForecastHours=Config.getint(iSite,'Forecast_Hours')
HindCastRange=range(1,2*ForecastHours*4,2)
HindCasts=len(HindCastRange)

ANLForecastIrr=np.zeros((ForecastHours*4))
ANLForecastTime=np.zeros((ForecastHours*4))

ANLHindCastIrr=np.zeros((ForecastHours*4,HindCasts+1))
ANLHindCastTime=np.zeros((ForecastHours*4,HindCasts+1))

ANLForecastIrr=np.minimum(np.maximum(ANLClimIrr[LenANLObs-1:LenANLObs-1+ForecastHours*4]-(ANLClimIrr[LenANLObs-1]-ANLObsIrr[LenANLObs-1]),0),ANLClimIrrMax[LenANLObs-1:LenANLObs-1+ForecastHours*4])
ANLForecastTime=ANLClimTime[LenANLObs-1:LenANLObs-1+ForecastHours*4]

##################


ANLTodayClimatology=arrirrm.ravel()
ANLTodayClimatologyTime=arrdate[0][:].ravel()
ANLTodayClimatologyMax=arrirrmax.ravel()
ANLTodayObservations=np.array(ANLIrradiance).ravel()
ANLTodayObservationsTime=np.array(ANLTime).ravel()
ANLTodayObservationLen=len(ANLTodayObservations)
ANLMAEClimatology=MeanThreshold(np.absolute(ANLTodayClimatology[0:ANLTodayObservationLen]-ANLTodayObservations),1.0)
ANLERRORForecastLen=np.minimum(6*4,ANLTodayObservationLen)

ANLCLIPERRange=range(ANLTodayObservationLen)
ANLCLIPERCasts=len(ANLCLIPERRange)
ANLCLIPERCastIrr=np.zeros((ANLERRORForecastLen,ANLCLIPERCasts))
ANLCLIPERCastTime=np.zeros((ANLERRORForecastLen,ANLCLIPERCasts))
ANLMAECLIPERCastIrr=np.zeros((ANLERRORForecastLen))
ANLMAECLIPERCastIrrNo=np.zeros((ANLERRORForecastLen))

for idx in range(ANLCLIPERCasts):
    i=ANLCLIPERRange[idx]
    idxset_first=i
    idxset_last=np.minimum(24*4-1,i+ANLERRORForecastLen)
    ANLCLIPERCastIrr[0:idxset_last-idxset_first,idx]=np.minimum(np.maximum(ANLTodayClimatology[idxset_first:idxset_last]-(ANLTodayClimatology[i]-ANLTodayObservations[i]),0),ANLTodayClimatologyMax[idxset_first:idxset_last])
    ANLCLIPERCastTime[0:idxset_last-idxset_first,idx]=ANLTodayClimatologyTime[idxset_first:idxset_last]

for k in range(ANLERRORForecastLen):
    for j in range(k,ANLTodayObservationLen):
        if(ANLTodayObservations[j]>1.0):
            ANLMAECLIPERCastIrr[k]+=np.absolute(ANLCLIPERCastIrr[k,j-k]-ANLTodayObservations[j])
            ANLMAECLIPERCastIrrNo[k]+=1

for k in range(ANLERRORForecastLen):
    if(ANLMAECLIPERCastIrrNo[k]==0):
        ANLMAECLIPERCastIrr[k]=0
    else:
        ANLMAECLIPERCastIrr[k]=ANLMAECLIPERCastIrr[k]/ANLMAECLIPERCastIrrNo[k]

##################


for idx in range(HindCasts):
    i=HindCastRange[idx]
    ANLHindCastIrr[:,idx]=np.minimum(np.maximum(ANLClimIrr[LenANLObs-1-i:LenANLObs-1+ForecastHours*4-i]-(ANLClimIrr[LenANLObs-1-i]-ANLObsIrr[LenANLObs-1-i]),0),ANLClimIrrMax[LenANLObs-1-i:LenANLObs-1+ForecastHours*4-i])
    ANLHindCastTime[:,idx]=ANLClimTime[LenANLObs-1-i:LenANLObs-1+ForecastHours*4-i]

#### CliPer Figure ####

rightnow=datetime.now()
YYYY=rightnow.year
current_year_shift=datetime(YYYY,1,1)
d = datetime(rightnow.year, rightnow.month, rightnow.day, 0, 0, 0)
JulD=DateToJulian(d)

plt.plot(list((arrdate[0][:]-arrdate[0][0])*24),list(arrirrm),'b-',label='Climatology')
plt.plot(list((arrdate[0][:]-arrdate[0][0])*24),list(arrirrmax),'b--',label='Max clim.')

for idx in range(HindCasts):
    if(idx==0):
        plt.plot(list((ANLHindCastTime[:,idx]-arrdate[0][0])*24),list(ANLHindCastIrr[:,idx]),'g-',label='CliPer (prev)')
    else:
        plt.plot(list((ANLHindCastTime[:,idx]-arrdate[0][0])*24),list(ANLHindCastIrr[:,idx]),'g-')

plt.plot(list((ANLForecastTime-arrdate[0][0])*24),list(ANLForecastIrr),'r-',label='CliPer (curr)', lw=2)

#print list((ANLForecastTime-arrdate[0][0])*24)

plt.plot(ANLTime,ANLIrradiance,'k*',label='Observations')

plt.xlim((0,24))
datestr=str(rightnow.day)+'-'+months[rightnow.month-1].title()+'-'+str(rightnow.year)+'\nArgonne Weather\nObservatory'
plt.text(1,max(arrirrmax)/1.6,datestr, fontsize=15)
plt.xlabel('Time of day [CST]\nArgonne National Laboratory - ecampos@anl.gov\n ')
plt.ylabel('Global solar irradiance [W/m^2]')
plt.legend()

#plt.text(.5,1.015,'Argonne National Laboratory - ecampos@anl.gov',horizontalalignment='center')

fname=Config.get(iSite,'SaveGraphFile')
print 'Saving CliPer output in '+fname
plt.savefig(fname)

#print 'Cliper error(MAE): ' + str(ANLMAECLIPERCastIrr)

#Code originally developed on 2014 May 22 by: Emil Constantinescu, emconsta@mcs.anl.gov
#Recent modifications made by: Edwin Campos, ecampos@anl.gov
#last modified: July 9, 2014

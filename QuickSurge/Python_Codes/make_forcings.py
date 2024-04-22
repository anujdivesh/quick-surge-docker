# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 11:01:07 2023

@author: moritzw
"""
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 13:27:04 2023
Script to generate fort.19 files
@author: moritzw
"""
import numpy as np
import numpy.matlib
import pandas as pd
import os
import pyTMD
from adcircpy import AdcircMesh
import datetime as dt
import matplotlib.dates as mdates
import shutil
from scipy import interpolate
import re
import fnmatch
from  get_periodic_tidal_cond import get_periodic_tidal_cond
import scipy.interpolate as interp
import datetime
from shapely import geometry
###

def define_run_directory(pathres):
    folders = fnmatch.filter(os.listdir(pathres), 'Fore_*')
    if not folders: 
        run = 0
        ramp = True
    else:   
        ramp = False
        runs = [] 
        for i in folders:
            runs.append(int(re.findall(r'\d+', i)[0])) 
        run = np.max(runs)+1     
    newrun = 'Run_' + str(run).zfill(3)     
    return(newrun,ramp)

def make_run_folder(pathres):
    if not os.path.exists(pathres):
           os.mkdir(pathres)
           print("Directory " , pathres ,  " Created ")
    else:    
           print("Directory " , pathres ,  " already exists")
    return()

def correct_x(x):
    for i in range(len(x)):
        if x[i][-1] == 'W':
            x = x.replace(x[i][:],' '+ str(1800 + (1800 -  int(x[i][0:-1]))) + 'E')
    return(x)

def get_lonstr_from_x(x):
    lonstr = [];
    for i in range(len(x)):
        if x[i] < 0:
           x[i] = round((x[i]+360)*10)
        else:
           x[i] = round(x[i]*10)     
        lonstr.append(str(int(x[i]))+'E')
    return(lonstr)

def get_latstr_from_y(y):
    latstr = [];
    for i in range(len(y)):
         if y[i] < 0:
            y[i] = round((-y[i])*10)
            latstr.append(str(int(y[i]))+'S')
         else:
            y[i] = round(y[i]*10) 
            latstr.append(str(int(y[i]))+'N')
    return(latstr)


def get_lat_from_y(y):
    #lat = y.copy()
    lat = np.zeros(np.shape(y))
    for i in range(len(y)):
        #lat = lat.replace(lat[i][:], int(lat[i][0:-1])/10)
        lat[i] = int(y[i][0:-1])/10
    return(lat)

def get_rmax_with_knaff(Wind,lat):
   
    #Knaff et al. (2016)
    r_max = 218.3784 - 1.2014*Wind + (Wind/10.9844)**2 - (Wind/35.3052)**3 - (145.5090*np.cos(np.deg2rad(lat)))
    return(r_max)


def create_wind_file(wind_out,time,y,x,Wind,Pmin,Rmax):
    if os.path.exists(wind_out):
        os.remove(wind_out)
    frmt='SH, 01,%11s, 00, BEST,   0,%5s,%6s,%4s,%5s,   ,    ,    ,     ,     ,     ,     , 1013,    ,%4s,    ,    ,    ,    ,    ,    ,    ,    Unnamed,  ,   ,    ,    ,    ,    ,   0\n'
    with open(wind_out,'a') as wind_file:
        for i in range(len(time)):
            wind_file.write(frmt % (time[i],y[i],x[i],int(Wind[i]),int(Pmin[i]),int(Rmax[i])))
    return()

def get_x_and_y(grid_in):
    pmesh = AdcircMesh.open(grid_in,crs=4326)
    #pmesh = read_mesh(folder_name+"/fort.14")
    px = np.array(pmesh.x)
    py = np.array(pmesh.y)
    pbnd = pmesh.boundaries.ocean
    pbnd_val = pbnd.indexes[:]
    pbnd_val_ = []
    for sublist in pbnd_val :
        for item in sublist:
            pbnd_val_.append(item)
    
    pbnd_x = np.squeeze(px[pbnd_val_])
    pbnd_y = np.squeeze(py[pbnd_val_])
    return(pbnd_x,pbnd_y)

def get_tidal_elevation_matrix(tide_dir,start_time,end_time,pbnd_x,pbnd_y):
    grid_file = os.path.join(tide_dir,'TPXO8','DATA','grid_tpxo8atlas_30')## check if we can go for a newer version, others models can be used:'GOT4.8','FES2014',...
    model_file = os.path.join(tide_dir,'TPXO8','DATA','hf.tpxo8_atlas_30')
       
    model_format = 'OTIS'
    EPSG = '4326'
    TYPE = 'z'
    
    time_tide = mdates.drange(start_time,end_time,dt.timedelta(minutes=10))
    time_tmd=time_tide-mdates.date2num(np.datetime64('1992-01-01T00:00:00')) 
    
    z = np.zeros((len(pbnd_x),len(time_tmd)))
  
    ##LON,LAT
    # LON = pbnd_x[i]
    # LAT = pbnd_y[i]
    # amp,ph,D,c = pyTMD.io.extract_constants(np.array([LON]), np.array([LAT]),
    #                            grid_file,model_file,EPSG,TYPE=TYPE,METHOD='spline',GRID=model_format)
    
    amp,ph,D,c = pyTMD.io.extract_constants(pbnd_x, pbnd_y,
                               grid_file,model_file,EPSG,TYPE=TYPE,METHOD='spline',GRID=model_format,extrapolate = True,cutoff=50)

    #-- calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    #-- calculate constituent oscillation
    hc = amp*np.exp(cph)
    #-- predict tidal elevations at time 1 and infer minor corrections
    for i in range(len(pbnd_x)):
        print(str(i+1) + ' / ' + str(len(pbnd_x)))
        kk=hc[i,:]
        kk = kk[np.newaxis]
        TIDE = pyTMD.predict.time_series(time_tmd, kk, c)
        MINOR = pyTMD.predict.infer_minor(time_tmd, kk, c, CORRECTIONS=model_format)
        TIDE.data[:] += MINOR.data
        z[i,:] = TIDE.data[:]
    return(z)

# from matplotlib import pyplot as plt

# plt.figure(figsize=(25,10)) 
  

# kk =z[:,100] 
# kk[kk>2]=np.nan

# plt.plot(z[:,100] ,'-k',linewidth=1,label='Storm Tide')
# plt.plot(time_tmd,z[100,:] ,'-k',linewidth=1,label='Storm Tide')
# plt.plot(time_tmd,z[200,:] ,'-k',linewidth=1,label='Storm Tide')
# plt.plot(time_tmd,z[300,:] ,'-k',linewidth=1,label='Storm Tide')


def create_tide_file(tide_out,start_time,end_time,tide,pbnd_x,pbnd_y):
    time_tide = mdates.drange(start_time,end_time,dt.timedelta(minutes=10))
    if os.path.exists(tide_out):
            os.remove(tide_out)
    with open(tide_out,'a') as tide_file:
        tide_file.write('%i\n' % 600)
        for i in range(len(time_tide)):
            for k in range(len(pbnd_x)):
                tide_file.write('%4.4f\n' % tide[k,i])
    return()

# def generate_fort22_file(wind_in,pathres):
# # wind_in = 'F:/Adcirc_SWAN/PARTneR2/BestTrack/Harold_tonga.csv'
# # pathres='F:/Adcirc_SWAN/PARTneR2/Test_Runs/Test_06_all_py/'
# # generate_fort22_file(wind_in,pathres)
#     wind_out = pathres +'fort.22'
#     make_run_folder(pathres)
#     besttrack = pd.read_csv(wind_in,header=None,delimiter=',')
    
#     time = besttrack[2][:]
#     x = besttrack[7][:]
#     x = correct_x(x)
#     y = besttrack[6]
#     lat = get_lat_from_y(y)
    
#     Pmin = besttrack[9]
#     Wind = besttrack[8]
#     r_max = get_rmax_with_knaff(Wind,lat)
#     Rmax = r_max
#     create_wind_file(wind_out,time,y,x,Wind,Pmin,Rmax)
#     return()

def get_start_and_end_time_from_besttrack(f22):
    
    # f22 = pathrun +'fort.22'
    besttrack = pd.read_csv(f22,header=None,delimiter=',')
    besttrack[2] = besttrack[2].astype(str)
    # Series of datetime values from Column
    times = pd.to_datetime(besttrack[2], format='%Y%m%d%H').dt.to_pydatetime()
    
    start_time = times[0]   #dt.datetime(2020,4,5,00)
    start_time = start_time.replace(minute=0, second=0, microsecond=0)
    end_time = times[-1]    #dt.datetime(2020,4,10,13)
    end_time = end_time.replace(minute=0, second=0, microsecond=0)
    
    return(start_time,end_time)
    

def generate_fort19_file(input_folder,pathres,tide_model_dir):
# input_folder = 'F:/Adcirc_SWAN/PARTneR2/Test_Runs/input_files/'
# pathres='F:/Adcirc_SWAN/PARTneR2/Test_Runs/Test_06_all_py/'
# tide_model_dir = 'F:/Adcirc_SWAN/PARTneR2'
# generate_fort19_file(input_folder,pathres,tide_model_dir)
    f22 = pathres +'fort.22'
    tide_out = pathres +'fort.19'
    grid_in = input_folder + 'fort.14'
    
    pbnd_x,pbnd_y = get_x_and_y(grid_in)
    start_time,end_time = get_start_and_end_time_from_besttrack(f22)
    #start_time = start_time-dt.timedelta(days= 3) #- dt.timedelta(hours = 1)
    end_time = end_time + dt.timedelta(hours = 1)
    tide = get_tidal_elevation_matrix(tide_model_dir,start_time,end_time,pbnd_x,pbnd_y)
    create_tide_file(tide_out,start_time,end_time,tide,pbnd_x,pbnd_y)
    return()



def fill_tide_fort15(pathres,tide_model_dir,start_time,end_time):
# input_folder = 'F:/Adcirc_SWAN/PARTneR2/Test_Runs/input_files/'
# pathres='F:/Adcirc_SWAN/PARTneR2/Test_Runs/Test_06_all_py/'
# tide_model_dir = 'F:/Adcirc_SWAN/PARTneR2'
# generate_fort19_file(input_folder,pathres,tide_model_dir)
    # f22 = pathres +'fort.22'
    grid_in = pathres + 'fort.14'
    
    pbnd_x,pbnd_y = get_x_and_y(grid_in)
    # start_time,end_time = get_start_and_end_time_from_besttrack(f22)
    lista = get_periodic_tidal_cond(tide_model_dir,start_time,end_time,pbnd_x,pbnd_y)
    
    # lista = get_periodic_tidal_cond(tide_model_dir,start_time,end_time,139,-10)
    
    with open(pathres+"fort.15", 'r') as f:
        filedata = f.read()
    filedata = filedata.replace('TIDE_DATA', lista)
    
    with open(pathres+"fort.15", 'w') as fout:
            fout.write(filedata)

def change_dates_and_copy_f26(input_folder,pathres,start_time,end_time):
    with open(input_folder+"fort.26", 'r') as f:
                        filedata = f.read()
                        
    t1 = start_time.strftime('%Y%m%d.%H%M%S')
    t2 = end_time.strftime('%Y%m%d.%H%M%S')
    filedata = filedata.replace('%%dateinis%%', t1)
    filedata = filedata.replace('%%dateends%%', t2)
    
    with open(pathres+"fort.26", 'w') as fout:
        fout.write(filedata)
    return()

def change_dates_and_copy_f15(pathres,start_time,end_time):

    with open(pathres+"fort.15", 'r') as f:
                            filedata = f.read()
                        
    YYYY = start_time.strftime('%Y')
    MM = start_time.strftime('%m')
    DD = start_time.strftime('%d')
    HH = start_time.strftime('%H')
    delta = end_time - start_time
    DELTA_IN_DAYS = delta.days + delta.seconds/3600/24
    filedata = filedata.replace('YYYY', YYYY)
    filedata = filedata.replace('MM', MM)
    filedata = filedata.replace('DD', DD)
    filedata = filedata.replace('HH', HH)
    filedata = filedata.replace('DELTA_IN_DAYS', str(DELTA_IN_DAYS))
    
    with open(pathres+"fort.15", 'w') as fout:
            fout.write(filedata)
    return()

def change_dates_and_copy_f15_swan(input_folder,pathres,start_time,end_time):

    with open(input_folder+"fort.15_swan", 'r') as f:
                            filedata = f.read()
                        

    delta = end_time - start_time
    DELTA_IN_DAYS = delta.days + delta.seconds/3600/24
    
    YYYY = start_time.strftime('%Y')
    MM = start_time.strftime('%m')
    DD = start_time.strftime('%d')
    HH = start_time.strftime('%H')
    filedata = filedata.replace('YYYY', YYYY)
    filedata = filedata.replace('MM', MM)
    filedata = filedata.replace('DD', DD)
    filedata = filedata.replace('HH', HH)
    filedata = filedata.replace('DELTA_IN_DAYS', str(DELTA_IN_DAYS))
    
    with open(pathres+"fort.15", 'w') as fout:
            fout.write(filedata)
    return()

def change_dates_and_copy_f15_ramp(input_folder,pathres,start_time,end_time):

    with open(input_folder+"fort.15_ramp", 'r') as f:
                            filedata = f.read()
                        
    delta = end_time - start_time
    DELTA_IN_DAYS = delta.days + delta.seconds/3600/24
    filedata = filedata.replace('DELTA_IN_DAYS', str(DELTA_IN_DAYS))
    with open(pathres+"fort.15", 'w') as fout:
            fout.write(filedata)
    return()

def copy_remaining_forcing_files_and_change_dates(input_folder,pathres, pathnewrun):
    #pathres='F:/Adcirc_SWAN/PARTneR2/Test_Runs/Test_06_all_py/'
    #input_folder = 'F:/Adcirc_SWAN/PARTneR2/Test_Runs/input_files/'
    # f22 = pathres +'fort.22'
    
    shutil.copyfile(input_folder+'fort.13', pathnewrun+'fort.13')
    shutil.copyfile(input_folder+'fort.14', pathnewrun+'fort.14')
    shutil.copyfile(input_folder+'swaninit', pathnewrun+'swaninit')
    try:
        shutil.copyfile(pathres+'/ramp/fort.67', pathnewrun+'fort.67')
    except:
        print('Hotstart fort.67 not found')
    # start_time,end_time = get_start_and_end_time_from_besttrack(f22)
    # change_dates_and_copy_f26(input_folder,pathres,start_time,end_time)
    # change_dates_and_copy_f15(input_folder,pathres,start_time,end_time)
    return()



def fill_nan(A):
    '''
    interpolate to fill nan values
    '''
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = interpolate.interp1d(inds[good], A[good],bounds_error=False)
    B = np.where(np.isfinite(A),A,f(inds))
    return B



def complete_storm_forecast(storm,forecasts):
    
        
    orig_times =   forecasts['fhr']  
        
    interp_fhr = np.arange(0, forecasts['fhr'][-1]+1, 6)
       
    
    f = interp.interp1d(orig_times, forecasts['vmax'], kind='linear',
                            fill_value='extrapolate')
    vmax = f(interp_fhr)
    
    lona=np.array(forecasts['lon'])
    lona[(lona<=0)] = lona[(lona<0)] + 360
    f = interp.interp1d(orig_times, lona, kind='linear',
                            fill_value='extrapolate')
    lon = f(interp_fhr)
    
    f = interp.interp1d(orig_times, forecasts['lat'], kind='linear',
                            fill_value='extrapolate')
    lat = f(interp_fhr)
    
    
    mslp =-(vmax/3.86)**(1/0.645)+1010
    
    lon_ = np.concatenate((storm.lon,np.round(lon[1:-1]*10)/10),axis = None)
    lat_ = np.concatenate((storm.lat,np.round(lat[1:-1]*10)/10),axis = None)
    vmax_ = np.concatenate((storm.vmax,np.round(vmax[1:-1])),axis = None)
    mslp_ = np.concatenate((storm.mslp,np.round(mslp[1:-1])),axis = None)
    
    
    
    data = [storm.time]
    dt = storm.time[-1]
    for i in range(1,len(interp_fhr)-1):
        delta = datetime.timedelta(hours = int(interp_fhr[i]))
        dtnew = dt + delta
        data = np.concatenate((data,dtnew), axis=None)
    
    class storm_():

       lon = lon_
       lat = lat_
       vmax = vmax_
       mslp = mslp_
       time = data
    
    
    return storm_


def trackinsidemesh(storm,pathres):
    
    grid_in = pathres + 'fort.14'
    pbnd_x,pbnd_y = get_x_and_y(grid_in)
    
    pol = np.vstack((pbnd_x,pbnd_y))
    polygon = geometry.Polygon(pol.T)
    
    track = np.vstack((storm.lon,storm.lat)).T
    
    points_in = []
    for p in range(0,len(track)):
        kk = geometry.Point(track[p,:])
        if kk.within(polygon) is True:
            points_in.append(p)
            
    if len(track)>points_in[-1]:
        points_in.append(points_in[-1]+1)
        
    if points_in[0] >= 1:
        points_in = np.insert(points_in, 0, points_in[0]-1, axis=0)
        
       
    storm.lon = storm.lon[points_in]
    storm.lat = storm.lat[points_in]
    storm.vmax = storm.vmax[points_in]
    storm.mslp = storm.mslp[points_in]
    storm.time = storm.time[points_in]
    
       
    return (storm)

def generate_fort22_from_tropycal(storm,pathres,rampdays):
    
    # storm = storm_
    # pathres = pathramp 
    # rampdays=4
    
    wind_out = pathres +'fort.22'
    make_run_folder(pathres)
    
    time = []
    for i in range(len(storm.time)):                                                       
       ts = storm.time[i]
       time.append(ts.strftime('%Y%m%d%H'))
       
    vmax = fill_nan(storm.vmax)
    mspl = fill_nan(storm.mslp)
       
    lon_ = np.array(storm.lon).copy()
    lat_ = np.array(storm.lat).copy()

    lon_ = get_lonstr_from_x(lon_)
    lat_ = get_latstr_from_y(lat_)

    r_max = get_rmax_with_knaff(vmax,storm.lat)
    Rmax = r_max

    ## add t    
    time0 = storm.time[0]-dt.timedelta(days=rampdays)
    timeramp_ = mdates.drange(time0,storm.time[0]-dt.timedelta(hours=6),dt.timedelta(hours=6))
    timeramp = []
    for i in range(len(timeramp_ )):                                                       
       ts = timeramp_[i]
       ts= pd.to_datetime(mdates.num2date(ts))
       timeramp.append(ts.strftime('%Y%m%d%H'))

    time = np.concatenate((np.array(timeramp),time), axis=0)

    vmax = np.concatenate((np.matlib.repmat(vmax[0],1,len(timeramp))[0,:],vmax), axis=0)
    mspl  = np.concatenate((np.matlib.repmat(mspl[0],1,len(timeramp))[0,:],mspl), axis=0)
    lon_  = np.concatenate((np.matlib.repmat(lon_[0],1,len(timeramp))[0,:],lon_), axis=0)
    lat_  = np.concatenate((np.matlib.repmat(lat_[0],1,len(timeramp))[0,:],lat_), axis=0)
    Rmax  = np.concatenate((np.matlib.repmat(Rmax[0],1,len(timeramp))[0,:],Rmax), axis=0)

    
    create_wind_file(wind_out,time,lat_,lon_,vmax,mspl,Rmax)
    return()



def read_RMSC_Fiji_TCtracks(trackfile):
    
    data = pd.read_csv(trackfile,skiprows=8,usecols=[0,1,2,4,5,6,10],names=['time','lat','lon','cat','pmin','poci','wind'])
    # 10-minute sustained wind speeds

    time = np.array(data['time'])    
    lat = np.array(data['lat'])
    lon = np.array(data['lon'])
    cat = np.array(data['cat'])
    pmin = np.array(data['pmin'])
    # pmin = fill_nan(pmin)
    wind = np.array(data['wind'])
    # wind = fill_nan(wind)
    
    class storm():
        pass

    storm.lon = lon
    storm.lat = lat
    storm.cat = cat
    storm.vmax = wind
    storm.mslp = pmin
    storm.time = np.array([datetime.datetime.strptime(n,"%Y-%m-%dT%H:%M:%S%z") for n in time])

    return (storm)

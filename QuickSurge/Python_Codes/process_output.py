# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 13:27:33 2023

@author: moritzw
"""
import os
import numpy as np
import pandas as pd
import netCDF4
import math
import matplotlib.dates as mdates
import datetime as dt
import scipy.interpolate as interp
from joblib import Parallel, delayed
from make_forcings import get_start_and_end_time_from_besttrack
import matplotlib.tri as tri
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
#plt.rcParams["figure.figsize"] = (30, 15)
import cartopy.crs as ccrs
from adcircpy import AdcircMesh
import datetime
# from matplotlib.colors import LinearSegmentedColormap
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
     new_cmap = colors.LinearSegmentedColormap.from_list(
           'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
     return new_cmap 


def create_times(start_time,nh):
    times = [str(start_time)]
    for i in range(nh):
        delta = datetime.timedelta(hours = i)
        dtnew = start_time + delta
        times.append(str(dtnew))
    return(times)
    

def make_results_folder(pathres):
    if not os.path.exists(pathres + '/results'):
           os.mkdir(pathres+ '/results')
           print("Directory " , pathres+ '/results' ,  " Created ")
    else:    
           print("Directory " , pathres+ '/results' ,  " already exists")
    return()

def load_zeta_and_hs_max(pathres):

    data = netCDF4.Dataset(pathres+ 'maxele.63.nc')
    zeta_max = np.array(data['zeta_max'])
    lon = np.array(data['x'])
    lat = np.array(data['y'])
    
    data = netCDF4.Dataset(pathres+'fort.63.nc')
    zeta = np.array(data['zeta'])
    

    try:
        data= netCDF4.Dataset(pathres+ 'swan_HS_max.63.nc')
        hs_max = np.array(data['swan_HS_max'])
    except:
        hs_max = []
        

    try:
        data= netCDF4.Dataset(pathres+'swan_HS.63.nc')
        hs = np.array(data['swan_HS'])
    except:
        hs = []
        

    try:
        data = netCDF4.Dataset(pathres+'swan_TPS.63.nc')
        tp = np.array(data['swan_TPS'])
    except:
        tp = []   
        
    
    try:
        data = netCDF4.Dataset(pathres+'swan_DIR.63.nc')
        dire = np.array(data['swan_DIR'])
    except:
        dire = []   
        
    return(lon,lat,zeta_max,hs_max,zeta,hs,tp,dire)

def get_target_lon_and_lat(output_locations_csv):    
    output_pts = pd.read_csv(output_locations_csv,delimiter=',')
    target_lat = output_pts['lat']
    target_lon = output_pts['lon']
    return(target_lon,target_lat)


def haversine_distance(lat1, lon1, lat2, lon2):
    p = 0.017453292519943295
    hav = 0.5 - math.cos((lat2-lat1)*p)/2 + math.cos(lat1*p)*math.cos(lat2*p) * (1-math.cos((lon2-lon1)*p)) / 2
    return(12742 * math.asin(math.sqrt(hav)))

def get_indices_for_closest_lat_lon(lat,lon,target_lat,target_lon):
    posix = np.zeros(np.size(target_lon))
    for jx in range(len(target_lon)):
        dis = np.zeros(np.size(lon))
        for ix in range(len(lon)):
            dis[ix] = haversine_distance(lat[ix],lon[ix],target_lat[jx],target_lon[jx])
        pos_ix = np.where(dis == np.min(dis))
        posix[jx] = pos_ix[0][0]
        print(jx)
    return(posix)

def get_ix_of_min_dist(lon,lat,trgt_ln,trgt_lt):
    dis = np.zeros(np.size(lon))
    for ix in range(len(lon)):
        dis[ix] = haversine_distance(lat[ix],lon[ix],trgt_lt,trgt_ln)
    pos_ix = np.where(dis == np.min(dis))
    posix = pos_ix[0][0]
    return(posix)

def get_zeta_and_hs_at_point_locations(zeta_max,hs_max,posix):
    zeta_output = np.zeros(np.size(posix))
    hs_output = np.zeros(np.size(posix))
    for i in range(len(posix)):
        zeta_output[i] = zeta_max[posix[i]]
        hs_output[i] = hs_max[posix[i]]
    return(zeta_output,hs_output)

def calc_max_TWL_nearshore_based_on_Merrifield(hs_ts,tp_ts,zeta_ts):
    tp_ts[tp_ts<=0] = 0
    hs_ts[hs_ts<=0] = 0
    zeta_ts[zeta_ts<=-5] = -5
    Hs = hs_ts
    Per = tp_ts
    MSLA = zeta_ts
    
    angles = 1# np.cos(np.deg2rad(Dir-theta_N))
    # ix = np.where(angles < 0)[0]
    # angles[ix] = 0
    Hb = np.zeros(np.shape(hs_ts))
    RU_2_perc = np.zeros(np.shape(hs_ts))
    for i in range(len(Hb)):
        Hb[i] = (Hs[i]**2.* Per[i] * (4*np.pi)**(-1)* angles* np.sqrt(1.0*9.81))**(2/5)
        RU_2_perc[i] = 0.33 * Hb[i] + (-0.1)
    TWL_guess = RU_2_perc + MSLA
    max_TWL_guess = np.nanmax(TWL_guess)
    return(max_TWL_guess)

def store_output_max_at_point_locations(pathres,output_locations_csv,outfilename):
# pathres='F:/Adcirc_SWAN/PARTneR2/Test_Runs/Harold_Test/'
# input_folder = 'F:/Adcirc_SWAN/PARTneR2/input_files/'
# output_locations_csv = input_folder + 'TO.csv'
    lon,lat,zeta_max,hs_max,zeta,hs,tp,dire = load_zeta_and_hs_max(pathres)
    target_lon,target_lat = get_target_lon_and_lat(output_locations_csv)
    #posix2 = get_indices_for_closest_lat_lon(lat,lon,target_lat,target_lon) # can do in serial but it's very slow.
    posix = Parallel(n_jobs=12)(delayed(get_ix_of_min_dist)(lon,lat,target_lon[i],target_lat[i]) for i in range(len(target_lon)))
    zeta_output,hs_output = get_zeta_and_hs_at_point_locations(zeta_max,hs_max,posix)
    max_TWL_guess = np.zeros(np.shape(posix))
    for i in range(len(posix)):
        hs_ts = hs[:,posix[i]]
        tp_ts = tp[:,posix[i]]
        zeta_ts = zeta[:,posix[i]]
        max_TWL_guess[i] = calc_max_TWL_nearshore_based_on_Merrifield(hs_ts,tp_ts,zeta_ts)
    
    # make_results_folder(pathres)
    df = pd.DataFrame({"lon": target_lon, "lat": target_lat, "zeta_max[m]": zeta_output, "hs_max[m]": hs_output, "max_TWL_nearshore[m]": max_TWL_guess})
    df.to_csv(pathres + '/results/'+outfilename,index=False)
    return()

def store_output_timser_at_point_locations(pathres,input_folder,output_locations_csv,outfilename):
# pathres='F:/Adcirc_SWAN/PARTneR2/Test_Runs/Harold_Test/'
# input_folder = 'F:/Adcirc_SWAN/PARTneR2/input_files/'
# output_locations_csv = input_folder + 'TO.csv'
    start_time,end_time = get_start_and_end_time_from_besttrack(pathres+'fort.22')
    lon,lat,zeta_max,hs_max,zeta,hs,tp,dire = load_zeta_and_hs_max(pathres)
    times = create_times(start_time,len(hs[:,0])-1)
    target_lon,target_lat = get_target_lon_and_lat(output_locations_csv)
    #posix2 = get_indices_for_closest_lat_lon(lat,lon,target_lat,target_lon) # can do in serial but it's very slow.
    posix = Parallel(n_jobs=12)(delayed(get_ix_of_min_dist)(lon,lat,target_lon[i],target_lat[i]) for i in range(len(target_lon)))
    # zeta_output,hs_output = get_zeta_and_hs_at_point_locations(zeta_max,hs_max,posix)
    # max_TWL_guess = np.zeros(np.shape(posix))
    for i in range(len(posix)):
        hs_ts = hs[:,posix[i]]
        hs_ts = np.where(hs_ts <0, 0, hs_ts)
        tp_ts = tp[:,posix[i]]
        tp_ts = np.where(tp_ts <0, 0, tp_ts)
        zeta_ts = zeta[:,posix[i]]
        dir_ts = dire[:,posix[i]]
        # max_TWL_guess = calc_max_TWL_nearshore_based_on_Merrifield(hs_ts,tp_ts,zeta_ts)
        
        make_results_folder(pathres)   
        df = pd.DataFrame({ "time": times, "hs[m]": hs_ts, "tp[s]": tp_ts, "dir[N=0,E=90]": dir_ts, "eta[m]": zeta_ts})
        df.to_csv(pathres + '/results/Pnt_' + str(i).zfill(4) + '.csv',index=False)
    return()


def plot_timser_at_point_locations(pathres,input_folder,target_locations,start_time,save):
   
    # start_time,end_time = get_start_and_end_time_from_besttrack(pathres+'fort.22')
    # start_time = start_time+dt.timedelta(days = 4)
    # end_time = end_time-dt.timedelta(days = 0.5)
    lon,lat,zeta_max,hs_max,zeta,hs,tp,dire = load_zeta_and_hs_max(pathres)
    times = create_times(start_time,len(zeta[:,0])-1)
    end_time = start_time+dt.timedelta(hours =len(times)) 

    # target_lon,target_lat = get_target_lon_and_lat(target_locations)    
    posix1 = get_indices_for_closest_lat_lon(lat,lon,target_locations[:,1],target_locations[:,0]) # can do in serial but it's very slow.
    # posix = Parallel(n_jobs=12)(delayed(get_ix_of_min_dist)(lon,lat,target_lon[i],target_lat[i]) for i in range(len(target_lon)))
   
    for i in range(len(posix1)):

        zeta_ts = zeta[:,int(posix1[i])]
       
        # np.datetime64("1990-01-01")
        # np.arange(times[0],times[-1],dtype="M8[h]")

        xori= mdates.drange(start_time,end_time,dt.timedelta(minutes=60))
        xnew= mdates.drange(start_time,end_time,dt.timedelta(minutes=10))
        # kk = xori_.strftime("%Y%m%d%H%M%S")
        
        # integer_list = [int(item) for item in kk]
        
        
        xori_ = pd.to_datetime(mdates.num2date(xori))
        xnew_ = pd.to_datetime(mdates.num2date(xnew))
        
      
        f = interp.interp1d(xori ,zeta_ts,  kind='quadratic',fill_value='extrapolate')
       
        twl_min=f(xnew)
        
        plt.figure(figsize=(25,10)) 
      
        plt.plot(xnew_,twl_min,'-k',linewidth=1,label='Storm Tide')
        plt.xlabel('UTC Time',fontsize=15)
        plt.ylabel('Water level [m]',fontsize=15)
        plt.yticks(fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        ylbl = plt.yticks()[0]
        for h in range(0,len(  xori_ ),24):
           plt.axvline(  xori_ [h].replace(hour=0),color='gray',lw=0.5,ls='dashed',dashes=(11,5))
        for h in range(0,len(ylbl)):
            plt.axhline(ylbl[h],color='gray',lw=0.5,ls='dashed',dashes=(11,5))
        plt.legend(loc='upper center', ncol = 3,fontsize=20)
        fileprint = pathres + '/results/Pnt_' + str(i).zfill(4) 
        
        if save == True:
            plt.savefig(fileprint)
            plt.close()
        else:
            print('Figure is not saved')
            
            
    return()



def plot_unswan(lon,lat,z,title,cbarlabel,levels,inter,cmap,pathres):
    f14 = pathres + 'fort.14'
    # open mesh file
    mesh = AdcircMesh.open(f14, crs=4326)
    trx = np.array(mesh.triangles)
    trx = trx-1
    triang = tri.Triangulation(lon, lat, trx)

    fig = plt.figure(figsize=(12, 8), dpi=300)
    img_extent = (np.min(lon),np.max(lon), np.min(lat), np.max(lat))
    #img_extent = (176, 181.5, -19.5, -15.5)
    proj = ccrs.PlateCarree(central_longitude=180)
    ax = plt.axes(projection=proj)

    ax.set_extent(img_extent, crs=ccrs.PlateCarree())

    plt.title(title)
    ax.use_sticky_edges = False
    # set a margin around the data
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)

    # -----------------------------------------------------------------------------
    # Refine data
    # -----------------------------------------------------------------------------
    refiner = tri.UniformTriRefiner(triang)
    tri_refi, z_test_refi = refiner.refine_field(z, subdiv=2)
    # levels = np.arange(0., 1., 0.025)

    im = ax.tricontourf(tri_refi, z_test_refi, levels=levels,
                        cmap=cmap, transform=ccrs.PlateCarree())
    # im = ax.tricontour(tri_refi, z_test_refi, levels=levels,
    #                colors=['0.25', '0.5', '0.5', '0.5', '0.5'],
    #                linewidths=np.arange(0.,max(levels),inter), transform=ccrs.PlateCarree())
        
    
    #im = ax.tricontourf(triang, z, cmap='turbo', shading='gouraud', transform=ccrs.PlateCarree())
    #im = ax.tripcolor(triang, z, cmap='turbo', transform=ccrs.PlateCarree())

    # ax.coastlines(resolution='50m', color='black', linewidth=1)
    ax.gridlines(crs=ccrs.PlateCarree(central_longitude=180), draw_labels=True,
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    cb = fig.colorbar(im,extend="both")
    cb.set_label(cbarlabel)
    return(fig,ax)

def plot_merrifield(pathres,outfilename,out_fig_name):
    output_pts = pd.read_csv( pathres + '/results/'+outfilename,delimiter=',')
    lat = output_pts['lat']
    lon = output_pts['lon']
    max_TWL_nearshore = output_pts['max_TWL_nearshore[m]']
    fig = plt.figure(figsize=(12, 8))
    img_extent = (np.min(lon-1),np.max(lon+1),np.min(lat-1),np.max(lat+1))
    
    proj = ccrs.PlateCarree(central_longitude=180)
    ax = plt.axes(projection=proj)
    
    ax.set_extent(img_extent,crs=ccrs.PlateCarree())
    ax.use_sticky_edges = False    
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)

    ax.coastlines(resolution='10m', color='black', linewidth=1)
    ax.gridlines(crs=ccrs.PlateCarree(central_longitude=180), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    im = ax.scatter(lon,lat,c=max_TWL_nearshore,s = 5*lon/lon,cmap='turbo', transform=ccrs.PlateCarree())
    cb = fig.colorbar(im, extend="both")
    cb.set_label('TWL nearshore (m) based on an empirical equation by Merrifield et al. (2014)')
    plt.savefig(pathres+ '/results/'+out_fig_name)
    plt.close(fig)
    return()

def plot_and_save_figures(pathres):
    if not os.path.exists(pathres + '/results'):
           os.mkdir(pathres+ '/results')
           print("Directory " , pathres+ '/results' ,  " Created ")
    else:    
           print("Directory " , pathres+ '/results' ,  " already exists")
    lon,lat,zeta_max,hs_max,zeta,hs,tp,dire = load_zeta_and_hs_max(pathres)
    levels = np.linspace(0, max(zeta_max), 25)
    #levels = np.linspace(min(zeta_max), 2, 91)
    #levels = np.linspace(0,4, 91)
    cmap = cm.get_cmap(name='gist_ncar', lut=None)
    new_cmap = truncate_colormap(cmap, 0.1, 1)
    inter = 0.25
    fig,ax = plot_unswan(lon,lat,zeta_max,'zeta_max','zeta_max',levels,inter,new_cmap ,pathres)
    plt.savefig(pathres+'results/zeta_max.png')
    # plt.close(fig)
    hs = hs_max
    hs = np.where(hs <0, 0, hs)
    levels = np.linspace(min(hs), max(hs), 25)
    #levels = np.linspace(0,8, 91)
    inter = 1
    fig,ax = plot_unswan(lon,lat,hs,'hs_max','hs_max',levels,inter,new_cmap ,pathres)
    plt.savefig(pathres+'results/hs_max.png')
    # plt.close(fig)
    return()


def load_zeta(pathres):
    nc_zeta_ts = netCDF4.Dataset(pathres+'fort.63.nc')
    lon = np.array(nc_zeta_ts['x'])
    lat = np.array(nc_zeta_ts['y'])
    zeta = np.array(nc_zeta_ts['zeta'])
    return(lon,lat,zeta)



def load_hs(pathres):
    nc_hs_ts = netCDF4.Dataset(pathres+'swan_HS.63.nc')
    lon = np.array(nc_hs_ts['x'])
    lat = np.array(nc_hs_ts['y'])
    hs = np.array(nc_hs_ts['swan_HS'])
    return(lon,lat,hs)

def plot_last_time_eta(pathres):
    lon,lat,zeta = load_zeta(pathres)
    levels = np.linspace(np.nanmin(zeta), np.nanmax(zeta), 91)
    #levels = np.linspace(0,4, 91)
    fig,ax = plot_unswan(lon,lat,zeta[-1,:],'zeta','zeta',levels,"gist_ncar",pathres)
    return()

def plot_last_time_hs(pathres):
    lon,lat,hs = load_hs(pathres)
    hs = hs[-1,:]
    hs = np.where(hs<0, 0, hs)
    levels = np.linspace(np.nanmin(hs), np.nanmax(hs), 91)
    #levels = np.linspace(0,4, 91)
    fig,ax = plot_unswan(lon,lat,hs,'hs','hs',levels,"gist_ncar",pathres)
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
    
    opbnd_x = np.squeeze(px[pbnd_val_])
    opbnd_y = np.squeeze(py[pbnd_val_])
    
    pbnd = pmesh.boundaries.interior
    pbnd_val = pbnd.indexes[:]
    pbnd_val_ = []
    
    pbnd_x = []
    pbnd_y = []
    for sublist in pbnd_val :     
        pbnd_x.append(np.append(px[sublist],np.nan))
        pbnd_y.append(np.append(py[sublist],np.nan))
        
    lpbnd_x = [item for sublist in pbnd_x for item in sublist]
    lpbnd_y = [item for sublist in pbnd_y for item in sublist]

    return(opbnd_x,opbnd_y,lpbnd_x,lpbnd_y)

def plot_track_domain(input_folder,pathres,storm,stormname,stormyear,newrun):
    
    grid_in = input_folder + 'fort.14'
    opbnd_x,opbnd_y,lpbnd_x,lpbnd_y = get_x_and_y(grid_in)

    # Define colors
    colors = ['lightblue', 'yellow', 'orange','red',(0.8, 0.6, 0.9), 'purple']    
    # Define corresponding values
    # values = np.linspace(0, 5, 6)

    # fig = plt.figure(figsize=(12, 8))
    
    img_extent = (np.min(opbnd_x-2),np.max(opbnd_x+2),np.min(opbnd_y-2),np.max(opbnd_y+2))
    
    proj = ccrs.PlateCarree(central_longitude=180)
    ax = plt.axes(projection=proj)
    
    ax.plot(opbnd_x,opbnd_y,color='k', transform=ccrs.PlateCarree())
    ax.plot(lpbnd_x,lpbnd_y,color='#00FF00', transform=ccrs.PlateCarree())
    
    ax.plot(storm.lon,storm.lat,color='k', transform=ccrs.PlateCarree(),zorder=1)
    
    for i in range(0,6):       
        indices = np.where(storm.cat == i)[0]
        if i == 0:
            ax.scatter(storm.lon[indices],storm.lat[indices], facecolors=colors[i], s=100, 
                   marker = "o", alpha = 1,edgecolor = 'k',transform=ccrs.PlateCarree(), label='Tropical Storm',zorder=2)
        else:
            ax.scatter(storm.lon[indices],storm.lat[indices], facecolors=colors[i], s=100,  
                       marker = "o", alpha = 1,edgecolor = 'k',transform=ccrs.PlateCarree(), label='Category {}'.format(i),zorder=2)
    
    
    ax.set_extent(img_extent,crs=ccrs.PlateCarree())
    ax.use_sticky_edges = False    
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)

    ax.legend()

    ax.gridlines(crs=ccrs.PlateCarree(central_longitude=180), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')

    ax.set_title('Tropical Cyclone '+ stormname + ' ' + str(stormyear))
    plt.savefig(pathres+ '/'+ newrun)
    # plt.close(fig)
    return()

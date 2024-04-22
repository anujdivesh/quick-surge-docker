#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 10:56:58 2023

@author: judith
"""


import numpy as np
import pandas as pd
import os
import pyTMD
import datetime as dt
import matplotlib.dates as mdates
import scipy.io
from datetime import date
from scipy.signal import savgol_filter
import helper as gpath

#Load directory paths
paths = gpath.load_paths()

mat = scipy.io.loadmat(paths.baseDir+"Tides/tide_fac_constants.mat")
shallowiname = mat['shallow']['iname'][0,0].copy()
shallowcoef = mat['shallow']['coef'][0,0].copy()
ishallow= mat['const']['ishallow'][0,0].copy()
nshallow= mat['const']['nshallow'][0,0].copy()
   
nshallow = np.ma.masked_invalid(nshallow).astype(int)
ishallow = np.ma.masked_invalid(ishallow).astype(int) - 1
not_shallow = ishallow.mask  # True where it was masked.
nshallow = nshallow.compressed()
ishallow = ishallow.compressed()
kshallow = np.nonzero(~not_shallow)[0]
# (comments based on t_tide)
# Coefficients of the formulas in the Explan. Suppl.
_sc = np.array([270.434164, 13.1763965268, -0.0000850, 0.000000039])
_hc = np.array([279.696678, 0.9856473354, 0.00002267, 0.000000000])
_pc = np.array([334.329556, 0.1114040803, -0.0007739, -0.00000026])
_npc = np.array([-259.183275, 0.0529539222, -0.0001557, -0.000000050])
# First coeff was 281.220833 in Foreman but Expl. Suppl. has 44.
_ppc = np.array([281.220844, 0.0000470684, 0.0000339, 0.000000070])
_coefs = np.vstack((_sc, _hc, _pc, _npc, _ppc))


def ut_astron(jd):
    """
    Compute the astronomical variables and their time derivatives.

    Parameters
    ----------
    jd : float, scalar or sequence
        Time (UTC) in days starting with 1 on 1 Jan. of the year 1
        in the proleptic Gregorian calendar as in
        `datetime.date.toordinal`.

    Returns
    -------
    astro : array, (6, nt)
        rows are tau, s, h, p, np, pp (cycles)
    ader : array, (6, nt)
        time derivatives of the above (cycles/day)

    Notes
    -----
    2-D arrays are always returned.

    Variables are:

    ===  ====================================================
    tau  lunar time
    s    mean longitude of the moon
    h    mean longitude of the sun
    p    mean longitude of the lunar perigee
    np   negative of the longitude of the mean ascending node
    pp   mean longitude of the perihelion (solar perigee)
    ===  ====================================================


    Based on UTide v1p0 9/2011 d.codiga@gso.uri.edu, which
    in turn came from t_tide's t_astron.m, Pawlowicz et al 2002

    For more background information from t_tide, see the t_tide_doc
    string variable in this module.

    """


    jd = np.atleast_1d(jd).flatten()

    # Shift epoch to 1899-12-31 at noon:
    # daten = 693961.500000000  Matlab datenum version

    daten = 693595.5  # Python epoch is 366 days later than Matlab's

    d = jd - daten
    D = d / 10000

    args = np.vstack((np.ones(jd.shape), d, D * D, D**3))

    astro = np.fmod((np.dot(_coefs, args) / 360), 1)

    # lunar time: fractional part of solar day
    #             plus hour angle to longitude of sun
    #             minus longitude of moon
    tau = jd % 1 + astro[1, :] - astro[0, :]
    astro = np.vstack((tau, astro))

    # derivatives (polynomial)
    dargs = np.vstack(
        (np.zeros(jd.shape), np.ones(jd.shape), 2.0e-4 * D, 3.0e-4 * D * D),
    )

    ader = np.dot(_coefs, dargs) / 360.0
    dtau = 1.0 + ader[1, :] - ader[0, :]
    ader = np.vstack((dtau, ader))

    return astro, ader



def FUV(t, tref, lind, lat, ngflgs):

    """
    UT_FUV()
    compute nodal/satellite correction factors and astronomical argument
    inputs
      t = times [datenum UTC] (nt x 1)
      tref = reference time [datenum UTC] (1 x 1)
      lind = list indices of constituents in ut_constants.mat (nc x 1)
      lat = latitude [deg N] (1 x 1)
      ngflgs = [NodsatLint NodsatNone GwchLint GwchNone] each 0/1
    output
      F = real nodsat correction to amplitude [unitless] (nt x nc)
      U = nodsat correction to phase [cycles] (nt x nc)
      V = astronomical argument [cycles] (nt x nc)
    UTide v1p0 9/2011 d.codiga@gso.uri.edu
    (uses parts of t_vuf.m from t_tide, Pawlowicz et al 2002)
    """

    t = np.atleast_1d(t).flatten()
    nt = len(t)
    nc = len(lind)



    if ngflgs[1]:
        F = np.ones((nt, nc))
        U = np.zeros((nt, nc))
    else:
        if ngflgs[0]:
            tt = np.array([tref])
        else:
            tt = t
        ntt = len(tt)

        astro, ader = ut_astron(tt)

        if abs(lat) < 5:
            lat = np.sign(lat) * 5

        slat = np.sin(np.deg2rad(lat))

        
        rr = mat['sat']['amprat'][0,0].copy()


        j = mat['sat']['ilatfac'][0,0] == 1
        rr[j] *= 0.36309 * (1.0 - 5.0 * slat**2) / slat

        

        j = mat['sat']['ilatfac'][0,0] == 2
        rr[j] *= 2.59808 * slat


        # sat.deldood is (162, 3); all other sat vars are (162,)
        uu = np.dot(mat['sat']['deldood'][0,0],astro[3:6, :]) + mat['sat']['phcorr'][0,0 ][:, ]
        
        
        np.fmod(uu, 1, out=uu)  # fmod is matlab rem; differs from % op
        matt = rr[:, ] * np.exp(1j * 2 * np.pi * uu)
       

        nfreq = len(mat['const']['isat'][0,0])  # 162
        F = np.ones((nfreq, ntt), dtype=complex)

        iconst = mat['sat']['iconst'][0,0][:,0] - 1
        ind = np.unique(iconst)
        
        # ii =7
        # k = 1 + np.sum(matt[iconst== ii], axis=0) 
        # print(str(kk[0,0]))
        
        for ii in ind:
            F[ii, :] = 1 + np.sum(matt[iconst== ii], axis=0)

        U = np.angle(F) / (2 * np.pi)  # cycles
        F = np.abs(F)
        

        
        for i0, nshal, k in zip(ishallow, nshallow, kshallow):
            ik = i0 + np.arange(nshal)
            j = shallowiname[ik] - 1
            exp1 = shallowcoef[ik, None]
            exp2 = np.abs(exp1)
            F[k, :] = np.prod(F[j, :] ** exp2, axis=0)
            U[k, :] = np.sum(U[j, :] * exp1, axis=0)

        F = F[lind, :].T
        U = U[lind, :].T

    # gwch (astron arg)
    if ngflgs[3]:  # None (raw phase lags not Greenwich phase lags).
        freq = linearized_freqs(tref)
        V = 24 * (t[:, np.newaxis] - tref) * freq[lind]

    else:
        if ngflgs[2]:  # Linearized times.
            tt = np.array([tref])
        else:
            tt = t  # Exact times.
        ntt = len(tt)

        astro, ader = ut_astron(tt)

        V = np.dot(mat['const']['doodson'][0,0], astro) + mat['const']['semi'][0,0][:, ]
        # V has nan values from both const.* arrays
        with np.errstate(invalid="ignore"):
            np.fmod(V, 1, out=V)

 
        for i0, nshal, k in zip(ishallow, nshallow, kshallow):
            ik = i0 + np.arange(nshal)
            j = shallowiname[ik] - 1
            exp1 = shallowcoef[ik, None]
            V[k, :] = np.sum(V[j, :] * exp1, axis=0)

        V = V[lind, :].T

        if ngflgs[2]:  # linearized times
            freq = linearized_freqs(tref)
            V = V + 24 * (t[:, None] - tref) * freq[None, lind]

    return F, U, V



def linearized_freqs(tref):
    astro, ader = ut_astron(tref)
    freq = mat['const']['freq'][0,0]
    selected = np.dot(mat['const']['doodson'][0,0][not_shallow, :], ader) / 24
    freq[not_shallow] = selected.squeeze()
    for i0, nshal, k in zip(ishallow, nshallow, kshallow):
        ik = i0 + np.arange(nshal)
        freq[k] = (freq[shallowiname[ik] - 1] * shallowcoef[ik]).sum()
    return freq

def get_periodic_tidal_cond(tide_dir,start_time,end_time,pbnd_x,pbnd_y):
    
  
    grid_file = os.path.join(tide_dir,'TPXO9','DATA','grid_tpxo9v2')## check if we can go for a newer version, others models can be used:'GOT4.8','FES2014',...
    model_file = os.path.join(tide_dir,'TPXO9','DATA','h_tpxo9.v2')
    
    
    model_format = 'OTIS'
    EPSG = '4326'
    TYPE = 'z'
    
    time_tide = mdates.drange(start_time,end_time,dt.timedelta(minutes=10))
    time_tmd=time_tide-mdates.date2num(np.datetime64('1992-01-01T00:00:00')) 
    
 
    amp,ph,D,constants = pyTMD.io.extract_constants(pbnd_x, pbnd_y,
                               grid_file,model_file,EPSG,TYPE=TYPE,METHOD='spline',GRID=model_format,extrapolate = True,cutoff=50)

    # Convert real and imaginary parts to amplitude 
    amp_r = np.abs(amp)
    phase_r = ph

    names, amp, red, frq, F, phs, indices = get_tidal_pot(start_time,end_time,pbnd_x,pbnd_y,constants)
    

    amp_r = amp_r [:,indices]
    amp_rr = np.vstack((amp_r,amp_r,amp_r))
    smoothed = savgol_filter(amp_rr.T, window_length=5, polyorder=2).T
    
    amp_r =  smoothed [len(amp_r[:,1]):2*len(amp_r[:,1]),:]
    
    phase_r = phase_r[:,indices]
    
    
    valid =np.where(np.isnan(amp)==False)
    valid = valid[0]
    
    lista = []

    
    list_harm  = ['Q1  ','O1  ','P1  ','K1  ','N2  ','M2  ','S2  ','K2  ']
    
    
    lista.append(str(len( list_harm ))+ '             !NTIF\n')
    for ip in range(len(valid)):
        if names.values[ip] in list_harm:       
            lista.append(names.values[ip]+'\n')
            lista.append(str(amp[ip,0]) +' '+ str(frq[ip,0])+ ' ' +str(red[ip,0]) +' ' +str(F[ip]) + ' ' +str(phs[ip]) + '   !  TPK, AMIGT, ETRF, FFT, FACET\n')
 
    lista.append(str(len(list_harm ))+ '             !NBFR\n') 
    for ip in range(len(valid)):
        if names.values[ip] in list_harm: 
            lista.append(names.values[ip]+ '\n')
            lista.append( str(frq[ip,0])+ ' ' +str(F[ip]) + ' ' +str(phs[ip]) +  '   \n')
    
    for ip in range(len(valid)):
        if names.values[ip] in list_harm: 
            lista.append(names.values[ip]+ '\n')
            for ib in range(len(amp_r)):
                # if ip==len(valid)-1 and ib==len(amp_r)-1:
                if names.values[ip]==list_harm[-1] and ib==len(amp_r)-1:
                    lista.append( str(amp_r[ib,ip])+ ' ' + str(phase_r[ib,ip]))
                else :
                    lista.append( str(amp_r[ib,ip])+ ' ' + str(phase_r[ib,ip])+ '\n')
            
    lista_ = str()     
    lista_ = lista_.join([str(a) for a in lista]) 
   
    return lista_

def get_tidal_pot(start_time,end_time,pbnd_x,pbnd_y,constants):
    
    # mat = scipy.io.loadmat('/media/judith/10TB/QuickSurge/Tides/tide_fac_constants.mat')
    cnames = pd.Series(mat['const']['name'][0,0])
      
    t0 = date.toordinal(start_time)#+1+365
    tend = date.toordinal(end_time)#+1+365 
    tref = np.mean([t0,tend]) 
    t = np.linspace(t0,tend,num = 1000)
    lat = np.mean(pbnd_y)

    ind= np.ones(len(constants),dtype =int)
    for c in range(len(constants)):
        
        cn = constants[c]
        i = np.where(cnames.str.match(cn.upper()))
        ind[c] = i[0]
     
    indices = np.argsort(ind)
    ind = np.sort(ind)
    
    freq_ = mat['const']['freq'][0,0]
    frq = freq_[ind]
    names = cnames[ind]
    
    # Doodson amp multiplied by Doodson constant/g:  
    # (See : https://www.whoi.edu/fileserver.do?id=21351&pt=10&p=17272)
    doodsonamp = mat['const']['doodsonamp'][0,0]
    earthreduc = mat['const']['earthreduc'][0,0]
    amp  = np.abs(doodsonamp[ind])*0.2675
    red  = earthreduc[ind] 
    
    # Now get the F, U, and V (using the U_tide function ut_FUV)
    
    F,U,V = FUV(t, tref, ind, lat,np.zeros([4,1]))


    #Final factors
    F = np.mean(F,axis = 0) #This is the average nodal factor
    phs = (np.mean(U,axis = 0)+V[0,:])*360  # This is the average nodal correction + astronomical argument at beginning ofsimulation
                               
    # Make sure phase between 0 and 360 for aesthetic purposes
    while any(phs < 0):
        phs[phs<0] = phs[phs<0] + 360;
    
    frq = frq *2*np.pi/3600
    
    
 
    
    # dictionary = {'name':names, 'amp':amp[:,0],'freq':freq[:,0],'red':red[:,0],'F':F,'phs':phs}
    
    # import pandas as pd 
    # company_df = pd.DataFrame(dictionary )
    
    # list = [names.T,amp.T]
    
    return names, amp, red, frq, F, phs, indices
    
    # obj.f15.ntif = length(cnstit.NR.name);
    # for k = 1:obj.f15.ntif
    #     obj.f15.tipotag(k).name   = cnstit.NR.name{k}; % name in frq order
    #     obj.f15.tipotag(k).val(1) = cnstit.NR.amp(k); % potential amplitude of species
    #     obj.f15.tipotag(k).val(2) = cnstit.NR.frq(k)*2*pi/3600; % frq in rad/s format
    #     obj.f15.tipotag(k).val(3) = cnstit.NR.red(k); % earth rigidity reduction factor
    #     obj.f15.tipotag(k).val(4) = F(k);   % average nodal factor
    #     obj.f15.tipotag(k).val(5) = phs(k); % average nodal correction + 
    #                          % astronomical argument at beginning of simulation
    # end
    
            
            
        
    
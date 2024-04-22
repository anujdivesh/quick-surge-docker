# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 15:32:25 2023

@author: moritzw
"""
import shutil
import glob
import os
import helper as gpath

#Load directory paths
paths = gpath.load_paths()

def run_model(nproc,pathres):
    # os.chdir(pathres)
    os.system('cd %s;  %s/adcprep --np %s --partmesh' %(pathres,paths.baseDir+"Model",str(nproc)))
    os.system('cd %s;  %s/adcprep  --np %s --prepall' %(pathres,paths.baseDir+"Model",str(nproc)))
    os.system('cd %s; mpirun -np %s %s/padcswan' %(pathres,str(nproc),paths.baseDir+"Model"))
    
    try: 
        os.system('taskkill /F /IM padcswan_SH.exe')
    except:
        print('ADCIRC  +  SWAN run finished succesfully')
    
    dir_list = glob.iglob(os.path.join(pathres, "PE0*"))
    for path in dir_list:
        if os.path.isdir(path):
            shutil.rmtree(path)
            
            
def run_ramp_model(nproc,pathres):
    # os.chdir(pathres)
    os.system('cd %s;  %s/adcprep --np %s --partmesh' %(pathres,paths.baseDir+"Model",str(nproc)))
    os.system('cd %s;  %s/adcprep  --np %s --prepall' %(pathres,paths.baseDir+"Model",str(nproc)))
    os.system('cd %s; mpirun -np %s %s/padcirc' %(pathres,str(nproc),paths.baseDir+"Model"))
    
    try: 
        os.system('taskkill /F /IM padcswan_SH.exe')
    except:
        print('ADCIRC  +  SWAN run finished succesfully')
        
    os.system('cp %sfort.67 %s' %(os.path.join(pathres + 'PE0000/'),pathres))
    
    dir_list = glob.iglob(os.path.join(pathres, "PE0*"))
    for path in dir_list:
        if os.path.isdir(path):
            shutil.rmtree(path)
            
            

def run_model_hotstart(nproc,pathres):
    # os.chdir(pathres)
    os.system('cd %s;  %s/adcprep --np %s --partmesh' %(pathres,paths.baseDir+"Model",str(nproc)))
    os.system('cd %s;  %s/adcprep  --np %s --prepall' %(pathres,paths.baseDir+"Model",str(nproc)))
    os.system('cd %s; %s/adcprep --np %s  --hotlocalize 67'%(pathres,paths.baseDir+"Model",str(nproc)))
    os.system('cd %s; %s/adcprep --np %s  --hotglobalize 67'%(pathres,paths.baseDir+"Model",str(nproc)))
    os.system('cd %s; mpirun -np %s %s/padcswan' %(pathres,str(nproc),paths.baseDir+"Model"))
    
    try: 
        os.system('taskkill /F /IM padcswan_SH.exe')
    except:
        print('ADCIRC  +  SWAN run finished succesfully')
    
    dir_list = glob.iglob(os.path.join(pathres, "PE0*"))
    for path in dir_list:
        if os.path.isdir(path):
            shutil.rmtree(path)
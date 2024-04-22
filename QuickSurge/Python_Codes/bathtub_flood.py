# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 10:26:36 2024

@author: antonioh
"""
import matplotlib.tri as mtri
import pandas as pd
import numpy as np
import rasterio
from rasterio.mask import mask
import geopandas as gpd
from shapely.geometry import mapping


def bathtub_flood(input_geotiff,subdomains,twl_points,aux_points,outputdir):
    
    print('Generating inundation extends')

    # input_geotiff = 'D:\\Projects_SPC\\PARTneR2\\Data\\DEM_VU\\merged_WGS84.tif'
    # subdomains = 'D:\\Projects_SPC\\PARTneR2\\GIS\\bathtub_domains_VU.shp'
    
    # twl_points = 'D:/Projects_SPC/PARTneR2/GIS/TWL_results.csv'
    # aux_points = 'D:\\Projects_SPC\\PARTneR2\\Data\\DEM_VU\\inner_vertices_VU.shp'
    
    # Load TWL results
    output_pts = pd.read_csv(twl_points,delimiter=',')
    lat = output_pts['lat']
    lon = output_pts['lon']
    max_TWL_nearshore = output_pts['max_TWL_nearshore[m]']
    
    # Load inner poits
    auxpoints = gpd.read_file(aux_points)
    
    
    # Convert the DataFrame to a GeoDataFrame
    points_gdf = gpd.GeoDataFrame(
        output_pts,
        geometry=gpd.points_from_xy(lon,lat,max_TWL_nearshore),
        crs='EPSG:4326'  # Assuming the coordinates are in WGS84 (longitude, latitude)
    )
    
    
    # Load the GeoTIFF file
    # with rasterio.open(input_geotiff) as src:
    # Load the shapefile
    src =  rasterio.open(input_geotiff)
    # Load subdomains shapefile
    shapefile = gpd.read_file(subdomains)
    # Reproject the shapefile to match the coordinate reference system (CRS) of the raster
    shapefile = shapefile.to_crs(src.crs)
    
    # loop over all the national subdomains
    for i in range(0,len(shapefile.geometry)):
        # Extract the geometry in GeoJSON format
        subd = mapping(shapefile.geometry[i])
        # Clip the raster with the polygon
        clipped_raster, transform = mask(src, [subd], crop=True)
        clipped_raster = clipped_raster.astype(np.float16)
        # mask SEM bellow 0 m
        clipped_raster[clipped_raster<=0] = np.nan
        
    
        # Get the grid coordinates (pixel coordinates)
        height = clipped_raster.shape[1]
        width = clipped_raster.shape[2]
        cols, rows = np.meshgrid(np.arange(width), np.arange(height))
        xs, ys = rasterio.transform.xy(transform, rows, cols)
        lons= np.array(xs)
        lats = np.array(ys)
        
        # Check TWL poimts inside the subdomain polygon
        twl_within = points_gdf[points_gdf.geometry.within(shapefile.geometry[i], align = True)]
        twl_within = np.vstack((np.array(twl_within['lon']),np.array(twl_within['lat']),np.array(twl_within['max_TWL_nearshore[m]']))).T
        
        # Fill the ausxiliar points for interpolation with the minimum TWL around the island
        aux_within = auxpoints[auxpoints.geometry.within(shapefile.geometry[i], align = True)] 
        aux_within = aux_within.geometry.apply(lambda geom: np.array(geom.coords))  
        aux_within = np.array([ item[1][0] for item in aux_within.items()])
        
        # thera are some ilands where the auxiliary points are not needed
        if len(aux_within) !=  0:
    
            x =  np.append(twl_within[:,0],aux_within[:,0],axis = 0 )
            y =  np.append(twl_within[:,1],aux_within[:,1],axis = 0 )
            z =  np.append(twl_within[:,2],aux_within[:,1]*0 + twl_within[:,2].min(),axis = 0 )
        
        else:
            
            x =  twl_within[:,0]
            y =  twl_within[:,1]
            z =  twl_within[:,2]
        
    
        # Create a Triangulation object
        triang = mtri.Triangulation(x, y)
        # Create a LinearTriInterpolator object
        interp = mtri.LinearTriInterpolator(triang, z)
        # Interpolate z values 
        zi = interp(np.squeeze(lons),np.squeeze(lats))
        
    
        # inundation depths
        inundation = clipped_raster[0,:,:] - zi
        inundation[inundation>=0] = np.nan
        inundation = -inundation
        clipped_raster[0,:,:] = inundation
        
        # Get the metadata of the clipped raster
        meta = src.meta.copy()
        meta.update({"driver": "GTiff",
                     "height": clipped_raster.shape[1],
                     "width": clipped_raster.shape[2],
                     "transform": transform})
    
        # Write the clipped raster to a new GeoTIFF file
        output_geotiff = outputdir + 'clipped_output_'+ str(i).zfill(2) + '.tif'
        with rasterio.open(output_geotiff, "w", **meta,compress='deflate') as dest:
            dest.write(clipped_raster)
    
    
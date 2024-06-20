#!/usr/bin/env python
# coding: utf-8

# NDVI anomalies

import geopandas as gpd
import shapely.geometry
from pystac_client import Client
from odc.stac import configure_rio, stac_load
import pathlib
import json
import rasterio
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
import rioxarray
import os
from osgeo import gdal
import shutil
import datetime
from dateutil.relativedelta import relativedelta


# Now we need to configure and start a session in the HDC, to do this, we need to log-in using a token


root = pathlib.Path().resolve().parent.parent
root

TOKEN_PATH = "C:/Users/oscar.bautista/OneDrive - World Food Programme/Scripts/tk.json"
HDC_STAC_URL= "https://api.earthobservation.vam.wfp.org/stac/"

import os
def _get_hdc_stac_param_from_env():
    
    if "JUPYTERHUB_USER" in os.environ:
    
        signer = None
        header = None
        aws={}   # Get credentials for accessing S3 bucket 

    else:

        def make_signer(fname="./tk.json"):
            """
            Loads token from file at fname, and returns a function patching request urls with said token
            """
            tk = ""
            with open(fname, "rt") as src:
                tk = json.load(src)["tk"]

            def sign(url, _tk=tk):
                signed = f"{url}?{_tk}"
                return signed

            return sign

        signer = make_signer(TOKEN_PATH)
        header = {"origin": "https://wfp.org"}
        aws = None
        
    # Instantiate an API client pointing to the WFP HDC STAC API
    hdc_stac_client = Client.open(HDC_STAC_URL, headers=header)
    
    # Set up GDAL/rasterio configuration.
    configure_rio(cloud_defaults=True, verbose=True, aws=aws)
        
    return hdc_stac_client, signer


#
# STAC CLIENTS
#

hdc_stac_client, signer = _get_hdc_stac_param_from_env() 


# Ask the user to select an option
print('Please select the pilot area:')
print('0. GLOBAL')
print('1. COL')
print('2. CHAD')
print('3. IRAQ Dahuk')
print('4. IRAQ Najaf')
print('5. IRAQ')
print('6. LBN')
print('7. VEN')
print('9. SOM')
print('10. BGD')

pilot = input()
assert pilot in ['0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10'], "Invalid pilot area selected."

if pilot == "0":
    #input_shp = "C:/Geotar/COL/geodata/Processed/Education/Education_facilities.shp"
    pilot_name = "GLOBAL"
    mask_shp = f"C:/Geotar/{pilot_name}/geodata/workspace/test_mask.shp"
    period = "2021-05-01/2022-01-31"
    #output = f"C:/Geotar/COL/geodata/Processed/{res_folder}/dist_"+ out_name
    print("You selected GLobal")
elif pilot == "1":
    #input_shp = "C:/Geotar/COL/geodata/Processed/Education/Education_facilities.shp"
    pilot_name = "COL"
    mask_shp = f"C:/Geotar/{pilot_name}/geodata/Processed/Mask/COL_mask.shp"
    period = "2021-05-01/2022-01-31"
    #output = f"C:/Geotar/COL/geodata/Processed/{res_folder}/dist_"+ out_name
    print("You selected Colombia")
elif pilot == "2":
    pilot_name = "CHAD"
    #input_shp = "zip://C:/Geotar/CHAD/geodata/Processed/Education/hotosm_chad_education_facilities_points_shp.zip/hotosm_chad_education_facilities_points.shp"
    mask_shp = f"C:/Geotar/{pilot_name}/geodata/Processed/Mask/Chad_mask.shp"
    period = "2022-05-01/2023-01-31"
    #periodlta = "1970-01-01/1970-12-31"
    #output = f"C:/Geotar/CHAD/geodata/Processed/{res_folder}/dist_{res_folder}_"+ out_name
    print("You selected CHAD")
elif pilot == "3":
    pilot_name = "IRAQ_D"
    #input_shp = "zip://C:/Geotar/IRAQ_D/geodata/Raw/Education/hotosm_irq_education_facilities_points_shp.zip/hotosm_irq_education_facilities_points.shp"
    mask_shp = f"C:/Geotar/{pilot_name}/geodata/Processed/Mask/Dahuk_mask.shp"
    period = "2021-11-01/2022-05-31"
    #output = f"C:/Geotar/IRAQ_D/geodata/Processed/{res_folder}/dist_"+ out_name
    print("You selected IRAQ Dahuk")
elif pilot == "4":
    pilot_name = "IRAQ_N"
    period = "2021-11-01/2022-05-31"
    #input_shp = "zip://C:/Geotar/IRAQ_N/geodata/Raw/Education/hotosm_irq_education_facilities_points_shp.zip/hotosm_irq_education_facilities_points.shp"
    mask_shp = f"C:/Geotar/{pilot_name}/geodata/Processed/Mask/Najaf_mask.shp"
    #output = f"C:/Geotar/IRAQ_N/geodata/Processed/{res_folder}/dist_"+ out_name
    print("You selected IRAQ Najaf")
elif pilot == "5":
    pilot_name = "IRAQ"
    period = "2021-11-01/2022-05-31"
    #input_shp = "zip://C:/Geotar/IRAQ_N/geodata/Raw/Education/hotosm_irq_education_facilities_points_shp.zip/hotosm_irq_education_facilities_points.shp"
    mask_shp = f"C:/Geotar/{pilot_name}/geodata/Processed/Mask/Iraq_mask.shp"
    #output = f"C:/Geotar/IRAQ_N/geodata/Processed/{res_folder}/dist_"+ out_name
    print("You selected IRAQ")
elif pilot == "6":
    pilot_name = "LBN"
    period = "2021-10-01/2022-04-30"
    #input_shp = "zip://C:/Geotar/IRAQ_N/geodata/Raw/Education/hotosm_irq_education_facilities_points_shp.zip/hotosm_irq_education_facilities_points.shp"
    mask_shp = f"C:/Geotar/{pilot_name}/geodata/Processed/Mask/LBN_mask.shp"
    #output = f"C:/Geotar/IRAQ_N/geodata/Processed/{res_folder}/dist_"+ out_name
    print("You selected Lebanon")
elif pilot == "7":
    pilot_name = "VEN"
    period = "2023-01-01/2023-07-30"
    #input_shp = "zip://C:/Geotar/IRAQ_N/geodata/Raw/Education/hotosm_irq_education_facilities_points_shp.zip/hotosm_irq_education_facilities_points.shp"
    mask_shp = f"C:/Geotar/{pilot_name}/geodata/Processed/Mask/VEN_mask.shp"
    #output = f"C:/Geotar/IRAQ_N/geodata/Processed/{res_folder}/dist_"+ out_name
    print("You selected Venezuela")
elif pilot == "8":
    pilot_name = "AFG"
    period = "2023-04-01/2023-07-30"
    #input_shp = "zip://C:/Geotar/IRAQ_N/geodata/Raw/Education/hotosm_irq_education_facilities_points_shp.zip/hotosm_irq_education_facilities_points.shp"
    mask_shp = f"C:/Geotar/{pilot_name}/geodata/Processed/Mask/{pilot_name}_mask.shp"
    #output = f"C:/Geotar/IRAQ_N/geodata/Processed/{res_folder}/dist_"+ out_name
    print("You selected Afghanistan")
elif pilot == "9":
    pilot_name = "SOM"
    period = "2023-04-01/2023-07-30"
    #input_shp = "zip://C:/Geotar/IRAQ_N/geodata/Raw/Education/hotosm_irq_education_facilities_points_shp.zip/hotosm_irq_education_facilities_points.shp"
    mask_shp = f"C:/Geotar/{pilot_name}/geodata/Processed/Mask/{pilot_name}_mask.shp"
    #output = f"C:/Geotar/IRAQ_N/geodata/Processed/{res_folder}/dist_"+ out_name
    print("You selected Somalia")
elif pilot == "10":
    pilot_name = "BGD"
    period = "2023-04-01/2023-07-30"
    #input_shp = "zip://C:/Geotar/IRAQ_N/geodata/Raw/Education/hotosm_irq_education_facilities_points_shp.zip/hotosm_irq_education_facilities_points.shp"
    mask_shp = f"C:/Geotar/{pilot_name}/geodata/Processed/Mask/{pilot_name}_mask.shp"
    #output = f"C:/Geotar/IRAQ_N/geodata/Processed/{res_folder}/dist_"+ out_name
    print("You selected Bangladesh")
    



def process_ndvi(bbox, period ):
    """
    Fetches and computes the NDVI for the specified period.
    Returns:
    - Raster file with mean NDVI for the selected period.
    """
    NDVI = hdc_stac_client.search(bbox=bbox,
                                  # collections=["mod13q1_vim_native"],
                                  collections=["mxd13q1_vim_dekad"],
                                  datetime=period,  # emulates the period of data cube files
                                  ).get_all_items()

    res = 0.0022457882102988 # 250 or 0.01 for 1km
    ndvi_stack = stac_load(NDVI,  output_crs='EPSG:4326', resolution= res, patch_url=signer, bbox=bbox)
    ndvi_stack
    # Aggregate the dekadal NDVI by month
    m_ndvi = ndvi_stack.groupby('time.month').mean('time')
    m_ndvi = m_ndvi * 0.0001
    m_ndvi

    # Mask out extreme values
    m_ndvi = xr.where(m_ndvi < -1, np.nan, m_ndvi)
    m_ndvi = xr.where(m_ndvi > 1, np.nan, m_ndvi)
    m_ndvi

    #Aggregate the data by season
    ndvi_m_s = m_ndvi.mean(dim=["month"])
    ndvi_m_s

    output_dir_s = f"C:/Geotar/{pilot_name}/geodata/Processed/vegetation/season"
    filename_m_s = f'{output_dir_s}/ndvi_m_s.tif'
    if not os.path.exists(output_dir_s):
    os.makedirs(output_dir_s)

    # write the data to a geotiff file
    #ndvi_m_s.rio.to_raster(filename_m_s, driver='GTiff')
    print(f"{filename_m_s} saved successfully")

    # Mask NDVI using land cover
    tiff_path = f"C:/Geotar/{pilot_name}/geodata/Processed/LandCover/Worldcover_{pilot_name}.tif"
    mask_array = rioxarray.open_rasterio(tiff_path)
    mask_array = mask_array.squeeze("band", drop=True)
    mask_array = mask_array.rename({'x': 'longitude','y': 'latitude'})
    mask_array
    mask_array = mask_array.transpose('latitude', 'longitude')
    mask_array
    latitude = ndvi_m_s['latitude'].values
    longitude = ndvi_m_s['longitude'].values
    # Convert mask_array to an xarray DataArray
    mask_dataarray = xr.DataArray(mask_array, coords={'latitude': latitude, 'longitude': longitude}, dims=['latitude', 'longitude'])
    xr.align(ndvi_m_s, mask_dataarray, join='exact')  # will raise a ValueError if not aligned
    # Create a mask where conditions are not met and set to 0 where conditions are met
    ndvi_masked = xr.where((mask_dataarray == 40) | (mask_dataarray == 50), ndvi_m_s, 0)
    ndvi_masked = ndvi_masked.drop_vars('spatial_ref')
    ndvi_masked

    # save the masked ndvi data
    ndvi_masked.rio.to_raster(f'{output_dir_s}/ndvi_m_s.tif', driver='GTiff')

    # Save the monthly ndvi files

    # create a directory to store the geotiff files
    # output_dir = f"C:/Geotar/{pilot_name}/geodata/Processed/Vegetation"
    if not os.path.exists(output_dir):
    os.makedirs(output_dir)

    # loop over all timesteps in the dataset
    for i in range(len(m_ndvi.month)):
        # extract a single timestep as a DataArray
        image_ndvi = m_ndvi.isel(month=i)
        # create a file path for the geotiff file
        filename = f'{output_dir}/ndvi_{m_ndvi.month.values[i]}.tif'

        # write the data to a geotiff file
        image_ndvi.rio.to_raster(filename, driver='GTiff')
        print(f"{filename} saved successfully")
    return()

process_ndvi()

# ## Fetch the dekadal NDVI anomaly data from HDC

def process_ndvi_anomaly():
    """
        Fetches and computes the NDVI anomaly for the specified period.
        Returns:
         - raster file with the mean NDVI anomaly for the period selected.
         - raster file with the maximum NDVI anomaly for the period selected.
        """
    #lt_dates = "2002-07-01/2018-07-01"
    query_ndvi_anom = hdc_stac_client.search(bbox=bbox,
    #collections=["mod13q1_vim_native"],
    collections=["mxd13q1_viq_dekad"], #mxd13a2_vim_dekad_lta
    datetime= period).get_all_items()#1970-01-01T00:00:00Z/1970-12-31T00:00:00Z
    res = 0.0022457882102988 # 250 or 0.01 for 1km
    ndvi_anom = stac_load(query_ndvi_anom, patch_url=signer, output_crs='EPSG:4326', resolution= res, bbox=bbox, chunks ={})
    ndvi_anom
    # create a directory to store the geotiff files
    output_dir_zarr = f"C:/Geotar/{pilot_name}/geodata/zarr"
    if not os.path.exists(output_dir_zarr):
    os.makedirs(output_dir_zarr)
    outfile = output_dir_zarr + "/ndvi_an_stac.zarr"

    # Delete the existing Zarr store if it exists
    if os.path.exists(outfile):
    shutil.rmtree(outfile)

    ndvi_anom.to_zarr(outfile)
    print(f"{outfile} saved")

    output_dir_zarr = f"C:/Geotar/{pilot_name}/geodata/zarr"
    outfile_zarr = output_dir_zarr + "/ndvi_an_stac.zarr"
    #Load zarr file containing NDVI anomaly dekadal data for the season
    ndvi_anom = xr.open_zarr(outfile_zarr)
    ndvi_anom


    # Mask out no data/extreme values

    ndvi_anom_masked = xr.where(ndvi_anom < -150, np.nan, ndvi_anom)
    ndvi_anom_masked = xr.where(ndvi_anom_masked > 150, np.nan, ndvi_anom_masked)


    # Group the data by month

    m_ndvi_anom = ndvi_anom_masked.groupby('time.month').mean('time')
    m_ndvi_lta_r = m_ndvi_anom #/100  this rescaling is not working
    m_ndvi_lta_r

    s_ndvi_anom = m_ndvi_anom.mean(dim=["month"])
    s_ndvi_anom
    image = s_ndvi_anom*0.01
    image


    # Mask NDVI anomaly using land cover

    # Create a mask where conditions are not met and set to 0 where conditions are met
    ndvi_anom_masked = xr.where((mask_dataarray == 40) | (mask_dataarray == 50), image, 0)
    ndvi_anom_masked = ndvi_anom_masked.drop_vars('spatial_ref')
    ndvi_anom_masked


    output_dir_s = f"C:/Geotar/{pilot_name}/geodata/Processed/Vegetation/season"

    filename_s = f'{output_dir_s}/ndvi_a_m.tif'
    if not os.path.exists(output_dir_s):
        os.makedirs(output_dir_s)

    # write the data to a geotiff file
    ndvi_anom_masked.rio.to_raster(filename_s, driver='GTiff', band_indices=[1] , mask_and_scale=True)
    print(f"{filename_s} saved successfully")

    s_ndvi_anom_max = m_ndvi_anom.max(dim=["month"])
    s_ndvi_anom_max
    image_max = s_ndvi_anom_max*0.01
    output_dir_s = f"C:/Geotar/{pilot_name}/geodata/Processed/Vegetation/season"
    filename_s_max = f'{output_dir_s}/ndvi_a_ma.tif'
    if not os.path.exists(output_dir_s):
        os.makedirs(output_dir_s)

    # write the data to a geotiff file
    image_max.rio.to_raster(filename_s_max, driver='GTiff')
    print(f"{filename_s_max} saved successfully")


    # Exports a geotiff file for each date

    # create a directory to store the geotiff files
    output_dir = f"C:/Geotar/{pilot_name}/geodata/Processed/Vegetation"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # loop over all timesteps in the dataset
    for i in range(len(m_ndvi_lta_r.month)):
        # extract a single timestep as a DataArray
        image = m_ndvi_lta_r.isel(month=i)
        # create a file path for the geotiff file
        filename = f'{output_dir}/ndvi_a_{m_ndvi_lta_r.month.values[i]}.tif'

        # write the data to a geotiff file
        image.rio.to_raster(filename, driver='GTiff')
        print(f"{filename} saved successfully")
    return()

process_ndvi_anomaly()


# # CHIRPS processing

# The dates for the precipitation need to be adjusted, the idea is to account for a lag effect of a change the precipitation. The date shift was arbitrarily modified by one month in advance compared with the NDVI dates.

def process_CHIRPS():

    # print the NDVI period:
    print("NDVI period:", period)
    # Convert the period string to datetime objects
    start_date, end_date = period.split('/')
    start_date = datetime.datetime.strptime(start_date, "%Y-%m-%d")
    end_date = datetime.datetime.strptime(end_date, "%Y-%m-%d")

    # Subtract one month from the start and end dates while preserving the day of the month
    modified_start_date = start_date - relativedelta(months=1)
    modified_end_date = end_date - relativedelta(months=1)

    # Format the modified dates back to the desired string format
    modified_period = modified_start_date.strftime("%Y-%m-%d") + '/' + modified_end_date.strftime("%Y-%m-%d")

    # Update the period variable with the modified period
    period = modified_period

    # Print the updated period
    print("CHIRPS period:", period)

    CHIRPS = hdc_stac_client.search(bbox=bbox,
        #collections=['mod13q1_vim_native'],
        collections=['rfh_dekad'],
        datetime= period, #'2022-01-01/2022-12-31'
    ).get_all_items()

    res = 0.0022457882102988  #50.045454545Km #0.0022457882102988 # 250 or 0.01 for 1km

    CHIRPS_stac = stac_load(CHIRPS, output_crs='EPSG:4326', resolution=res , patch_url=signer, bbox=bbox)
    CHIRPS_stac


    # Paths to Save the zarr file

    # create a directory to store the geotiff files
    output_dir_zarr = f'C:/Geotar/{pilot_name}/geodata/zarr'
    if not os.path.exists(output_dir_zarr):
        os.makedirs(output_dir_zarr)
    outfile = output_dir_zarr + '/CHIRPS_stac.zarr'


    # Save the zarr file


    # Delete the existing Zarr store if it exists
    if os.path.exists(outfile):
        shutil.rmtree(outfile)
    # save the zarr file
    CHIRPS_stac.to_zarr(outfile)
    print(f'{outfile} saved')


    # Load xarray from zarr file (optional)

    #CHIRPS_stac = xr.open_zarr(outfile)


    # Mask out no data values


    CHIRPS_stac = xr.where(CHIRPS_stac == -9999, np.nan, CHIRPS_stac)
    CHIRPS_stac


    # Group CHIRPS data by month


    CHIRPS_m = CHIRPS_stac.groupby('time.month').sum('time', skipna=False)
    CHIRPS_m


    # Creates the seasonal monthly mean precipitation

    CHIRPS_m_s = CHIRPS_m.mean(dim=['month'])
    CHIRPS_m_s

    output_dir_s = f'C:/Geotar/{pilot_name}/geodata/Processed/precipitation/season'
    filename_m_s = f'{output_dir_s}/rain_m_s.tif'
    if not os.path.exists(output_dir_s):
        os.makedirs(output_dir_s)

    # write the data to a geotiff file
    CHIRPS_m_s.rio.to_raster(filename_m_s, driver='GTiff')
    print(f'{filename_m_s} saved successfully')


    # Creates the seasonal monthly sum precipitation

    # In[ ]:


    CHIRPS_s_s = CHIRPS_m.sum(dim=['month'])
    CHIRPS_s_s

    output_dir_s = f'C:/Geotar/{pilot_name}/geodata/Processed/precipitation/season'
    filename_s_s = f'{output_dir_s}/rain_s_s.tif'
    if not os.path.exists(output_dir_s):
        os.makedirs(output_dir_s)

    # write the data to a geotiff file
    CHIRPS_s_s.rio.to_raster(filename_s_s, driver='GTiff')
    print(f'{filename_s_s} saved successfully')
    return()

process_CHIRPS()

# Query and access the dekadal **CHIRPS precipitation anomaly** stored in the datacube

def process_CHIRPS_Anomaly():

    CHIRPS_an = hdc_stac_client.search(bbox=bbox,
        #collections=['mod13q1_vim_native'],
        collections=['rfq_dekad'],
        datetime= period #'2022-01-01/2022-12-31',
    ).get_all_items()
    #print(stac_items)
    res = 0.0022457882102988

    CHIRPS_an_stac = stac_load(CHIRPS_an, output_crs='EPSG:4326', resolution=res , patch_url=signer, bbox=bbox)
    CHIRPS_an_stac

    # create a directory to store the geotiff files
    output_dir_zarr = f'C:/Geotar/{pilot_name}/geodata/zarr'
    if not os.path.exists(output_dir_zarr):
        os.makedirs(output_dir_zarr)
    outfile = output_dir_zarr + '/CHIRPS_an_stac.zarr'

    # Delete the existing Zarr store if it exists
    if os.path.exists(outfile):
        shutil.rmtree(outfile)
    CHIRPS_an_stac.to_zarr(outfile)
    print(f'{outfile} saved')


    # Load the xarray anomay from a zarr file (optional)

    #CHIRPS_an_stac = xr.open_zarr(outfile)
    #CHIRPS_an_stac


    # Mask out/set no data values

    CHIRPS_an_stac = xr.where(CHIRPS_an_stac == -9999, np.nan, CHIRPS_an_stac)
    CHIRPS_an_stac


    # Aggregate the dekadal Chirps anomaly data by month

    CHIRPS_an_m = CHIRPS_an_stac.groupby('time.month').mean('time')
    #rescale to have values from 0 to 1
    CHIRPS_an_m = CHIRPS_an_m/100
    CHIRPS_an_m


    # Aggregate mean anomaly data by season

    CHIRPS_an_s = CHIRPS_an_m.mean(dim=['month'])
    CHIRPS_an_s

    image_an = CHIRPS_an_s
    output_dir_s = f'C:/Geotar/{pilot_name}/geodata/Processed/precipitation/season'
    filename_s = f'{output_dir_s}/rain_an_m.tif'
    if not os.path.exists(output_dir_s):
        os.makedirs(output_dir_s)

    # write the data to a geotiff file
    image_an.rio.to_raster(filename_s, driver='GTiff')
    print(f'{filename_s} saved successfully')


    # Agreggate anomaly max data season

    CHIRPS_an_s_max = CHIRPS_an_m.max(dim=['month'])
    CHIRPS_an_s_max
    image_an_m = CHIRPS_an_s_max # Rescaling applied here
    filename_s = f'{output_dir_s}/rain_an_ma.tif'
    # write the data to a geotiff file
    image_an_m.rio.to_raster(filename_s, driver='GTiff')
    print(f'{filename_s} saved successfully')


    # Export CHIRPS Rainfall monthly data

    # create a directory to store the geotiff files
    output_dir = f'C:/Geotar/{pilot_name}/geodata/Processed/Precipitation'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # loop over all timesteps in the dataset
    for i in range(len(CHIRPS_m.month)):
        # extract a single timestep as a DataArray
        da = CHIRPS_m.isel(month=i)
        # Set NaN as NoData
        #print("The no data value is:",da.rio.nodata)

        # create a file path for the geotiff file
        filename = f'{output_dir}/rain_{CHIRPS_m.month.values[i]}.tif'

        # write the data to a geotiff file
        da.rio.to_raster(filename)
        print(f'{filename} saved successfully')


    # Export the monthly anomaly data to GeoTiff files

    # create a directory to store the geotiff files
    output_dir = f'C:/Geotar/{pilot_name}/geodata/Processed/Precipitation'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # loop over all timesteps in the dataset
    for i in range(len(CHIRPS_an_m.month)):
        # extract a single timestep as a DataArray
        da = CHIRPS_an_m.isel(month=i)

        # create a file path for the geotiff file
        filename = f'{output_dir}/rain_a_{CHIRPS_an_m.month.values[i]}.tif'

        # write the data to a geotiff file
        da.rio.to_raster(filename, driver='GTiff')
        print(f'{filename} saved successfully')
    return()

process_CHIRPS_Anomaly()


# # Land Surface temperature processing

# Need to readjust the dates of the period to have same dates as the NDVI
def process_LST():

    # print the NDVI period:
    print("CHIRPS period:", period)
    # Convert the period string to datetime objects
    start_date, end_date = period.split('/')
    start_date = datetime.datetime.strptime(start_date, "%Y-%m-%d")
    end_date = datetime.datetime.strptime(end_date, "%Y-%m-%d")

    # add one month from the start and end dates while preserving the day of the month
    modified_start_date = start_date + relativedelta(months=1)
    modified_end_date = end_date + relativedelta(months=1)

    # Format the modified dates back to the desired string format
    modified_period = modified_start_date.strftime("%Y-%m-%d") + '/' + modified_end_date.strftime("%Y-%m-%d")

    # Update the period variable with the modified period
    period = modified_period

    # Print the updated period
    print("LST period:", period)

    LST_query = hdc_stac_client.search(bbox=bbox,
        #collections=['mod13q1_vim_native'],
        collections=['myd11a2_txa_dekad'],
        datetime= period #'2022-01-01/2022-12-31'
                                    ).get_all_items()
    res = 0.0022457882102988 # 250 or 0.01 for 1km
    LST = stac_load(LST_query, output_crs='EPSG:4326', resolution= res , patch_url=signer, bbox=bbox)
    LST

    output_dir_zarr = f'C:/Geotar/{pilot_name}/geodata/zarr'
    outfile = output_dir_zarr + '/LST_stac.zarr'


    # create a directory to store the geotiff files
    if not os.path.exists(output_dir_zarr):
        os.makedirs(output_dir_zarr)

    # Delete the existing Zarr file if it exists
    if os.path.exists(outfile):
        shutil.rmtree(outfile)

    LST.to_zarr(outfile)
    print(f'{outfile} saved')

    # create a directory to store the geotiff files
    if not os.path.exists(output_dir_zarr):
        os.makedirs(output_dir_zarr)

    # Delete the existing Zarr file if it exists
    if os.path.exists(outfile):
        shutil.rmtree(outfile)


    # Group LST and LST anomalies by month

    # group the lST data by month
    LST_m = LST.drop('tna')
    LST_m = LST_m.drop('spatial_ref')
    LST_m = LST_m.groupby('time.month').mean('time')
    LST_m = (LST_m * 0.02)- 273.15
    LST_m

    return()

process_LST()


# Seasonal mean anomaly

def process_LST_anomaly():
    LST_anom_query = hdc_stac_client.search(bbox=bbox,
                                            # collections=['mod13q1_vim_native'],
                                            collections=['myd11a2_txd_dekad'],
                                            datetime=period  # '2022-01-01/2022-12-31'
                                            ).get_all_items()
    res = 0.0022457882102988  # 250 or 0.01 for 1km
    LST_anom = stac_load(LST_anom_query, output_crs='EPSG:4326', resolution=res, patch_url=signer, bbox=bbox)
    LST_anom
    # group the lST anomaly data by month
    LST_an_m = LST_anom.drop('tnd')
    LST_an_m = LST_an_m.drop('spatial_ref')
    LST_an_m = LST_an_m.groupby('time.month').mean('time')
    LST_an_m = LST_an_m  * 0.02
    #LST_an_m
    LST_an_m

    LST_an_s = LST_an_m.mean(dim=['month'])

    image_an = LST_an_s/100 # Rescaling applied here
    output_dir_s = f'C:/Geotar/{pilot_name}/geodata/Processed/Temperature/season'
    filename_s = f'{output_dir_s}/LST_an_m.tif'
    if not os.path.exists(output_dir_s):
        os.makedirs(output_dir_s)

    # write the data to a geotiff file
    image_an.rio.to_raster(filename_s, driver='GTiff')
    print(f'{filename_s} saved successfully')
    #LST_an_s


    # Seasonal max anomaly

    # In[ ]:


    LST_an_s_max = LST_an_m.max(dim=['month'])
    LST_an_s_max
    image_an_max = LST_an_s_max/100 # Rescaling applied here
    filename_s_max = f'{output_dir_s}/LST_an_ma.tif'
    # write the data to a geotiff file
    image_an_max.rio.to_raster(filename_s_max, driver='GTiff')
    print(f'{filename_s_max} saved successfully')


    # Export the mean seasonal temperature

    # In[ ]:


    LST_s = LST_m.mean(dim=['month'])
    LST_s
    image_m = LST_s # Rescaling applied here
    output_dir_s = f'C:/Geotar/{pilot_name}/geodata/Processed/Temperature/season'
    filename_s = f'{output_dir_s}/LST_m.tif'
    if not os.path.exists(output_dir_s):
        os.makedirs(output_dir_s)

    # write the data to a geotiff file
    image_m.rio.to_raster(filename_s, driver='GTiff')
    print(f'{filename_s} saved successfully')


    # Export the max seasonal temperature

    # In[ ]:


    LST_s_max = LST_m.max(dim=['month'])
    LST_s_max
    image_max = LST_s_max # Rescaling applied here
    output_dir_s_max = f'C:/Geotar/{pilot_name}/geodata/Processed/Temperature/season'
    filename_s_max = f'{output_dir_s_max}/LST_ma.tif'
    if not os.path.exists(output_dir_s_max):
        os.makedirs(output_dir_s_max)

    # write the data to a geotiff file
    image_max.rio.to_raster(filename_s_max, driver='GTiff')
    print(f'{filename_s_max} saved successfully')


    # Export monthly LST files

    # In[ ]:


    #LST_m = LST_m.drop('tnd')
    # create a directory to store the geotiff files
    output_dir = f'C:/Geotar/{pilot_name}/geodata/Processed/Temperature/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # loop over all timesteps in the dataset
    for i in range(len(LST_m.month)):
        # extract a single timestep as a DataArray
        da = LST_m.isel(month=i)

        # create a file path for the geotiff file
        filename = f'{output_dir}/LST_{LST_m.month.values[i]}.tif'
        # write the data to a geotiff file
        da.tda.rio.to_raster(filename, driver='GTiff')
        print(f'{filename} saved successfully')


    # Export the monthly LST anomaly data to GeoTiff files

    # In[ ]:


    # create a directory to store the geotiff files
    output_dir = f'C:/Geotar/{pilot_name}/geodata/Processed/Temperature'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # loop over all timesteps in the dataset
    for i in range(len(LST_an_m.month)):
        # extract a single timestep as a DataArray
        da = LST_an_m.isel(month=i)

        # create a file path for the geotiff file
        filename = f'{output_dir}/LST_an_{LST_an_m.month.values[i]}.tif'

        # write the data to a geotiff file
        da.rio.to_raster(filename, driver='GTiff')
        print(f'{filename} saved successfully')
    return()

process_LST_anomaly()
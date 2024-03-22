#!/usr/bin/env python
# coding: utf-8

# NDVI anomalies

import geopandas as gpd
import shapely.geometry
from pystac_client import Client
from odc.stac import configure_rio, stac_load
import pathlib
import json

from process_functions import *


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

hdc_stac_client, signer = _get_hdc_stac_param_from_env() 

masks = {
    'GLOBAL': {  "pilot_name" : "GLOBAL'",
            "mask_shp" : f"C:/Geotar/GLOBAL/geodata/workspace/test_mask.shp",
            "period" : "2021-05-01/2022-01-31", 'country_name':'All Countries'},
    'COL': {"pilot_name" : "COL'",
            "mask_shp" : f"C:/Geotar/COL/geodata/workspace/COL_mask.shp",
            "period" : "2021-05-01/2022-01-31"},
    'CHAD': {"pilot_name": "CHAD",
            "mask_shp": f"C:/Geotar/CHAD/geodata/workspace/COL_mask.shp",
            "period": "2021-05-01/2022-01-31"},
}



if __name__ =='__main__':

    #TODO: 1. requirement files 2. fill in the dictionary 3. define inputs of the functions

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

    area_shp = gpd.read_file(masks[pilot]['mask_shp'])
    period = masks[pilot]['period']
    pilot_name = masks[pilot]['pilot_name']


    # Get the bounding box of the shapefile
    bbox = area_shp.total_bounds

    process_ndvi(bbox=bbox, period=period, pilot_name=pilot_name)
    process_CHIRPS_Anomaly(bbox=bbox, period=period)
    process_ndvi_anomaly(bbox=bbox, period=period)
    process_LST_anomaly(bbox=bbox, period=period)
    process_LST(bbox=bbox, period=period)
    process_CHIRPS(bbox=bbox, period=period)
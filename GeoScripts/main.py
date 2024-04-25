#!/usr/bin/env python
# coding: utf-8

# NDVI anomalies

import geopandas as gpd
import shapely.geometry
from pystac_client import Client
from odc.stac import configure_rio, stac_load
import pathlib
import json
import os
from process_functions import *
#from get_polygon_data import *
#from prox_GDAL import *


root = pathlib.Path().resolve().parent.parent
root

TOKEN_PATH = "C:/Users/oscar.bautista/OneDrive - World Food Programme/Scripts/tk.json"
HDC_STAC_URL= "https://api.earthobservation.vam.wfp.org/stac/"


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

masks = {
    'GLOBAL': { "pilot_name" : "GLOBAL'",
            "mask_shp" : f"C:/Geotar/GLOBAL/geodata/workspace/test_mask.shp",
            "period" : "2021-05-01/2022-01-31", 'country_name':'All Countries'},
    'COL': {"pilot_name" : "COL",
            "mask_shp" : f"C:/Geotar/COL/geodata/processed/mask/COL_mask.shp",
            "period" : "2021-05-01/2022-01-31"},
    'CHAD': {"pilot_name": "CHAD",
            "mask_shp": f"C:/Geotar/CHAD/geodata/processed/mask/CHAD_mask.shp",
            "period": "2021-05-01/2022-01-31"},
    'IRAQ': {"pilot_name": "IRAQ",
            "mask_shp": f"C:/Geotar/IRAQ/geodata/processed/mask/IRAQ_mask.shp",
            "period": "2021-05-01/2022-01-31"},
    'LBN': {"pilot_name": "LBN",
            "mask_shp": f"C:/Geotar/LBN/geodata/processed/mask/LBN_mask.shp",
            "period": "2023-10-01/2024-03-21"},
    'BGD': {"pilot_name": "BGD",
            "mask_shp": f"C:/Geotar/BGD/geodata/processed/mask/BGD_mask.shp",
            "period": "2021-05-01/2022-01-31"},
    'ETH': {"pilot_name": "ETH",
            "mask_shp": f"C:/Geotar/ETH/geodata/processed/mask/ETH_mask.shp",
            "period": "2023-01-01/2023-12-31"},
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
    print('11. ETH')

    pilot = input()
    hdc_stac_client, signer = _get_hdc_stac_param_from_env()

    area_shp = gpd.read_file(masks[pilot]['mask_shp'])
    period = masks[pilot]['period']
    pilot_name = masks[pilot]['pilot_name']


    # Get the bounding box of the shapefile
    bbox = area_shp.total_bounds

    # Instantiate the Process class
    process_obj = Process(bbox=bbox,
                  period=period,
                  pilot_name=pilot_name,
                  hdc_stac_client=hdc_stac_client,
                  signer=signer)
    process_obj.process_ndvi()
    process_obj.process_ndvi_anomaly()
    process_obj.process_CHIRPS()
    process_obj.process_CHIRPS_Anomaly()
    process_obj.process_LST()
    process_obj.process_LST_anomaly()


    # process_ndvi(bbox=bbox,
    #              period=period,
    #              pilot_name=pilot_name,
    #              hdc_stac_client=hdc_stac_client,
    #              signer=signer)
    # process_ndvi_anomaly(bbox=bbox, period=period, pilot_name=pilot_name)
    # process_CHIRPS(bbox=bbox, period=period, pilot_name=pilot_name)
    # process_CHIRPS_Anomaly(bbox=bbox, period=period, pilot_name=pilot_name)
    # process_LST(bbox=bbox, period=period, pilot_name=pilot_name)
    # process_LST_anomaly(bbox=bbox, period=period, pilot_name=pilot_name)

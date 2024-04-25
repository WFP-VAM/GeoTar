import os
import geopandas as gpd
import shapely.geometry
from pystac_client import Client
from odc.stac import configure_rio, stac_load

import json
import rasterio
import xarray as xr
import numpy as np
import pandas as pd
import zarr
import rioxarray
from osgeo import gdal
import shutil
import os
from typing import List
import os
import datetime
from dateutil.relativedelta import relativedelta


class Vector:
    def __init__(self, bbox: List, period: str, pilot_name: str):

        self.bbox=bbox
        self.period=period
        self.pilot_name=pilot_name
        self.hdc_stac_client=hdc_stac_client
        self.signer=signer

   def
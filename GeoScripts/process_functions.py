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


class Process:
    def __init__(self, bbox: List, period: str, pilot_name: str, hdc_stac_client, signer):

        self.bbox=bbox
        self.period=period
        self.pilot_name=pilot_name
        self.hdc_stac_client=hdc_stac_client
        self.signer=signer

    def process_ndvi(self):
        """
        Fetches and computes the NDVI for the specified period.
        Returns:
        - Raster file with mean NDVI for the selected period.
        """
        print("Processing NDVI")
        NDVI = self.hdc_stac_client.search(bbox=self.bbox,
                                      # collections=["mod13q1_vim_native"],
                                      collections=["mxd13q1_vim_dekad"],
                                      datetime=self.period,  # emulates the period of data cube files
                                      ).get_all_items()

        res = 0.0022457882102988 # 250 or 0.01 for 1km
        ndvi_stack = stac_load(NDVI,  output_crs='EPSG:4326', resolution= res, patch_url=self.signer, bbox=self.bbox)
        # Aggregate the dekadal NDVI by month
        m_ndvi = ndvi_stack.groupby('time.month').mean('time')
        m_ndvi = m_ndvi * 0.0001

        # Mask out extreme values
        m_ndvi = xr.where(m_ndvi < -1, np.nan, m_ndvi)
        m_ndvi = xr.where(m_ndvi > 1, np.nan, m_ndvi)

        #Aggregate the data by season
        ndvi_m_s = m_ndvi.mean(dim=["month"])

        output_dir_s = f"C:/Geotar/{self.pilot_name}/geodata/Processed/vegetation/season"
        filename_m_s = f'{output_dir_s}/ndvi_m_s.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        ndvi_m_s.rio.to_raster(filename_m_s, driver='GTiff')
        print(f"{filename_m_s} saved successfully")

        # Mask NDVI using land cover
        tiff_path = f"C:/Geotar/{self.pilot_name}/geodata/Processed/LandCover/Worldcover_{self.pilot_name}.tif"
        mask_array = rioxarray.open_rasterio(tiff_path)
        mask_array = mask_array.squeeze("band", drop=True)
        mask_array = mask_array.rename({'x': 'longitude','y': 'latitude'})

        mask_array = mask_array.transpose('latitude', 'longitude')

        latitude = ndvi_m_s['latitude'].values
        longitude = ndvi_m_s['longitude'].values
        # Convert mask_array to an xarray DataArray
        mask_dataarray = xr.DataArray(mask_array, coords={'latitude': latitude, 'longitude': longitude}, dims=['latitude', 'longitude'])
        xr.align(ndvi_m_s, mask_dataarray, join='exact')  # will raise a ValueError if not aligned
        # Create a mask where conditions are not met and set to 0 where conditions are met
        ndvi_masked = xr.where((mask_dataarray == 40) | (mask_dataarray == 50), ndvi_m_s, 0)
        ndvi_masked = ndvi_masked.drop_vars('spatial_ref')

        # save the masked ndvi data
        #ndvi_masked.rio.to_raster(f'{output_dir_s}/ndvi_m_s.tif', driver='GTiff')

        # Save the monthly ndvi files

        # create a directory to store the geotiff files
        output_dir = f"C:/Geotar/{self.pilot_name}/geodata/Processed/Vegetation"
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
        return

#
#
    def process_ndvi_anomaly(self):
        """
            Fetches and computes the NDVI anomaly for the specified period.
            Returns:
             - raster file with the mean NDVI anomaly for the period selected.
             - raster file with the maximum NDVI anomaly for the period selected.
            """
        print("Processing NDVI anomalies")
        #lt_dates = "2002-07-01/2018-07-01"
        query_ndvi_anom = self.hdc_stac_client.search(bbox=self.bbox,
        #collections=["mod13q1_vim_native"],
        collections=["mxd13q1_viq_dekad"], #mxd13a2_vim_dekad_lta
        datetime= self.period).get_all_items()#1970-01-01T00:00:00Z/1970-12-31T00:00:00Z
        res = 0.0022457882102988 # 250 or 0.01 for 1km
        ndvi_anom = stac_load(query_ndvi_anom, patch_url=self.signer, output_crs='EPSG:4326', resolution= res, bbox=self.bbox, chunks ={})
        ndvi_anom
        # create a directory to store the geotiff files
        # output_dir_zarr = f"C:/Geotar/{self.pilot_name}/geodata/zarr"
        # if not os.path.exists(output_dir_zarr):
        #     os.makedirs(output_dir_zarr)
        #     outfile = output_dir_zarr + "/ndvi_an_stac.zarr"

        # Delete the existing Zarr store if it exists
        # if os.path.exists(outfile):
        #     shutil.rmtree(outfile)

        # ndvi_anom.to_zarr(outfile)
        # print(f"{outfile} saved")

        # output_dir_zarr = f"C:/Geotar/{self.pilot_name}/geodata/zarr"
        # outfile_zarr = output_dir_zarr + "/ndvi_an_stac.zarr"
        # #Load zarr file containing NDVI anomaly dekadal data for the season
        # ndvi_anom = xr.open_zarr(outfile_zarr)
        # ndvi_anom


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



        # Mask NDVI anomaly using land cover
        cover_path = f'C:/Geotar/{self.pilot_name}/geodata/Processed/LandCover/Worldcover_{self.pilot_name}.tif'
        mask_array = rioxarray.open_rasterio(cover_path)
        mask_array = mask_array.squeeze("band", drop=True)
        mask_array = mask_array.rename({'x': 'longitude', 'y': 'latitude'})
        mask_array = mask_array.transpose('latitude', 'longitude')
        latitude = image['latitude'].values
        longitude = image['longitude'].values
        # Convert mask_array to an xarray DataArray
        mask_dataarray = xr.DataArray(mask_array, coords={'latitude': latitude, 'longitude': longitude},
                                      dims=['latitude', 'longitude'])
        xr.align(image, mask_dataarray, join='exact')

        # Create a mask where conditions are not met and set to 0 where conditions are met
        ndvi_anom_masked = xr.where((mask_dataarray == 40) | (mask_dataarray == 50), image, 0)
        ndvi_anom_masked = ndvi_anom_masked.drop_vars('spatial_ref')
        ndvi_anom_masked


        output_dir_s = f"C:/Geotar/{self.pilot_name}/geodata/Processed/Vegetation/season"

        filename_s = f'{output_dir_s}/ndvi_a_m.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        ndvi_anom_masked.rio.to_raster(filename_s, driver='GTiff', band_indices=[1] , mask_and_scale=True)
        print(f"{filename_s} saved successfully")

        s_ndvi_anom_max = m_ndvi_anom.max(dim=["month"])
        s_ndvi_anom_max
        image_max = s_ndvi_anom_max*0.01
        output_dir_s = f"C:/Geotar/{self.pilot_name}/geodata/Processed/Vegetation/season"
        filename_s_max = f'{output_dir_s}/ndvi_a_ma.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        image_max.rio.to_raster(filename_s_max, driver='GTiff')
        print(f"{filename_s_max} saved successfully")


        # Exports a geotiff file for each date

        # create a directory to store the geotiff files
        output_dir = f"C:/Geotar/{self.pilot_name}/geodata/Processed/Vegetation"
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
        return


    def process_CHIRPS(self):
        """
        :
        return:
        """
        print("Processing CHIRPS")

        # print the NDVI period:
        print("NDVI period:", self.period)
        # Convert the period string to datetime objects
        start_date, end_date = self.period.split('/')
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

        CHIRPS = self.hdc_stac_client.search(bbox=self.bbox,
            #collections=['mod13q1_vim_native'],
            collections=['rfh_dekad'],
            datetime= period, #'2022-01-01/2022-12-31'
        ).get_all_items()

        res = 0.0022457882102988  #50.045454545Km #0.0022457882102988 # 250 or 0.01 for 1km

        CHIRPS_stac = stac_load(CHIRPS, output_crs='EPSG:4326', resolution= res , patch_url=self.signer, bbox=self.bbox)
        CHIRPS_stac


        # Paths to Save the zarr file

        # create a directory to store the geotiff files
        #output_dir_zarr = f'C:/Geotar/{self.pilot_name}/geodata/zarr'
        #if not os.path.exists(output_dir_zarr):
        #    os.makedirs(output_dir_zarr)
        #outfile = output_dir_zarr + '/CHIRPS_stac.zarr'


        # Save the zarr file


        # Delete the existing Zarr store if it exists
        #if os.path.exists(outfile):
        #    shutil.rmtree(outfile)
        # save the zarr file
        #CHIRPS_stac.to_zarr(outfile)
        #print(f'{outfile} saved')


        # Load xarray from zarr file (optional)

        #CHIRPS_stac = xr.open_zarr(outfile)


        # Mask out no data values


        CHIRPS_stac = xr.where(CHIRPS_stac == -9999, np.nan, CHIRPS_stac)
        #CHIRPS_stac


        # Group CHIRPS data by month


        CHIRPS_m = CHIRPS_stac.groupby('time.month').sum('time', skipna=False)
        #CHIRPS_m


        # Creates the seasonal monthly mean precipitation

        CHIRPS_m_s = CHIRPS_m.mean(dim=['month'])
        #CHIRPS_m_s

        output_dir_s = f'C:/Geotar/{self.pilot_name}/geodata/Processed/precipitation/season'
        filename_m_s = f'{output_dir_s}/rain_m_s.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        CHIRPS_m_s.rio.to_raster(filename_m_s, driver='GTiff')
        print(f'{filename_m_s} saved successfully')

        # Export CHIRPS Rainfall monthly data

        # create a directory to store the geotiff files
        output_dir = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Precipitation'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # loop over all timesteps in the dataset
        for i in range(len(CHIRPS_m.month)):
            # extract a single timestep as a DataArray
            da = CHIRPS_m.isel(month=i)
            # Set NaN as NoData
            # print("The no data value is:",da.rio.nodata)

            # create a file path for the geotiff file
            filename = f'{output_dir}/rain_{CHIRPS_m.month.values[i]}.tif'

            # write the data to a geotiff file
            da.rio.to_raster(filename)
            print(f'{filename} saved successfully')


        # Creates the seasonal monthly sum precipitation

        CHIRPS_s_s = CHIRPS_m.sum(dim=['month'])
        #CHIRPS_s_s

        output_dir_s = f'C:/Geotar/{self.pilot_name}/geodata/Processed/precipitation/season'
        filename_s_s = f'{output_dir_s}/rain_s_s.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        CHIRPS_s_s.rio.to_raster(filename_s_s, driver='GTiff')
        print(f'{filename_s_s} saved successfully')
        return
#
#
#
#
    def process_CHIRPS_Anomaly(self):
        """

        Returns:

        """

        print("Processing CHIRPS anomalies")

        CHIRPS_an = self.hdc_stac_client.search(bbox=self.bbox,
            #collections=['mod13q1_vim_native'],
            collections=['rfq_dekad'],
            datetime= self.period #'2022-01-01/2022-12-31',
        ).get_all_items()
        #print(stac_items)

        res = 0.0022457882102988
        CHIRPS_an_stac = stac_load(CHIRPS_an, output_crs='EPSG:4326', resolution= res , patch_url=self.signer, bbox=self.bbox)
        CHIRPS_an_stac

        # # create a directory to store the geotiff files
        # output_dir_zarr = f'C:/Geotar/{pilot_name}/geodata/zarr'
        # if not os.path.exists(output_dir_zarr):
        #     os.makedirs(output_dir_zarr)
        # outfile = output_dir_zarr + '/CHIRPS_an_stac.zarr'
        #
        # # Delete the existing Zarr store if it exists
        # if os.path.exists(outfile):
        #     shutil.rmtree(outfile)
        # CHIRPS_an_stac.to_zarr(outfile)
        # print(f'{outfile} saved')



        # Mask out/set no data values

        CHIRPS_an_stac = xr.where(CHIRPS_an_stac == -9999, np.nan, CHIRPS_an_stac)
        #CHIRPS_an_stac


        # Aggregate the dekadal Chirps anomaly data by month

        CHIRPS_an_m = CHIRPS_an_stac.groupby('time.month').mean('time')
        #rescale to have values from 0 to 1
        CHIRPS_an_m = CHIRPS_an_m/100
        #CHIRPS_an_m


        # Aggregate mean anomaly data by season

        CHIRPS_an_s = CHIRPS_an_m.mean(dim=['month'])
        #CHIRPS_an_s

        image_an = CHIRPS_an_s
        output_dir_s = f'C:/Geotar/{self.pilot_name}/geodata/Processed/precipitation/season'
        filename_s = f'{output_dir_s}/rain_an_m.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        image_an.rio.to_raster(filename_s, driver='GTiff')
        print(f'{filename_s} saved successfully')


        # Agreggate anomaly max data season

        CHIRPS_an_s_max = CHIRPS_an_m.max(dim=['month'])
        #CHIRPS_an_s_max
        image_an_m = CHIRPS_an_s_max # Rescaling applied here
        filename_s = f'{output_dir_s}/rain_an_ma.tif'
        # write the data to a geotiff file
        image_an_m.rio.to_raster(filename_s, driver='GTiff')
        print(f'{filename_s} saved successfully')


        # Export the monthly anomaly data to GeoTiff files

        # create a directory to store the geotiff files
        output_dir = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Precipitation'
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
        return
#

# # # Land Surface temperature processing
#
    # Need to readjust the dates of the period to have same dates as the NDVI
    def process_LST(self):

        print("Processing Land surface temperature")

        # print the NDVI period:
        print("CHIRPS period:", self.period)
        # Convert the period string to datetime objects
        start_date, end_date = self.period.split('/')
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

        LST_query = self.hdc_stac_client.search(bbox=self.bbox,
            #collections=['mod13q1_vim_native'],
            collections=['myd11a2_txa_dekad'],
            datetime= self.period #'2022-01-01/2022-12-31'
                                        ).get_all_items()
        res = 0.0022457882102988 # 250 or 0.01 for 1km
        LST = stac_load(LST_query, output_crs='EPSG:4326', resolution= res , patch_url=self.signer, bbox=self.bbox)
        #LST

        # output_dir_zarr = f'C:/Geotar/{pilot_name}/geodata/zarr'
        # outfile = output_dir_zarr + '/LST_stac.zarr'
        #
        #
        # # create a directory to store the geotiff files
        # if not os.path.exists(output_dir_zarr):
        #     os.makedirs(output_dir_zarr)
        #
        # # Delete the existing Zarr file if it exists
        # if os.path.exists(outfile):
        #     shutil.rmtree(outfile)
        #
        # LST.to_zarr(outfile)
        # print(f'{outfile} saved')
        # output_dir_zarr = f'C:/Geotar/{self.pilot_name}/geodata/zarr'
        # outfile = output_dir_zarr + '/LST_an_stac.zarr'

        # # create a directory to store the geotiff files
        # if not os.path.exists(output_dir_zarr):
        #     os.makedirs(output_dir_zarr)
        #
        # # Delete the existing Zarr file if it exists
        # if os.path.exists(outfile):
        #     shutil.rmtree(outfile)
        #
        # LST_anom.to_zarr(outfile)
        # print(f'{outfile} saved')


        # # Load xarray from zarr file (optional)
        #
        # LST_anom= xr.open_zarr(outfile)
        # LST_anom




        # group the lST data by month
        LST_m = LST.drop('tna')
        LST_m = LST_m.drop('spatial_ref')
        LST_m = LST_m.groupby('time.month').mean('time')
        LST_m = (LST_m * 0.02)- 273.15
        #LST_m
        LST_s = LST_m.mean(dim=['month'])
        LST_s
        image_m = LST_s # Rescaling applied here
        output_dir_s = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Temperature/season'
        filename_s = f'{output_dir_s}/LST_m.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        image_m.rio.to_raster(filename_s, driver='GTiff')
        print(f'{filename_s} saved successfully')


        # Export the max seasonal temperature
        LST_s_max = LST_m.max(dim=['month'])
        LST_s_max
        image_max = LST_s_max # Rescaling applied here
        output_dir_s_max = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Temperature/season'
        filename_s_max = f'{output_dir_s_max}/LST_ma.tif'
        if not os.path.exists(output_dir_s_max):
            os.makedirs(output_dir_s_max)

        # write the data to a geotiff file
        image_max.rio.to_raster(filename_s_max, driver='GTiff')
        print(f'{filename_s_max} saved successfully')

        # create a directory to store the geotiff files
        output_dir = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Temperature/'
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

        return
#
#
#
#
    def process_LST_anomaly(self):
        print("Processing Land surface temperature anomalies")
        LST_anom_query = self.hdc_stac_client.search(bbox=self.bbox,
            #collections=['mod13q1_vim_native'],
            collections=['myd11a2_txd_dekad'],
            datetime= self.period #'2022-01-01/2022-12-31'
                                        ).get_all_items()
        res = 0.0022457882102988 # 250 or 0.01 for 1km
        LST_anom = stac_load(LST_anom_query, output_crs='EPSG:4326', resolution= res , patch_url=self.signer, bbox=self.bbox)
        LST_anom




    # group the lST anomaly data by month
        LST_an_m = LST_anom.drop('tnd')
        LST_an_m = LST_an_m.drop('spatial_ref')
        LST_an_m = LST_an_m.groupby('time.month').mean('time')
        LST_an_m = LST_an_m  * 0.02
        #LST_an_m
        #LST_an_m

        LST_an_s = LST_an_m.mean(dim=['month'])

        image_an = LST_an_s / 100  # Rescaling applied here
        output_dir_s = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Temperature/season'
        filename_s = f'{output_dir_s}/LST_an_m.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        image_an.rio.to_raster(filename_s, driver='GTiff')
        print(f'{filename_s} saved successfully')
        # LST_an_s

        LST_an_s_max = LST_an_m.max(dim=['month'])
        LST_an_s_max
        image_an_max = LST_an_s_max/100 # Rescaling applied here
        filename_s_max = f'{output_dir_s}/LST_an_ma.tif'
        # write the data to a geotiff file
        image_an_max.rio.to_raster(filename_s_max, driver='GTiff')
        print(f'{filename_s_max} saved successfully')


        # Export the monthly LST anomaly data to GeoTiff files

        # create a directory to store the geotiff files
        output_dir = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Temperature'
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
        return

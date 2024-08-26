import rasterio.crs
from odc.stac import configure_rio, stac_load
import xarray as xr
import numpy as np
import zarr
import rioxarray
from typing import List
import os
import datetime
from dateutil.relativedelta import relativedelta


class Process:
    def __init__(self, bbox: List, period: str, pilot_name: str, hdc_stac_client, signer):

        self.bbox = bbox
        self.period = period
        self.pilot_name = pilot_name
        self.hdc_stac_client = hdc_stac_client
        self.signer = signer

    def process_ndvi(self):
        """
        Fetches and computes the NDVI for the specified period.
        Returns:
        - Raster file with mean NDVI for the selected period.
        """
        print("Processing NDVI")
        ndvi = self.hdc_stac_client.search(bbox=self.bbox,
                                      collections=["mxd13a2_vim_dekad"],
                                      # collections=["mxd13q1_vim_dekad"],
                                      datetime=self.period,  # emulates the period of data cube files
                                      ).get_all_items()

        res = 0.0022457882102988 # 250 or 0.01 for 1km
        crs = rasterio.crs.CRS.from_epsg(4326)
        ndvi_stack = stac_load(ndvi,  output_crs='EPSG:4326', resolution=res, patch_url=self.signer, bbox=self.bbox)
        # Aggregate the dekadal NDVI by month
        m_ndvi = ndvi_stack.resample(time='1M').mean()
        m_ndvi = m_ndvi * 0.0001

        # Mask out extreme values
        m_ndvi = xr.where(m_ndvi < -1, np.nan, m_ndvi)
        m_ndvi = xr.where(m_ndvi > 1, np.nan, m_ndvi)

        #Aggregate the data by season
        ndvi_m_s = m_ndvi.mean(dim="time")
        # ndvi_m_s = ndvi_m_s.set_crs(crs)
        year = m_ndvi.time.dt.year.item(0)

        output_dir_s = f"C:/Geotar/{self.pilot_name}/geodata/Processed/vegetation/season"
        filename_m_s = f'{output_dir_s}/ndvi_m_s_{year}.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        ndvi_m_s.rio.to_raster(filename_m_s, driver='GTiff', crs='EPSG:4326')
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
        mask_dataarray = xr.DataArray(mask_array,
                                      coords={'latitude': latitude, 'longitude': longitude},
                                      dims=['latitude', 'longitude'])
        xr.align(ndvi_m_s, mask_dataarray, join='exact')  # will raise a ValueError if not aligned
        # Create a mask where conditions are not met and set to 0 where conditions are met
        ndvi_masked = xr.where((mask_dataarray == 40) | (mask_dataarray == 50), ndvi_m_s, 0)
        # ndvi_masked = ndvi_masked.drop_vars('spatial_ref')
        # ndvi_masked.set_crs(crs)

        # save the masked ndvi data
        # ndvi_masked.rio.to_raster({filename_m_s}, driver='GTiff')

        # Save the monthly ndvi files

        # create a directory to store the geotiff files
        output_dir = f"C:/Geotar/{self.pilot_name}/geodata/Processed/Vegetation"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # loop over all timesteps in the dataset
        for time_val in m_ndvi.time:
            year = time_val.dt.year.item()
            month = time_val.dt.month.item()

            # extract a single timestep as a DataArray
            image_ndvi = m_ndvi.sel(time=time_val)

            # create a file path for the geotiff file with both year and month in filename
            filename = f'{output_dir}/ndvi_{year}_{month:02d}.tif'

            # write the data to a geotiff file
            image_ndvi.rio.to_raster(filename, driver='GTiff', crs='EPSG:4326')
            print(f"{filename} saved successfully")
        return


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
        collections=["mxd13a2_viq_dekad"], #mxd13q1_viq_dekad
        datetime=self.period).get_all_items()#1970-01-01T00:00:00Z/1970-12-31T00:00:00Z
        res = 0.0022457882102988 # 250 or 0.01 for 1km
        crs = rasterio.crs.CRS.from_epsg(4326)
        ndvi_anom = stac_load(query_ndvi_anom, patch_url=self.signer, output_crs='EPSG:4326',
                              resolution=res, bbox=self.bbox, chunks={})
        #ndvi_anom
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

        m_ndvi_anom = ndvi_anom_masked.resample(time='1M').mean()
        m_ndvi_lta_r = m_ndvi_anom #/100  this rescaling is not working
        #m_ndvi_lta_r

        s_ndvi_anom = m_ndvi_anom.mean(dim="time")
        year = m_ndvi_anom.time.dt.year.item(0)
        #s_ndvi_anom
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
        # ndvi_anom_masked.set_crs(crs)
        #ndvi_anom_masked

        output_dir_s = f"C:/Geotar/{self.pilot_name}/geodata/Processed/Vegetation/season"

        filename_s = f'{output_dir_s}/ndvi_a_m_{year}.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        ndvi_anom_masked.rio.to_raster(filename_s, driver='GTiff', band_indices=[1], mask_and_scale=True)
        print(f"{filename_s} saved successfully")

        s_ndvi_anom_max = m_ndvi_anom.max(dim="time")
        #s_ndvi_anom_max
        image_max = s_ndvi_anom_max*0.01
        image_max.rio.set_crs(crs)
        output_dir_s = f"C:/Geotar/{self.pilot_name}/geodata/Processed/Vegetation/season"
        filename_s_max = f'{output_dir_s}/ndvi_a_ma_{year}.tif'
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
        for time_val in m_ndvi_lta_r.time:
            month = time_val.dt.month.item()  # Extract month from tuple
            year = time_val.dt.year.item()  # Extract year from tuple
            # Extract a single timestep as a DataArray
            image = m_ndvi_lta_r.sel(time=time_val)  # Squeeze to remove unnecessary dimensions
            # create a file path for the geotiff file
            filename = f'{output_dir}/ndvi_a_{year}_{month:02d}.tif'

            # write the data to a geotiff file
            image.rio.to_raster(filename, driver='GTiff')
            print(f"{filename} saved successfully")
        return

    def process_CHIRPS(self):
        """
        Fetches CHIRPS data from HDC
        return:
        A xarray object with the timeseries and the exported tif files
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

        chirps = self.hdc_stac_client.search(bbox=self.bbox,
            # collections=['mod13q1_vim_native'],
            collections=['rfh_dekad'],
            datetime= period,# '2022-01-01/2022-12-31'
        ).get_all_items()

        res= 0.0022457882102988  #50.045454545Km #0.0022457882102988 # 250 or 0.01 for 1km

        chirps_stac = stac_load(chirps, output_crs='EPSG:4326', resolution=res, patch_url=self.signer, bbox=self.bbox)
        #chirps_stac

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

        chirps_stac = xr.where(chirps_stac == -9999, np.nan, chirps_stac)
        #CHIRPS_stac

        # Group CHIRPS data by month
        chirps_m = chirps_stac.resample(time='1M').sum()
        #CHIRPS_m
        # Creates the seasonal monthly mean precipitation
        chirps_m_s = chirps_m.mean(dim='time')
        #CHIRPS_m_s

        output_dir_s = f'C:/Geotar/{self.pilot_name}/geodata/Processed/precipitation/season'

        filename_m_s = f'{output_dir_s}/rain_m_s.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        chirps_m_s.rio.to_raster(filename_m_s, driver='GTiff')
        print(f'{filename_m_s} saved successfully')

        # Export CHIRPS Rainfall monthly data

        # create a directory to store the geotiff files
        output_dir = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Precipitation'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # loop over all timesteps in the dataset
        for time_val in chirps_m.time:
            month = time_val.dt.month.item()
            year = time_val.dt.year.item()
            # extract a single timestep as a DataArray
            da = chirps_m.sel(time=time_val)
            # Set NaN as NoData
            # print("The no data value is:",da.rio.nodata)

            # create a file path for the geotiff file
            filename = f'{output_dir}/rain_{year}_{month:02d}.tif'

            # write the data to a geotiff file
            da.rio.to_raster(filename)
            print(f'{filename} saved successfully')

        # Creates the seasonal monthly sum precipitation
        chirps_s_s = chirps_m.sum(dim='time')
        #CHIRPS_s_s

        output_dir_s = f'C:/Geotar/{self.pilot_name}/geodata/Processed/precipitation/season'
        filename_s_s = f'{output_dir_s}/rain_s_s_{year}.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        chirps_s_s.rio.to_raster(filename_s_s, driver='GTiff')
        print(f'{filename_s_s} saved successfully')
        return

    def process_CHIRPS_Anomaly(self):
        """
        Fetches CHIRPS anomaly data from HDC
        Returns:
        A xarray object with the CHIRPS anomaly timeseries and the exported tif files
        """

        print("Processing CHIRPS anomalies")

        start_date, end_date = self.period.split('/')
        start_date = datetime.datetime.strptime(start_date, "%Y-%m-%d")
        end_date = datetime.datetime.strptime(end_date, "%Y-%m-%d")

        # Subtract one month from the start and end dates while preserving the day of the month
        modified_start_date = start_date - relativedelta(months=1)
        modified_end_date = end_date - relativedelta(months=1)

        # Format the modified dates back to the desired string format
        modified_period = modified_start_date.strftime("%Y-%m-%d") + '/' + modified_end_date.strftime("%Y-%m-%d")

        chirps_an = self.hdc_stac_client.search(bbox=self.bbox,
            #collections=['mod13q1_vim_native'],
            collections=['rfq_dekad'],
            datetime=modified_period #'2022-01-01/2022-12-31',
        ).get_all_items()
        #print(stac_items)

        res = 0.0022457882102988
        crs = rasterio.crs.CRS.from_epsg(4326)
        chirps_an_stac = stac_load(chirps_an, output_crs='EPSG:4326', resolution= res,
                                   patch_url=self.signer, bbox=self.bbox)
        #CHIRPS_an_stac

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
        chirps_an_stac = xr.where(chirps_an_stac == -9999, np.nan, chirps_an_stac)
        #CHIRPS_an_stac

        # Aggregate the dekadal Chirps anomaly data by month
        chirps_an_m = chirps_an_stac.resample(time='1M').mean()
        #rescale to have values from 0 to 1
        chirps_an_m = chirps_an_m/100
        #CHIRPS_an_m

        # Aggregate mean anomaly data by season

        chirps_an_s = chirps_an_m.mean(dim='time')
        #CHIRPS_an_s

        image_an = chirps_an_s
        image_an.rio.set_crs(crs)
        output_dir_s = f'C:/Geotar/{self.pilot_name}/geodata/Processed/precipitation/season'
        filename_s = f'{output_dir_s}/rain_an_m.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        image_an.rio.to_raster(filename_s, driver='GTiff', crs='EPSG:4326')
        print(f'{filename_s} saved successfully')

        # Agreggate anomaly max data season
        chirps_an_s_max = chirps_an_m.max(dim='time')
        #CHIRPS_an_s_max
        image_an_m = chirps_an_s_max
        image_an_m.rio.set_crs(crs)


        # Export the monthly anomaly data to GeoTiff files

        # create a directory to store the geotiff files
        output_dir = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Precipitation'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # loop over all timesteps in the dataset
        for time_val in chirps_an_m.time:
            month = time_val.dt.month.item()
            year = time_val.dt.year.item()
            da = chirps_an_m.sel(time=time_val)
            da.rio.set_crs(crs)

            # create a file path for the geotiff file
            filename = f'{output_dir}/rain_a_{year}_{month:02d}.tif'

            # write the data to a geotiff file
            da.rio.to_raster(filename, driver='GTiff', crs='EPSG:4326')
            print(f'{filename} saved successfully')

        filename_s = f'{output_dir_s}/rain_an_ma_{year}.tif'
        # write the data to a geotiff file
        image_an_m.rio.to_raster(filename_s, driver='GTiff', crs='EPSG:4326')
        print(f'{filename_s} saved successfully')
        return

# # # Land Surface temperature processing
#
    # Need to readjust the dates of the period to have same dates as the NDVI
    def process_LST(self):
        '''
        Fetches data Land surface temperature data from HDC
        Returns:
        A xarray object with the lST timeseries and the exported tif files
        '''

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

        lst_query = self.hdc_stac_client.search(bbox=self.bbox,
            # collections=['mod13q1_vim_native'],
            collections=['myd11a2_txa_dekad'],
            datetime= self.period #'2022-01-01/2022-12-31'
                                        ).get_all_items()
        res = 0.0022457882102988 # 250 or 0.01 for 1km
        lst = stac_load(lst_query, output_crs='EPSG:4326', resolution=res, patch_url=self.signer, bbox=self.bbox)
        #lst

        # output_dir_zarr = f'C:/Geotar/{pilot_name}/geodata/zarr'
        # outfile = output_dir_zarr + '/lst_stac.zarr'
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
        # lst.to_zarr(outfile)
        # print(f'{outfile} saved')
        # output_dir_zarr = f'C:/Geotar/{self.pilot_name}/geodata/zarr'
        # outfile = output_dir_zarr + '/lst_an_stac.zarr'

        # # create a directory to store the geotiff files
        # if not os.path.exists(output_dir_zarr):
        #     os.makedirs(output_dir_zarr)
        #
        # # Delete the existing Zarr file if it exists
        # if os.path.exists(outfile):
        #     shutil.rmtree(outfile)
        #
        # lst_anom.to_zarr(outfile)
        # print(f'{outfile} saved')

        # # Load xarray from zarr file (optional)
        #
        # lst_anom= xr.open_zarr(outfile)
        # lst_anom

        # group the lST data by month
        lst_m = lst.drop('tna')
        lst_m = lst_m.drop('spatial_ref')
        lst_m = lst_m.resample(time='1M').mean()
        lst_m = (lst_m * 0.02) - 273.15
        #lst_m
        lst_s = lst_m.mean(dim='time')
        #lst_s
        image_m = lst_s # Rescaling applied here
        output_dir_s = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Temperature/season'

        # Export the max seasonal temperature
        lst_s_max = lst_m.max(dim='time')
        #lst_s_max
        image_max = lst_s_max # Rescaling applied here
        output_dir_s_max = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Temperature/season'

        # create a directory to store the geotiff files
        output_dir = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Temperature'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # loop over all timesteps in the dataset
        for time_val in lst_m.time:
            month = time_val.dt.month.item()
            year = time_val.dt.year.item()
            # extract a single timestep as a DataArray
            da = lst_m.sel(time=time_val)

            # create a file path for the geotiff file
            filename = f'{output_dir}/LST_{year}_{month:02d}.tif'
            # write the data to a geotiff file
            da.tda.rio.to_raster(filename, driver='GTiff')
            print(f'{filename} saved successfully')

        filename_s_max = f'{output_dir_s_max}/LST_ma_{year}.tif'
        if not os.path.exists(output_dir_s_max):
            os.makedirs(output_dir_s_max)

        # write the data to a geotiff file
        image_max.rio.to_raster(filename_s_max, driver='GTiff')
        print(f'{filename_s_max} saved successfully')

        filename_s = f'{output_dir_s}/LST_m_{year}.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        image_m.rio.to_raster(filename_s, driver='GTiff')
        print(f'{filename_s} saved successfully')

        return

    def process_LST_anomaly(self):
        """
        Fetches Lands surface temperature anomaly data from HDC
        Returns:
        A xarray object with the LST anomaly timeseries and the exported tif files
        """
        print("Processing Land surface temperature anomalies")
        lst_anom_query = self.hdc_stac_client.search(bbox=self.bbox,
            # collections=['mod13q1_vim_native'],
            collections=['myd11a2_txd_dekad'],
            datetime= self.period #'2022-01-01/2022-12-31'
                                        ).get_all_items()
        res = 0.0022457882102988 # 250 or 0.01 for 1km
        lst_anom = stac_load(lst_anom_query, output_crs='EPSG:4326',
                             resolution=res, patch_url=self.signer, bbox=self.bbox)
        #lst_anom

    # group the lST anomaly data by month
        lst_an_m = lst_anom.drop('tnd')
        lst_an_m = lst_an_m.drop('spatial_ref')
        lst_an_m = lst_an_m.resample(time='1M').mean()
        lst_an_m = lst_an_m*0.02
        #lst_an_m

        lst_an_s = lst_an_m.mean(dim='time')

        image_an = lst_an_s / 100  # Rescaling applied here
        output_dir_s = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Temperature/season'
        filename_s = f'{output_dir_s}/LST_an_m.tif'
        if not os.path.exists(output_dir_s):
            os.makedirs(output_dir_s)

        # write the data to a geotiff file
        image_an.rio.to_raster(filename_s, driver='GTiff')
        print(f'{filename_s} saved successfully')
        # lst_an_s

        lst_an_s_max = lst_an_m.max(dim='time')



        # Export the monthly LST anomaly data to GeoTiff files

        # create a directory to store the geotiff files
        output_dir = f'C:/Geotar/{self.pilot_name}/geodata/Processed/Temperature'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # loop over all timesteps in the dataset
        for time_val in lst_an_m.time:
            month = time_val.dt.month.item()
            year = time_val.dt.year.item()
            # extract a single timestep as a DataArray
            da = lst_an_m.sel(time=time_val)
            # create a file path for the geotiff file
            filename = f'{output_dir}/LST_an_{year}_{month:02d}.tif'

            # write the data to a geotiff file
            da.rio.to_raster(filename, driver='GTiff')
            print(f'{filename} saved successfully')

            # lst_an_s_max
            image_an_max = lst_an_s_max / 100  # Rescaling applied here
            filename_s_max = f'{output_dir_s}/LST_an_ma_{year}.tif'
            # write the data to a geotiff file
            image_an_max.rio.to_raster(filename_s_max, driver='GTiff')
            print(f'{filename_s_max} saved successfully')
        return

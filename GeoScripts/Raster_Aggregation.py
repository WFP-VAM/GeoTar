import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
import os
from rasterio.features import geometry_mask
import boto3
from rasterio.io import MemoryFile
from rasterio.windows import from_bounds
from rasterio.crs import CRS
import warnings


def aggregate_tiffs_mean(boundary_file, s3_dir, bucket_name):
    # Prepare the point files by creating a buffer
    boundaries = gpd.read_file(boundary_file)

    # Select the polygons in poly_shp that are completely contained in mask_shp
    bbox = boundaries.total_bounds

    # Define the center of the map
    center = [(bbox[1] + bbox[3]) / 2, (bbox[0] + bbox[2]) / 2]

    s3 = boto3.client('s3', region_name='eu-central-1')

    # List files in the specified S3 directory
    response = s3.list_objects_v2(Bucket=bucket_name, Prefix=s3_dir)
    # create an empty dictionary to hold the raster data
    rasters = {}

    if 'Contents' not in response:
        print(f"No files found in the directory {s3_dir}.")

    for obj in response['Contents']:
        file_key = obj['Key']
        filename = os.path.basename(file_key)

        # Check if the file is a GeoTIFF
        if file_key.endswith('.tif'):
            # Read the raster data from S3 into memory
            s3_object = s3.get_object(Bucket=bucket_name, Key=file_key)
            data_bytes = s3_object['Body'].read()

            # Use MemoryFile to open the raster data from in-memory bytes
            with MemoryFile(data_bytes) as memfile:
                with memfile.open() as src:
                    data = src.read(1)
            rasters[filename] = data
    # create a list to hold the mean values for each raster
    mean_values = []

    multicol_df = pd.DataFrame()

    # loop over each raster in the rasters dictionary
    for name, data in rasters.items():
        print(f'Processing: {name}')

        raster = src
        # create an empty list to hold the mean values for this raster
        raster_mean_values = []

        # loop over each polygon in the GeoDataFrame
        for i, poly in enumerate(boundaries['geometry']):
            # create a mask from the polygon geometry
            mask = geometry_mask([poly], out_shape=data.shape, transform=raster.transform, invert=True)

            # read in the raster data using the mask
            masked_data = data[mask]

            # get the nodata value for the raster
            # nodata = np.isnan(masked_data)
            # print(nodata)
            # Check if masked_data has nodata values (assuming nodata is represented by a specific value like -9999)
            # if np.any(nodata > 0):  # Change -9999 to your actual nodata value
            #    print(f"Raster '{name}' has nodata values in polygon {i}")
            # To avoid No data values in the resulting table this part catches warnings and creates a buffer on the polygon
            # the buffer created increases the polygon area in a third of the pixel size that way the center of the pixel
            # is recorded
            with warnings.catch_warnings():
                warnings.filterwarnings("error")  # Raise warnings as errors

                try:
                    # print(f"{buffered_points.loc[i,'ADM3_ES']}")
                    mean_value = np.nanmean(masked_data)
                    # print("Mean Value:", mean_value)
                except Warning as w:
                    # A warning occurred when calculating the mean
                    # print(f"Warning: {w}")
                    poly = poly.buffer(0.0022457882102988 * 4)
                    mask = geometry_mask([poly], out_shape=data.shape, transform=raster.transform, invert=True)
                    masked_data = data[mask]
                    mean_value = np.nanmean(masked_data)
                except Warning as w:
                    mean_value = np.nan
                    # print("Mean Value:", mean_value)
            # calculate the mean value of the masked raster data
            # mean_value = np.mean(masked_data)
            # add the mean value to the list for this raster
            raster_mean_values.append(mean_value)
            # append the raster_mean_values to mean_values list
        mean_values.append(raster_mean_values)

        # Create a DataFrame to hold the mean values for each raster
        result_df = pd.DataFrame()
        # print the mean value for this polygon
        # print(f"mean value of {name} for {buffered_points.loc[i,'event_id_c']} polygon: {mean_value} \n")
        result_df[name] = raster_mean_values
        multicol_df = pd.concat([multicol_df, result_df], axis=1, join='outer')
    return multicol_df

def aggregate_floods(boundary_file, s3_dir, bucket_name):
    s3 = boto3.client('s3', region_name='eu-central-1')

    # List files in the specified S3 directory
    response = s3.list_objects_v2(Bucket=bucket_name, Prefix=s3_dir)
    # create an empty dictionary to hold the raster data
    rasters_floods = {}

    # Prepare the point files by creating a buffer
    boundaries = gpd.read_file(boundary_file)

    # Select the polygons in poly_shp that are completely contained in mask_shp
    bbox = boundaries.total_bounds

    if 'Contents' not in response:
        print(f"No files found in the directory {s3_dir}.")

    for obj in response['Contents']:
        file_key = obj['Key']
        filename = os.path.basename(file_key)

        # Check if the file is a GeoTIFF
        if file_key.endswith('.tif'):
            # Read the raster data from S3 into memory
            s3_object = s3.get_object(Bucket=bucket_name, Key=file_key)
            data_bytes = s3_object['Body'].read()

            # Use MemoryFile to open the raster data from in-memory bytes
            with MemoryFile(data_bytes) as memfile:
                with memfile.open() as src:
                    data_floods = src.read(1)
            rasters_floods[filename] = data_floods
    # create a list to hold the mean values for each raster
    mean_values = []

    flood_df = pd.DataFrame()

    raster_floods = src

    # loop over each raster in the rasters dictionary
    for name_floods, data_floods in rasters_floods.items():
        print(f'processing {name_floods}')
        # create an empty list to hold the mean values for this raster
    raster_mean_values = []

    # loop over each polygon in the GeoDataFrame
    for i, poly_floods in enumerate(boundaries['geometry']):
        # create a mask from the polygon geometry
        mask_floods = geometry_mask([poly_floods], out_shape=data_floods.shape, transform=raster_floods.transform,
                                    invert=True)

        # read in the raster data using the mask
        masked_data_floods = data_floods[mask_floods]
        total_pixels = masked_data_floods.size

        with warnings.catch_warnings():
            warnings.filterwarnings("error")  # Raise warnings as errors

            try:
                mean_value_floods = np.sum(masked_data_floods > 0) / total_pixels
                # print the mean value for each polygon with name for debugging
                # print(f"mean value of {name_floods} for {i}) {buffered_points.loc[i,'ADM3_ES']} polygon: {mean_value_floods} \n")
            except Exception as e:
                # A warning occurred when calculating the mean
                # print(f"Warning: {e}")
                print(f'buffering {i}')
                try:
                    poly_floods = poly_floods.buffer(0.0022457882102988)
                    mask_floods = geometry_mask([poly_floods], out_shape=data_floods.shape,
                                                transform=raster_floods.transform, invert=True)
                    masked_data_floods = data_floods[mask_floods]
                    mean_value_floods = np.sum(masked_data_floods > 0) / total_pixels
                except Exception as e2:
                    mean_value_floods = 0
                    print(f'value 0 assigned {i}')
        raster_mean_values.append(mean_value_floods)

        # Create a DataFrame to hold the mean values for each raster
    result_df = pd.DataFrame()
    # print the mean value for this polygon
    # print(f"mean value of {name} for {buffered_points.loc[i,'event_id_c']} polygon: {mean_value} \n")
    result_df[name_floods] = raster_mean_values
    flood_df = pd.concat([flood_df, result_df], axis=1, join='outer')
    return flood_df








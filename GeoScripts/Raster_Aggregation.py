import geopandas as gpd
import pandas as pd
import numpy as np
import os
from rasterio.features import geometry_mask
import boto3
from rasterio.io import MemoryFile
import warnings

def select_unique_column_index(df):
    """
    This function iterates through the columns of the given DataFrame and prints only those with unique values.
    It then prompts the user to select a column index from the printed indices.

    Args:
        df (pd.DataFrame): The DataFrame to process.

    Returns:
        int: The selected column index with unique values.
    """
    # Collect valid indices and column names
    valid_indices = []

    # Iterate through the columns and print only those with unique values
    for idx, column_name in enumerate(df.columns):
        if df[column_name].is_unique:
            print(f"{idx}: {column_name}")
            valid_indices.append((idx, column_name))

    # Ensure valid_indices is not empty
    if not valid_indices:
        print("the selected boundary file does not have a unique column, try another boundary file.")
        return None
    else:
        print("Type the column index number to track the admin name processed:")
        while True:
            try:
                col_index = int(input())
                if any(col_index == idx for idx, _ in valid_indices):
                    selected_column = next(name for idx, name in valid_indices if idx == col_index)
                    return col_index
                else:
                    print(f"Invalid index. Please choose from {[idx for idx, _ in valid_indices]}.")
            except ValueError:
                print("Invalid input. Please enter a number.")
def aggregate_tiffs_mean(boundary_file, s3_dir, bucket_name):
    # Prepare the point files by creating a buffer
    boundaries = gpd.read_file(boundary_file)
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
        base_name = os.path.splitext(filename)[0]

        # Check if the file is a GeoTIFF
        if file_key.endswith('.tif'):
            # Read the raster data from S3 into memory
            s3_object = s3.get_object(Bucket=bucket_name, Key=file_key)
            data_bytes = s3_object['Body'].read()

            # Use MemoryFile to open the raster data from in-memory bytes
            with MemoryFile(data_bytes) as memfile:
                with memfile.open() as src:
                    data = src.read(1)
            rasters[base_name] = data
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
            with warnings.catch_warnings():
                warnings.filterwarnings("error")  # Raise warnings as errors

                try:
                    # print(f"{buffered_points.loc[i,'ADM3_ES']}")
                    mean_value = np.nanmean(masked_data)
                    # print("Mean Value:", mean_value)
                except Warning:
                    # A warning occurred when calculating the mean
                    # print(f"Warning: {w}")
                    poly = poly.buffer(0.0022457882102988 * 4)
                    mask = geometry_mask([poly], out_shape=data.shape, transform=raster.transform, invert=True)
                    masked_data = data[mask]
                    mean_value = np.nanmean(masked_data)
                except Warning:
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

    if 'Contents' not in response:
        print(f"No files found in the directory {s3_dir}.")

    for obj in response['Contents']:
        file_key = obj['Key']
        filename = os.path.basename(file_key)
        base_name = os.path.splitext(filename)[0]

        # Check if the file is a GeoTIFF
        if file_key.endswith('.tif'):
            # Read the raster data from S3 into memory
            s3_object = s3.get_object(Bucket=bucket_name, Key=file_key)
            data_bytes = s3_object['Body'].read()

            # Use MemoryFile to open the raster data from in-memory bytes
            with MemoryFile(data_bytes) as memfile:
                with memfile.open() as src:
                    data_floods = src.read(1)
            rasters_floods[base_name] = data_floods
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
def aggregate_fatalities(boundary_file, pilot_name, column_name):
    boundaries = gpd.read_file(boundary_file)
    # Load the point dataset
    df_period = gpd.read_file(fr's3://geotar.s3.hq/Geotar/{pilot_name}/geodata/Raw/Conflict/filtered_acled.geojson')

    # Perform a spatial join to associate points with polygons
    spatial_join = gpd.sjoin(df_period, boundaries, how='inner', predicate='intersects')
    # spatial_join.head()
    # print(spatial_join.columns)
    # Group by the 'ADM3_PCODE' and calculate the sum of fatalities
    aggregated_data = spatial_join.groupby(column_name)['fatalities'].sum().reset_index()
    # print(aggregated_data)
    # # Merge the aggregated data back into the boundary polygon dataset
    boundaries = boundaries.merge(aggregated_data, on=column_name, how='left')
    boundaries['fatalities'] = boundaries['fatalities'].fillna(0)

    # Replace NaN values with 0 in the "fatalities" column for polygons with no points
    result_df = boundaries[['fatalities']]
    return result_df
def aggregate_nightlights(boundary_file, s3_ntl_dir, bucket_name):
    # Prepare the point files by creating a buffer
    boundaries = gpd.read_file(boundary_file)
    s3 = boto3.client('s3', region_name='eu-central-1')

    # List files in the specified S3 directory
    response = s3.list_objects_v2(Bucket=bucket_name, Prefix=s3_ntl_dir)
    # create an empty dictionary to hold the raster data
    rasters = {}

    if 'Contents' not in response:
        print(f"No files found in the directory {s3_ntl_dir}.")

    for obj in response['Contents']:
        file_key = obj['Key']
        filename = os.path.basename(file_key)
        base_name = os.path.splitext(filename)[0]

        # Check if the file is a GeoTIFF
        if file_key.endswith('.tif'):
            # Read the raster data from S3 into memory
            s3_object = s3.get_object(Bucket=bucket_name, Key=file_key)
            data_bytes = s3_object['Body'].read()

            # Use MemoryFile to open the raster data from in-memory bytes
            with MemoryFile(data_bytes) as memfile:
                with memfile.open() as src:
                    data = src.read(1)
            rasters[base_name] = data
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
            with warnings.catch_warnings():
                warnings.filterwarnings("error")  # Raise warnings as errors

                try:
                    # print(f"{buffered_points.loc[i,'ADM3_ES']}")
                    mean_value = np.max(masked_data)
                    # print("Mean Value:", mean_value)
                except Warning:
                    # A warning occurred when calculating the mean
                    # print(f"Warning: {w}")
                    poly = poly.buffer(0.0022457882102988 * 4)
                    mask = geometry_mask([poly], out_shape=data.shape, transform=raster.transform, invert=True)
                    masked_data = data[mask]
                    mean_value = np.max(masked_data)
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
def consolidate_results(boundary_file, pilot_name, bucket_name):
    # output_file = f'Geotar/{pilot}/geodata/Outputs/{pilot_name}_output_file.geojson'
    output_file = f'{pilot_name}_output_file.geojson'
    # Load the boundary file as a GeoDataFrame
    boundaries = gpd.read_file(boundary_file)
    col_index = select_unique_column_index(boundaries)
    column_name = boundaries.columns[col_index]
    print(f'column {column_name} selected')
    all_results = pd.DataFrame()
    s3_dirs = [
        f'Geotar/{pilot_name}/geodata/Processed/Vegetation/season',
        f'Geotar/{pilot_name}/geodata/Processed/Precipitation/season',
        f'Geotar/{pilot_name}/geodata/Processed/Temperature/season',
        f'Geotar/{pilot_name}/geodata/Processed/250m',
        f'Geotar/{pilot_name}/geodata/Processed/Travel_time'
    ]

    s3_dirs_floods = [
        f'Geotar/{pilot_name}/geodata/Processed/Floods'
    ]

    s3_ntl_dir = [f'Geotar/{pilot_name}/geodata/Processed/NTL']

    for file in s3_dirs:
        print(f"Processing files in: {file}")
        result = aggregate_tiffs_mean(boundary_file, file, bucket_name)

        if result is not None:
            if all_results.empty:
                all_results = result
            else:
                all_results = pd.concat([all_results, result], axis=1, join='outer')

    for file in s3_dirs_floods:
        print(f"Processing files in: {file}")
        result1 = aggregate_floods(boundary_file, file, bucket_name)

        if result1 is not None:
            if all_results.empty:
                all_results = result1
            else:
                all_results = pd.concat([all_results, result1], axis=1, join='outer')

    print(f"Processing files fatalities")
    fatalities = aggregate_fatalities(boundary_file, pilot_name, column_name)
    if fatalities is not None:
        if all_results.empty:
            all_results = fatalities
        else:
            all_results = pd.concat([all_results, fatalities], axis=1, join='outer')

    for file in s3_ntl_dir:
        print(f"Processing files in: {file}")
        ntl_result = aggregate_nightlights(boundary_file, file, bucket_name)

        if ntl_result is not None:
            if all_results.empty:
                all_results = ntl_result
            else:
                all_results = pd.concat([all_results, ntl_result], axis=1, join='outer')

    # Merge the consolidated results with the boundaries GeoDataFrame
    if not all_results.empty:

        # Combine the results DataFrame with the boundaries GeoDataFrame
        all_results_gdf = boundaries.join(all_results)
        all_results_gdf.to_file(f"output_file.geojson")
        # all_results_gdf.to_file(f"s3://geotar.s3.hq/{output_file}")

        # Also save a CSV without the geometry
        csv_output_file = os.path.splitext(output_file)[0] + '.csv'
        csv = all_results_gdf.drop('geometry', axis=1)
        # csv.to_csv(f"s3://geotar.s3.hq/{csv_output_file}", index=False)
        csv.to_csv(f"output_file.csv", index=False)
        print(f"Consolidated results saved to {output_file} and {csv_output_file}")
    else:
        print("No data was processed.")











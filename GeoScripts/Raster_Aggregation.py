#!/usr/bin/env python
# coding: utf-8

# # Processing data from raster to polygon

# This script can be run in the environment AutoGIS


import geopandas as gpd
import rasterio
import numpy as np
import os
from rasterio.features import geometry_mask



# Select the pilot area

# Ask the user to select an option
print("Please select the pilot area:")
print("0. GLOBAL")
print("1. COL")
print("2. CHAD")
print("3. IRAQ Dahuk")
print("4. IRAQ Najaf")
print("5. IRAQ")
print("6. LBN")
print("7. VEN")
print("8. AFG")
pilot_num = input()
assert pilot_num in ["0","1", "2", "3", "4", "5", "6", "7","8"], "Invalid pilot area selected."

if pilot_num == "0":
    pilot_name = "GLOBAL"
    mask = gpd.read_file(f"C:/Geotar/VEN/geodata/Processed/Mask/VEN_mask.shp")
    period = "2021-05-01/2022-01-31"
    polygons_path = f"C:/Geotar/{pilot_name}/geodata/workspace/Test_boundaries.shp"
    print("You selected Global")
elif pilot_num == "1":
    input_shp = "C:/Geotar/COL/geodata/Processed/Education/Education_facilities.shp"
    pilot = "COL"
    mask = gpd.read_file(f"C:/Geotar/{pilot}/geodata/Processed/Mask/{pilot}_mask.shp")
    polygons_path = f"C:/Geotar/{pilot}/geodata/workspace/veredas_sel.shp"
    print("You selected Colombia")
elif pilot_num == "2":
    pilot = "CHAD"
    point_shp = f"C:\Geotar\{pilot}\geodata\workspace\Village_Baseline_FCS_rCSI_06_02.shp"
    mask = gpd.read_file(f"C:/Geotar/{pilot}/geodata/Processed/Mask/{pilot}_mask.shp")
    polygons_path = f"C:/Geotar/{pilot}/geodata/workspace/Village_bound_manual.shp"
    print("You selected CHAD")
elif pilot_num == "3":
    pilot = "IRAQ_D"
    point_shp = ""
    mask = gpd.read_file(f"C:/Geotar/{pilot}/geodata/Processed/Mask/dahuk_mask.shp")
    polygons_path = f"C:/Geotar/{pilot}/geodata/Processed/boundaries/Dahuk_adm3.shp"
    print("You selected IRAQ Dahuk")
elif pilot_num == "4":
    pilot = "IRAQ_N"
    point_shp = ""
    mask = gpd.read_file(f"C:/Geotar/{pilot}/geodata/Processed/Mask/najaf_mask.shp")
    polygons_path = f"C:/Geotar/{pilot}/geodata/Processed/boundaries/Najaf_adm3.shp"
    print("You selected IRAQ Najaf")
elif pilot_num == "5":
    pilot = "IRAQ"
    point_shp = ""
    mask = gpd.read_file(f"C:/Geotar/{pilot}/geodata/Processed/Mask/Iraq_mask.shp")
    polygons_path = f"C:/Geotar/{pilot}/geodata/Processed/boundaries/Iraq_mask.shp"
    print("You selected IRAQ")
elif pilot_num == "6":
    pilot = "LBN"
    point_shp = ""
    mask = gpd.read_file(f"C:/Geotar/{pilot}/geodata/Processed/Mask/{pilot}_mask.shp")
    #polygons_path = f"C:/Geotar/{pilot}/geodata/Raw/boundaries/lbn_bnd_cadastral_unhcr_wfpco.shp"
    polygons_path = 'zip://C:/Geotar/LBN/geodata/Raw/Co_data/lbn_bnd_cadastral_unhcr_wfpco.zip/lbn_bnd_cadastral_unhcr_wfpco.shp'
    print("You selected Lebanon")
elif pilot_num == "7":
    pilot = "VEN"
    point_shp = ""
    mask = gpd.read_file(f"C:/Geotar/{pilot}/geodata/Processed/Mask/{pilot}_mask.shp")
    polygons_path = f"C:/Geotar/{pilot}/geodata/Processed/boundaries/adm3.shp"
    print(f"You selected {pilot}")
elif pilot_num == "8":
    pilot = "AFG"
    point_shp = ""
    mask = gpd.read_file(f"C:/Geotar/{pilot}/geodata/Processed/Mask/{pilot}_mask.shp")
    polygons_path = f"C:/Geotar/{pilot}/geodata/Processed/boundaries/afg_admbnda_adm2_agcho_20211117.shp"
    print(f"You selected {pilot}")


def Raster_aggregation():
    '''
    Function to aggregate all the indicators at admin polygon level

    Returns:
        shapefile with aggregated values on each polygon of the boundary file
        csv file with aggregated values on each polygon of the boundary file

    '''
    # Prepare the point files by creating a buffer

    buffered_points = gpd.read_file(polygons_path)

    # Select the polygons in poly_shp that are completely contained in mask_shp

    buffered_points.head()

    bbox = buffered_points.total_bounds
    bbox
    # Define the center of the map
    center = [(bbox[1]+bbox[3])/2, (bbox[0]+bbox[2])/2]


    # This part of the script reads the NDVI anomalies data from each geotiff file, then it uses the polygons from the previous cell to compute the mean NDVI anomaly on each of the resulting polygons. Finally it writes the mean NDVI anomaly values computed to a column on the polygon file.

    import warnings

    # set the path to the folder containing the raster files
    path = f'C:/Geotar/{pilot}/geodata/Processed/Vegetation/season'

    # create an empty dictionary to hold the raster data
    rasters = {}

    # loop over each file in the folder
    for file in os.listdir(path):

        # check if the file is a GeoTIFF
        if file.endswith('.tif'):

            # read in the raster data
            with rasterio.open(os.path.join(path, file)) as src:
                data = src.read(1)
                #print("CRS of the raster:", src.crs)

            # add the raster data to the dictionary
            rasters[file] = data


    # create a list to hold the mean values for each raster
    mean_values = []
    raster = src

    # loop over each raster in the rasters dictionary
    for name, data in rasters.items():
        print (f'Processing: {name}')
        # create an empty list to hold the mean values for this raster
        raster_mean_values = []

        # loop over each polygon in the GeoDataFrame
        for i, poly in enumerate (buffered_points['geometry']):
            # create a mask from the polygon geometry
            mask = geometry_mask([poly], out_shape=data.shape, transform=raster.transform, invert=True)

            # read in the raster data using the mask
            masked_data = data[mask]

            # get the nodata value for the raster
            nodata = np.isnan(masked_data)
            #print(nodata)
            # Check if masked_data has nodata values (assuming nodata is represented by a specific value like -9999)
            #if np.any(nodata > 0):  # Change -9999 to your actual nodata value
            #    print(f"Raster '{name}' has nodata values in polygon {i}")
            #To avoid No data values in the resulting table this part catches warnings and creates a buffer on the polygon
            #the buffer created increases the polygon area in a third of the pixel size that way the center of the pixel
            #is recorded
            with warnings.catch_warnings():
                warnings.filterwarnings("error")  # Raise warnings as errors

                try:
                    #print(f"{buffered_points.loc[i,'ADM3_ES']}")
                    mean_value = np.nanmean(masked_data)
                    #print("Mean Value:", mean_value)
                except Warning as w:
                    # A warning occurred when calculating the mean
                    #print(f"Warning: {w}")
                    poly = poly.buffer(0.0022457882102988/3)
                    mask = geometry_mask([poly], out_shape=data.shape, transform=raster.transform, invert=True)
                    masked_data = data[mask]
                    mean_value = np.nanmean(masked_data)
                except Warning as w:
                    mean_value =np.nan
                    #print("Mean Value:", mean_value)
            #calculate the mean value of the masked raster data
            # mean_value = np.mean(masked_data)
            # add the mean value to the list for this raster
            raster_mean_values.append(mean_value)
            # append the raster_mean_values to mean_values list
        mean_values.append(raster_mean_values)
            # print the mean value for this polygon
            #print(f"mean value of {name} for {buffered_points.loc[i,'event_id_c']} polygon: {mean_value} \n")
    buffered_points_copy = buffered_points.copy()
    # loop over the mean values for each raster and add them to the GeoDataFrame as new columns
    for i, (name, data) in enumerate(rasters.items()):
        name_fix = name[:-4]
        #name_fix = name_fix.replace("chirps_anomaly", "rain")
        #name_fix = name_fix.replace("anomaly", "a")
        column_name = f'{name_fix}'
        #print(column_name)
    #    buffered_points[column_name] = mean_values[i]
    #    buffered_points.loc[:, column_name] = mean_values[i]
        buffered_points_copy[column_name] = mean_values[i]


    buffered_points_copy.head()


    # The next part of the script will execute an identical computation of the mean CHIRPS

    # set the path to the folder containing the raster files
    path_p = f'C:/Geotar/{pilot}/geodata/Processed/Precipitation/season'

    # create an empty dictionary to hold the raster data
    rasters_p = {}

    # loop over each file in the folder
    for file_p in os.listdir(path_p):

        # check if the file is a GeoTIFF
        if file_p.endswith('.tif'):

            # read in the raster data
            with rasterio.open(os.path.join(path_p, file_p)) as src_p:
                data_p = src_p.read(1)

            # add the raster data to the dictionary
            rasters_p[file_p] = data_p
            #print(file_p)
    # read in the buffered polygon shapefile
    #buffered_p = gpd.read_file('buffered.shp')

    # create a list to hold the mean values for each raster
    mean_values_p = []
    raster_p = src_p

    # loop over each raster in the rasters dictionary
    for name_p, data_p in rasters_p.items():
        print (f'Processing: {name_p}')

        # create an empty list to hold the mean values for this raster
        raster_mean_values_p = []

        # loop over each polygon in the GeoDataFrame
        for i, poly_p in enumerate (buffered_points['geometry']):
            # create a mask from the polygon geometry
            mask_p = geometry_mask([poly_p], out_shape=data_p.shape, transform=raster_p.transform, invert=True)

            # read in the raster data using the mask
            masked_data_p = data_p[mask_p]
            #print(masked_data_p)
            # get the nodata value for the raster
            nodata_p = raster_p.nodata

            # check if any values in the masked array are nodata
            #if nodata_p is not None:
            #    masked_data_p = masked_data_p[masked_data_p != nodata_p]
            with warnings.catch_warnings():
                warnings.filterwarnings("error")  # Raise warnings as errors

                try:
                    mean_value_p = np.nanmean(masked_data_p)
                except Warning as w:
                    # A warning occurred when calculating the mean
                    #print(f"Warning: {w}")
                    poly_p = poly_p.buffer(0.0022457882102988/3)
                    mask_p = geometry_mask([poly_p], out_shape=data_p.shape, transform=raster_p.transform, invert=True)
                    masked_data_p = data_p[mask_p]
                    mean_value_p = np.nanmean(masked_data_p)
                except Warning as w:
                    mean_value =np.nan
            # calculate the mean value of the masked raster data
            #mean_value_p = np.mean(masked_data_p)
            # add the mean value to the list for this raster
            raster_mean_values_p.append(mean_value_p)
            # append the raster_mean_values to mean_values list
        mean_values_p.append(raster_mean_values_p)
            # print the mean value for this polygon
            #print(f"mean value of {name_p} for {buffered_points_p.loc[i_p,'event_id_c']} polygon: {mean_value_p} \n")

    # loop over the mean values for each raster and add them to the GeoDataFrame as new columns
    for i, (name_p, data_p) in enumerate(rasters_p.items()):
        name_fix_p = name_p[:-4]
        #name_fix_p = name_fix_p.replace("chirps_anomaly", "rain")
        #name_fix_p = name_fix_p.replace("anomaly", "an")
        column_name_p = f'{name_fix_p}'
        #print(column_name_p)
        buffered_points_copy[column_name_p] = mean_values_p[i]
        #buffered_point.loc[:, column_name_p] = mean_values_p[i]

    buffered_points_copy.head()


    # set the path to the folder containing the raster files
    path_lst = f'C:/Geotar/{pilot}/geodata/Processed/Temperature/season'

    # create an empty dictionary to hold the raster data
    rasters_lst = {}

    # loop over each file in the folder
    for file_lst in os.listdir(path_lst):

        # check if the file is a GeoTIFF
        if file_lst.endswith('.tif'):

            # read in the raster data
            with rasterio.open(os.path.join(path_lst, file_lst)) as src_lst:
                data_lst = src_lst.read(1)

            # add the raster data to the dictionary
            rasters_lst[file_lst] = data_lst
            #print(file_lst)
    # read in the buffered polygon shapefile
    #buffered_lst = gpd.read_file('buffered.shp')

    # create a list to hold the mean values for each raster
    mean_values_lst = []
    raster_lst = src_lst

    # loop over each raster in the rasters dictionary
    for name_lst, data_lst in rasters_lst.items():
        print(f'processing: {name_lst}')

        # create an empty list to hold the mean values for this raster
        raster_mean_values_lst = []

        # loop over each polygon in the GeoDataFrame
        for i, poly_lst in enumerate (buffered_points['geometry']):
            # create a mask from the polygon geometry
            mask_lst = geometry_mask([poly_lst], out_shape=data_lst.shape, transform=raster_lst.transform, invert=True)

            # read in the raster data using the mask
            masked_data_lst = data_lst[mask_lst]

            # get the nodata value for the raster
            nodata_lst = raster_lst.nodata

            # check if any values in the masked array are nodata
            #if nodata_lst is not None:
            #    masked_data_lst = masked_data_lst[masked_data_lst != nodata_lst]

            # calculate the mean value of the masked raster data
            #mean_value_lst = np.mean(masked_data_lst)
            with warnings.catch_warnings():
                warnings.filterwarnings("error")  # Raise warnings as errors

                try:
                    mean_value_lst = np.nanmean(masked_data_lst)
                except Warning as w:
                    # A warning occurred when calculating the mean
                    #print(f"Warning: {w}")
                    poly_lst = poly_lst.buffer(0.0022457882102988/3)
                    mask_lst = geometry_mask([poly_lst], out_shape=data.shape, transform=raster.transform, invert=True)
                    masked_data_lst = data_lst[mask_lst]
                    mean_value_lst = np.nanmean(masked_data_lst)
                except Warning as w:
                    mean_value =np.nan

            # add the mean value to the list for this raster
            raster_mean_values_lst.append(mean_value_lst)

        # append the raster_mean_values to mean_values list
        mean_values_lst.append(raster_mean_values_lst)

        # print the mean value for this polygon
        #print(f"mean value of {name_lst} for {i} polygon: {mean_value_lst} \n")

    # loop over the mean values for each raster and add them to the GeoDataFrame as new columns
    for i, (name_lst, data_lst) in enumerate(rasters_lst.items()):
        name_fix_lst = name_lst[:-4]
        #name_fix_lst = name_fix_lst.replace("anomaly", "an")
        column_name_lst = f'{name_fix_lst}'
        #print(column_name_lst)
        buffered_points_copy[column_name_lst] = mean_values_lst[i]

    buffered_points_copy.head()


    # Extract the values of Nighttime lights

    # set the path to the folder containing the raster files
    path_ntl = f'C:/Geotar/{pilot}/geodata/Processed/NTL'

    # create an empty dictionary to hold the raster data
    rasters_ntl = {}

    # loop over each file in the folder
    for file_ntl in os.listdir(path_ntl):

        # check if the file is a GeoTIFF
        if file_ntl.endswith('.tif'):
            #print(file_ntl)
            # read in the raster data
            with rasterio.open(os.path.join(path_ntl, file_ntl)) as src_ntl:
                data_ntl = src_ntl.read(1)

            # add the raster data to the dictionary
            rasters_ntl[file_ntl] = data_ntl
            #print(data_ntl)
        else:
            print('Tiff files not available')
    # read in the buffered polygon shapefile
    #buffered_ntl = gpd.read_file('buffered.shp')

    # create a list to hold the mean values for each raster
    max_values_ntl = []
    raster_ntl = src_ntl

    # loop over each raster in the rasters dictionary
    for name_ntl, data_ntl in rasters_ntl.items():
        print(f'processing {name_ntl}')
        # create an empty list to hold the mean values for this raster
        raster_max_values_ntl = []

        # loop over each polygon in the GeoDataFrame
        for i, poly_ntl in enumerate (buffered_points['geometry']):
            # create a mask from the polygon geometry
            mask_ntl = geometry_mask([poly_ntl], out_shape=data_ntl.shape, transform=raster_ntl.transform, invert=True)

            # read in the raster data using the mask
            masked_data_ntl = data_ntl[mask_ntl]

            # get the nodata value for the raster
            #nodata_ntl = raster_ntl.nodata
            #print("NTL no data values:", nodata_ntl)
            # check if any values in the masked array are nodata
            #if nodata_ntl is not None:
            #    masked_data_ntl = masked_data_ntl[masked_data_ntl != nodata_ntl]

            # calculate the mean value of the masked raster data
            #mean_value_ntl = np.max(masked_data_ntl)
            with warnings.catch_warnings():
                warnings.filterwarnings("error")  # Raise warnings as errors

                try:
                    max_value_ntl = np.max(masked_data_ntl)
                except Exception as e:
                    # A warning occurred when calculating the mean
                    #print(f"Warning: {e}")
                    poly_ntl = poly_ntl.buffer(0.0022457882102988/3)
                    mask_ntl = geometry_mask([poly_ntl], out_shape=data_ntl.shape, transform=raster_ntl.transform, invert=True)
                    masked_data_ntl = data_ntl[mask_ntl]
                    max_value_ntl = np.max(masked_data_ntl)
            raster_max_values_ntl.append(max_value_ntl)
        # append the raster_mean_values to mean_values list
        max_values_ntl.append(raster_max_values_ntl)

            # print the mean value for each polygon
            # print(f"mean value of NLT for {i}) {buffered_points.loc[i,'C_Ref_name']} polygon: {mean_value_ntl} \n")

    # loop over the mean values for each raster and add them to the GeoDataFrame as new columns
    for i, (name_ntl, data_ntl) in enumerate(rasters_ntl.items()):
        # name_fix_ntl = name_ntl[:-4]
        column_name_ntl = 'NTL_med_ma'
        #print(column_name_ntl)
        max_values_ntl.append(raster_max_values_ntl)
        buffered_points_copy[column_name_ntl] = max_values_ntl[i]

    buffered_points_copy.head()


    # Extract the mean proximity values for each polygon

    # set the path to the folder containing the raster files
    path_dst = f'C:/Geotar/{pilot}/geodata/Processed/250m'
    single_polygon = buffered_points['geometry'].unary_union

    # create an empty dictionary to hold the raster data
    rasters_dst = {}

    # loop over each file in the folder
    for file_dst in os.listdir(path_dst):

        # check if the file is a GeoTIFF
        if file_dst.endswith('.tif'):

            # read in the raster data
            with rasterio.open(os.path.join(path_dst, file_dst)) as src_dst:
                data_dst = src_dst.read(1)

            # add the raster data to the dictionary
            rasters_dst[file_dst] = data_dst

    # read in the buffered polygon shapefile
    #buffered_dst = gpd.read_file('buffered.shp')

    # create a list to hold the mean values for each raster
    mean_values_dst = []
    raster_dst = src_dst

    # loop over each raster in the rasters dictionary
    for name_dst, data_dst in rasters_dst.items():
        print (f'Processing: {name_dst}')
        #print(f'data dst {data_dst.shape[0]}')
        # create an empty list to hold the mean values for this raster
        raster_mean_values_dst = []

        # loop over each polygon in the GeoDataFrame
        for i, poly_dst in enumerate (buffered_points['geometry']):

            # create a mask from the polygon geometry
            mask_dst = geometry_mask([poly_dst], out_shape=data_dst.shape, transform=raster_dst.transform, invert=True)

            # read in the raster data using the mask
            masked_data_dst = data_dst[mask_dst]
            # get the nodata value for the raster
            #nodata_dst = raster_dst.nodata

            # check if any values in the masked array are nodata
            #if nodata_dst is not None:
    #        masked_data_dst = masked_data_dst[masked_data_dst != nodata_dst]

            # calculate the mean value of the masked raster data
            mean_value_dst = np.nanmean(masked_data_dst)
            with warnings.catch_warnings():
                warnings.filterwarnings("error")  # Raise warnings as errors

                try:
                    mean_value = np.nanmean(masked_data_dst)
                except Warning as w:
                    # A warning occurred when calculating the mean
                    #print(f"Warning: {w}")
                    poly_dst_bu = poly_dst.buffer(0.0022457882102988/3)
                    #print(type(poly_dst_bu))
                    poly_gdf = gpd.GeoDataFrame({'geometry': [poly_dst_bu]})
                    mask_dst = geometry_mask([poly_dst_bu], out_shape=data_dst.shape, transform=raster.transform, invert=True)
                    masked_data_dst = data_dst[mask_dst]
                    mean_value_dst = np.nanmean(masked_data_dst)
                    # print the mean value for this polygon
                    print(f"mean value of {name_dst} for {i} polygon: {mean_value_dst} \n")
                # add the mean value to the list for this raster
                raster_mean_values_dst.append(mean_value_dst)

        # append the raster_mean_values to mean_values list
        mean_values_dst.append(raster_mean_values_dst)



    # loop over the mean values for each raster and add them to the GeoDataFrame as new columns
    for i, (name_dst, data_dst) in enumerate(rasters_dst.items()):
        name_fix_dst = name_dst[5:-4]
        column_name_dst = f'{name_fix_dst}'
        #print(column_name_dst)
        buffered_points_copy[column_name_dst] = mean_values_dst[i]

    buffered_points_copy.head()


    # Compute the number of fatalities for the selected period

    # Load the point dataset
    df_period = gpd.read_file(fr'C:/Geotar/{pilot}/geodata/Raw/conflict/filtered_acled.shp')

    # Load the boundary polygon dataset
    boundary_polygons = buffered_points

    # Perform a spatial join to associate points with polygons
    spatial_join = gpd.sjoin(df_period, boundary_polygons, how='inner', predicate='intersects')

    # # Group by the 'ADM3_PCODE' and calculate the sum of fatalities
    aggregated_data = spatial_join.groupby('ADM3_PCODE')['fatalities'].sum().reset_index()
    #print(aggregated_data)
    # # Merge the aggregated data back into the boundary polygon dataset
    boundary_polygons = boundary_polygons.merge(aggregated_data, on='ADM3_PCODE', how='left')
    boundary_polygons['fatalities']=boundary_polygons['fatalities'].fillna(0)
    #boundary_polygons.head()
    # Replace NaN values with 0 in the "fatalities" column for polygons with no points
    buffered_points_copy['fatalities'] = boundary_polygons['fatalities']
    buffered_points_copy.head()


    # Compute the mean travel time to multiple city sizes

    # set the path to the folder containing the raster files
    path_tt = f'C:/Geotar/GLOBAL/geodata/raw/Travel_time/selected'
    # create an empty dictionary to hold the raster data
    rasters_tt = {}

    # loop over each file in the folder
    print(f"available files: ")
    for file_tt in os.listdir(path_tt):

        # check if the file is a GeoTIFF
        if file_tt.endswith(".tif"):
            print(file_tt)
            # read in the raster data
            with rasterio.open(os.path.join(path_tt, file_tt)) as src_tt:
                data_tt = src_tt.read(1)

            # add the raster data to the dictionary
            rasters_tt[file_tt] = data_tt
            #print(data_tt)
    # read in the buffered polygon shapefile
    #buffered_ntl = gpd.read_file('buffered.shp')

    # create a list to hold the mean values for each raster
    mean_values_tt = []
    raster_tt = src_tt

    # loop over each raster in the rasters dictionary
    for name_tt, data_tt in rasters_tt.items():
        print(f'processing {name_tt}')
        # create an empty list to hold the mean values for this raster
        raster_mean_values_tt = []

        # loop over each polygon in the GeoDataFrame
        for i, poly_tt in enumerate (buffered_points['geometry']):
            # create a mask from the polygon geometry
            mask_tt = geometry_mask([poly_tt], out_shape=data_tt.shape, transform=raster_tt.transform, invert=True)

            # read in the raster data using the mask
            masked_data_tt = data_tt[mask_tt]

            # get the nodata value for the raster
            #nodata_ntl = raster_ntl.nodata
            #print("NTL no data values:", nodata_ntl)
            # check if any values in the masked array are nodata
            #if nodata_ntl is not None:
            #    masked_data_ntl = masked_data_ntl[masked_data_ntl != nodata_ntl]

            # calculate the mean value of the masked raster data
            #mean_value_ntl = np.max(masked_data_ntl)
            with warnings.catch_warnings():
                warnings.filterwarnings("error")  # Raise warnings as errors

                try:
                    mean_value_tt = np.nanmean(masked_data_tt)
                    # print the mean value for each polygon with name for debugging
                    #print(f"mean value of {name_tt} for {i}) {buffered_points.loc[i,'ADM3_ES']} polygon: {mean_value_tt} \n")
                except Exception as e:
                    # A warning occurred when calculating the mean
                    #print(f"Warning: {e}")
                    poly_tt = poly_tt.buffer(0.0022457882102988)
                    mask_tt = geometry_mask([poly_tt], out_shape=data_tt.shape, transform=raster_tt.transform, invert=True)
                    masked_data_tt = data_tt[mask_tt]
                    mean_value_tt = np.nanmean(masked_data_tt)
                    # print the mean value for each polygon with name for debugging
                    #print(f"mean value of Travel time for {i}) {buffered_points.loc[i,'ADM3_ES']} polygon: {mean_value_tt} \n")
            raster_mean_values_tt.append(mean_value_tt)
        # append the raster_mean_values to mean_values list
        mean_values_tt.append(raster_mean_values_tt)
        # print the mean value for each polygon
        #print(f"mean value of Travel time for {i}) {buffered_points.loc[i,'ADM3_ES']} polygon: {mean_value_tt} \n")

    # loop over the mean values for each raster and add them to the GeoDataFrame as new columns
    for i, (name_tt, data_tt) in enumerate(rasters_tt.items()):
        #name_fix_ntl = name_ntl[:-4]
        num= str(i+1)
        column_name_tt = "t_time"+num
        #print(column_name_tt)
        mean_values_tt.append(raster_mean_values_tt)
        buffered_points_copy[column_name_tt] = mean_values_tt[i]


    buffered_points_copy.head()


    import datetime

    current_time = datetime.datetime.now()
    print("Current System Time:", current_time)

    print(buffered_points_copy.columns)


    # Get the centroid of each polygon and extract the latitude and longitude
    buffered_points_copy['lat'] = buffered_points_copy.centroid.y
    buffered_points_copy['lon'] = buffered_points_copy.centroid.x
    buffered_points_copy.head()


    # create a directory to store the geotiff files
    output_path = f'C:/Geotar/{pilot}/geodata/outputs/'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    buffered_points_copy.to_file(output_path + f'{pilot}_output_file.shp')
    print(f'{pilot}_output_file.shp saved in {output_path}' )


    from dbfread import DBF
    import pandas as pd

    # Specify the path to the DBF file
    dbf_file_path = f'C:/Geotar/{pilot}/geodata/outputs/{pilot}_output_file.dbf'

    # Specify the path where you want to save the CSV file
    csv_file_path = f'C:/Geotar/{pilot}/geodata/outputs/{pilot}_output_file.csv'

    # Read the DBF file using dbfread
    dbf = DBF(dbf_file_path, encoding= 'Latin1')

    # Convert the DBF data to a list of dictionaries
    dbf_data = [record for record in dbf]

    # Create a DataFrame from the list of dictionaries
    dbf_df = pd.DataFrame(dbf_data)

    # Save the DataFrame as a CSV file
    dbf_df.to_csv(csv_file_path, index=False)
    return()

Raster_aggregation()


#buffered_points_copy[buffered_points_copy['ADM3_ES'].str.contains('El Rosario', case=False)]


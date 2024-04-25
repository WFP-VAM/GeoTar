#!/usr/bin/env python
# coding: utf-8

# # Proximity rasters processing 

# This is script can be run using geo or cloudserver environment.
# This script uses GDAL to rasterize line, polygon and point datasets.

import geopandas as gpd
from osgeo import gdal, ogr
import os

# Ask the user to select an option
print("Please select the pilot area:")
print("1. COL")
print("2. CHAD")
print("3. IRAQ Dahuk")
print("4. IRAQ Najaf")
print("5. IRAQ")
print("6. LBN")
print("7. VEN")
print("8. AFG")
pilot = input()
assert pilot in ["1", "2", "3", "4", "5","6","7","8"], "Invalid pilot area selected."

print("Please select data to process:")
print("1. Education facilities")
print("2. Healthcare facilities")
print("3. Roads")
p_data = input()
assert p_data in ['1', '2', '3'], "Invalid option selected."

print("Please select the spatial resolution:")
print("1. 1 kilometer")
print("2. 250 meters")
resolution = input()
assert resolution in ["1", "2"], "Invalid resolution selected."

if p_data == "1":
    out_name = "education"
elif p_data == "2":
    out_name = "health"
elif p_data == "3":
    out_name = "roads"
else:
    print("Invalid selection")
    
if resolution == "1":
    res_folder = "1K"
elif resolution == "2":
    res_folder = "250m"
else:
    print("Invalid spatial resolution")


# Check the user"s selection and print a message
if pilot == "1" and p_data == "1":
    pilot_nam = "COL"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Education/Education_facilities.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name + ".tif"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    print(f"You selected {pilot_nam} Education facilities.")
elif pilot == "1" and p_data == "2":
    pilot_nam = "COL"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Raw/Health/Colombia-node.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/COL_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Healthcare facilities.")
elif pilot == "1" and p_data == "3":
    pilot_nam = "COL"
    input_shp = ""
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/COL_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Roads.")
elif pilot == "2" and p_data == "1":
    pilot_nam = "CHAD"
    input_shp = "zip://C:/Geotar/CHAD/geodata/Processed/Education/hotosm_chad_education_facilities_points_shp.zip/hotosm_chad_education_facilities_points.shp"
    mask_shp = "C:/Geotar/CHAD/geodata/Processed/Mask/Chad_mask.shp"
    output = f"C:/Geotar/CHAD/geodata/Processed/{res_folder}/dist"+ out_name+ ".tif"
    print("You selected CHAD Education facilities.")
elif pilot == "2" and p_data == "2":
    pilot_nam = "CHAD"
    input_shp = "C:/Geotar/CHAD/geodata/Processed/Health/Chad_Health_facilities.shp"
    mask_shp = "C:/Geotar/CHAD/geodata/Processed/Mask/Chad_mask.shp"
    output = f"C:/Geotar/CHAD/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print("You selected CHAD Healthcare facilities.")
elif pilot == "2" and p_data == "3":
    pilot_nam = "CHAD"
    input_shp = "zip://C:/Geotar/CHAD/geodata/Raw/Roads/tcd_trs_roads_ocha.zip/tcd_trs_roads_ocha/tcd_trs_roads_ocha.shp"
    mask_shp = "C:/Geotar/CHAD/geodata/Processed/Mask/Chad_mask.shp"
    output = f"C:/Geotar/CHAD/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print("You selected CHAD Roads.")
elif pilot == "3" and p_data == "1":
    pilot_nam = "IRAQ_D"
    input_shp = "zip://C:/Geotar/IRAQ_D/geodata/Raw/Education/hotosm_irq_education_facilities_points_shp.zip/hotosm_irq_education_facilities_points.shp"
    mask_shp = "C:/Geotar/IRAQ_D/geodata/Processed/Mask/Dahuk_mask.shp"
    output = f"C:/Geotar/IRAQ_D/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print("You selected IRAQ Dahuk Education facilities.")
elif pilot == "3" and p_data == "2":
    pilot_nam = "IRAQ_D"
    input_shp = "zip://C:/Geotar/IRAQ_D/geodata/Raw/Health/iraq-shapefiles.zip/shapefiles/healthsites.shp"
    mask_shp = "C:/Geotar/IRAQ_D/geodata/Processed/Mask/Dahuk_mask.shp"
    output = f"C:/Geotar/IRAQ_D/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print("You selected IRAQ Dahuk Healthcare facilities.")
elif pilot == "3" and p_data == "3":
    pilot_nam = "IRAQ_D"
    input_shp = ""
    mask_shp = "C:/Geotar/IRAQ_D/geodata/Processed/Mask/Dahuk_mask.shp"
    output = f"C:/Geotar/IRAQ_D/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print("You selected IRAQ Dahuk Roads.")
elif pilot == "4" and p_data == "1":
    pilot_nam = "IRAQ_N"
    input_shp = "zip://C:/Geotar/IRAQ_N/geodata/Raw/Education/hotosm_irq_education_facilities_points_shp.zip/hotosm_irq_education_facilities_points.shp"
    mask_shp = "C:/Geotar/IRAQ_N/geodata/Processed/Mask/Najaf_mask.shp"
    output = f"C:/Geotar/IRAQ_N/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print("You selected IRAQ Najaf Education facilities.")
elif pilot == "4" and p_data == "2":
    pilot_nam = "IRAQ_N"
    input_shp = "zip://C:/Geotar/IRAQ_N/geodata/Raw/Health/iraq-shapefiles.zip/shapefiles/healthsites.shp"
    mask_shp = "C:/Geotar/IRAQ_N/geodata/Processed/Mask/Najaf_mask.shp"
    output = f"C:/Geotar/IRAQ_N/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print("You selected IRAQ Najaf Healthcare facilities.")
elif pilot == "4" and p_data == "3":
    pilot_nam = "IRAQ_N"
    input_shp = ""
    mask_shp = "C:/Geotar/IRAQ_N/geodata/Processed/Mask/Najaf_mask.shp"
    output = f"C:/Geotar/IRAQ_N/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print("You selected IRAQ Najaf Roads.")
elif pilot == "5" and p_data == "1":
    pilot_nam = "IRAQ"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Education/education.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Education facilities.")
elif pilot == "5" and p_data == "2":
    pilot_nam = "IRAQ"
    input_shp = "zip://C:/Geotar/IRAQ/geodata/Raw/Health/iraq-shapefiles.zip/shapefiles/healthsites.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Healthcare facilities.")
elif pilot == "5" and p_data == "3":
    pilot_nam = "IRAQ"
    input_shp = "zip://C:/Geotar/IRAQ/geodata/Raw/Roads/irq_roads.zip/IRQ_roads.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Roads.")
elif pilot == "6" and p_data == "1":
    pilot_nam = "LBN"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Education/education.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Education facilities.")
elif pilot == "6" and p_data == "2":
    pilot_nam = "LBN"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Health/healthsites.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Healthcare facilities.")
elif pilot == "6" and p_data == "3":
    pilot_nam = "LBN"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Roads/lbn_trs_roads_osm.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Roads.")
elif pilot == "7" and p_data == "1":
    pilot_nam = "VEN"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Education/education.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Education facilities.")
elif pilot == "7" and p_data == "2":
    pilot_nam = "VEN"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Health/healthsites.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} health facilities.")
elif pilot == "7" and p_data == "3":
    pilot_nam = "VEN"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Roads/roads.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Roads.")
elif pilot == "8" and p_data == "1":
    pilot_nam = "AFG"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Education/education.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Education facilities.")
elif pilot == "8" and p_data == "2":
    pilot_nam = "AFG"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Health/healthsites.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Healthcare facilities.")
elif pilot == "8" and p_data == "3":
    pilot_nam = "AFG"
    input_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Roads/roads.shp"
    mask_shp = f"C:/Geotar/{pilot_nam}/geodata/Processed/Mask/{pilot_nam}_mask.shp"
    output = f"C:/Geotar/{pilot_nam}/geodata/Processed/{res_folder}/dist_"+ out_name+ ".tif"
    print(f"You selected {pilot_nam} Roads.")



else:
    print("Invalid selection.")

def proximity_rasters():
    """
    Process vector data into proximity raster data

    Returns:
        Raster proximity files

    """
    print('starting process')
    # Read the GeoPandas object
    gdf = gpd.read_file(input_shp)
    mask = gpd.read_file(mask_shp)
    # Clip the func_gdf using the mask
    gdf = gpd.clip(gdf, mask)

    #clipped shapefile
    clipped_shape = f"C:/Geotar/{pilot_nam}/geodata/workspace/{out_name}.shp"
    #export shapefile
    gdf.to_file(clipped_shape)


    # Define NoData value of new raster
    NoData_value = -9999

    # set the name and location of the output raster file
    dst_filename = f"C:/Geotar/{pilot_nam}/geodata/workspace/{out_name}.tif"

    # Open the data source and read in the extent
    pixel_size = 0.0022457882102988 #250m

    vector_ds  = ogr.Open(mask_shp)
    shp_layer = vector_ds.GetLayer()

    xmin, xmax, ymin, ymax = shp_layer.GetExtent()

    # check if the output file already exists, and delete it if it does
    if os.path.exists(dst_filename):
        print(f'{dst_filename} exists, deleting...')
        drv = gdal.GetDriverByName('GTiff')
        drv.Delete(dst_filename)

    # rasterize the vectori file with the spatial resolution defined
    ds = gdal.Rasterize(dst_filename, clipped_shape, xRes=pixel_size, yRes=pixel_size,
                        burnValues=1,outputBounds=[xmin, ymin, xmax, ymax],
                        outputType=gdal.GDT_Byte, allTouched=True)
    ds = None
    source_ds = None


    src_ds = gdal.Open(dst_filename)

    # get the first band of the source raster file
    srcband = src_ds.GetRasterBand(1)

    # check if the output file already exists, and delete it if it does
    if os.path.exists(output):
        print(f'{output} exists, deleting...')
        drv = gdal.GetDriverByName('GTiff')
        drv.Delete(dst_filename)

    # create a new raster file with the same dimensions and data type as the source raster file
    # but with only one band of Float32 data type
    empty_raster = gdal.GetDriverByName('GTiff')
    dst_ds = empty_raster.Create(output,
                        src_ds.RasterXSize,
                        src_ds.RasterYSize, 1,
                        gdal.GetDataTypeByName('Float32'))

    # set the geotransform and projection of the output raster file
    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjectionRef())

    #get the first band of the output raster file
    dstband = dst_ds.GetRasterBand(1)

    # Compute the proximity of the input raster values to the raster value of 1
    # and write the resulting distances to the output raster file
    prox = gdal.ComputeProximity(srcband, dstband, ["VALUES=1", "DISTUNITS=GEO"])
    print(f'{output} file processed')
    # close the input and output raster files and bands to free up memory
    srcband = None
    dstband = None
    src_ds = None
    dst_ds = None
    prox = None
    return


proximity_rasters()

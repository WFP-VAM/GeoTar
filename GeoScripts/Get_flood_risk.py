import geopandas as gpd
from osgeo import gdal
import os
def get_original_resolution(dataset):
    # Assuming the resolution is stored in the GeoTransform (affine transformation parameters)
    x_res = abs(dataset.GetGeoTransform()[1])  # Get pixel width (x resolution)
    y_res = abs(dataset.GetGeoTransform()[5])  # Get pixel height (y resolution)
    return x_res, y_res
# Filter raster files based on bounding box intersection
def get_flood_risk(mask_shp, pilot):
    # Function to check if raster's bounding box intersects with the given bounding box
    def bbox_intersects(bbox1, bbox2):
        return not (bbox1[2] < bbox2[0] or bbox1[0] > bbox2[2] or bbox1[3] < bbox2[1] or bbox1[1] > bbox2[3])

    # Read the shapefile
    area_shp = gpd.read_file(mask_shp)

    # Get the bounding box of the shapefile
    bbox = area_shp.total_bounds

    # Define the FTP server URL and directory
    ftp_host = 'jeodpp.jrc.ec.europa.eu'
    ftp_directory = '/ftp/jrc-opendata/CEMS-GLOFAS/flood_hazard/RP10/'

    scenario = ftp_directory.rstrip('/').split('/')[-1]

    # Construct the FTP URL
    ftp_url = f'https://{ftp_host}{ftp_directory}'
    ftp_vsi_url = f'/vsicurl/{ftp_url}'

    # List raster files in the FTP directory
    try:
        file_list = gdal.ReadDir(ftp_vsi_url)
    except Exception as e:
        print(f"Error occurred while listing files from FTP: {e}")
        file_list = []

    print('Fetching flood data from the web...')
    selected_files = []
    for filename in file_list:
        if filename.lower().endswith('.tif'):
            file_path = f"{ftp_vsi_url}{filename}"
            try:
                # Open raster file using GDAL
                ds = gdal.Open(file_path)
                if ds:
                    # Get raster's bounding box
                    raster_bbox = ds.GetGeoTransform()  # (minX, pixel_width, rotation, maxY, rotation, pixel_height)
                    # Calculate maxX and minY using raster size
                    width = ds.RasterXSize
                    height = ds.RasterYSize
                    maxX = raster_bbox[0] + (width * raster_bbox[1])
                    minY = raster_bbox[3] + (height * raster_bbox[5])
                    raster_bbox = (raster_bbox[0], minY, maxX, raster_bbox[3])  # (minX, minY, maxX, maxY)
                    x_res, y_res = get_original_resolution(ds)

                    # Check if raster's bounding box intersects with the shapefile's bounding box
                    if bbox_intersects(raster_bbox, bbox):
                        selected_files.append(filename)
                else:
                    print(f"Failed to open raster file: {filename}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")
            finally:
                # Close the GDAL dataset
                ds = None

    # Print the selected files
    print('Selected Files within Bounding Box:')
    for file in selected_files:
        print(file)

    # # GDAL Warp options
    kwargs = {'format': 'GTiff',
              'xRes': x_res,
              'yRes': y_res,
              'dstSRS': 'EPSG:4326'}

    print('creating mosaic...')
    base_path = f'C://Geotar/{pilot}/geodata/processed/floods/'
    mosaic_path = f'{base_path}flood_risk_{scenario}.tif'

    if not os.path.exists(base_path):
        os.makedirs(base_path)

    gdal.Warp(mosaic_path, [ftp_vsi_url + file for file in selected_files], **kwargs)
    print(f'mosaic saved as:{mosaic_path}')
    return

get_flood_risk(f'C:/Geotar/SOM/geodata/Processed/Mask/SOM_mask.shp', 'SOM')

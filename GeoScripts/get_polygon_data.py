import osmnx as ox
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.geometry import Point
import pathlib
import Paths
import importlib
import requests
import json
import requests
import pandas as pd
from osgeo import gdal, ogr
import os


masks = {
    'SOM': { "pilot_name" : "SOM",
            "mask_shp" : f"C:/Geotar/SOM/geodata/Processed/Mask/SOM_mask.shp",
            "period" : "2021-05-01/2022-01-31", 'country_name':'Somalia'}}


def proximity_rasters(input_shp, mask_shp, pilot_name, out_name, output):
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
    clipped_shape = f"C:/Geotar/{pilot_name}/geodata/workspace/{out_name}.shp"
    #export shapefile
    gdf.to_file(clipped_shape)


    # Define NoData value of new raster
    NoData_value = -9999

    # set the name and location of the output raster file
    dst_filename = f"C:/Geotar/{pilot_name}/geodata/workspace/{out_name}.tif"

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


def get_roads(output_file):
    """
    Function to get road network from the OSM API
    :return: 
    returns the roads dataset for the mask area and saves it in the data folder
    """

    if not os.path.exists(output_file):

        # Load the GeoDataFrame from the shapefile
        area_shp = gpd.read_file(masks["SOM"]['mask_shp'])

        # Extract the bounding box coordinates
        bbox = area_shp.total_bounds

        # Create a Polygon object from the bounding box coordinates
        polygon = Polygon([(bbox[0], bbox[1]), (bbox[2], bbox[1]), (bbox[2], bbox[3]), (bbox[0], bbox[3])])

        # Download road data intersecting the polygon
        graph = ox.graph_from_polygon(polygon, network_type='all', simplify=True)

        # Convert the graph to a GeoDataFrame
        gdf = ox.graph_to_gdfs(graph, nodes=False, edges=True)

        # Filter roads intersecting the polygon
        roads_intersecting_polygon = gdf[gdf.geometry.intersects(polygon)]

        # Print the resulting GeoDataFrame
        print(roads_intersecting_polygon)
        # Save the file
        roads_intersecting_polygon.to_file(output_file)
        print(f"roads file saved for", masks["SOM"])
    else:
        print("Roads file already exists")
    return

def get_conflict(pilot_name, period):
    path = r'C:/Users/oscar.bautista/OneDrive - World Food Programme/GLOBAL'
    api_path = path / 'Geodata' / 'Raw' / 'ACLED'
    try:
        with open(api_path / 'acled_key.txt', 'r') as file:
            contents = file.read()
            print('Key loaded')  # Instead of return
    except FileNotFoundError:
        print('The file acled_key.txt does not exist.')
    except Exception as e:
        print(f'An error occurred: {e}')

    # set the period for the API call
    # Split the string using the '/' character
    start_date_str, end_date_str = period.split('/')

    # Convert the start date and end date strings to datetime objects
    start_date = pd.to_datetime(start_date_str)
    start_date = start_date.strftime('%Y-%m-%d')
    end_date = pd.to_datetime(end_date_str)
    end_date = end_date.strftime('%Y-%m-%d')
    print('Start Date:', start_date)
    print('End Date:', end_date)

    response = requests.get(f'https://api.acleddata.com/acled/read?'
                            f'key={contents}&email=oscar.bautista@wfp.org&'
                            f'country={country}&'
                            f'&limit=15000')
    # f'first_event_date={start_date}'
    print(response)

    # Check if the request was successful and the response exists
    if response and response.status_code == 200:
        try:
            # Convert the response to a JSON object
            data = response.json()
            # Check if the 'data' key exists in the JSON response
            if 'data' in data:
                formatted_data = []
                for record in data['data']:
                    formatted_record = {
                        'event_id_cnty': record.get('event_id_cnty', None),
                        'event_date': record.get('event_date', None),
                        'year': record.get('event_date', '').split('-')[0],
                        'month': record.get('event_date', '').split('-')[1],
                        'day': record.get('event_date', '').split('-')[2],
                        'actor1': record.get('actor1', None),
                        'actor2': record.get('actor2', None),
                        'interaction': record.get('interaction', None),
                        # 'iso3': record.get('iso3', None),
                        'country': record.get('country', None),
                        'admin1': record.get('admin1', None),
                        'admin2': record.get('admin2', None),
                        'admin3': record.get('admin3', None),
                        'location': record.get('location', None),
                        'fatalities': record.get('fatalities', None),
                        'latitude': record.get('latitude', None),
                        'longitude': record.get('longitude', None),
                        'geo_precision': record.get('geo_precision', None)

                    }
                    formatted_data.append(formatted_record)

                # Convert the formatted data to a pandas data frame
                df = pd.DataFrame(formatted_data)
            else:
                print('The data key is missing in the JSON response.')
        except ValueError as e:
            print('Error while parsing JSON:', e)
    else:
        print('Request failed with status code:', response.status_code)
    # Print the first 5 rows of the data frame
    # print(type(df))

    print('number of events: ', len(df['fatalities']))

    # chage data type of coordinate columns
    df['latitude'] = df['latitude'].astype(float)
    df['longitude'] = df['longitude'].astype(float)

    # chage data type of date columns
    df['year'] = df['year'].astype(int)
    df['month'] = df['month'].astype(int)
    df['day'] = df['day'].astype(int)
    df['fatalities'] = df['fatalities'].astype(int)
    df.loc[:, 'event_date'] = pd.to_datetime(df['event_date']).dt.strftime('%Y-%m-%d')
    # df['event_date'] = df['event_date'].strftime('%Y-%m-%d')

    df_period = df[df['event_date'].between(start_date, end_date)].copy()
    print(f'date filtered: {len(df_period["event_date"])}\ndate unfiltered: {len(df["event_date"])}')

    min_date_df = df['event_date'].min()
    max_date_df = df['event_date'].max()
    print(f'unfiltered dataset: \nfirst date: {min_date_df} last year: {max_date_df}')
    min_date = df_period['event_date'].min()
    max_date = df_period['event_date'].max()
    print(f'filtered dataset: \nfirst date: {min_date} last year: {max_date}')

    # Create a Point geometry column from the latitude and longitude columns
    geometry = [Point(xy) for xy in zip(df_period['longitude'], df_period['latitude'])]

    # Create a GeoDataFrame from the pandas DataFrame and the geometry column
    gdf = gpd.GeoDataFrame(df_period, geometry=geometry, crs='EPSG:4326')
    print(len(gdf))

    area_shp = gpd.read_file(mask_shp)

    geo_acled = gpd.clip(gdf, area_shp)
    print(len(geo_acled))

    min_year = geo_acled['year'].min()
    max_year = geo_acled['year'].max()
    print(f'first date: {min_year} last year: {max_year}')

    fatalities = geo_acled['fatalities'].sum()
    print(f'fatalities: {fatalities}')

    geo_acled.loc[:, 'event_date'] = geo_acled.loc[:, 'event_date'].astype(str)

    print(f'records in data {len(geo_acled["event_date"])}')
    conf_shp = fr'C:/Geotar/{pilot_name}/geodata/Raw/conflict/filtered_acled.shp'
    geo_acled.to_file(conf_shp)
    print(f'Conflict shapefile saved as {conf_shp}')

    # Define NoData value of new raster
    NoData_value = -9999

    # set the name and location of the output raster file
    # dst_filename = f'C:/Geotar/{pilot}/geodata/Processed/250m/dist_roads.tif'
    dst_filename = f'C:/Geotar/{pilot_name}/geodata/workspace/conflict_ras.tif'

    vector_ds = ogr.Open(mask_shp)
    shp_layer = vector_ds.GetLayer()

    # Open the data source and read in the extent
    source_ds = vector_ds
    pixel_size = 0.0022457882102988  # 250m

    xmin, xmax, ymin, ymax = shp_layer.GetExtent()

    # check if the output file already exists, and delete it if it does
    if os.path.exists(dst_filename):
        drv = gdal.GetDriverByName('GTiff')
        drv.Delete(dst_filename)

    # rasterize the vectori file suing the spatial resolution defined
    ds = gdal.Rasterize(dst_filename, mask_shp, xRes=pixel_size, yRes=pixel_size,
                        burnValues=1, outputBounds=[xmin, ymin, xmax, ymax],
                        outputType=gdal.GDT_Byte, allTouched=True)
    ds = None
    source_ds = None

    # Define NoData value of the new raster
    NoData_value = -9999

    # Set the name and location of the output raster file
    dst_filename = f'C:/Geotar/{pilot_name}/geodata/workspace/conflict_points.tif'

    # Define the path to your shapefile (conf_shp should be defined)
    conf_shp = f'C:/Geotar/{pilot_name}/geodata/Raw/conflict/filtered_acled.shp'

    # Read the extent of the shapefile layer
    pixel_size = 0.0022457882102988  # 250m
    xmin, xmax, ymin, ymax = shp_layer.GetExtent()

    # Check if the output file already exists and delete it if it does
    if os.path.exists(dst_filename):
        drv = gdal.GetDriverByName('GTiff')
        drv.Delete(dst_filename)

    # Rasterize the vector file using the spatial resolution defined
    ds = gdal.Rasterize(dst_filename, conf_shp, xRes=pixel_size, yRes=pixel_size,
                        burnValues=1, outputBounds=[xmin, ymin, xmax, ymax],
                        outputType=gdal.GDT_Byte, allTouched=True)
    ds = None
    print(f'Rasterized conflict data saved as: {dst_filename}')

    output = f'C:/Geotar/{pilot_name}/geodata/Processed/250m/dist_conflict.tif'

    src_ds = gdal.Open(dst_filename)

    # get the first band of the source raster file
    srcband = src_ds.GetRasterBand(1)

    # set the name and location of the output raster file
    # dst_filename = f'C:/Geotar/{pilot}/geodata/Processed/250m/dist_roads.tif'

    # check if the output file already exists, and delete it if it does
    if os.path.exists(dst_filename):
        drv = gdal.GetDriverByName('GTiff')
        drv.Delete(dst_filename)

    # create a new raster file with the same dimensions and data type as the source raster file
    # but with only one band of Float32 data type
    drv = gdal.GetDriverByName('GTiff')
    dst_ds = drv.Create(output,
                        src_ds.RasterXSize,
                        src_ds.RasterYSize, 1,
                        gdal.GetDataTypeByName('Float32'))

    # set the geotransform and projection of the output raster file
    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjectionRef())

    # get the first band of the output raster file
    dstband = dst_ds.GetRasterBand(1)

    # Compute the proximity of the input raster values to the raster value of 1
    # and write the resulting distances to the output raster file
    prox = gdal.ComputeProximity(srcband, dstband, ['VALUES=1', 'DISTUNITS=GEO'])

    # close the input and output raster files and bands to free up memory
    srcband = None
    dstband = None
    src_ds = None
    dst_ds = None
    prox = None
    return

def get_schools():
    '''
    Function to get school locations from OSM using the API
    Returns:


    '''
    print()
    tags = {"amenity": "school"}
    schools_query = ox.features_from_place(country_name, tags)
    schools_query.shape
    schools_query.reset_index(inplace=True)
    # get just the schools
    export_schools = schools_query[schools_query["amenity"] == "school"]
    export_schools = export_schools[['osmid', 'amenity', 'name', 'geometry']]
    # Filter out geometries that are not points
    export_schools = export_schools[export_schools['geometry'].apply(lambda geom: geom_type == 'Point')]
    school_file = f"C:/Geotar/{pilot_name}/geodata/Processed/Education/{pilot_name}_education.shp"
    export_schools.to_file(school_file)
    print(school_file, " successfully processed")

    return

def get_healthsites():


    return

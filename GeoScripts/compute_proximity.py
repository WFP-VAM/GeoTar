import geopandas as gpd
from osgeo import gdal, ogr
from S3_functions import put_tif_to_s3

def proximity_rasters(pilot_name: str, input_shp: str, mask_shp: str, out_name: str):
    """
    Process vector data into proximity raster data
    Args:
        pilot_name:
        input_shp:
        mask_shp:
        out_name:
        output:
    Returns:
        Raster proximity files
    """
    print(f"Starting to create {out_name} proximity data")
    # print(f"mask: {mask_shp}")
    # Read the GeoPandas object
    gdf = gpd.read_file(input_shp)
    mask = gpd.read_file(mask_shp)
    # Clip the func_gdf using the mask
    gdf = gpd.clip(gdf, mask)

    root = "s3://geotar.s3.hq/"
    dir_s3 = out_name[0].upper() + out_name[1:]
    # clipped shapefile
    clipped_shape = f"{root}Geotar/{pilot_name}/geodata/Processed/{dir_s3}/{out_name}.geojson"
    # export file
    gdf.to_file(clipped_shape)

    # set the name and location of the output raster file
    dst_filename = f"/vsimem/{out_name}.tif"

    # Open the data source and read in the extent
    pixel_size = 0.0022457882102988  # 250m

    # load the mask
    vector_ds = ogr.Open(f'/vsis3/geotar.s3.hq/Geotar/{pilot_name}/geodata/Processed/Mask/{pilot_name}_mask.geojson')
    print(vector_ds)
    shp_layer = vector_ds.GetLayer()

    xmin, xmax, ymin, ymax = shp_layer.GetExtent()


    points = f'/vsis3/geotar.s3.hq/Geotar/{pilot_name}/geodata/Processed/{dir_s3}/{out_name}.geojson'
    gdal.Rasterize(dst_filename, points,
                   xRes=pixel_size, yRes=pixel_size,
                   burnValues=1, outputBounds=[xmin, ymin, xmax, ymax],
                   outputType=gdal.GDT_Byte, allTouched=True)
    # ds = None
    # source_ds = None
    # print('gdal rasterize ended')

    src_ds = gdal.Open(dst_filename)

    # print(gdal.Info(src_ds))

    # get the first band of the source raster file
    srcband = src_ds.GetRasterBand(1)

    # output = f"{rootgdal}Geotar/{pilot_name}/geodata/Processed/250m/dist_{out_name}.tif"
    output = f"/vsimem/dist_{out_name}.tif"

    # # check if the output file already exists, and delete it if it does
    # if os.path.exists(output):
    #     print(f'{output} exists, deleting...')
    #     drv = gdal.GetDriverByName('GTiff')
    #     drv.Delete(output)

    # create a new raster file with the same dimensions and data type as the source raster file
    # but with only one band of Float32 data type
    empty_raster = gdal.GetDriverByName('GTiff')
    dst_ds = empty_raster.Create(output,
                                 src_ds.RasterXSize,
                                 src_ds.RasterYSize, 1,
                                 gdal.GetDataTypeByName('Float32'))
    # print(gdal.Info(dst_ds))

    # set the geotransform and projection of the output raster file
    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjectionRef())

    # get the first band of the output raster file
    dstband = dst_ds.GetRasterBand(1)

    # Compute the proximity of the input raster values to the raster value of 1
    # and write the resulting distances to the output raster file
    gdal.ComputeProximity(srcband, dstband, ["VALUES=1", "DISTUNITS=GEO"])

    put_tif_to_s3(output, f'Geotar/{pilot_name}/geodata/Processed/250m/dist_{out_name}.tif', 'geotar.s3.hq', retries=3)

    gdal.Unlink(dst_filename)
    gdal.Unlink(output)

    print(f'{output} file processed')
    # close the input and output raster files and bands to free up memory
    srcband = None
    dstband = None
    src_ds = None
    dst_ds = None
    return

def get_travel_time(iso3:str, bbox:list):
    in_key = 'Geotar/GLOBAL/Geodata/Raw/Travel_time/travel_time_to_cities_11.tif'
    bucket_name = 'geotar.s3.hq'
    # Create S3 URLs
    input_tif = f'/vsis3/{bucket_name}/{in_key}'
    clipped_tif = '/vsimem/clipped_output.tif'

    # Open the input file
    ds = gdal.Open(input_tif)

    # Get the spatial reference of the input file
    input_srs = osr.SpatialReference()
    input_srs.ImportFromWkt(ds.GetProjection())

    # Define the target spatial reference (EPSG:4326)
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(4326)

    # Create a coordinate transformation
    transform = osr.CoordinateTransformation(input_srs, target_srs)

    # Transform the bounding box coordinates to the target CRS
    (xmin, ymin, _) = transform.TransformPoint(bbox[0], bbox[1])
    (xmax, ymax, _) = transform.TransformPoint(bbox[2], bbox[3])

    translate_options = gdal.TranslateOptions(
        projWin=[xmin, ymax, xmax, ymin],
        projWinSRS='EPSG:4326',
        outputSRS='EPSG:4326',
        noData=-9999)

    # Clip the raster using gdal.Translate
    gdal.Translate(
        clipped_tif,
        ds,
        options=translate_options
    )

    # Upload the clipped file to S3
    out_key = f'Geotar/{iso3}/geodata/Processed/Travel_time/travel_time_to_cities_11.tif'
    # Upload the file
    put_tif_to_s3(clipped_tif, out_key, bucket_name, retries=3)
    print(f'travel time uploaded to s3://{bucket_name}/{out_key}')
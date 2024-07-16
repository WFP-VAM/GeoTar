import geopandas as gpd
from osgeo import gdal, ogr
import boto3
import os


def check_and_delete_s3_object(bucket_name, key):
    # Initialize a session using Amazon S3
    s3_client = boto3.client('s3')

    # Check if the object exists in S3
    try:
        s3_client.head_object(Bucket=bucket_name, Key=key)
        print(f'{key} exists in bucket {bucket_name}, deleting...')

        # Delete the object
        s3_client.delete_object(Bucket=bucket_name, Key=key)


    except s3_client.exceptions.ClientError as e:
        # If a 404 error is raised, the object does not exist
        if e.response['Error']['Code'] == '404':
            print(f'{key} does not exist in bucket {bucket_name}')
        else:
            # Something else has gone wrong.
            raise


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
    print(f"mask: {mask_shp}")
    # Read the GeoPandas object
    gdf = gpd.read_file(input_shp)
    mask = gpd.read_file(mask_shp)
    # Clip the func_gdf using the mask
    gdf = gpd.clip(gdf, mask)

    root = "s3://geotar.s3.hq/"
    rootgdal = "/vsis3/geotar.s3.hq/"

    # clipped shapefile
    clipped_shape = f"{root}Geotar/{pilot_name}/geodata/workspace/{out_name}.geojson"
    # export shapefile
    gdf.to_file(clipped_shape)

    # set the name and location of the output raster file
    dst_filename = f"/vsimem/{out_name}.tif"

    # Open the data source and read in the extent
    pixel_size = 0.0022457882102988  # 250m

    vector_ds = ogr.Open(f'{rootgdal}Geotar/{pilot_name}/geodata/Processed/Mask/{pilot_name}_mask.shp')
    print(vector_ds)
    shp_layer = vector_ds.GetLayer()

    xmin, xmax, ymin, ymax = shp_layer.GetExtent()

    # check if the output file already exists, and delete it if it does
    if os.path.exists(dst_filename):
        print(f'{dst_filename} exists, deleting...')
        drv = gdal.GetDriverByName('GTiff')
        drv.Delete(dst_filename)
    else:
        print(f'{dst_filename} not created, starting to write...')

    # rasterize the vectori file with the spatial resolution defined
    # print()
    # print('rasterize conflict file')
    ras = gdal.Rasterize(dst_filename, f"{rootgdal}Geotar/{pilot_name}/geodata/workspace/{out_name}.geojson",
                         xRes=pixel_size, yRes=pixel_size,
                         burnValues=1, outputBounds=[xmin, ymin, xmax, ymax],
                         outputType=gdal.GDT_Byte, allTouched=True)
    ds = None
    # source_ds = None
    # print('gdal rasterize ended')

    src_ds = gdal.Open(dst_filename)
    gdal.Unlink(dst_filename)
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
    check_and_delete_s3_object('geotar.s3.hq', f'Geotar/{pilot_name}/geodata/Processed/250m/dist_{out_name}.tif')

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
    prox = gdal.ComputeProximity(srcband, dstband, ["VALUES=1", "DISTUNITS=GEO"])

    def put_tif_from_gdal_mem_dataset(key, bucket_name):
        # Load the dataset into the virtual filesystem
        temp_name = f"/vsimem/dist_{out_name}.tif"
        # Read the raw data from the virtual filesystem
        f = gdal.VSIFOpenL(temp_name, 'rb')
        gdal.VSIFSeekL(f, 0, 2)  # seek to end
        size = gdal.VSIFTellL(f)
        gdal.VSIFSeekL(f, 0, 0)  # seek to beginning
        data = gdal.VSIFReadL(1, size, f)
        gdal.VSIFCloseL(f)
        # Upload the raw data to s3
        boto3.client('s3').put_object(Key=key, Bucket=bucket_name, Body=data, ContentLength=size)
        gdal.Unlink(temp_name)

    gdal.Unlink(output)

    put_tif_from_gdal_mem_dataset(f'Geotar/{pilot_name}/geodata/Processed/250m/dist_{out_name}.tif', 'geotar.s3.hq')

    print(f'{output} file processed')
    # close the input and output raster files and bands to free up memory
    srcband = None
    dstband = None
    src_ds = None
    dst_ds = None
    prox = None
    return

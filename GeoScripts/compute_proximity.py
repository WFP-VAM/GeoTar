import geopandas as gpd
from osgeo import gdal, ogr
import os


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
    # Read the GeoPandas object
    gdf = gpd.read_file(input_shp)
    mask = gpd.read_file(mask_shp)
    # Clip the func_gdf using the mask
    gdf = gpd.clip(gdf, mask)

    # clipped shapefile
    clipped_shape = f"C:/Geotar/{pilot_name}/geodata/workspace/{out_name}.shp"
    # export shapefile
    gdf.to_file(clipped_shape)

    # set the name and location of the output raster file
    dst_filename = f"C:/Geotar/{pilot_name}/geodata/workspace/{out_name}.tif"

    # Open the data source and read in the extent
    pixel_size = 0.0022457882102988  # 250m

    vector_ds = ogr.Open(mask_shp)
    shp_layer = vector_ds.GetLayer()

    xmin, xmax, ymin, ymax = shp_layer.GetExtent()

    # check if the output file already exists, and delete it if it does
    if os.path.exists(dst_filename):
        print(f'{dst_filename} exists, deleting...')
        drv = gdal.GetDriverByName('GTiff')
        drv.Delete(dst_filename)

    # rasterize the vectori file with the spatial resolution defined
    _ = gdal.Rasterize(dst_filename, clipped_shape, xRes=pixel_size, yRes=pixel_size,
                        burnValues=1, outputBounds=[xmin, ymin, xmax, ymax],
                        outputType=gdal.GDT_Byte, allTouched=True)
    #ds = None
    #source_ds = None

    src_ds = gdal.Open(dst_filename)

    # get the first band of the source raster file
    srcband = src_ds.GetRasterBand(1)

    output = f"C:/Geotar/{pilot_name}/geodata/Processed/250m/dist_{out_name}.tif"

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

    # get the first band of the output raster file
    dstband = dst_ds.GetRasterBand(1)

    # Compute the proximity of the input raster values to the raster value of 1
    # and write the resulting distances to the output raster file
    prox = gdal.ComputeProximity(srcband, dstband, ["VALUES=1", "DISTUNITS=GEO"])
    print(f'{output} file processed')
    # close the input and output raster files and bands to free up memory
    # srcband = None
    # dstband = None
    # src_ds = None
    # dst_ds = None
    # prox = None
    return

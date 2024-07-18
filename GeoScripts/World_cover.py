from osgeo import gdal
import geopandas as gpd
from shapely.geometry import box
from tqdm.auto import tqdm
from pystac_client import Client
from odc.stac import configure_rio, stac_load
import os
import json
from S3_functions import put_tif_to_s3


def WorldcovertoMODIS(ref_tif, bucket_name, iso3):
    # area_shp = gpd.read_file(mask_shp)

    # Get the bounding box of the shapefile
    # bbox = area_shp.total_bounds
    # bbox
    # ## Create reference raster from MODIS stored in HDC

    TOKEN_PATH = "C:/Users/oscar.bautista/OneDrive - World Food Programme/Scripts/tk.json"
    HDC_STAC_URL = "https://api.earthobservation.vam.wfp.org/stac/"

    def _get_hdc_stac_param_from_env():

        if "JUPYTERHUB_USER" in os.environ:

            signer = None
            header = None
            aws = {}  # Get credentials for accessing S3 bucket

        else:

            def make_signer(fname="./tk.json"):
                """
                Loads token from file at fname, and returns a function patching request urls with said token
                """
                tk = ""
                with open(fname, "rt") as src:
                    tk = json.load(src)["tk"]

                def sign(url, _tk=tk):
                    signed = f"{url}?{_tk}"
                    return signed

                return sign

            signer = make_signer(TOKEN_PATH)
            header = {"origin": "https://wfp.org"}
            aws = None

        # Instantiate an API client pointing to the WFP HDC STAC API
        hdc_stac_client = Client.open(HDC_STAC_URL, headers=header)

        # Set up GDAL/rasterio configuration.
        configure_rio(cloud_defaults=True, verbose=False, aws=aws)

        return hdc_stac_client, signer

    # STAC CLIENTS
    hdc_stac_client, signer = _get_hdc_stac_param_from_env()

    NDVI = hdc_stac_client.search(bbox=bbox,
                                  # collections=["mod13q1_vim_native"],
                                  collections=["mxd13q1_vim_dekad"],
                                  datetime="2018-04-01/2018-04-30",  # emulates the period of data cube files
                                  ).get_all_items()

    res = 0.0022457882102988  # 250 or 0.01 for 1km
    ndvi_stack = stac_load(NDVI, output_crs='EPSG:4326', resolution=res, patch_url=signer, bbox=bbox)
    ndvi = ndvi_stack.resample(time='1M').mean()
    ndvi = ndvi.mean(dim="time")

    # ref_tif = f'/vsimem/MODIS_mask.tif'
    # write the data to a geotiff file
    # ndvi.rio.to_raster(ref_tif, driver='GTiff')
    # print(f"{ref_tif} created successfully")

    # # define the output folder path
    # s3_dir = f"/Geotar/{iso3}/geodata/Processed/LandCover/tiles"
    # create_s3_folder(bucket_name, s3_dir)

    # get AOI geometry (select a country name)
    # country = 'Somalia'

    s3_url_prefix = "https://esa-worldcover.s3.eu-central-1.amazonaws.com"

    # Create a GeoDataFrame with a single geometry representing the bounding box
    geom1 = gpd.GeoDataFrame(geometry=[box(bbox[0], bbox[1], bbox[2], bbox[3])], crs='EPSG:4326')
    geom1 = geom1.geometry.iloc[0]

    # load worldcover grid
    url = f'{s3_url_prefix}/esa_worldcover_grid.geojson'
    # print(url)
    grid = gpd.read_file(url)

    # get grid tiles intersecting AOI
    tiles = grid[grid.intersects(geom1)]
    print(tiles)
    # tiles = gpd.overlay(grid, geom1, how='intersection')

    year = 2021  # setting this to 2020 will download the v100 product instead

    # select version tag, based on the year
    version = {2020: 'v100',
               2021: 'v200'}[year]

    tiffslist = []

    for tile in tqdm(tiles.ll_tile):
        url = f"{s3_url_prefix}/{version}/{year}/map/ESA_WorldCover_10m_{year}_{version}_{tile}_Map.tif"
        tiffslist.append(f"/vsicurl/{url}")

    '''Function to mosaic the tiff files downloaded from the worldcover dataset AWS bucket
    Inputs are:
    1. the output tiff path
    2. path to the tiled tif files
    3. NDVI MODIS processed file
    '''
    # List all TIFF files in the working directory with complete file paths
    # tiffslist = [os.path.join(tiffs_path, f) for f in os.listdir(tiffs_path) if f.endswith('.tif')]
    # tiffs_list = list_tiff_files_s3(bucket_name, )
    print("Files available: ")

    for i in tiffslist:
        print(i)

    def get_extent(ref_tif):
        # print(ref_tif)
        dataset = gdal.Open(ref_tif)
        if dataset is None:
            print("Failed to open the raster file.")
            return None

        # Get geotransform information
        geotransform = dataset.GetGeoTransform()
        if geotransform is None:
            print("Failed to get geotransform information.")
            return None

        # Extract extent
        minX = geotransform[0]
        maxY = geotransform[3]
        maxX = minX + geotransform[1] * dataset.RasterXSize
        minY = maxY + geotransform[5] * dataset.RasterYSize

        # Extract x and y resolutions
        xRes = abs(geotransform[1])
        yRes = abs(geotransform[5])

        return minX, minY, maxX, maxY, xRes, yRes

    minX, minY, maxX, maxY, xRes, yRes = get_extent(ref_tif)
    # keyword arguments that define extent and pixel size, these match the footprint og the modis ndvi data
    kwargs = {'format': 'GTiff',
              'outputBounds': [minX, minY, maxX, maxY],
              'outputBoundsSRS': 'EPSG:4326',
              'xRes': xRes,
              'yRes': yRes,
              # 'width':4984,
              # 'height':6932,
              }

    dst_file = f'/vsimem/Worldcover_{iso3}.tif'
    # warp opetation
    gdal.Warp(dst_file, tiffslist, **kwargs)
    key = f'Geotar/{iso3}/geodata/Processed/LandCover/Worldcover_{iso3}.tif'
    put_tif_to_s3(dst_file, key, bucket_name)
    gdal.Unlink(dst_file)
    gdal.Unlink(ref_tif)
    print(key, "processed successfully")

    return



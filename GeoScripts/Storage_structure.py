import boto3
from botocore.exceptions import ClientError
import httpx
from typing import Any, Dict, Tuple, Union
from shapely.geometry import Polygon, box
import geopandas as gpd


def create_s3_folder(bucket_name, folder_path):
    s3 = boto3.client('s3')
    key = 'Geotar/' + folder_path + '/'

    # Check if the folder exists
    try:
        s3.head_object(Bucket=bucket_name, Key=key)
        print(f"Folder {folder_path} already exists in bucket {bucket_name}")
    except ClientError as e:
        if e.response['Error']['Code'] == '404':
            # Create the folder if it doesn't exist
            try:
                s3.put_object(Bucket=bucket_name, Key=key)
                print(f"Folder {folder_path} created successfully in bucket {bucket_name}")
            except ClientError as error:
                print(f"Creation of the folder {folder_path} failed: {error}")
        else:
            raise  # If the error is not '404 Not Found', re-raise the exception


def storage_structure_s3(bucket_name:str, iso3:str):
    s3 = boto3.client('s3')

    # Define the folder structure
    subfolders = [iso3]

    for folder in subfolders:
        folder_path = folder
        create_s3_folder(bucket_name, folder_path)

        geodata_path = f"{folder}/geodata"
        create_s3_folder(bucket_name, geodata_path)

        subfolders_level_2 = ["Raw", "Processed", "workspace", "arcGIS", "outputs", "Docs"]
        for subfolder in subfolders_level_2:
            subfolder_path = f"{geodata_path}/{subfolder}"
            create_s3_folder(bucket_name, subfolder_path)

            # Add subfolders to Raw and Processed folders
            if subfolder == "Raw" or subfolder == "Processed":
                data_path = f"{geodata_path}/{subfolder}"
                subfolders_level_3 = ["Boundaries", "Roads", "Health", "Education", "Vegetation",
                                      "Precipitation", "Conflict", "Population", "Elevation",
                                      "Mask", "1K", "250m", "NTL", "Temperature", "Buildings", 'LandCover']
                for subsubfolder in subfolders_level_3:
                    subsubfolder_path = f"{data_path}/{subsubfolder}"
                    create_s3_folder(bucket_name, subsubfolder_path)

    return

def get_admin_shapes(iso3: str, admin_level: int = 1):
    GeoJSON = Dict[str, Any]
    # An async version of this function, whose design would be preferable, can be found here:
    # https://github.com/WFP-VAM/cleo-applied-research/blob/climate-analysis/topics/climate-analysis/adm_rainfall_analysis.ipynb
    # However, it is very difficult to support code that runs in the same way in both notebooks and
    # python programs because nested event loops are not supported. See https://blog.jupyter.org/ipython-7-0-async-repl-a35ce050f7f7
    # Additionally, rate limiting from the VAM GeoAPI means that end-to-end times of the two variants of this function are similar.
    # For these reasons, here we implement synchronous calls to the API.
    def _fetch(adm0: str, admCode: str, client: httpx.AsyncClient) -> GeoJSON:
        while True:
            try:
                r = client.get(
                    "https://api.vam.wfp.org/geodata/GetGeoAdmins",
                    params=dict(adm0=adm0, admCode=admCode),
                )
                r.raise_for_status()
            except httpx.HTTPStatusError as err:
                if r.status_code == 429:
                    _ = time.sleep(4)
                    continue
                raise err

            break

        return r.json()

    def _fetch_country(adm0: str) -> GeoJSON:
        while True:
            try:
                r = client.get(
                    "https://api.vam.wfp.org/geodata/Country",
                    params=dict(countryCode=adm0),
                )
                r.raise_for_status()
            except httpx.HTTPStatusError as err:
                if r.status_code == 429:
                    _ = time.sleep(4)
                    continue
                raise err

            break

        return r.json()

    def _get_adm0_info(iso3: str) -> Tuple[str, str]:
        r = httpx.get(
            "https://api.vam.wfp.org/geodata/CountriesAndSubdivisions",
            params=dict(iso3=iso3),
        )
        r.raise_for_status()
        return r.json()[0]["adm0Code"], r.json()[0]["name"]

    if admin_level not in [0, 1, 2]:
        raise ValueError(
            f"admin_level must be either 0, 1, or 2. Received {admin_level}"
        )

    adm0_code, adm0_name = _get_adm0_info(iso3)

    client = httpx.Client(limits=httpx.Limits(max_connections=1))

    adm1_geo = _fetch(adm0_code, adm0_code, client)
    adm1 = gpd.GeoDataFrame.from_features(adm1_geo)
    adm1["adm0_Code"] = adm0_code

    if admin_level == 0:
        adm0_geo = _fetch_country(adm0_code)
        return gpd.GeoDataFrame.from_features(adm0_geo["features"])

    if admin_level == 1:
        return adm1

    if admin_level == 2:
        adm_sub_codes_list = [x["properties"]["Code"] for x in adm1_geo["features"]]
        sub_adm_vec = [_fetch(adm0_code, x, client) for x in adm_sub_codes_list]
        for i, x in enumerate(sub_adm_vec):
            for y in x["features"]:
                y["properties"]["adm1_Code"] = adm_sub_codes_list[i]

        adm2 = gpd.GeoDataFrame.from_features(
            [y for x in sub_adm_vec for y in x["features"]]
        )

        adm2["adm0_Code"] = adm0_code
        return adm2


# Example usage
# bucket_name = "geotar.s3.hq"
# iso3 = "LBN"
# storage_structure_s3(bucket_name, iso3)

polygon = get_admin_shapes(iso3, 0)

def generate_bbox(polygon, iso3):
    """
    Args:
        polygon:
        iso3:

    Returns:

    """
    # Extract the bounding box
    bounding_box = polygon.total_bounds  # returns (minx, miny, maxx, maxy)

    # Create a Shapely box from the bounding box coordinates
    bbox_polygon = box(*bounding_box)
    # Create a GeoDataFrame with the bounding box
    bbox_gdf = gpd.GeoDataFrame(geometry=[bbox_polygon], crs=polygon.crs)

    # Save the GeoDataFrame to a GeoJSON file
    bbox_gdf.to_file(f"s3://geotar.s3.hq/Geotar/{iso3}/geodata/Processed/Mask/{iso3}_Mask.geojson", driver="GeoJSON")
    return

generate_bbox(polygon, iso3)
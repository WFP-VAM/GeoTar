import boto3
from urllib.parse import urlparse
import geopandas as gpd
import os

def read_geospatial_file_from_s3(s3_path, local_dir='/tmp'):
        """
        Downloads a GeoJSON or Shapefile from S3 and reads it into a GeoDataFrame.

        Parameters:
        - s3_path: str, complete S3 path (e.g., 's3://bucket_name/object_key.geojson' or 's3://bucket_name/object_key.shp').
        - local_dir: str, local directory to temporarily save the downloaded file. Defaults to '/tmp'.

        Returns:
        - gdf: GeoDataFrame, the data read from the file.
        """
        # Parse the S3 path to get the bucket name and object key
        parsed_url = urlparse(s3_path)
        bucket_name = parsed_url.netloc
        object_key = parsed_url.path.lstrip('/')  # Remove leading slash

        # Determine the file type and set the local path
        file_extension = os.path.splitext(object_key)[1].lower()
        local_path = os.path.join(local_dir, os.path.basename(object_key))

        # Initialize a session using Amazon S3
        s3 = boto3.client('s3')

        # Download the file
        s3.download_file(bucket_name, object_key, local_path)

        base_name = os.path.splitext(local_path)[0]  # Base name for cleanup
        if file_extension == '.shp':
            # If the file is a Shapefile, also download associated files (.shx, .dbf, etc.)
            for ext in ['.shx', '.dbf', '.prj', '.cpg']:
                associated_key = f"{os.path.splitext(object_key)[0]}{ext}"
                associated_local_path = os.path.join(local_dir, os.path.basename(associated_key))
                try:
                    s3.download_file(bucket_name, associated_key, associated_local_path)
                except:
                    print(f"Associated file {associated_key} not found, continuing without it.")

        # Read the file using geopandas
        gdf = gpd.read_file(local_path)

        # Optionally, clean up the local files if you don't need them
        for file in os.listdir(local_dir):
            if file.startswith(os.path.basename(base_name)):
                os.remove(os.path.join(local_dir, file))

        return gdf
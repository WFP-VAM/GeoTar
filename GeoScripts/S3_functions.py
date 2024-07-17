import boto3
from osgeo import gdal
import logging
from botocore.config import Config
from botocore.exceptions import SSLError, EndpointConnectionError
import time

# Enable detailed logging
# logging.basicConfig(level=logging.DEBUG)

# Specify the S3 client with a specific configuration
s3_client = boto3.client('s3', region_name='eu-central-1', config=Config(retries={'max_attempts': 10, 'mode': 'standard'}))

def put_tif_to_S3(in_file, key, bucket_name, retries=3):
    # Load the dataset into the virtual filesystem
    temp_name = in_file
    # Read the raw data from the virtual filesystem
    f = gdal.VSIFOpenL(temp_name, 'rb')
    gdal.VSIFSeekL(f, 0, 2)  # seek to end
    size = gdal.VSIFTellL(f)
    gdal.VSIFSeekL(f, 0, 0)  # seek to beginning
    data = gdal.VSIFReadL(1, size, f)
    gdal.VSIFCloseL(f)
    # Upload the raw data to s3
    for attempt in range(retries):
        try:
            s3_client.put_object(Key=key, Bucket=bucket_name, Body=data, ContentLength=size)
            print("Upload successful")
            return
        except (SSLError, EndpointConnectionError) as e:
                print(f"Attempt {attempt + 1} failed: {e}")
                if attempt < retries - 1:
                    time.sleep(2 ** attempt)  # Exponential backoff
                else:
                    raise e
    gdal.Unlink(temp_name)
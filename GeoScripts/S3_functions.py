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

def put_tif_to_s3(in_file, key, bucket_name, retries=3):
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

def read_s3_acled_key(bucket_name, s3_key):
    try:
        # Read the file from S3
        obj = s3_client.get_object(Bucket=bucket_name, Key=s3_key)
        contents = obj['Body'].read().decode('utf-8')
        print('Key loaded')
        return contents  # Return the contents of the file if needed
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == 'NoSuchKey':
            print('The file does not exist.')
        else:
            print(f'An error occurred: {e}')
    return(contents)
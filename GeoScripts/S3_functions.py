import boto3
from osgeo import gdal
def put_tif_to_S3(in_file, key, bucket_name):
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
    boto3.client('s3').put_object(Key=key, Bucket=bucket_name, Body=data, ContentLength=size)
    gdal.Unlink(temp_name)
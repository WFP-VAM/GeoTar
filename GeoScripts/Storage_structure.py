import boto3
from botocore.exceptions import ClientError


def create_s3_folder(bucket_name, folder_path):
    s3 = boto3.client('s3')
    try:
        s3.put_object(Bucket=bucket_name, Key=('Geotar/' + folder_path + '/'))
        print(f"Folder {folder_path} created successfully in bucket {bucket_name}")
    except ClientError as error:
        print(f"Creation of the folder {folder_path} failed: {error}")


def stodarage_structure_s3(bucket_name, iso3):
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


# Example usage
bucket_name = 'geotar.s3.hq'
iso3 = 'LBN'
stodarage_structure_s3(bucket_name, iso3)

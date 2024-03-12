# specify the path for the new folders
parent_path = "C:\\Geotar"
#subfolders = ["GLOBAL","COL", "CHAD", "IRAQ_N", "IRAQ_D", 'IRAQ']




def stodarage_sctructure(iso3):
    subfolders = [iso3]
    # create the parent directory if it does not exist
    if not os.path.exists(parent_path):
        os.makedirs(parent_path)

    # create the subfolders
    for folder in subfolders:
        folder_path = os.path.join(parent_path, folder)
        if not os.path.exists(folder_path):
            try:
                os.mkdir(folder_path)
                print(f"Folder {folder} created successfully")
            except OSError as error:
                print(f"Creation of the folder {folder} failed: {error}")
        else:
            print(f"Folder {folder} already exists.")

        # create subfolders under each geodata folder
        geodata_path = os.path.join(folder_path, "geodata")
        if not os.path.exists(geodata_path):
            try:
                os.mkdir(geodata_path)
                print(f"Folder geodata created successfully under {folder_path}")
            except OSError as error:
                print(f"Creation of the folder geodata under {folder_path} failed: {error}")
        else:
            print(f"Folder geodata under {folder_path} already exists.")

        subfolders = ["Raw", "Processed", "workspace", "arcGIS", "outputs", "Docs"]
        for subfolder in subfolders:
            subfolder_path = os.path.join(geodata_path, subfolder)
            if not os.path.exists(subfolder_path):
                try:
                    os.mkdir(subfolder_path)
                    print(f"Folder {subfolder} created successfully under {geodata_path}")
                except OSError as error:
                    print(f"Creation of the folder {subfolder} under {geodata_path} failed: {error}")
            else:
                print(f"Folder {subfolder} under {geodata_path} already exists.")

            # add subfolders to Raw and Processed folders
            if subfolder == "Raw" or subfolder == "Processed":
                data_path = os.path.join(geodata_path, subfolder)
                subfolders = ["Boundaries", "Roads", "Health", "Education", "Vegetation",
                              "Precipitation", "Conflict", "Population", "Elevation",
                              "Mask", "1K","250m", "NTL", "Temperature", "Buildings", 'LandCover']
                for subfolder in subfolders:
                    subsubfolder_path = os.path.join(data_path, subfolder)
                    if not os.path.exists(subsubfolder_path):
                        try:
                            os.mkdir(subsubfolder_path)
                            print(f"Folder {subfolder} created successfully under {data_path}")
                        except OSError as error:
                            print(f"Creation of the folder {subfolder} under {data_path} failed: {error}")
                    else:
                        print(f"Folder {subfolder} under {geodata_path} already exists.")
    return()

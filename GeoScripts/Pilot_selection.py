def select_pilot_area():
    print('Current pilot options')
    pilot_options = {
        '0': {'name': 'GLOBAL', 'mask_shp': 's3://geotar.s3.hq//Geotar/GLOBAL/geodata/workspace/test_mask.shp',
              'period': '2021-05-01/2022-01-31'},
        '1': {'name': 'COL', 'mask_shp': 's3://geotar.s3.hq//Geotar/COL/geodata/Processed/Mask/COL_mask.shp',
              'period': '2021-05-01/2022-01-31'},
        '2': {'name': 'CHAD', 'mask_shp': 's3://geotar.s3.hq//Geotar/CHAD/geodata/Processed/Mask/Chad_mask.shp',
              'period': '2022-05-01/2023-01-31'},
        '3': {'name': 'IRAQ_D', 'mask_shp': 's3://geotar.s3.hq//Geotar/IRAQ_D/geodata/Processed/Mask/Dahuk_mask.shp',
              'period': '2021-11-01/2022-05-31'},
        '4': {'name': 'IRAQ_N', 'mask_shp': 's3://geotar.s3.hq//Geotar/IRAQ_N/geodata/Processed/Mask/Najaf_mask.shp',
              'period': '2021-11-01/2022-05-31'},
        '5': {'name': 'IRAQ', 'mask_shp': 's3://geotar.s3.hq//Geotar/IRAQ/geodata/Processed/Mask/Iraq_mask.shp',
              'period': '2021-11-01/2022-05-31'},
        '6': {'name': 'LBN', 'mask_shp': 's3://geotar.s3.hq//Geotar/LBN/geodata/Processed/Mask/LBN_mask.shp',
              'period': '2021-10-01/2022-04-30'},
        '7': {'name': 'VEN', 'mask_shp': 's3://geotar.s3.hq//Geotar/VEN/geodata/Processed/Mask/VEN_mask.shp',
              'period': '2023-01-01/2023-07-30'},
        '8': {'name': 'AFG', 'mask_shp': 's3://geotar.s3.hq//Geotar/AFG/geodata/Processed/Mask/AFG_mask.shp',
              'period': '2023-04-01/2023-07-30'},
        '9': {'name': 'SOM', 'mask_shp': 's3://geotar.s3.hq//Geotar/SOM/geodata/Processed/Mask/SOM_mask.shp',
              'period': '2023-04-01/2023-07-30'},
        '10': {'name': 'BGD', 'mask_shp': 's3://geotar.s3.hq//Geotar/BGD/geodata/Processed/Mask/BGD_mask.shp',
               'period': '2023-04-01/2023-07-30'}
    }

    for num, details in pilot_options.items():
        print(f"{num}. {details['name']}")

    print("Please select the pilot area:")
    while True:
        pilot = input()
        if pilot in pilot_options:
            return pilot_options[pilot]
        else:
            print("Invalid pilot area selected. Try again")



def print_selected_pilot(pilot_name):
    print(f"You selected {pilot_name}")


def main():
    selected_pilot = select_pilot_area()
    print_selected_pilot(selected_pilot['name'])

    # Use selected_pilot['mask_shp'] and selected_pilot['period'] as needed


if __name__ == "__main__":
    main()
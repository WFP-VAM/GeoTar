import osmnx as ox
from shapely.geometry import Polygon
from shapely.geometry import Point
import requests
import geopandas as gpd
import pandas as pd
from typing import List
import os
from compute_proximity import proximity_rasters
from S3_functions import read_s3_acled_key

class get_vectors:
    def __init__(self, bbox: List, period: str, pilot_name: str, country_name: str, mask_shp: str, root: str):

        self.bbox=bbox
        self.period=period
        self.pilot_name=pilot_name
        self.country_name = country_name
        self.mask_shp = mask_shp
        self.root = root
        #self.hdc_stac_client=hdc_stac_client
        #self.signer=signer

    def get_roads(self):
        """
        Function to get road network from the OSM API
        :return:
        returns the roads dataset for the mask area and saves it in the data folder
        """
        base_dir = f'{root}Geotar/{self.pilot_name}/geodata/Processed/Roads'
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
        output_file = os.path.normpath(os.path.join(base_dir, 'roads.geojson'))
        if not os.path.exists(output_file):
            print('Fetching Roads')

            # Load the GeoDataFrame from the shapefile
            #area_shp = self.area_shp

            # Extract the bounding box coordinates
            #bbox = self.bbox #area_shp.total_bounds

            # Create a Polygon object from the bounding box coordinates
            polygon = Polygon([(self.bbox[0], self.bbox[1]), (self.bbox[2], self.bbox[1]), (self.bbox[2], self.bbox[3]), (self.bbox[0], self.bbox[3])])

            # Download road data intersecting the polygon
            graph = ox.graph_from_polygon(polygon, network_type='all', simplify=True)

            # Convert the graph to a GeoDataFrame
            gdf = ox.graph_to_gdfs(graph, nodes=False, edges=True)

            # Filter roads intersecting the polygon
            roads_intersecting_polygon = gdf[gdf.geometry.intersects(polygon)]
            roads_intersecting_polygon.reset_index(inplace=True)
            roads_to_convert = roads_intersecting_polygon[['key', 'geometry', 'length']]

            roads_to_convert.to_file(output_file)

            proximity_rasters(self.pilot_name, output_file, self.mask_shp, "roads")

            print(f"roads file saved for", self.pilot_name)
        else:
            print("Roads file already exists")
        return

    def get_conflict(self):
        print('fetching ACLED data')
        # path = r'C:/Users/oscar.bautista/OneDrive - World Food Programme/GLOBAL/Geodata/Raw/ACLED'
        contents = read_s3_acled_key('geotar.s3.hq', 'Geotar/GLOBAL/Geodata/Raw/ACLED/acled_key.txt')

        # set the period for the API call
        # Split the string using the '/' character
        start_date_str, end_date_str = self.period.split('/')

        # Convert the start date and end date strings to datetime objects
        start_date = pd.to_datetime(start_date_str)
        start_date = start_date.strftime('%Y-%m-%d')
        end_date = pd.to_datetime(end_date_str)
        end_date = end_date.strftime('%Y-%m-%d')
        print('Start Date:', start_date)
        print('End Date:', end_date)

        response = requests.get(f'https://api.acleddata.com/acled/read?'
                                f'key={contents}&email=oscar.bautista@wfp.org&'
                                f'country={self.country_name}&'
                                f'&limit=15000')
        # f'first_event_date={start_date}'
        # print(response)

        # Check if the request was successful and the response exists
        if response and response.status_code == 200:
            try:
                # Convert the response to a JSON object
                data = response.json()
                # Check if the 'data' key exists in the JSON response
                if 'data' in data:
                    formatted_data = []
                    for record in data['data']:
                        formatted_record = {
                            'event_id_cnty': record.get('event_id_cnty', None),
                            'event_date': record.get('event_date', None),
                            'year': record.get('event_date', '').split('-')[0],
                            'month': record.get('event_date', '').split('-')[1],
                            'day': record.get('event_date', '').split('-')[2],
                            'actor1': record.get('actor1', None),
                            'actor2': record.get('actor2', None),
                            'interaction': record.get('interaction', None),
                            # 'iso3': record.get('iso3', None),
                            'country': record.get('country', None),
                            'admin1': record.get('admin1', None),
                            'admin2': record.get('admin2', None),
                            'admin3': record.get('admin3', None),
                            'location': record.get('location', None),
                            'fatalities': record.get('fatalities', None),
                            'latitude': record.get('latitude', None),
                            'longitude': record.get('longitude', None),
                            'geo_precision': record.get('geo_precision', None)

                        }
                        formatted_data.append(formatted_record)

                    # Convert the formatted data to a pandas data frame
                    df = pd.DataFrame(formatted_data)
                else:
                    print('The data key is missing in the JSON response.')
            except ValueError as e:
                print('Error while parsing JSON:', e)
        else:
            print('Request failed with status code:', response.status_code)
        # Print the first 5 rows of the data frame
        # print(type(df))

        print()
        print('number of events: ', len(df['fatalities']))

        # chage data type of coordinate columns
        df['latitude'] = df['latitude'].astype(float)
        df['longitude'] = df['longitude'].astype(float)

        # chage data type of date columns
        df['year'] = df['year'].astype(int)
        df['month'] = df['month'].astype(int)
        df['day'] = df['day'].astype(int)
        df['fatalities'] = df['fatalities'].astype(int)
        df.loc[:, 'event_date'] = pd.to_datetime(df['event_date']).dt.strftime('%Y-%m-%d')
        # df['event_date'] = df['event_date'].strftime('%Y-%m-%d')

        df_period = df[df['event_date'].between(start_date, end_date)].copy()
        print(f'date filtered: {len(df_period["event_date"])}\ndate unfiltered: {len(df["event_date"])}')

        min_date_df = df['event_date'].min()
        max_date_df = df['event_date'].max()
        print(f'unfiltered dataset: \nfirst date: {min_date_df} last year: {max_date_df}')
        min_date = df_period['event_date'].min()
        max_date = df_period['event_date'].max()
        print()
        print(f'filtered dataset: \nfirst date: {min_date} last year: {max_date}')

        # Create a Point geometry column from the latitude and longitude columns
        geometry = [Point(xy) for xy in zip(df_period['longitude'], df_period['latitude'])]

        # Create a GeoDataFrame from the pandas DataFrame and the geometry column
        gdf = gpd.GeoDataFrame(df_period, geometry=geometry, crs='EPSG:4326')
        # print(len(gdf))

        area_shp = gpd.read_file(self.mask_shp)
        geo_acled = gpd.clip(gdf, area_shp)
        # print(len(geo_acled))

        min_year = geo_acled['year'].min()
        max_year = geo_acled['year'].max()
        print(f'first date: {min_year} last year: {max_year}')

        fatalities = geo_acled['fatalities'].sum()
        print(f'fatalities: {fatalities}')

        geo_acled['event_date'] = geo_acled['event_date'].astype(str)

        print(f'records in data {len(geo_acled["event_date"])}')
        conf_shp = fr's3://geotar.s3.hq/Geotar/{self.pilot_name}/geodata/Raw/Conflict/filtered_acled.geojson'
        geo_acled.to_file(conf_shp)
        print(f'Conflict file saved as {conf_shp}')
        print()
        proximity_rasters(self.pilot_name, conf_shp, self.mask_shp, "conflict")

        return

    def get_schools(self):
        '''
        Function to get school locations from OSM using the API
        Returns:
        Point shapefile with school locations
        '''
        print('fetching schools...')
        tags = {"amenity": "school"}
        schools_query = ox.features_from_place(self.country_name, tags)
        schools_query.shape
        schools_query.reset_index(inplace=True)
        # get just the schools
        export_schools = schools_query[schools_query["amenity"] == "school"]
        export_schools = export_schools[['osmid', 'amenity', 'name', 'geometry']]
        # Filter out geometries that are not points
        export_schools = export_schools[export_schools['geometry'].apply(lambda geom: geom.type == 'Point')]
        school_file = f"C:/Geotar/{self.pilot_name}/geodata/Processed/Education/{self.pilot_name}_education.shp"
        export_schools.to_file(school_file)
        proximity_rasters(self.pilot_name, school_file, self.mask_shp, "education")
        return

    def get_healthsites(self):

        print('Fetching health sites...')
        input_shp = r"C:/Geotar/Global/Geodata/Raw/Health/health_facilities.shp"

        # Read the GeoPandas object
        gdf = gpd.read_file(input_shp)
        mask = gpd.read_file(self.mask_shp)

        # Check if gdf intersects with mask
        if gdf.geometry.intersects(mask.unary_union).any():
            # If there is an intersection, proceed with clipping
            gdf_clipped = gpd.clip(gdf, mask)

            output_shape = f"C:/Geotar/{self.pilot_name}/geodata/Processed/health/healthsites.shp"
            gdf_clipped.to_file(output_shape)
            print(f'health sites saved at: {output_shape}')
            proximity_rasters(self.pilot_name, output_shape, self.mask_shp, "health")
            return
        else:
            # No intersection found
            print("No intersection found between the health sites database and the mask.")
            return None



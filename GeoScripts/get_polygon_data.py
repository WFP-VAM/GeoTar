import osmnx as ox
import geopandas as gpd
from shapely.geometry import Polygon

masks = {
    'SOM': { "pilot_name" : "SOM",
            "mask_shp" : f"C:/Geotar/SOM/geodata/Processed/Mask/SOM_mask.shp",
            "period" : "2021-05-01/2022-01-31", 'country_name':'Somalia'}}


def get_roads():
    """
    Function to get road network from the OSM API
    :return: 
    returns the roads dataset for the mask area and saves it in the data folder
    """
    output_file = f"C:/Geotar/SOM/geodata/Processed/Roads/SOM_roads.gpkg"
    if not os.path.exists(output_file):

        # Load the GeoDataFrame from the shapefile
        area_shp = gpd.read_file(masks["SOM"]['mask_shp'])

        # Extract the bounding box coordinates
        bbox = area_shp.total_bounds

        # Create a Polygon object from the bounding box coordinates
        polygon = Polygon([(bbox[0], bbox[1]), (bbox[2], bbox[1]), (bbox[2], bbox[3]), (bbox[0], bbox[3])])

        # Download road data intersecting the polygon
        graph = ox.graph_from_polygon(polygon, network_type='all', simplify=True)

        # Convert the graph to a GeoDataFrame
        gdf = ox.graph_to_gdfs(graph, nodes=False, edges=True)

        # Filter roads intersecting the polygon
        roads_intersecting_polygon = gdf[gdf.geometry.intersects(polygon)]

        # Print the resulting GeoDataFrame
        print(roads_intersecting_polygon)
        # Save the file
        roads_intersecting_polygon.to_file(output_file)
        print(f"roads file saved for", masks["SOM"])
    else:
        print("Roads file already exists")
    return()

def get_conflict():

    return()

def get_schools():

    return()

def get_healthsites():

    return

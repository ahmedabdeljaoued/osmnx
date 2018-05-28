################################################################################
# Module: leisures.py
# Description: Download and plot leisure (leisure) from OpenStreetMap
# License: MIT, see full license in LICENSE.txt

################################################################################

import geopandas as gpd
from shapely.geometry import box
from shapely.geometry import LineString
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import Polygon

from .core import bbox_from_point
from .core import gdf_from_place
from .core import overpass_request
from .utils import bbox_to_poly
from .utils import geocode
from .utils import log


def parse_leisure_query(north, south, east, west, leisures=None, timeout=180, maxsize=''):
    """
    Parse the Overpass QL query based on the list of leisures.

    Parameters
    ----------

    north : float
        Northernmost coordinate from bounding box of the search area.
    south : float
        Southernmost coordinate from bounding box of the search area.
    east : float
        Easternmost coordinate from bounding box of the search area.
    west : float
        Westernmost coordinate of the bounding box of the search area.
    leisures : list
        List of leisures that will be used for finding the leisure(S) from the selected area.
    timeout : int
        Timeout for the API request.
    """
    if leisures:
        # Overpass QL template
        query_template = ('[out:json][timeout:{timeout}]{maxsize};((node["leisure"~"{leisures}"]({south:.6f},'
                          '{west:.6f},{north:.6f},{east:.6f});(._;>;););(way["leisure"~"{leisures}"]({south:.6f},'
                          '{west:.6f},{north:.6f},{east:.6f});(._;>;););(relation["leisure"~"{leisures}"]'
                          '({south:.6f},{west:.6f},{north:.6f},{east:.6f});(._;>;);););out;')

        # Parse leisures
        query_str = query_template.format(leisures="|".join(leisures), north=north, south=south, east=east, west=west,
                                          timeout=timeout, maxsize=maxsize)
    else:
        # Overpass QL template
        query_template = ('[out:json][timeout:{timeout}]{maxsize};((node["leisure"]({south:.6f},'
                          '{west:.6f},{north:.6f},{east:.6f});(._;>;););(way["leisure"]({south:.6f},'
                          '{west:.6f},{north:.6f},{east:.6f});(._;>;););(relation["leisure"]'
                          '({south:.6f},{west:.6f},{north:.6f},{east:.6f});(._;>;);););out;')

        # Parse amenties
        query_str = query_template.format(north=north, south=south, east=east, west=west,
                                          timeout=timeout, maxsize=maxsize)

    return query_str


def osm_leisure_download(polygon=None, leisures=None, north=None, south=None, east=None, west=None,
                     timeout=180, max_query_area_size=50*1000*50*1000):
    """
    Get points of interests (leisures) from OpenStreetMap based on selected leisure types.

    Parameters
    ----------
    poly : shapely.geometry.Polygon
        Polygon that will be used to limit the leisure search.
    leisures : list
        List of leisures that will be used for finding the leisure(s) from the selected area.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        Points of interest and the tags associated with them as geopandas GeoDataFrame.
    """

    if polygon:
        # Bounds
        west, south, east, north = polygon.bounds

        # Parse the Overpass QL query
        query = parse_leisure_query(leisures=leisures, west=west, south=south, east=east, north=north)

    elif not (north is None or south is None or east is None or west is None):
        # TODO: Add functionality for subdividing search area geometry based on max_query_area_size
        # Parse Polygon from bbox
        #polygon = bbox_to_poly(north=north, south=south, east=east, west=west)

        # Parse the Overpass QL query
        query = parse_leisure_query(leisures=leisures, west=west, south=south, east=east, north=north)

    else:
        raise ValueError('You must pass a polygon or north, south, east, and west')

    # Get the leisure(s)
    responses = overpass_request(data={'data': query}, timeout=timeout)

    return responses


def parse_vertice_nodes(osm_response):
    """
    Parse node vertices from OSM response.

    Parameters
    ----------
    osm_response : JSON
        OSM response JSON

    Returns
    -------
    Dict of vertex IDs and their lat, lon coordinates.
    """

    vertices = {}
    for result in osm_response['elements']:
        if 'type' in result and result['type'] == 'node':
            vertices[result['id']] = {'lat': result['lat'],
                                      'lon': result['lon']}
    return vertices


def parse_osm_way(vertices, response):
    """
    Parse ways (areas) from OSM node vertices.

    Parameters
    ----------
    vertices : Python dict
        Node vertices parsed from OSM response.

    Returns
    -------
    Dict of vertex IDs and their lat, lon coordinates.
    """

    if 'type' in response and response['type'] == 'way':
        nodes = response['nodes']
        try:
            polygon = Polygon([(vertices[node]['lon'], vertices[node]['lat']) for node in nodes])

            leisure = {'nodes': nodes,
                   'geometry': polygon,
                   'osmid': response['id']}

            if 'tags' in response:
                for tag in response['tags']:
                    leisure[tag] = response['tags'][tag]
            return leisure

        except Exception:
            log('Polygon has invalid geometry: {}'.format(nodes))
    return None


def parse_osm_node(response):
    """
    Parse points from OSM nodes.

    Parameters
    ----------
    response : JSON
        Nodes from OSM response.

    Returns
    -------
    Dict of vertex IDs and their lat, lon coordinates.
    """

    try:
        point = Point(response['lon'], response['lat'])

        leisure = {
            'osmid': response['id'],
            'geometry': point
        }
        if 'tags' in response:
            for tag in response['tags']:
                leisure[tag] = response['tags'][tag]

    except Exception:
        log('Point has invalid geometry: {}'.format(response['id']))

    return leisure


def invalid_multipoly_handler(gdf, relation, way_ids):
    """
    Handles invalid multipolygon geometries when there exists e.g. a feature without
    geometry (geometry == NaN)

    Parameters
    ----------

    gdf : gpd.GeoDataFrame
        GeoDataFrame with Polygon geometries that should be converted into a MultiPolygon object.
    relation : dict
        OSM 'relation' dictionary
    way_ids : list
        A list of 'way' ids that should be converted into a MultiPolygon object.
    """

    try:
        gdf_clean = gdf.dropna(subset=['geometry'])
        multipoly = MultiPolygon(list(gdf_clean['geometry']))
        return multipoly

    except Exception:
        log("Invalid geometry at relation id %s.\nWay-ids of the invalid MultiPolygon:" % (
        relation['id'], str(way_ids)))
        return None

def parse_osm_relations(relations, osm_way_df):
    """
    Parses the osm relations from osm ways and nodes.
    See more information about relations from OSM documentation: http://wiki.openstreetmap.org/wiki/Relation

    Parameters
    ----------
    relations : list
        OSM 'relation' items (dictionaries) in a list.
    osm_way_df : gpd.GeoDataFrame
        OSM 'way' features as a GeoDataFrame that contains all the 'way' features that will constitute the multipolygon relations.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame with MultiPolygon representations of the relations and the attributes associated with them.
    """

    gdf_relations = gpd.GeoDataFrame()

    # Iterate over relations and extract the items
    for relation in relations:
        if relation['tags']['type'] == 'multipolygon':
            try:
                # Parse member 'way' ids
                member_way_ids = [member['ref'] for member in relation['members'] if member['type'] == 'way']
                # Extract the ways
                member_ways = osm_way_df.reindex(member_way_ids)
                # Extract the nodes of those ways
                member_nodes = list(member_ways['nodes'].values)
                try:
                    # Create MultiPolygon from geometries (exclude NaNs)
                    multipoly = MultiPolygon(list(member_ways['geometry']))
                except Exception:
                    multipoly = invalid_multipoly_handler(gdf=member_ways, relation=relation, way_ids=member_way_ids)

                if multipoly:
                    # Create GeoDataFrame with the tags and the MultiPolygon and its 'ways' (ids), and the 'nodes' of those ways
                    geo = gpd.GeoDataFrame(relation['tags'], index=[relation['id']])
                    # Initialize columns (needed for .loc inserts)
                    geo = geo.assign(geometry=None, ways=None, nodes=None, element_type=None, osmid=None)
                    # Add attributes
                    geo.loc[relation['id'], 'geometry'] = multipoly
                    geo.loc[relation['id'], 'ways'] = member_way_ids
                    geo.loc[relation['id'], 'nodes'] = member_nodes
                    geo.loc[relation['id'], 'element_type'] = 'relation'
                    geo.loc[relation['id'], 'osmid'] = relation['id']

                    # Append to relation GeoDataFrame
                    gdf_relations = gdf_relations.append(geo, sort=False)
                    # Remove such 'ways' from 'osm_way_df' that are part of the 'relation'
                    osm_way_df = osm_way_df.drop(member_way_ids)
            except Exception:
                log("Could not handle OSM 'relation': {}".format(relation['id']))

    # Merge 'osm_way_df' and the 'gdf_relations'
    osm_way_df = osm_way_df.append(gdf_relations, sort=False)
    return osm_way_df


def create_leisure_gdf(polygon=None, leisures=None, north=None, south=None, east=None, west=None):
    """
    Parse GeoDataFrames from leisure json that was returned by Overpass API.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        geographic shape to fetch the building footprints within
    leisures: list
        List of leisures that will be used for finding the leisure(s) from the selected area.
        See available leisures from: http://wiki.openstreetmap.org/wiki/Key:leisure
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box

    Returns
    -------
    Geopandas GeoDataFrame with leisure(s) and the associated attributes.
    """

    responses = osm_leisure_download(polygon=polygon, leisures=leisures, north=north, south=south, east=east, west=west)

    # Parse vertices
    vertices = parse_vertice_nodes(responses)

    # leisure nodes
    leisure_nodes = {}

    # leisure ways
    leisure_ways = {}

    # A list of leisure relations
    relations = []

    for result in responses['elements']:
        if result['type'] == 'node' and 'tags' in result:
            leisure = parse_osm_node(response=result)
            # Add element_type
            leisure['element_type'] = 'node'
            # Add to 'leisures'
            leisure_nodes[result['id']] = leisure
        elif result['type'] == 'way':
            # Parse leisure area Polygon
            leisure_area = parse_osm_way(vertices=vertices, response=result)
            if leisure_area:
                # Add element_type
                leisure_area['element_type'] = 'way'
                # Add to 'leisure_ways'
                leisure_ways[result['id']] = leisure_area

        elif result['type'] == 'relation':
            # Add relation to a relation list (needs to be parsed after all nodes and ways have been parsed)
            relations.append(result)

    # Create GeoDataFrames
    gdf_nodes = gpd.GeoDataFrame(leisure_nodes).T
    gdf_nodes.crs = {'init': 'epsg:4326'}

    gdf_ways = gpd.GeoDataFrame(leisure_ways).T
    gdf_ways.crs = {'init': 'epsg:4326'}

    # Parse relations (MultiPolygons) from 'ways'
    gdf_ways = parse_osm_relations(relations=relations, osm_way_df=gdf_ways)

    # Combine GeoDataFrames
    gdf = gdf_nodes.append(gdf_ways, sort=False)

    return gdf


def leisures_from_point(point, distance=None, leisures=None):
    """
    Get point of interests (leisures) within some distance north, south, east, and west of
    a lat-long point.

    Parameters
    ----------
    point : tuple
        a lat-long point
    distance : numeric
        distance in meters
    leisures : list
        List of leisures that will be used for finding the Leisure(s) from the selected area.
        See available leisures from: http://wiki.openstreetmap.org/wiki/Key:leisure

    Returns
    -------
    GeoDataFrame
    """

    bbox = bbox_from_point(point=point, distance=distance)
    north, south, east, west = bbox
    return create_leisure_gdf(leisures=leisures, north=north, south=south, east=east, west=west)


def leisures_from_address(address, distance, leisures=None):
    """
    Get OSM points of Interests within some distance north, south, east, and west of
    an address.

    Parameters
    ----------
    address : string
        the address to geocode to a lat-long point
    distance : numeric
        distance in meters
    leisures : list
        List of leisures that will be used for finding the leisure(s) from the selected area. See available
        leisures from: http://wiki.openstreetmap.org/wiki/Key:leisure

    Returns
    -------
    GeoDataFrame
    """

    # geocode the address string to a (lat, lon) point
    point = geocode(query=address)

    # get buildings within distance of this point
    return leisures_from_point(point=point, leisures=leisures, distance=distance)


def leisures_from_polygon(polygon, leisures=None):
    """
    Get OSM points of interest within some polygon.

    Parameters
    ----------
    polygon : Polygon
        Polygon where the leisure(s) are search from.
    leisures : list
        List of leisures that will be used for finding the leisure(s) from the selected area.
        See available leisures from: http://wiki.openstreetmap.org/wiki/Key:leisure

    Returns
    -------
    GeoDataFrame
    """

    return create_leisure_gdf(polygon=polygon, leisures=leisures)


def leisures_from_place(place, leisures=None):
    """
    Get points of interest (leisures) within the boundaries of some place.

    Parameters
    ----------
    place : string
        the query to geocode to get geojson boundary polygon.
    leisures : list
        List of leisures that will be used for finding the leisure(s) from the selected area.
        See available leisures from: http://wiki.openstreetmap.org/wiki/Key:leisure

    Returns
    -------
    GeoDataFrame
    """

    city = gdf_from_place(place)
    polygon = city['geometry'].iloc[0]
    return create_leisure_gdf(polygon=polygon, leisures=leisures)

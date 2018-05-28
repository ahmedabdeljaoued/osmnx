################################################################################
# Module: public_transport.py
# Description: Download and plot public_transport (public_transport) from OpenStreetMap
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


def parse_public_transport_query(north, south, east, west, public_transport=None, timeout=180, maxsize=''):
    """
    Parse the Overpass QL query based on the list of public_transport.

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
    public_transport : list
        List of public_transport that will be used for finding the public_transport(S) from the selected area.
    timeout : int
        Timeout for the API request.
    """
    if public_transport:
        # Overpass QL template
        query_template = ('[out:json][timeout:{timeout}]{maxsize};((node["public_transport"~"{public_transport}"]({south:.6f},'
                          '{west:.6f},{north:.6f},{east:.6f});(._;>;););(way["public_transport"~"{public_transport}"]({south:.6f},'
                          '{west:.6f},{north:.6f},{east:.6f});(._;>;););(relation["public_transport"~"{public_transport}"]'
                          '({south:.6f},{west:.6f},{north:.6f},{east:.6f});(._;>;);););out;')

        # Parse public_transport
        query_str = query_template.format(public_transport="|".join(public_transport), north=north, south=south, east=east, west=west,
                                          timeout=timeout, maxsize=maxsize)
    else:
        # Overpass QL template
        query_template = ('[out:json][timeout:{timeout}]{maxsize};((node["public_transport"]({south:.6f},'
                          '{west:.6f},{north:.6f},{east:.6f});(._;>;););(way["public_transport"]({south:.6f},'
                          '{west:.6f},{north:.6f},{east:.6f});(._;>;););(relation["public_transport"]'
                          '({south:.6f},{west:.6f},{north:.6f},{east:.6f});(._;>;);););out;')

        # Parse public_trnsport
        query_str = query_template.format(north=north, south=south, east=east, west=west,
                                          timeout=timeout, maxsize=maxsize)

    return query_str


def osm_public_transport_download(polygon=None, public_transport=None, north=None, south=None, east=None, west=None,
                     timeout=180, max_query_area_size=50*1000*50*1000):
    """
    Get points of interests (public_transport) from OpenStreetMap based on selected public_transport types.

    Parameters
    ----------
    poly : shapely.geometry.Polygon
        Polygon that will be used to limit the public_transport search.
    public_transport : list
        List of public_transport that will be used for finding the public_transport(s) from the selected area.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        Points of interest and the tags associated with them as geopandas GeoDataFrame.
    """

    if polygon:
        # Bounds
        west, south, east, north = polygon.bounds

        # Parse the Overpass QL query
        query = parse_public_transport_query(public_transport=public_transport, west=west, south=south, east=east, north=north)

    elif not (north is None or south is None or east is None or west is None):
        # TODO: Add functionality for subdividing search area geometry based on max_query_area_size
        # Parse Polygon from bbox
        #polygon = bbox_to_poly(north=north, south=south, east=east, west=west)

        # Parse the Overpass QL query
        query = parse_public_transport_query(public_transport=public_transport, west=west, south=south, east=east, north=north)

    else:
        raise ValueError('You must pass a polygon or north, south, east, and west')

    # Get the public_transport(s)
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

            public_transport = {'nodes': nodes,
                   'geometry': polygon,
                   'osmid': response['id']}

            if 'tags' in response:
                for tag in response['tags']:
                    public_transport[tag] = response['tags'][tag]
            return public_transport

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

        public_transport = {
            'osmid': response['id'],
            'geometry': point
        }
        if 'tags' in response:
            for tag in response['tags']:
                public_transport[tag] = response['tags'][tag]

    except Exception:
        log('Point has invalid geometry: {}'.format(response['id']))

    return public_transport


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
                    gdf_relations = gdf_relations.append(geo)
                    # Remove such 'ways' from 'osm_way_df' that are part of the 'relation'
                    osm_way_df = osm_way_df.drop(member_way_ids)
            except Exception:
                log("Could not handle OSM 'relation': {}".format(relation['id']))

    # Merge 'osm_way_df' and the 'gdf_relations'
    osm_way_df = osm_way_df.append(gdf_relations)
    return osm_way_df


def create_public_transport_gdf(polygon=None, public_transport=None, north=None, south=None, east=None, west=None):
    """
    Parse GeoDataFrames from public_transport json that was returned by Overpass API.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        geographic shape to fetch the building footprints within
    public_transport: list
        List of public_transport that will be used for finding the public_transport(s) from the selected area.
        See available public_transport from: http://wiki.openstreetmap.org/wiki/Key:public_transport
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
    Geopandas GeoDataFrame with public_transport(s) and the associated attributes.
    """

    responses = osm_public_transport_download(polygon=polygon, public_transport=public_transport, north=north, south=south, east=east, west=west)

    # Parse vertices
    vertices = parse_vertice_nodes(responses)

    # public_transport nodes
    public_transport_nodes = {}

    # public_transport ways
    public_transport_ways = {}

    # A list of public_transport relations
    relations = []

    for result in responses['elements']:
        if result['type'] == 'node' and 'tags' in result:
            public_transport = parse_osm_node(response=result)
            # Add element_type
            public_transport['element_type'] = 'node'
            # Add to 'public_transport'
            public_transport_nodes[result['id']] = public_transport
        elif result['type'] == 'way':
            # Parse public_transport area Polygon
            public_transport_area = parse_osm_way(vertices=vertices, response=result)
            if public_transport_area:
                # Add element_type
                public_transport_area['element_type'] = 'way'
                # Add to 'public_transport_ways'
                public_transport_ways[result['id']] = public_transport_area

        elif result['type'] == 'relation':
            # Add relation to a relation list (needs to be parsed after all nodes and ways have been parsed)
            relations.append(result)

    # Create GeoDataFrames
    gdf_nodes = gpd.GeoDataFrame(public_transport_nodes).T
    gdf_nodes.crs = {'init': 'epsg:4326'}

    gdf_ways = gpd.GeoDataFrame(public_transport_ways).T
    gdf_ways.crs = {'init': 'epsg:4326'}

    # Parse relations (MultiPolygons) from 'ways'
    gdf_ways = parse_osm_relations(relations=relations, osm_way_df=gdf_ways)

    # Combine GeoDataFrames
    gdf = gdf_nodes.append(gdf_ways)

    return gdf


def public_transport_from_point(point, distance=None, public_transport=None):
    """
    Get point of interests (public_transport) within some distance north, south, east, and west of
    a lat-long point.

    Parameters
    ----------
    point : tuple
        a lat-long point
    distance : numeric
        distance in meters
    public_transport : list
        List of public_transport that will be used for finding the public_transport(s) from the selected area.
        See available public_transport from: http://wiki.openstreetmap.org/wiki/Key:public_transport

    Returns
    -------
    GeoDataFrame
    """

    bbox = bbox_from_point(point=point, distance=distance)
    north, south, east, west = bbox
    return create_public_transport_gdf(public_transport=public_transport, north=north, south=south, east=east, west=west)


def public_transport_from_address(address, distance, public_transport=None):
    """
    Get OSM points of Interests within some distance north, south, east, and west of
    an address.

    Parameters
    ----------
    address : string
        the address to geocode to a lat-long point
    distance : numeric
        distance in meters
    public_transport : list
        List of public_transport that will be used for finding the POIs from the selected area. See available
        public_transport from: http://wiki.openstreetmap.org/wiki/Key:public_transport

    Returns
    -------
    GeoDataFrame
    """

    # geocode the address string to a (lat, lon) point
    point = geocode(query=address)

    # get buildings within distance of this point
    return public_transport_from_point(point=point, public_transport=public_transport, distance=distance)


def public_transport_from_polygon(polygon, public_transport=None):
    """
    Get OSM points of interest within some polygon.

    Parameters
    ----------
    polygon : Polygon
        Polygon where the public_transport(s) are search from.
    public_transport : list
        List of public_transport that will be used for finding the public_transport(s) from the selected area.
        See available public_transport from: http://wiki.openstreetmap.org/wiki/Key:public_transport

    Returns
    -------
    GeoDataFrame
    """

    return create_public_transport_gdf(polygon=polygon, public_transport=public_transport)


def public_transport_from_place(place, public_transport=None):
    """
    Get points of interest (public_transport) within the boundaries of some place.

    Parameters
    ----------
    place : string
        the query to geocode to get geojson boundary polygon.
    public_transport : list
        List of public_transport that will be used for finding the public_transport(s) from the selected area.
        See available public_transport from: http://wiki.openstreetmap.org/wiki/Key:public_transport

    Returns
    -------
    GeoDataFrame
    """

    city = gdf_from_place(place)
    polygon = city['geometry'].iloc[0]
    return create_public_transport_gdf(polygon=polygon, public_transport=public_transport)

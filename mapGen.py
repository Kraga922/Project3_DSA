import os
import math
import pandas as pd
from pyrosm import OSM
import networkx as nx
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import osmnx as ox
# from shapely.geometry import LineString, MultiLineString, Point
from pyproj import CRS
# import random
from pyproj import Transformer
#used this website for pyrosm features such as custom filter
#https://pyrosm.readthedocs.io/en/latest/basics.html


# Path to PBF file
PBF_FILE = "data/small_Alachua.osm.pbf"
GRAPHML_PATH = "maps/AlachuaTest.graphml"

def fix_boolean_fields(G):
    #iterates over edges and makes sure the oneway attribute is boolean instead of yes/no
    #converts the "yes" or "true" to True, and "no"/"false"/"none" to False
    for u, v, key, data in G.edges(keys=True, data=True):
        if "oneway" in data:
            val = data["oneway"]
            # If the value is None or is a string representing None or "no", then set to False.
            if val is None or (isinstance(val, str) and val.lower() in ["none", "no", "false"]):
                data["oneway"] = False
            elif isinstance(val, str) and val.lower() in ["yes", "true"]:
                data["oneway"] = True
            elif isinstance(val, bool):
                # Already a boolean, no action needed.
                continue
            else:
                data["oneway"] = False
    return G

def load_and_convert_public_transport_graph(filepath: str = PBF_FILE, walk_penalty: float = 5.0) -> nx.MultiDiGraph:
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} does not exist")
    osm = OSM(filepath)

    #Pull the pure walking network
    nodes_walk, edges_walk = osm.get_network(
        nodes=True,
        network_type="walking"
    )

    #Pull the drivable network
    nodes_drive, edges_drive = osm.get_network(
        nodes=True,
        network_type="driving+service"
    )

    # Merge edges
    edges = pd.concat([edges_walk, edges_drive], ignore_index=True)
    if {"u", "v", "key"}.issubset(edges.columns):
        edges = edges.drop_duplicates(subset=["u", "v", "key"])
    else:
        edges = edges.drop_duplicates(subset=["u", "v"])

    # Merge nodes
    nodes = pd.concat([nodes_walk, nodes_drive], ignore_index=True)
    nodes = nodes.drop_duplicates(subset=["id"])

    # get bus-stop Points and append them as transfer nodes
    stops = osm.get_data_by_custom_criteria(
        custom_filter={"highway": ["bus_stop"]},
        filter_type="keep",
        keep_nodes=True,
        keep_ways=False,
        keep_relations=False
    )
    stops = stops.loc[stops.geometry.geom_type == "Point"] \
                 .drop_duplicates(subset=["id"])
    if not stops.empty:
        nodes = pd.concat([nodes, stops], ignore_index=True) \
                 .drop_duplicates(subset=["id"])

    # Build the NetworkX graph
    G = osm.to_graph(
        nodes=nodes,
        edges=edges,
        graph_type="networkx",
        osmnx_compatible=True
    )

    # Project + lengths +  weights
    G.graph["crs"] = CRS.from_epsg(4326)
    G = ox.project_graph(G)
    G = add_length_attribute(G)
    G = assign_multimodal_weights(G, walk_penalty=walk_penalty)

    return G


def map_stat(G: nx.MultiDiGraph) -> None:
    #prints out the nodes and edges of graph
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")


def show_map(G, node_size=1000, edge_linewidth=0.5):
    fig, ax = ox.plot_graph(G, node_size=node_size, edge_linewidth=edge_linewidth, show=True)
    # OSMnx calls plt.show() internally


def add_length_attribute(G):
    #makes sure each edge has a length, if there's a length already use that value, otherwise
    # if length is missing we calulate it using distance in meters in graph
    # returns the graph with every edge having a length
    for u, v, key, data in G.edges(keys=True, data=True):
        if data.get("length") is None:
            x1, y1 = G.nodes[u]['x'], G.nodes[u]['y']
            x2, y2 = G.nodes[v]['x'], G.nodes[v]['y']
            data['length'] = math.hypot(x2 - x1, y2 - y1)
    return G

def get_final_graph():
    # If a cached GraphML already exists, try loading it
    if os.path.exists(GRAPHML_PATH):
        try:
            G = ox.load_graphml(GRAPHML_PATH)
            G = ox.project_graph(G)
            G = add_length_attribute(G)
        except Exception:
            os.remove(GRAPHML_PATH)
            G = load_and_convert_public_transport_graph(PBF_FILE)
    else:
        G = load_and_convert_public_transport_graph(PBF_FILE)

    G = fix_boolean_fields(G)
    ox.save_graphml(G, filepath=GRAPHML_PATH)
    G = assign_multimodal_weights(G, walk_penalty=5.0)
    return G

def find_closest_node_from_latlon(G, lon, lat):
    #Given a graph G (with projected x/y coords) and a (lon, lat) point,
    #returns the closest node ID to that point.

    # Convert (lon, lat) to projected x/y using the graph's CRS
    crs_graph = G.graph.get("crs", None)
    if crs_graph is None:
        raise ValueError("Graph does not have a defined CRS.")

    transformer = Transformer.from_crs("epsg:4326", crs_graph, always_xy=True)
    x_target, y_target = transformer.transform(lon, lat)

    # Search for the closest node in projected space
    closest_node = None
    min_dist = float('inf')
    for node, data in G.nodes(data=True):
        x, y = data.get("x"), data.get("y")
        if x is not None and y is not None:
            dist = (x - x_target) ** 2 + (y - y_target) ** 2
            if dist < min_dist:
                min_dist = dist
                closest_node = node
    return closest_node

def assign_multimodal_weights(G, walk_penalty: float = 5.0) -> nx.MultiDiGraph:
    #if it's a route, weight = length, else it's a walkable highway, so add a penalty so that we prefer bus route

    for u, v, key, data in G.edges(keys=True, data=True):
        length = data.get("length", 1.0)
        if data.get("route") is not None:
            # transit edge: no extra penalty
            data["weight"] = length
        else:
            # walking edge: apply penalty
            data["weight"] = length * walk_penalty
    return G
if __name__ == "__main__":
    try:
        G = get_final_graph()
        #Theory Gville
        start_lat = 29.655880
        start_lon = -82.337473

        # Publix at University Village (Publix #1103, 3195 SW Archer Rd)
        target_lat = 29.641389  # 29°38′28.9″N
        target_lon = -82.344722  # 82°20′41.0″W

        start = find_closest_node_from_latlon(G, start_lon, start_lat)
        target = find_closest_node_from_latlon(G, target_lon, target_lat)

        if nx.has_path(G, start, target):
            print("Hooray")
        else:
            print("Not Hooray")

        print("Start node:", start)
        print("Target node:", target)

        map_stat(G)
        node_attr_count = sum(len(data) for _, data in G.nodes(data=True))
        edge_attr_count = sum(len(data) for _, _, _, data in G.edges(keys=True, data=True))
        total_data_values = node_attr_count + edge_attr_count
        print("Total data values:", total_data_values)

        show_map(G)
    except Exception as e:
        print(f"Failed to generate graph: {e}")
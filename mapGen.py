import os
import math
import warnings
from pyrosm import OSM
import networkx as nx
import matplotlib
matplotlib.use("TkAgg")
#import matplotlib.pyplot as plt
import osmnx as ox
from shapely.geometry import LineString, MultiLineString
from pyproj import CRS
import random
#used this website for pyrosm features such as custom filter
#https://pyrosm.readthedocs.io/en/latest/basics.html


# warnings.filterwarnings("ignore", category=FutureWarning)

# Path to PBF file
PBF_FILE = "data/north_central_fl.osm.pbf"
GRAPHML_PATH = "maps/north_central_fl.graphml"

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

def load_and_convert_public_transport_graph(filepath: str = PBF_FILE) -> nx.MultiDiGraph:
    #loads pbf and creates graph with public transportation and highway info
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} does not exist")
    # Initialize the OSM
    osm = OSM(filepath)

    # custom filter for highway and public transportation
    highway_values = ["motorway", "primary", "secondary"]
    public_transport_filter = {
        "route": ["bus"],
        "highway": highway_values
    }

    transit_data = osm.get_data_by_custom_criteria(
        custom_filter=public_transport_filter,
        filter_type="keep",
        keep_nodes=False,
        keep_ways=True,
        keep_relations=True
    )

    if transit_data.empty:
        raise ValueError("No public transport data found in the file.")
    # Build a custom graph from the transit data by extracting start and end points of each geometry.
    G = nx.MultiDiGraph()
    for idx, row in transit_data.iterrows():
        geom_obj = row.geometry
        if geom_obj is None:
            continue

        # If geometry is a LineString, use its first and last coordinates.
        if isinstance(geom_obj, LineString):
            start = tuple(geom_obj.coords[0])
            end = tuple(geom_obj.coords[-1])
            G.add_node(start, pos=start)
            G.add_node(end, pos=end)
            # Extract a limited set of attributes
            edge_attrs = {k: v for k, v in row.to_dict().items() if k in ["highway", "name", "route", "oneway", "tags"]}
            G.add_edge(start, end, **edge_attrs)

        # If geometry is a MultiLineString, iterate over its parts.
        elif isinstance(geom_obj, MultiLineString):
            for line in geom_obj.geoms:
                if isinstance(line, LineString):
                    start = tuple(line.coords[0])
                    end = tuple(line.coords[-1])
                    G.add_node(start, pos=start)
                    G.add_node(end, pos=end)
                    edge_attrs = {k: v for k, v in row.to_dict().items() if
                                  k in ["highway", "name", "route", "oneway", "tags"]}
                    G.add_edge(start, end, **edge_attrs)
        else:
            continue

    #print(f"Constructed graph w {G.number_of_nodes()} nodes and {G.number_of_edges()} edges from data.")

    # Set the CRS attribute in the graph (using WGS84) coordinate ref
    G.graph["crs"] = CRS.from_epsg(4326)

    # Relabel the graph nodes from tuples to integers, but preserve the original coordinates in 'old_label'
    G = nx.convert_node_labels_to_integers(G, label_attribute='old_label')

    # For each node, assign "x" and "y" based on the original coordinate
    for node, data in G.nodes(data=True):
        old_label = data.get('old_label', None)
        if old_label and isinstance(old_label, tuple) and len(old_label) == 2:
            lon, lat = old_label
            data["x"] = lon
            data["y"] = lat
        else:
            data["x"] = None
            data["y"] = None

    G = fix_boolean_fields(G)
    tmp_filepath = "temp_graph.graphml"
    ox.save_graphml(G, filepath=tmp_filepath)
    G_standard = ox.load_graphml(tmp_filepath)
    G_proj = ox.project_graph(G_standard)
    G_standard = add_length_attribute(G_proj)
    os.remove(tmp_filepath)

    return G_standard


def map_stat(G: nx.MultiDiGraph) -> None:
    #prints out the nodes and edges of graph
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")


def show_map(G, node_size=10, edge_linewidth=0.5):
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
    # Check if the final graph already exists:
    if os.path.exists(GRAPHML_PATH):
        G = ox.load_graphml(GRAPHML_PATH)
        G = ox.project_graph(G)
        G = add_length_attribute(G)
    else:
        G = load_and_convert_public_transport_graph(PBF_FILE)
        ox.save_graphml(G, filepath=GRAPHML_PATH)
    return G

#to get random nodes for algs
def get_random_node_from_bbox(G, min_x, min_y, max_x, max_y):
    candidates = [
        node for node, data in G.nodes(data=True)
        if data.get("x") is not None and data.get("y") is not None and
           min_x <= data["x"] <= max_x and
           min_y <= data["y"] <= max_y
    ]
    if not candidates:
        raise ValueError("No nodes found in this bounding box.")
    return random.choice(candidates)

"""
Bounding boxes for gainesville and Orlando areas for testing 
gain_x_min, gain_x_max = 357000, 385000 
gain_y_min, gain_y_max = 3263290, 3290360

orlando_x_min, orlando_x_max = 434600, 479620
orlando_y_min, orlando_y_max = 3134170, 3166610
    """
def get_gainesville_to_orlando_nodes(G):
    gain_x_min, gain_x_max = 357000, 385000
    gain_y_min, gain_y_max = 3263290, 3290360

    orlando_x_min, orlando_x_max = 434600, 479620
    orlando_y_min, orlando_y_max = 3134170, 3166610

    start = get_random_node_from_bbox(G, gain_x_min, gain_y_min, gain_x_max, gain_y_max)
    target = get_random_node_from_bbox(G, orlando_x_min, orlando_y_min, orlando_x_max, orlando_y_max)
    return start, target

if __name__ == "__main__":
    try:
        G = get_final_graph()
        start, target = get_gainesville_to_orlando_nodes(G)

        print(f"Start node: {start}")
        print(f"Location: ({G.nodes[start]['x']:.2f}, {G.nodes[start]['y']:.2f})")

        print(f"Target node: {target}")
        print(f"Location: ({G.nodes[target]['x']:.2f}, {G.nodes[target]['y']:.2f})")
        map_stat(G)
        ox.plot_graph(G, node_color=["red" if n == start else "blue" if n == target else "gray" for n in G.nodes],
                      node_size=10, show=True)

        # Use get_final_graph to load or generate the graph
        # final_graph = get_final_graph()
        # map_stat(final_graph)
        # show_map(final_graph)

    except Exception as e:
        print(f"Failed to generate graph: {e}")
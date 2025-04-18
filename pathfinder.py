import os
import math
import argparse
import warnings
from pyrosm import OSM
import networkx as nx
import matplotlib
matplotlib.use("TkAgg")
import osmnx as ox
from shapely.geometry import LineString, MultiLineString
from pyproj import CRS

# Suppress warnings
warnings.filterwarnings("ignore", category=UserWarning)
ox.settings.log_console = False

# Path configurations
PBF_FILE = "north_central_fl.osm.pbf"
GRAPHML_PATH = "north_central_fl.graphml"

def fix_boolean_fields(G):
    """Convert oneway values to proper booleans"""
    for u, v, key, data in G.edges(keys=True, data=True):
        if "oneway" in data:
            val = data["oneway"]
            if val in ["yes", "true", True]:
                data["oneway"] = True
            else:
                data["oneway"] = False
    return G

def load_and_convert_public_transport_graph(filepath: str = PBF_FILE) -> nx.MultiDiGraph:
    """Load OSM data and create transportation graph"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"OSM file {filepath} not found")
    
    osm = OSM(filepath)
    # custom_filter = {"route": ["bus"], "highway": ["motorway", "primary", "secondary"]}
    custom_filter = {"highway": True}

    
    # transit_data = osm.get_data_by_custom_criteria(
    #     custom_filter=custom_filter,
    #     filter_type="keep",
    #     keep_nodes=False,
    #     keep_ways=True,
    #     keep_relations=True
    # )

    transit_data = osm.get_network(network_type="all")

    G = nx.MultiDiGraph()
    for idx, row in transit_data.iterrows():
        geom = row.geometry
        if isinstance(geom, LineString):
            coords = list(geom.coords)
        elif isinstance(geom, MultiLineString):
            coords = [list(line.coords) for line in geom.geoms]
        else:
            continue

        for segment in coords if isinstance(coords[0], list) else [coords]:
            start = segment[0]
            end = segment[-1]
            attrs = {k: v for k, v in row.items() if k in ["highway", "name", "route", "oneway"]}
            G.add_node(start, x=start[0], y=start[1])
            G.add_node(end, x=end[0], y=end[1])
            G.add_edge(start, end, **attrs)

    G = fix_boolean_fields(G)
    G = nx.convert_node_labels_to_integers(G, label_attribute="coords")
    G = ox.add_edge_lengths(G)
    return G

def get_final_graph():
    """Load or generate the transportation graph"""
    if os.path.exists(GRAPHML_PATH):
        return ox.load_graphml(GRAPHML_PATH)
    else:
        G = load_and_convert_public_transport_graph(PBF_FILE)
        ox.save_graphml(G, GRAPHML_PATH)
        return G

def astar_shortest_path(G, source_node_id, target_node_id):
    """A* pathfinding implementation"""
    def heuristic(u, v):
        return math.hypot(G.nodes[v]['x'] - G.nodes[u]['x'], 
                         G.nodes[v]['y'] - G.nodes[u]['y'])
    
    try:
        path = nx.astar_path(G, source_node_id, target_node_id, heuristic, 'length')
        length = sum(G[u][v][0]['length'] for u, v in zip(path[:-1], path[1:]))
        return path, length
    except nx.NetworkXNoPath:
        raise ValueError(f"No path between nodes {source_node_id} and {target_node_id}")

def visualize_path(G, path):
    """Generate interactive map using Folium"""
    import folium
    
    # Get center coordinates
    start_node = path[0]
    lat, lon = G.nodes[start_node]['y'], G.nodes[start_node]['x']
    
    # Create Folium map
    m = folium.Map(location=[lat, lon], zoom_start=14)
    
    # Add route as red line
    route_coords = []
    for node in path:
        route_coords.append((
            G.nodes[node]['y'],  # Latitude
            G.nodes[node]['x']   # Longitude
        ))
    
    folium.PolyLine(
        locations=route_coords,
        color='blue',
        weight=5,
        opacity=0.7
    ).add_to(m)
    
    m.save('optimal_path.html')
    return m

def print_sample_nodes(G, num_nodes=5):
    """Display sample nodes for testing"""
    print("\nSample nodes in graph:")
    for i, node in enumerate(list(G.nodes())[:num_nodes]):
        print(f"Node {node}: ({G.nodes[node]['x']:.5f}, {G.nodes[node]['y']:.5f})")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='A* Pathfinding for Transportation Network')
    parser.add_argument('-s', '--source', type=int, help='Source node ID')
    parser.add_argument('-t', '--target', type=int, help='Target node ID')
    parser.add_argument('-l', '--list', action='store_true', help='List sample nodes')
    args = parser.parse_args()

    try:
        G = get_final_graph()

        
        if args.list:
            print_sample_nodes(G)
            exit()

        # Sample nodes from initial XML data
        if not args.source or not args.target:
            print("Using default sample nodes from XML header:")
            args.source = 99  # (-81.386314, 28.625401)
            args.target = 110  # (-81.386360, 28.626255)
            print(f"Source: {args.source}, Target: {args.target}")

        path, length = astar_shortest_path(G, args.source, args.target)

        print(G)

        print(f"First node coordinates: {G.nodes[path[0]]['x']}, {G.nodes[path[0]]['y']}")
        print(f"Last node coordinates: {G.nodes[path[-1]]['x']}, {G.nodes[path[-1]]['y']}")

        print(f"\nOptimal path found ({path})")
        print(f"\nOptimal path found ({len(path)} nodes)")
        print(f"Total length: {length:.2f} meters")
        visualize_path(G, path)
        # print("Map saved to optimal_path.html")

    except Exception as e:
        print(f"\nError: {str(e)}")
        if "No path" in str(e):
            print("Possible causes: Disconnected nodes or invalid transportation routes")

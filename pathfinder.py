import os
import math
import argparse
import warnings
import folium
from pyrosm import OSM
import networkx as nx
import matplotlib
matplotlib.use("TkAgg")
import osmnx as ox
from shapely.geometry import LineString, MultiLineString
from pyproj import CRS, Transformer
from mapGen import *

# Suppress warnings
warnings.filterwarnings("ignore", category=UserWarning)
ox.settings.log_console = False



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
    """Generate interactive map using Folium with corrected coordinate projection"""

    # Convert projected coordinates to lat/lon
    transformer = Transformer.from_crs(G.graph["crs"], "epsg:4326", always_xy=True)

    route_coords = []
    for node in path:
        x, y = G.nodes[node]["x"], G.nodes[node]["y"]
        lon, lat = transformer.transform(x, y)
        route_coords.append((lat, lon))  # Folium uses (lat, lon)

    print("Start of path:", route_coords[0])
    print("End of path:", route_coords[-1])

    # Center map on first coordinate
    m = folium.Map(location=route_coords[0], zoom_start=15)

    folium.PolyLine(
        locations=route_coords,
        color='blue',
        weight=5,
        opacity=0.8
    ).add_to(m)

    # Optional: add markers for start/end
    folium.Marker(route_coords[0], popup="Start", icon=folium.Icon(color="green")).add_to(m)
    folium.Marker(route_coords[-1], popup="End", icon=folium.Icon(color="red")).add_to(m)

    m.save("optimal_path.html")
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

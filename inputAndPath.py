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
import heapq
import time

import tkinter as tk
from tkinter import messagebox

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
    """A* pathfinding implementation from scratch"""
    def heuristic(u, v):
        return math.hypot(G.nodes[v]['x'] - G.nodes[u]['x'],
                         G.nodes[v]['y'] - G.nodes[u]['y'])
    
    open_set = {source_node_id}
    came_from = {}
    g_score = {node: float('inf') for node in G.nodes}
    g_score[source_node_id] = 0
    f_score = {node: float('inf') for node in G.nodes}
    f_score[source_node_id] = heuristic(source_node_id, target_node_id)

    while open_set:
        # Get node with lowest f_score
        current = min(open_set, key=lambda node: f_score[node])
        
        if current == target_node_id:
            # Reconstruct path
            path = [current]
            while current in came_from:
                current = came_from[current]
                path.append(current)
            path.reverse()
            
            # Calculate total length
            length = 0
            for u, v in zip(path[:-1], path[1:]):
                # Get first edge between nodes (same as original code)
                length += G[u][v][0]['length']
            return path, length
        
        open_set.remove(current)
        
        # Iterate through neighbors
        for neighbor in G.neighbors(current):
            # Get first edge's length (same as original [0] index)
            edge_length = G[current][neighbor][0]['length']
            tentative_g = g_score[current] + edge_length
            
            if tentative_g < g_score[neighbor]:
                came_from[neighbor] = current
                g_score[neighbor] = tentative_g
                f_score[neighbor] = tentative_g + heuristic(neighbor, target_node_id)
                if neighbor not in open_set:
                    open_set.add(neighbor)
    
    raise ValueError(f"No path between nodes {source_node_id} and {target_node_id}")

def dijkstra_shortest_path(G, source_node_id, target_node_id):
    """Dijkstra's shortest path implementation from scratch"""
    distances = {node: float('inf') for node in G.nodes}
    distances[source_node_id] = 0
    prev = {}
    heap = [(0, source_node_id)]
    
    while heap:
        current_dist, current_node = heapq.heappop(heap)
        
        # Early exit if target reached
        if current_node == target_node_id:
            break
            
        # Skip if better path already found
        if current_dist > distances[current_node]:
            continue
            
        for neighbor in G.neighbors(current_node):
            # Get first edge's length (same as A* implementation)
            edge_length = G[current_node][neighbor][0]['length']
            new_dist = current_dist + edge_length
            
            if new_dist < distances[neighbor]:
                distances[neighbor] = new_dist
                prev[neighbor] = current_node
                heapq.heappush(heap, (new_dist, neighbor))
    
    # Path reconstruction
    if current_node != target_node_id:
        raise ValueError(f"No path between nodes {source_node_id} and {target_node_id}")
    
    path = []
    current = target_node_id
    while current in prev:
        path.append(current)
        current = prev[current]
    path.append(source_node_id)
    path.reverse()
    
    return path, distances[target_node_id]


def visualize_path(G, a_path, dij_path):
    """Generate interactive map using Folium with corrected coordinate projection"""

    # Convert projected coordinates to lat/lon
    transformer = Transformer.from_crs(G.graph["crs"], "epsg:4326", always_xy=True)

    a_route_coords = []
    for node in a_path:
        x, y = G.nodes[node]["x"], G.nodes[node]["y"]
        lon, lat = transformer.transform(x, y)
        a_route_coords.append((lat, lon))  # Folium uses (lat, lon)

    dij_route_coords = []
    for node in dij_path:
        x, y = G.nodes[node]["x"], G.nodes[node]["y"]
        lon, lat = transformer.transform(x, y)
        dij_route_coords.append((lat, lon))  # Folium uses (lat, lon)

    print("Start of path:", a_route_coords[0])
    print("End of path:", a_route_coords[-1])

    # Center map on first coordinate
    m = folium.Map(location=a_route_coords[0], zoom_start=15)

    folium.PolyLine(
        locations=a_route_coords,
        color='blue',
        weight=5,
        opacity=0.8
    ).add_to(m)
    
    folium.PolyLine(
        locations=dij_route_coords,
        color='red',
        weight=5,
        opacity=0.8,
        dashArray='5,12'  
    ).add_to(m)

    # Optional: add markers for start/end
    folium.Marker(a_route_coords[0], popup="Start", icon=folium.Icon(color="green")).add_to(m)
    folium.Marker(a_route_coords[-1], popup="End", icon=folium.Icon(color="red")).add_to(m)

    m.save("combined_path.html")
    return m

def print_sample_nodes(G, num_nodes=5):
    """Display sample nodes for testing"""
    print("\nSample nodes in graph:")
    for i, node in enumerate(list(G.nodes())[:num_nodes]):
        print(f"Node {node}: ({G.nodes[node]['x']:.5f}, {G.nodes[node]['y']:.5f})")

def compareAlgorithms():
    #Compare time for each algorithm to execute
    parser = argparse.ArgumentParser(description='A* Pathfinding for Transportation Network')

    parser.add_argument('-s', '--source', type=int, default=2168047757, help='Source node ID (default: 4418047440)')
    parser.add_argument('-t', '--target', type=int, default=8077888404, help='Target node ID (default: 8001930196)')

    parser.add_argument('-l', '--list', action='store_true', help='List sample nodes')
    args = parser.parse_args()

    try:
        G = get_final_graph()

        
        if args.list:
            print_sample_nodes(G)
            exit()

        # Sample nodes from initial XML data
        if not args.source or not args.target:
            print(f"Source: {args.source}, Target: {args.target}")

        
        #Timing A* algorithm
        begin_time_astar = float(time.time())
        a_path, a_length = astar_shortest_path(G, args.source, args.target)
        astar_time = float(time.time()) - begin_time_astar
        print("A* algorithm time:" + astar_time)

        #Timing Dijkstra's algorithm
        begin_time_dijkstra = float(time.time())
        dij_path, dij_length = dijkstra_shortest_path(G, args.source, args.target)
        dijkstra_time = float(time.time()) - begin_time_dijkstra
        print("Dijkstra algorithm time: " + dijkstra_time)

        #print(G)

        #Check output
        print(f"First node coordinates: {G.nodes[a_path[0]]['x']}, {G.nodes[a_path[0]]['y']}")
        print(f"Last node coordinates: {G.nodes[a_path[-1]]['x']}, {G.nodes[a_path[-1]]['y']}")

        # print(f"\nA* path found ({a_path})")
        print(f"\nA* path found ({len(a_path)} nodes)")
        print(f"Total length of A*: {a_length:.2f} meters")

        # print(f"\nDijkstra path found ({dij_path})")
        print(f"\nDijkstra path found ({len(dij_path)} nodes)")
        print(f"Total length of Dijkstra: {dij_length:.2f} meters")
        visualize_path(G, a_path, dij_path)
        # print("Map saved to optimal_path.html")

    except Exception as e:
        print(f"\nError: {str(e)}")
        if "No path" in str(e):
            print("Possible causes: Disconnected nodes or invalid transportation routes")


def getNodeFromLongLat(start_address, end_address, G):
    geolocator = Nominatim(user_agent="myApp", timeout=10)

    startLoc = geolocator.geocode(start_address)
    if startLoc is None:
        raise ValueError(f"Could not geocode the start address: {start_address}")

    endLoc = geolocator.geocode(end_address)
    if endLoc is None:
        raise ValueError(f"Could not geocode the end address: {end_address}")

    start_lat, start_lon = startLoc.latitude, startLoc.longitude
    end_lat, end_lon = endLoc.latitude, endLoc.longitude

    start_node = find_closest_node_from_latlon(G, start_lon, start_lat)
    end_node = find_closest_node_from_latlon(G, end_lon, end_lat)

    return start_node, end_node


def main(start, destination, dij, astar, testing):

    if (testing):
        compareAlgorithms()
        return

    parser = argparse.ArgumentParser(description='Pathfinding for Transportation Network')

    parser.add_argument('-s', '--source', type=int, default=2168047757, help='Source node ID (default: 4418047440)')
    parser.add_argument('-t', '--target', type=int, default=8077888404, help='Target node ID (default: 8001930196)')

    parser.add_argument('-l', '--list', action='store_true', help='List sample nodes')
    args = parser.parse_args()

    try:
        G = get_final_graph()

        
        if args.list:
            print_sample_nodes(G)
            exit()

        # Sample nodes from initial XML data
        if not args.source or not args.target:
            print(f"Source: {args.source}, Target: {args.target}")

        startNode, endNode = getNodeFromLongLat(start, destination, G)

        if (astar):
            a_path, a_length = astar_shortest_path(G, startNode, endNode)
        if (dij):
            dij_path, dij_length = dijkstra_shortest_path(G, startNode, endNode)


        #print(G)

        if (astar and dij):
            visualize_path(G, a_path, dij_path)
        elif (astar):
            visualize_path(G, a_path)
        elif (dij):
            visualize_path(G, dij_path)


        # print("Map saved to optimal_path.html")

    except Exception as e:
        print(f"\nError: {str(e)}")
        if "No path" in str(e):
            print("Possible causes: Disconnected nodes or invalid transportation routes")

#Start GUI

#Set up window
window = tk.Tk()
window.title("Best Bus Route")
window.geometry("300x215")

#Starting location input box
startLabel = tk.Label(window, text="Enter starting location:")
startLabel.grid(row=0, column=0, padx=5, pady=5)
startBox = tk.Entry(window)
startBox.grid(row=0, column=1, padx=5, pady=5)

#Destination input box
endLabel = tk.Label(window, text="Enter destination:")
endLabel.grid(row=1, column=0, padx=5, pady=5)
endBox = tk.Entry(window)
endBox.grid(row=1, column=1, padx=5, pady=5)

#Use Dijkstra's toggle
toggleDij = tk.BooleanVar()
toggleDij_checkbox = tk.Checkbutton(window, text="Use Dijkstra's", variable=toggleDij)
toggleDij_checkbox.grid(row=2, columnspan=2, padx=5, pady=5)

#Use A* toggle
toggleA = tk.BooleanVar()
toggleA_checkbox = tk.Checkbutton(window, text="Use A*", variable=toggleA)
toggleA_checkbox.grid(row=3, columnspan=2, padx=5, pady=5)

#Toggle testing mode
test = tk.BooleanVar()
testBox = tk.Checkbutton(window, text="Testing Mode", variable=test)
testBox.grid(row=4, columnspan=2, padx=5, pady=5)

#Call main with given parameters
button = tk.Button(window, text="Find Route", command=lambda: main(startBox, endBox, toggleDij, toggleA, testBox))
button.grid(row=5, columnspan=2, padx=5, pady=5)

window.mainloop()
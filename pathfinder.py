from geopy.geocoders import Nominatim
import argparse
import warnings
import folium
import matplotlib
matplotlib.use("TkAgg")
from mapGen import *
import heapq
import time

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
                length += G[u][v][0]['length']
            return path, length
        
        open_set.remove(current)
        
        # Iterate through neighbors
        for neighbor in G.neighbors(current):
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

        if current_node == target_node_id:
            break
            
        if current_dist > distances[current_node]:
            continue
            
        for neighbor in G.neighbors(current_node):
            edge_length = G[current_node][neighbor][0]['length']
            new_dist = current_dist + edge_length
            
            if new_dist < distances[neighbor]:
                distances[neighbor] = new_dist
                prev[neighbor] = current_node
                heapq.heappush(heap, (new_dist, neighbor))
    
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


def visualize_path(G, a_path=None, dij_path=None):
    """Draw only the selected path(s)"""

    transformer = Transformer.from_crs(G.graph["crs"], "epsg:4326", always_xy=True)
    m = None

    if dij_path:
        d_coords = [
            transformer.transform(G.nodes[node]["x"], G.nodes[node]["y"])[::-1]
            for node in dij_path
        ]
        m = folium.Map(location=d_coords[0], zoom_start=15)
        folium.PolyLine(
            locations=d_coords,
            color='blue',
            weight=5,
            opacity=0.8,
            tooltip="Dijkstra Path"
        ).add_to(m)
        folium.Marker(d_coords[0], popup="Start", icon=folium.Icon(color="green")).add_to(m)
        folium.Marker(d_coords[-1], popup="End", icon=folium.Icon(color="red")).add_to(m)

    if a_path:
        a_coords = [
            transformer.transform(G.nodes[node]["x"], G.nodes[node]["y"])[::-1]
            for node in a_path
        ]
        if not m:
            m = folium.Map(location=a_coords[0], zoom_start=15)
        folium.PolyLine(
            locations=a_coords,
            color='red',
            weight=5,
            opacity=0.8,
            dash_array='5,12',
            tooltip="A* Path"
        ).add_to(m)
        folium.Marker(a_coords[0], popup="Start", icon=folium.Icon(color="green")).add_to(m)
        folium.Marker(a_coords[-1], popup="End", icon=folium.Icon(color="red")).add_to(m)
    #saves based on algorithm chosen
    if m:
        filename = "combined_path.html"
        if a_path and not dij_path:
            filename = "astar_path.html"
        elif dij_path and not a_path:
            filename = "dijkstra_path.html"

        m.save(filename)
        return filename


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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='A* Pathfinding for Transportation Network')



    parser.add_argument('-s', '--start_address', type=str, 
                        default= "2088 Museum Rd, Gainesville, Florida 32603", help='Source node ID (default: 4418047440)')
    parser.add_argument('-e', '--end_address', type=str,  
                        default= "3019 SW 12th St, Gainesville, Florida 32608", help='Source node ID (default: 4418047440)')


    # parser.add_argument('-s', '--source', type=int, default=2168047757, help='Source node ID (default: 4418047440)')
    # parser.add_argument('-t', '--target', type=int, default=8077888404, help='Target node ID (default: 8001930196)')

    parser.add_argument('-l', '--list', action='store_true', help='List sample nodes')
    args = parser.parse_args()

    try:
        G = get_final_graph()

        start_node, end_node = getNodeFromLongLat(args.start_address, args.end_address, G)
        
        if args.list:
            print_sample_nodes(G)
            exit()

        # Sample nodes from initial XML data
        # if not args.source or not args.target:
        #     print(f"Source: {start_node}, Target: {end_node}")

        #     print(f"Source: {start_node}, Target: {end_node}")
        
        a_path, a_length = astar_shortest_path(G, start_node, end_node)
        dij_path, dij_length = dijkstra_shortest_path(G, start_node, end_node)


        print(G)

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

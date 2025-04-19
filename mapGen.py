import os
import math
from pyrosm import OSM
import networkx as nx
import osmnx as ox
import matplotlib
matplotlib.use("TkAgg")
from shapely.geometry import LineString
from shapely.geometry.base import BaseGeometry
from pyproj import CRS, Transformer
import pandas as pd
#used this website for pyrosm features such as custom filter
#https://pyrosm.readthedocs.io/en/latest/basics.html


# Path to PBF file
PBF_FILE = "data/small_Alachua.osm.pbf"
GRAPHML_PATH = "maps/AlachuaBusRoutes2.graphml"

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
        raise FileNotFoundError(filepath)
    osm = OSM(filepath)

    # Load walk and drive networks
    nodes_walk, edges_walk = osm.get_network(nodes=True, network_type="walking")
    nodes_drive, edges_drive = osm.get_network(nodes=True, network_type="driving+service")

    # Load bus stops
    stops = osm.get_data_by_custom_criteria(
        custom_filter={"highway": ["bus_stop"]},
        filter_type="keep",
        keep_nodes=True,
        keep_ways=False,
        keep_relations=False
    )
    stops = stops[stops.geometry.type == "Point"]

    # Load bus route geometries
    bus_routes = osm.get_data_by_custom_criteria(
        custom_filter={"route": ["bus"]},
        filter_type="keep",
        keep_nodes=False,
        keep_ways=True,
        keep_relations=True
    )
    bus_edges = bus_routes[bus_routes.geometry.type.isin(["LineString", "MultiLineString"])].copy()
    bus_edges["route"] = "bus"

    # Merge all nodes
    all_nodes = (
        pd.concat([nodes_walk, nodes_drive, stops], ignore_index=True)
        .drop_duplicates(subset=["id"])
    )

    # Create coordinate -> node ID map
    node_lookup = {}
    for _, row in all_nodes.iterrows():
        if hasattr(row.geometry, "x") and hasattr(row.geometry, "y"):
            coord = (round(row.geometry.x, 6), round(row.geometry.y, 6))
            node_lookup[coord] = row["id"]

    # Convert bus route geometries to edges
    def extract_uv_from_geometry(df, node_lookup):
        new_rows = []
        for _, row in df.iterrows():
            geom = row.geometry
            if geom is None or geom.is_empty:
                continue

            if geom.geom_type == "LineString":
                coords = list(geom.coords)
            elif geom.geom_type == "MultiLineString":
                coords = list(geom.geoms[0].coords)
            else:
                continue

            for i in range(len(coords) - 1):
                coord_u = (round(coords[i][0], 6), round(coords[i][1], 6))
                coord_v = (round(coords[i + 1][0], 6), round(coords[i + 1][1], 6))
                u = node_lookup.get(coord_u)
                v = node_lookup.get(coord_v)

                if u is not None and v is not None:
                    new_row = row.copy()
                    new_row["u"] = u
                    new_row["v"] = v
                    new_row["geometry"] = LineString([coords[i], coords[i + 1]])
                    new_rows.append(new_row)
        return pd.DataFrame(new_rows)

    bus_edges = extract_uv_from_geometry(bus_edges, node_lookup)

    # Make sure all u or v are valid node IDs
    valid_node_ids = set(all_nodes["id"])
    bus_edges = bus_edges[bus_edges["u"].isin(valid_node_ids) & bus_edges["v"].isin(valid_node_ids)]

    # Merge all edges
    all_edges = pd.concat([edges_walk, edges_drive, bus_edges], ignore_index=True)

    # Final cleanup before graph construction
    required_edge_cols = ["u", "v", "geometry"]
    for col in required_edge_cols:
        if col not in all_edges.columns:
            raise ValueError(f"Missing required column in edges: {col}")
        if all_edges[col].isnull().any():
            raise ValueError(f"Column {col} has null values")

    all_edges = all_edges[all_edges["geometry"].apply(lambda g: isinstance(g, BaseGeometry) and not g.is_empty)]

    # networkx to build graph
    G = nx.MultiDiGraph()

    for _, row in all_nodes.iterrows():
        node_id = row["id"]
        G.add_node(node_id)
        if hasattr(row.geometry, "x") and hasattr(row.geometry, "y"):
            G.nodes[node_id]["x"] = row.geometry.x
            G.nodes[node_id]["y"] = row.geometry.y

    for _, row in all_edges.iterrows():
        u = row["u"]
        v = row["v"]
        data = row.drop(["u", "v"]).to_dict()
        G.add_edge(u, v, **data)

    G.graph["crs"] = CRS.from_epsg(4326)
    G = ox.project_graph(G)
    G = fix_boolean_fields(G)
    G = add_length_attribute(G)
    G = assign_multimodal_weights(G, walk_penalty)

    #bus_edge_count = sum(1 for *_ , d in G.edges(data=True) if d.get("route") == "bus")
    # print(
    #     f"Graph built: {G.number_of_nodes():,} nodes ; "
    #     f"{G.number_of_edges():,} edges ; "
    #     f"{bus_edge_count:,} bus edges"
    # )

    return G


def map_stat(G: nx.MultiDiGraph) -> None:
    #prints out the nodes and edges of graph
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")


def show_map(G, node_size=0, edge_linewidth=0.5, bus_edge_color="blue", bus_edge_width=1.5, show_bus_edges=True):
   # Plots the full graph, with all edges in light gray and bus edges highlighted in color.
    import matplotlib.pyplot as plt
    nodes_gdf, edges_gdf = ox.graph_to_gdfs(G, nodes=True, edges=True)

    # Plot background edges (all)
    fig, ax = plt.subplots(figsize=(8, 8))
    edges_gdf.plot(ax=ax, linewidth=edge_linewidth, edgecolor="lightgray", zorder=1)

    # plot bus edges
    if show_bus_edges and "route" in edges_gdf.columns:
        bus_edges = edges_gdf[edges_gdf["route"] == "bus"]
        if not bus_edges.empty:
            bus_edges.plot(
                ax=ax,
                linewidth=bus_edge_width,
                edgecolor=bus_edge_color,
                label="bus routes",
                zorder=2
            )
            print(f"Highlighted {len(bus_edges)} bus edges.")
        else:
            print("No bus edges found to highlight.")

    # plot nodes
    if node_size > 0:
        nodes_gdf.plot(ax=ax, markersize=node_size, color="black", zorder=3)

    ax.set_title("Full Graph with Bus Routes Highlighted")
    ax.set_axis_off()
    ax.legend()
    plt.tight_layout()
    plt.show()



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
        # stadium:
        # target_lat = 29.6500368
        # target_lon = -82.3487666
        #walmart:
        target_lat = 29.6642231
        target_lon = -82.3010636

        start = find_closest_node_from_latlon(G, start_lon, start_lat)
        target = find_closest_node_from_latlon(G, target_lon, target_lat)
        #node id's for diff places:
        #theory 4418047440
        #stadium 8001930196
        #walmart 10086965103

        print(f"Start node: {start}")
        print(f"Target node: {target}")
        if nx.has_path(G, start, target):
            print("Hooray")
        else:
            print("Not Hooray")

        map_stat(G)
        show_map(G)
    except Exception as e:
        print(f"Failed to generate graph: {e}")
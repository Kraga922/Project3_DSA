import os
import math

from pyrosm import OSM
import networkx as nx
import osmnx as ox
import matplotlib
matplotlib.use("TkAgg")

from pyproj import CRS, Transformer
import pandas as pd
#used this website for pyrosm features such as custom filter
#https://pyrosm.readthedocs.io/en/latest/basics.html


# Path to PBF file
PBF_FILE = "data/small_Alachua.osm.pbf"
GRAPHML_PATH = "maps/AlachuaBusRoutes.graphml"

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
    #Loads a multimodal graph (walking, driving, bus stops, and bus routes)
    #from a PBF file using Pyrosm.
    if not os.path.exists(filepath):
        raise FileNotFoundError(filepath)
    osm = OSM(filepath)

    # walking and driving networks
    nodes_walk, edges_walk = osm.get_network(nodes=True, network_type="walking")
    nodes_drive, edges_drive = osm.get_network(nodes=True, network_type="driving+service")

    # all bus stops
    stops = osm.get_data_by_custom_criteria(
        custom_filter={"highway": ["bus_stop"]},
        filter_type="keep",
        keep_nodes=True,
        keep_ways=False,
        keep_relations=False
    )
    stops = stops[stops.geometry.type == "Point"]
    print(f"Found {len(stops)} bus stops in the PBF file")

    # bus edges(LineStrings)
    bus_routes = osm.get_data_by_custom_criteria(
        custom_filter={"route": ["bus"]},
        filter_type="keep",
        keep_nodes=False,
        keep_ways=True,
        keep_relations=True
    )
    bus_edges = bus_routes[bus_routes.geometry.type.isin(["LineString", "MultiLineString"])].copy()
    bus_edges["route"] = "bus"
    print(f"Found {len(bus_edges)} bus route edges in the PBF file")

    # combine all nodes
    all_nodes = (
        pd.concat([nodes_walk, nodes_drive, stops], ignore_index=True)
        .drop_duplicates(subset=["id"])
    )

    # Build coordinate → ID lookup for node matching
    node_lookup = {}
    for _, row in all_nodes.iterrows():
        if hasattr(row.geometry, "x") and hasattr(row.geometry, "y"):
            coord = (round(row.geometry.x, 6), round(row.geometry.y, 6))
            node_lookup[coord] = row["id"]

    # Match u, v for bus edges using node_lookup
    def extract_uv_from_geometry(df, node_lookup):
        def get_uv_ids(geom):
            if geom is None or geom.is_empty:
                return None, None
            if geom.geom_type == "LineString":
                coords = list(geom.coords)
            elif geom.geom_type == "MultiLineString":
                coords = list(geom.geoms[0].coords)
            else:
                return None, None
            start = (round(coords[0][0], 6), round(coords[0][1], 6))
            end   = (round(coords[-1][0], 6), round(coords[-1][1], 6))
            return node_lookup.get(start), node_lookup.get(end)

        uv_ids = df["geometry"].apply(get_uv_ids)
        df["u"] = uv_ids.apply(lambda x: x[0])
        df["v"] = uv_ids.apply(lambda x: x[1])
        return df

    bus_edges = extract_uv_from_geometry(bus_edges, node_lookup)
    bus_edges = bus_edges.dropna(subset=["u", "v"])

    # combine all edges
    all_edges = pd.concat([edges_walk, edges_drive, bus_edges], ignore_index=True)

    # build graph
    G = osm.to_graph(
        nodes=all_nodes,
        edges=all_edges,
        graph_type="networkx",
        osmnx_compatible=True
    )

    # processing
    G.graph["crs"] = CRS.from_epsg(4326)
    G = ox.project_graph(G)
    G = fix_boolean_fields(G)
    G = add_length_attribute(G)
    G = assign_multimodal_weights(G, walk_penalty)

    # print statements
    bus_edge_count = sum(1 for *_ , d in G.edges(data=True) if d.get("route") == "bus")
    print(
        f"Graph built: {G.number_of_nodes():,} nodes ; "
        f"{G.number_of_edges():,} edges ; "
        f"{bus_edge_count:,} bus edges"
    )

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

        # stadium
        target_lat = 29.6500368
        target_lon = -82.3487666

        start = find_closest_node_from_latlon(G, start_lon, start_lat)
        target = find_closest_node_from_latlon(G, target_lon, target_lat)
        if nx.has_path(G, start, target):
            print("Hooray")
        else:
            print("Not Hooray")

        map_stat(G)
        show_map(G)
    except Exception as e:
        print(f"Failed to generate graph: {e}")
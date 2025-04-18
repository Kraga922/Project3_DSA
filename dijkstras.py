import math
import queue
import networkx as nx

def dijkstrasAlgorithm(G, start):
    #Set up list of nodes, visited set, and priority queue
    result = {}
    visited = {start}
    result[start] = [0, -1]
    pq = queue.PriorityQueue()

    #Add source to priority queue and create an entry in result for each vertex
    pq.put(0, start)
    for (vertex in G.nodes):
        if vertex not in visited:
            result[vertex] = [math.inf, -1]
            pq.put(math.inf, vertex)

    #Update distances and paths of connected nodes in result
    while not pq.empty():
        current = pq.get()
        visited.add(current)
        distance = result[current][0]
        for (v in G.neighbors(current)):
            if v not in visited:
                pq.put(v)
                if (distance + G[current][v]["distance"] < result[v][0]):
                    result[v][0] = distance + G[current][v]["distance"]
                    result[v][1] = current

    return result








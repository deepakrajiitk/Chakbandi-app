
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point, LineString
import networkx as nx

class ChakmargGraph:
    def __init__(self):
        self.graph = nx.Graph()
        self.nodes_by_coord = {}
    
    def createNode(self, point_coord):
        if point_coord in self.nodes_by_coord:
            return self.nodes_by_coord[point_coord]
        node_name = f"Node_{len(self.graph.nodes) + 1}"  # Generate a unique node name
        self.graph.add_node(node_name, point_coord=point_coord)
        self.nodes_by_coord[point_coord] = node_name
        return node_name
    
    def addEdge(self, node1, node2, distance, owners_covered):
        self.graph.add_edge(node1, node2, distance=distance, owners_covered=owners_covered)
    
    def getNodeByCoord(coord):
        return self.nodes_by_coord.get(coord)   
import sys
import pandas as pd
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point, LineString
import networkx as nx
import geopandas as gpd
sys.path.append('..') 
# from chakmarg_module.helpers import *
from classes.chakmarg_graph_class import *
from consolidation_module.helpers import *

def convertPointToLatLong(x, y, cell_size, xmin, ymin):
    return x * cell_size + xmin, y * cell_size + ymin


def getRectangleIslandWise(owners_dict):
    island_rectangles = {}
    island_owners = {}
    for owner, owner_obj in owners_dict.items():
        rectangle = owner_obj.rectangle
        # print(rectangle.owner, owner)
        if rectangle:
            island_id = rectangle.island_id
            if island_id in island_rectangles:
                island_owners[island_id].append(owner)
                island_rectangles[island_id].append(rectangle)
            else:
                island_owners[island_id] = [owner]
                island_rectangles[island_id] = [rectangle]
    return island_rectangles, island_owners


def addItemToEdgeOwnerDict(owners_nodes_dict, owner, node1, node2):
    key = str(node1) + "_" + str(node2)
    if key in owners_nodes_dict:
        owners_nodes_dict[key].add(owner)
    else:
        owners_nodes_dict[key] = set()
        owners_nodes_dict[key].add(owner)


def updateBorderMatrix(border_matrix, x1, x2, y1, y2, min_row, min_col, max_row, max_col):
    for i in range(x1, x2+1):
        border_matrix[i-min_row][y1-min_col] = 1
        border_matrix[i-min_row][y2-min_col] = 1
    
    for j in range(y1, y2+1):
        border_matrix[x1-min_row][j-min_col] = 1
        border_matrix[x2-min_row][j-min_col] = 1
    

def getCornerPoints(graph, island_rectangles, border_matrix, min_row, min_col, max_row, max_col, island_id):
    corner_points = set()
    owners_nodes_dict = {}
    edge_owner_dict = {}
    for rectangle_obj in island_rectangles[island_id]:
        rectangle = rectangle_obj.coordinates
        x1, y1, x2, y2 = rectangle[0], rectangle[1], rectangle[2], rectangle[3]
        # addition is to align borders of two neighboring parcels
        x2+=1
        y2+=1
        corner_points.add((x1, y1))
        corner_points.add((x1, y2))
        corner_points.add((x2, y1))
        corner_points.add((x2, y2))

        owner = rectangle_obj.owner

        node1 = graph.createNode((x1, y1))
        node2 = graph.createNode((x1, y2))
        node3 = graph.createNode((x2, y1))
        node4 = graph.createNode((x2, y2))
        
        owners_nodes_dict[owner] = [node1]
        owners_nodes_dict[owner].append(node2)
        owners_nodes_dict[owner].append(node3)
        owners_nodes_dict[owner].append(node4)

        addItemToEdgeOwnerDict(edge_owner_dict, owner, node1, node2)
        addItemToEdgeOwnerDict(edge_owner_dict, owner, node2, node4)
        addItemToEdgeOwnerDict(edge_owner_dict, owner, node3, node4)
        addItemToEdgeOwnerDict(edge_owner_dict, owner, node1, node3)

        updateBorderMatrix(border_matrix, x1, x2, y1, y2, min_row, min_col, max_row, max_col)
    return corner_points, owners_nodes_dict, edge_owner_dict


# def printCornerPoints(corner_points):
#     points = []
#     for point in corner_points:
#         points.append(Point(convertPointToLatLong(point[0], point[1])))
#     # Create a GeoDataFrame from the list of points
#     point_gdf = gpd.GeoDataFrame(geometry=points, columns=['geometry'])
#     point_gdf.to_file("../preference_outputs/points.shp")


def find_covered_owner(x1, y1, x2, y2, p, q, square_matrix, island_id):
    while(x1!=x2 or y1!=y2):
        if x1>=square_matrix.shape[0] or y1>=square_matrix.shape[1] or square_matrix[x1][y1].island_id!=island_id:
            x1 += p
            y1 += q
            continue
        owner = square_matrix[x1][y1].owner
        if owner!=-1:
            return owner
        x1 += p
        y1 += q
    return -1


def isValid(x, y, border_matrix, min_row, min_col, max_row, max_col):
    if x<min_row or x>max_row or y<min_col or y>max_col or border_matrix[x-min_row][y-min_col]==0:
        return False
    return True

# NOTE: didn't use visited matrix here because a point can be marked by line in different direction which might affect traversal of other direction lines
def findEndPoint(x1, y1, border_matrix, corner_points, min_row, min_col, max_row, max_col, p, q):
    x2, y2 = x1+p, y1+q
    while (x2, y2) not in corner_points:
        if isValid(x2, y2, border_matrix, min_row, min_col, max_row, max_col):
            x2, y2 = x2+p, y2+q
        else:
            return None
    return (x2, y2)

def findEdgePoints(direc, start_point, border_matrix, square_matrix, corner_points, min_row, min_col, max_row, max_col, graph_obj, all_island_owners, island_id):
    x1, y1 = start_point
    if direc=="left":
        end_point = findEndPoint(x1, y1, border_matrix, corner_points, min_row, min_col, max_row, max_col, 0, -1)
        if end_point:
            owner1 = find_covered_owner(start_point[0], start_point[1]-1, end_point[0], end_point[1]-1, 0, -1, square_matrix, island_id)
            owner2 = find_covered_owner(start_point[0]-1, start_point[1]-1, end_point[0]-1, end_point[1]-1, 0, -1, square_matrix, island_id)
            # parcel1 = (end_point[0], end_point[1])
            # parcel2 = (end_point[0]-1, end_point[1])
    if direc=="right":
        end_point = findEndPoint(x1, y1, border_matrix, corner_points, min_row, min_col, max_row, max_col, 0, 1)
        if end_point:
            owner1 = find_covered_owner(start_point[0], start_point[1], end_point[0], end_point[1], 0, 1, square_matrix, island_id)
            owner2 = find_covered_owner(start_point[0]-1, start_point[1], end_point[0]-1, end_point[1], 0, 1, square_matrix, island_id)
            # parcel1 = (end_point[0], end_point[1])
            # parcel1 = (end_point[0], end_point[1]-1)
            # parcel2 = (end_point[0]-1, end_point[1]-1)
    if direc=="top":
        end_point = findEndPoint(x1, y1, border_matrix, corner_points, min_row, min_col, max_row, max_col, -1, 0)
        if end_point:
            owner1 = find_covered_owner(start_point[0]-1, start_point[1], end_point[0]-1, end_point[1], -1, 0, square_matrix, island_id)
            owner2 = find_covered_owner(start_point[0]-1, start_point[1]-1, end_point[0]-1, end_point[1]-1, -1, 0, square_matrix, island_id)
            # parcel1 = (end_point[0], end_point[1])
            # parcel1 = (end_point[0], end_point[1])
            # parcel2 = (end_point[0], end_point[1]-1)
    if direc=="bottom":
        end_point = findEndPoint(x1, y1, border_matrix, corner_points, min_row, min_col, max_row, max_col, 1, 0)
        if end_point:
            owner1 = find_covered_owner(start_point[0], start_point[1], end_point[0], end_point[1], 1, 0, square_matrix, island_id)
            owner2 = find_covered_owner(start_point[0], start_point[1]-1, end_point[0], end_point[1]-1, 1, 0, square_matrix, island_id)
            # parcel1 = (end_point[0], end_point[1])
            # parcel1 = (end_point[0]-1, end_point[1])
            # parcel2 = (end_point[0]-1, end_point[1]-1)
    if end_point:
        node1 = graph_obj.createNode(start_point)
        node2 = graph_obj.createNode(end_point)
        # owner1 = square_matrix[parcel1[0]][parcel1[1]].owner
        # owner2 = square_matrix[parcel2[0]][parcel2[1]].owner
        # print(owner1, owner2)
        if owner1 not in all_island_owners:
            owner1 = -1
        if owner2 not in all_island_owners:
            owner2 = -1
        edge_distance = abs(start_point[0]-end_point[0])+abs(start_point[1]-end_point[1])
        graph_obj.addEdge(node1, node2, edge_distance, [owner1, owner2])



# FIXME: Getting some extra lines in output
def createChakmargGraph(square_matrix, island_id, islands_rect_cover_dict, island_owners, island_rectangles):
    min_row, min_col, max_row, max_col = islands_rect_cover_dict[island_id]
    all_island_owners = island_owners[island_id]
    max_row+=1
    max_col+=1
    num_rows = max_row-min_row+1
    num_cols = max_col-min_col+1
    border_matrix = np.zeros((num_rows, num_cols))
    graph_obj = ChakmargGraph()
    corner_points, owners_nodes_dict, edge_owner_dict = getCornerPoints(graph_obj, island_rectangles, border_matrix, min_row, min_col, max_row, max_col, island_id)
    # printCornerPoints(corner_points)
    for corner_point in corner_points:
        findEdgePoints("left", corner_point, border_matrix, square_matrix, corner_points, min_row, min_col, max_row, max_col, graph_obj, all_island_owners, island_id)
        findEdgePoints("right", corner_point, border_matrix, square_matrix, corner_points, min_row, min_col, max_row, max_col, graph_obj, all_island_owners, island_id)
        findEdgePoints("top", corner_point, border_matrix, square_matrix, corner_points, min_row, min_col, max_row, max_col, graph_obj, all_island_owners, island_id)
        findEdgePoints("bottom", corner_point, border_matrix, square_matrix, corner_points, min_row, min_col, max_row, max_col, graph_obj, all_island_owners, island_id)
    return graph_obj, owners_nodes_dict, edge_owner_dict


def findEntryEdges(graph):
    entryEdges = SortedSet(key = lambda edge: (-sum(1 for owner in graph.edges[edge]['owners_covered'] if owner!=-1), graph.edges[edge]['distance']))

    for edge in graph.edges:
        owners_covered = graph.edges[edge]["owners_covered"]
        if owners_covered is not None and owners_covered.count(-1) == 1:
            entryEdges.add(edge)

    return entryEdges


def addNodeNeighboringEdges(node, graph, entryEdges, edge_seen):
    for edge in graph.edges(node):
        if edge not in entryEdges and edge not in edge_seen:
            entryEdges.add(edge)
    return entryEdges


def createIslandChakmarg(graph, entryEdges):
    visited_owners = {-1}
    edge_seen = set()
    marg = []
    while True and len(entryEdges):
        candidate_edge = entryEdges[0]
        # print(candidate_edge, "is the candidate")
        # print(visited_owners, "visited owner status")
        node1, node2 = candidate_edge
        owners_covered = graph.edges[candidate_edge]['owners_covered']
        first_owner = owners_covered[0]
        second_owner = owners_covered[1]
        # print("first owner", first_owner, "second owner", second_owner)
        if first_owner==-1 and second_owner==-1:
            edge_seen.add(candidate_edge)
            entryEdges = addNodeNeighboringEdges(node1, graph, entryEdges, edge_seen)
            entryEdges = addNodeNeighboringEdges(node2, graph, entryEdges, edge_seen)
            entryEdges.discard(candidate_edge)
            continue
        if first_owner in visited_owners and second_owner in visited_owners:
            # print("continuing because both owner have been seen before")
            entryEdges.discard(candidate_edge)
            # print("_________________________________________________________")
            continue
        visited_owners.add(first_owner)
        visited_owners.add(second_owner)
        entryEdges.discard(candidate_edge)
        marg.append(candidate_edge)
        edge_seen.add(candidate_edge)
        # print(candidate_edge, "added")
        entryEdges = addNodeNeighboringEdges(node1, graph, entryEdges, edge_seen)
        entryEdges = addNodeNeighboringEdges(node2, graph, entryEdges, edge_seen)
        # print(len(entryEdges))
        # print("------------------------------------------------")
    return marg


def create_marg_gdf(graph, marg, island_id, cell_size, xmin, ymin, crs):
    lines = []
    for edge in marg:
        node1, node2 = edge[0], edge[1]
        x1, y1 = graph.nodes[node1]['point_coord']
        x2, y2 = graph.nodes[node2]['point_coord']
        start_point = Point(convertPointToLatLong(x1, y1, cell_size, xmin, ymin))
        end_point = Point(convertPointToLatLong(x2, y2, cell_size, xmin, ymin))
        owners_covered = graph.edges[node1, node2]['owners_covered']
        line_geometry = LineString([start_point, end_point])
        lines.append([line_geometry, owners_covered[0], owners_covered[1]])

    marg_gdf = gpd.GeoDataFrame(lines, columns = ['geometry', 'owner1', 'owner2'], crs=crs)
    marg_gdf['island_id'] = island_id
    return marg_gdf


def create_points_gdf(graph, island_id, cell_size, xmin, ymin):
    points = []

    for node in graph.nodes:
        x1, y1 = graph.nodes[node]['point_coord']
        points.append([Point(convertPointToLatLong(x1, y1, cell_size, xmin, ymin)), node])

    point_gdf = gpd.GeoDataFrame(points, columns=['geometry', 'Node name'])
    point_gdf['island_id'] = island_id
    return point_gdf


def create_line_gdf(graph, island_id, cell_size, xmin, ymin):
    lines = []

    for edge in graph.edges:
        node1, node2 = edge
        x1, y1 = graph.nodes[node1]['point_coord']
        x2, y2 = graph.nodes[node2]['point_coord']
        start_point = Point(convertPointToLatLong(x1, y1, cell_size, xmin, ymin))
        end_point = Point(convertPointToLatLong(x2, y2, cell_size, xmin, ymin))
        owners_covered = graph.edges[node1, node2]['owners_covered']
        line_geometry = LineString([start_point, end_point])
        lines.append([line_geometry, owners_covered[0], owners_covered[1]])

    line_gdf = gpd.GeoDataFrame(lines, columns=['geometry', 'owner1', 'owner2'])
    line_gdf['island_id'] = island_id
    return line_gdf


def is_valid_index(matrix, row, col):
    # Check if the row index is within bounds
    if 0 <= row < len(matrix):
        # Check if the column index is within bounds
        if 0 <= col < len(matrix[0]):
            return True
    return False

def has_neighboring_minus_one(matrix, row, col, island_id):
    # Check if the current cell is within bounds
    if is_valid_index(matrix, row, col):
        # Check the neighboring cells in four directions
        neighbors = [
            (row - 1, col),  # Up
            (row + 1, col),  # Down
            (row, col - 1),  # Left
            (row, col + 1),  # Right
        ]

        if matrix[row][col].island_id!=island_id:
            return True

        # Check if each neighboring cell is within bounds and has a value of -1
        for neighbor_row, neighbor_col in neighbors:
            if is_valid_index(matrix, neighbor_row, neighbor_col) and matrix[neighbor_row][neighbor_col].island_id!=island_id:
                return True  # At least one neighboring cell is -1

    return False  # None of the neighboring cells are -1 or out of bounds



def getActiveNodes(graph, square_matrix, island_id, graph_obj):
    active_nodes = set()
    for node in graph.nodes:
        x, y = graph.nodes[node]['point_coord']
        if has_neighboring_minus_one(square_matrix, x, y, island_id):
            for neighbor in graph.neighbors(node):
                owners_covered = graph.edges[node, neighbor]['owners_covered']
                if owners_covered[0]!=-1 or owners_covered[1]!=-1:
                    active_nodes.add(node)
                    break
    return active_nodes


def findMinDistActiveNodes(active_nodes, shortest_path_lengths, owners_nodes_dict, min_dist_active_nodes_dict):
    for owner, corner_nodes in owners_nodes_dict.items():
        min_dist = float('inf')
        min_dist_active_node = None
        min_corner_node = None
        flag = False
        for corner_node in corner_nodes:
            for active_node in active_nodes:
                # if current corner node is the active node
                if corner_node==active_node:
                    min_dist_active_nodes_dict[owner] = (corner_node, corner_node, 0)
                    flag = True
                    break
                # if current corner node is not the active node
                dist = shortest_path_lengths[corner_node][active_node]
                if min_dist>dist:
                    min_dist = dist
                    min_corner_node = corner_node
                    min_dist_active_node = active_node
            if flag:
                break
        if not flag:
            min_dist_active_nodes_dict[owner] = (min_corner_node, min_dist_active_node, min_dist)


def getEdgeOwners(node1, node2, edge_owner_dict):
    key1 = str(node1) + "_" + str(node2)
    key2 = str(node2) + "_" + str(node1)
    if key1 in edge_owner_dict:
        return edge_owner_dict[key1]
    else:
        return edge_owner_dict[key2]


def markPath(path, min_dist_active_nodes_dict, graph, marg, edge_owner_dict):
    new_active_nodes = set()
    for i in range(len(path)-1):
        node1 = path[i]
        node2 = path[i+1]
        new_active_nodes.add(node1)
        new_active_nodes.add(node2)
        print("added new active nodes", node1, node2)
        edge = graph.edges[node1, node2]
        owners_covered = edge['owners_covered']
        for owner in owners_covered:
            print(owner, "covered")
            if owner in min_dist_active_nodes_dict:
                print("deleting", owner, "covered by", node1, node2)
                del min_dist_active_nodes_dict[owner]
        marg.append([node1, node2])
    return new_active_nodes


def updateMinDistActiveNodesDict(min_dist_active_nodes_dict, new_active_nodes, owners_nodes_dict, shortest_path_lengths):
    for owner, corner_nodes in owners_nodes_dict.items():
        if owner in min_dist_active_nodes_dict:
            min_corner_node, min_dist_active_node, min_dist = min_dist_active_nodes_dict[owner]
            for corner_node in corner_nodes:
                for active_node in new_active_nodes:
                    dist = shortest_path_lengths[corner_node][active_node]
                    if min_dist>dist:
                        min_dist = dist
                        min_dist_active_node = active_node
                        min_corner_node = corner_node
            min_dist_active_nodes_dict[owner] = (min_corner_node, min_dist_active_node, min_dist)


def printDict(mydict):
    for key, value in mydict.items():
        print(key, value)


def find_key_with_min_value(my_dict):
    if not my_dict:
        return None

    min_key = min(my_dict, key=lambda key: my_dict[key][2])
    return min_key


def findOptimalNode(graph, owners_nodes_dict, corner_node, min_dist_active_node, min_dist_active_nodes_dict, min_dist_owner, shortest_path_lengths):
    neighbor_nodes = graph.neighbors(corner_node)
    optimal_node = None
    print("we are searching for owner", min_dist_owner)
    for node in neighbor_nodes:
        print("neighbor", node)
        owners_covered = graph.edges[corner_node, node]['owners_covered']
        print("owners covered by edge", owners_covered)
        # FIXME: there is some error with owner covering of lines, there is an error on removal of second condition
        if min_dist_owner in owners_covered and node!=min_dist_active_node:
            owner1, owner2 = owners_covered[0], owners_covered[1]
            optimal_node = node
            # select the edge which covers maximum no of owners parcels
            if owner1 in min_dist_active_nodes_dict and owner2 in min_dist_active_nodes_dict:
                return optimal_node
    return optimal_node


def createAllChakMargs(island_ids, square_matrix, islands_rect_cover_dict, island_owners, island_rectangles, cell_size, xmin, ymin, crs):
    chakmarg_gdf = gpd.GeoDataFrame()
    lines_gdf = gpd.GeoDataFrame()
    for island_id in island_ids:
        print("working on island id", island_id)
        graph_obj, owners_nodes_dict, edge_owner_dict = createChakmargGraph(square_matrix, island_id, islands_rect_cover_dict, island_owners, island_rectangles)
        graph = graph_obj.graph

        # island_line_gdf = create_line_gdf(graph, island_id)
        # island_point_gdf = create_points_gdf(graph, island_id, cell_size, xmin, ymin)
        # island_line_gdf.to_file("../preference_outputs/lines.shp")
        # island_point_gdf.to_file("../preference_outputs/points.shp")

        shortest_path_lengths = dict(nx.all_pairs_shortest_path_length(graph))
        shortest_paths = dict(nx.all_pairs_shortest_path(graph))
        # print(island_id)

        active_nodes = getActiveNodes(graph, square_matrix, island_id, graph_obj)
        print("active nodes are", active_nodes)
        
        min_dist_active_nodes_dict = {}
        findMinDistActiveNodes(active_nodes, shortest_path_lengths, owners_nodes_dict, min_dist_active_nodes_dict)
        # print("here size of dict is", len(min_dist_active_nodes_dict))
        # printDict(min_dist_active_nodes_dict)
        
        island_marg = []
        while len(min_dist_active_nodes_dict)>0:
            min_dist_owner = find_key_with_min_value(min_dist_active_nodes_dict)
            corner_node, min_dist_active_node, min_dist = min_dist_active_nodes_dict[min_dist_owner]
            print("corner node ", corner_node, "active node nearest", min_dist_active_node, "distance", min_dist)
            
            optimal_node = findOptimalNode(graph, owners_nodes_dict, corner_node, min_dist_active_node, min_dist_active_nodes_dict, min_dist_owner, shortest_path_lengths)
            print("optimal node is", optimal_node)

            # FIXME: fix this, handle it properly
            if optimal_node is None:
                del min_dist_active_nodes_dict[min_dist_owner]
                continue
            
            # NOTE: path should be taken from corner to active node, not from optimal node to active(it might be different)
            min_path = shortest_paths[corner_node][min_dist_active_node]
            # FIXME: getting wrong path if same nodes are given 
            if corner_node==min_dist_active_node:
                min_path = [corner_node]
            min_path.insert(0, optimal_node)
            print("min path is", min_path)

            active_nodes = markPath(min_path, min_dist_active_nodes_dict, graph, island_marg, edge_owner_dict)
            print("new active nodes added are", active_nodes)

            print("remaining size of dict is", len(min_dist_active_nodes_dict))
            
            updateMinDistActiveNodesDict(min_dist_active_nodes_dict, active_nodes, owners_nodes_dict, shortest_path_lengths)
            print("-----------------------------------------------------")
        
        island_marg_gdf = create_marg_gdf(graph, island_marg, island_id, cell_size, xmin, ymin, crs)
        # island_marg_gdf.to_file("../preference_outputs/marg.shp")

        
        # # entryEdges = findEntryEdges(graph)
        # # islang_marg = createIslandChakmarg(graph, entryEdges)
        # island_marg_gdf = get_marg_gdf(graph, islang_marg, island_id)

        chakmarg_gdf = gpd.GeoDataFrame(pd.concat([chakmarg_gdf, island_marg_gdf], ignore_index=True), crs=island_marg_gdf.crs)
        # lines_gdf = gpd.GeoDataFrame(pd.concat([lines_gdf, island_line_gdf], ignore_index=True), crs=island_line_gdf.crs)
        chakmarg_gdf['geometry'] = chakmarg_gdf['geometry'].apply(replace_inf_with_zero)

    return chakmarg_gdf
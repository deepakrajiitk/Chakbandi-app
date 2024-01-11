# import libraries
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point, LineString
import networkx as nx
import math
import os
import pickle
import sys
sys.path.append('..') 

# import other modules
from classes.owner_class import *
from classes.rectangle_class import *
from classes.valuation_rectangle_sum_class import ValuationRectangleSum
from consolidation_module.initial_assignment import *
from consolidation_module.consolidation import *
from consolidation_module.helpers import *
from chakmarg_module.chakmarg_algo import *
# from chakmarg_module.helpers import *

# reading input parcels shp file
print("reading file...")
parcels_file_path = input("Enter the path to the parcels Shapefile: ")

folder_path = os.path.dirname(parcels_file_path)

try:
    parcels = gpd.read_file(parcels_file_path)
except FileNotFoundError:
    print(f"The file at '{file_path}' does not exist.")

parcels.to_file(folder_path+"/temp.shp")

# finding representative point of each parcel
parcels['representative_point'] = parcels['geometry'].representative_point()

# dividing parcels island wise
## graph to represent connectivity between 
print("creating graph...")
G = nx.Graph() 

## add all parcel edges as graph edges
for index, row in parcels.iterrows():
    parcel = row['geometry']
    neighbors = parcels[parcels.intersects(parcel)]['geometry']
    for neighbor in neighbors:
        if not parcel.equals(neighbor):
            G.add_edge(parcel, neighbor)
print("-------------graph created successfully-------------")

## identify islands
print("identifying all the islands")
island_id = 0
island_ids = []
for island in nx.connected_components(G):
    island_ids.append(island_id)
    for parcel in island:
        parcels.loc[parcels['geometry'] == parcel, 'island_id'] = island_id
    island_id += 1
parcels['island_id'] = parcels['island_id'].fillna(-1)
print("-----------------successfully identified all the islands--------------------------")

# create required square matrix
print("creating square matrix (it will take some time).....")
xmin, ymin, xmax, ymax = parcels.total_bounds
cell_size = 3 
xmin, ymin, xmax, ymax = int(xmin), int(ymin), int(xmax), int(ymax)

## Calculate the number of rows and columns in square grid
num_rows = math.ceil((xmax - xmin) / cell_size)
num_cols = math.ceil((ymax - ymin) / cell_size)
num_rows, num_cols

## Create a grid of squares
grid = []
for x in range(int(xmin), int(xmax), cell_size):
    for y in range(int(ymin), int(ymax), cell_size):
        square = Polygon([
            (x, y),
            (x + cell_size, y),
            (x + cell_size, y + cell_size),
            (x, y + cell_size)
        ])
        grid.append(square)
    
## Create a GeoDataFrame from the grid
grid_gdf = gpd.GeoDataFrame({'geometry': grid}, crs=parcels.crs)

## determine the v_factor of each squre in square grid
### Perform a spatial join with your original shapefile
result_gdf = gpd.sjoin(grid_gdf, parcels, how='left', predicate='intersects')

### calculate the intersection area for each square with intersecting parcel(a row might be intersecting with multiple parcels)
result_gdf['intersection_area'] = result_gdf.apply(calculate_intersection_area, parcels=parcels, axis=1)

### sort the result_gdf by the intersection area in descending order
result_gdf = result_gdf.sort_values(by='intersection_area', ascending=False)

### keep only that matching rows whose intersection area is maximum(this help to know the v_factor of square intersecting with
### multiple parcels, we take the v_factor of parcel with the maximum intersecting area)
result_gdf = result_gdf.drop_duplicates(subset='geometry', keep='first')

result_gdf['v_factor'] = result_gdf['v_factor'].fillna(0) 
result_gdf['Revenue_No'] = result_gdf['Revenue_No'].fillna(-1)
result_gdf['Revenue_No'] = result_gdf['Revenue_No'].astype(str)
result_gdf = result_gdf.sort_index()

## create final square matrix from the gdf created
square_matrix = []
index = 0
for _, row in result_gdf.iterrows():
    row_no, col_no = index/num_rows, index%num_rows
    rectangle = [row_no, col_no, row_no, col_no]
    square_matrix.append(create_rectangle(rectangle, row['island_id'], row['v_factor'], old_owner=row['Revenue_No'], owner=-1))
square_matrix = np.array(square_matrix).reshape(num_rows, num_cols)
print("----------------successfully created square matrix----------------------")

# create all owners objects and assign preference
print("creating owners object and assiging their preferences....")
all_owners = parcels['Revenue_No'].astype("str").unique()
all_owners = all_owners[all_owners!='nan']

# NOTE: preference is assigned randomly here to all the owners
owners_dict = {}
for owner in all_owners:
    total_valuation = parcels[parcels['Revenue_No']==owner]['v_factor'].sum()
    # randomly assign some preferences to the owner
    # NOTE: here preference point are float but while creating grid square we take ceil and int value of them, there they might not
    # align a bit in shp file
    preference_list = np.array(parcels.sample(n=5)['representative_point'])
    owners_dict[owner] = create_owner(owner, 0, preference_list)
print("---------------successfully created all owners object and assigned their preferences-------------------")

# find minimal rectangular cover of each island
print("finding minimal cover of each island")
islands_rect_cover_dict = find_minimal_rectangular_cover(square_matrix, island_ids)
print("-------------------successfully found minimal cover of each island--------------------------")

# create valuation rectangle class object of each to get sum of any give rectangle in that island
islands_rect_sum_dict = {}
for island_id in island_ids:
    islands_rect_sum_dict[island_id] = ValuationRectangleSum(square_matrix, island_id, islands_rect_cover_dict) 

# find total valuation of each owner according to old parcel structure
owners_dict = calculate_owners_total_valuation(all_owners, owners_dict, square_matrix)

# initial assignement according to preference of each owner
## Sort the owners based on valuation in descending order
print("doing initial assigment of rectangular parcels of the owner according to their preferences")
sorted_owners = sorted(owners_dict.items(), key=lambda x: x[1].old_valuation, reverse=True)
threshold = 0.1
for owner, owner_obj in sorted_owners:
    owner_valuation = owner_obj.old_valuation
    for preference in owner_obj.preference_list:
        x, y = preference.x, preference.y
        row = math.ceil((int(x) - int(xmin)) / cell_size)
        col = math.ceil((int(y) - int(ymin)) / cell_size)
        island_id = square_matrix[row][col].island_id
        if island_id == -1 or square_matrix[row][col].isOccupied:
            continue
        rect_sum = islands_rect_sum_dict[island_id]
        # allowed direction flag -> top, bottom, left, right
        allowed_dir = [1, 1, 1, 1]
        rectangle_coord, valuation, diff = find_optimal_rectangle(square_matrix, owner_valuation, row, col, rect_sum, island_id, islands_rect_cover_dict, allowed_dir)
        # if assigned valuation is greater than given threshold, skip it
        if diff/owner_valuation > threshold:
            continue
        owner_obj.rectangle = create_rectangle(rectangle_coord, island_id, valuation, owner=owner, isOccupied=True)
        owner_obj.curr_valuation = valuation
        markRectangleOccupied(square_matrix, rectangle_coord, owner, island_id)
        break
print("----------------------successfully completed initial assignement-------------------------")

# performing consolidation
print("performing consolidation")
owners_dict = perform_consolidation(square_matrix, island_ids, islands_rect_cover_dict, islands_rect_sum_dict, sorted_owners, owners_dict, rect_sum)
print("----------------------consolidation complete-------------------")

# coverting rectangles to geodataframe and saving them
print("getting result dataframe")
gdf2 = rectangles_to_gdf(owners_dict, cell_size, xmin, ymin, parcels.crs)
# gdf2.to_file("../output/consolidation_output_without_cut.shp")

# removing unnecessary part from gdf
result_gdf2 = overlay_gdf(gdf2, parcels, island_ids)
result_gdf2 = result_gdf2[~result_gdf2['owner'].str.startswith('no_owner')]
result_gdf2.to_file(folder_path + "/consolidation_output.shp")
print("---------------output shp file is saved successfully-----------------------")

print("finding optimal chakmarg......")
island_rectangles, island_owners = getRectangleIslandWise(owners_dict)
chakmarg_gdf = createAllChakMargs(island_ids, square_matrix, islands_rect_cover_dict, island_owners, island_rectangles, cell_size, xmin, ymin, parcels.crs)
chakmarg_gdf.to_file(folder_path + "/chakmarg.shp")
print("-----------------chakmarg shp file saved successfully-----------------------------")




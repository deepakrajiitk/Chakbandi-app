from chakmarg_algo import *
import pickle
import sys
sys.path.append('..') 
from owner_class import *
from rectangle_class import *

# Load variables
with open('../../data/pickles/variables.pkl', 'rb') as file:
    square_matrix = pickle.load(file)
    island_ids = pickle.load(file)
    owners_dict = pickle.load(file)
    islands_rect_cover_dict = pickle.load(file)

for i in owners_dict:
    print(type(owners_dict[i]))

print(square_matrix.shape)

print("finding optimal chakmarg......")
island_rectangles, island_owners = getRectangleIslandWise(owners_dict)
chakmarg_gdf = createAllChakMargs(island_ids, square_matrix, islands_rect_cover_dict, island_owners)
chakmarg_gdf.to_file("../output/chakmarg.shp")
print("-----------------chakmarg shp file saved successfully-----------------------------")
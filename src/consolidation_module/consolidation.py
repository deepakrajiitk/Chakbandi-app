import numpy as np
import sys
sys.path.append('..') 
from classes.rectangle_class import *
from classes.owner_class import *
from consolidation_module.initial_assignment import *
from shapely.geometry import Polygon, MultiPolygon, Point, LineString

# find empty space in the given square matrix and return all the rectangle with empty space
def get_empty_rectangles(square_matrix, island_id, islands_rect_cover_dict, islands_rect_sum_dict):
    min_row, min_col, max_row, max_col = islands_rect_cover_dict[island_id]
    H = max_row-min_row+1
    W = max_col-min_col+1
    visited_matrix = np.zeros((H, W))
    total_valuation = 0
    for i in range(H):
        for j in range(W):
            total_valuation += square_matrix[i+min_row][j+min_col].valuation if ((not square_matrix[i+min_row][j+min_col].isOccupied) and square_matrix[i+min_row][j+min_col].island_id == island_id) else 0
            if square_matrix[i+min_row][j+min_col].isOccupied and square_matrix[i+min_row][j+min_col].island_id==island_id:
                visited_matrix[i][j] = 1

    # print("island_id is ", island_id)
    # print("total valuation of the island is", total_valuation)
    rects = []
    curr_rectangles_val = 0

    def find_next_rect():
        found_corner = False
        rect = {'x1': 0, 'x2': H-1, 'y1': 0, 'y2': W-1}
        for i in range(H):
            for j in range(W):
                if visited_matrix[i][j] == 0:
                    rect['x1'] = i
                    rect['y1'] = j
                    found_corner = True
                    break
            if found_corner:
                break
        
        for i in range(rect['x1'], rect['x2'] + 1):
            if visited_matrix[i][rect['y1']] != 0:
                rect['x2'] = i - 1
                break
            for j in range(rect['y1'], rect['y2'] + 1):
                if visited_matrix[i][j] != 0:
                    rect['y2'] = j - 1
                    break
        return rect

    def mark_rect(rect):
        for i in range(rect['x1'], rect['x2'] + 1):
            for j in range(rect['y1'], rect['y2'] + 1):
                visited_matrix[i][j] = 1

    while round(curr_rectangles_val, 2) < round(total_valuation, 2):
        rect = find_next_rect()
        if rect == None:
            break
        rectangle_coord = (rect['x1']+min_row, rect['y1']+min_col, rect['x2']+min_row, rect['y2']+min_col)
        rect_val = islands_rect_sum_dict[island_id].get_sum(rectangle_coord)
        if round(rect_val, 2)>0:
            rectangle = create_rectangle(rectangle_coord, island_id, rect_val, owner=None, isOccupied=False)
            rects.append(rectangle)
        mark_rect(rect)
        curr_rectangles_val += rect_val
        # print(rectangle_coord, curr_rectangles_val)

    return rects
                
# find the most optimal rectangle from empty rectangle and assign it to owner
def assign_rectangles(square_matrix, sorted_owners, rectangles, threshold, rect_sum, islands_rect_cover_dict, islands_rect_sum_dict):
    assigned = 0
    for owner, owner_obj in sorted_owners:
        # check if rectangle is already assigned to owner or not
        if owner_obj.rectangle==None:
            min_diff = float('inf')
            closest_rectangle = None
            for rectangle in rectangles:
                # make sure candidate rectangle is not occupied
                if (not rectangle.isOccupied):
                    diff = abs(rectangle.valuation-owner_obj.old_valuation)
                    if diff<min_diff:
                        closest_rectangle = rectangle
                        min_diff = diff
            if min_diff/owner_obj.old_valuation < threshold:
                # if threshold is greater than normal, we get the create another rectangle inside that rectangle to fit it
                if threshold>0.1:
                    # FIXME: this will give error if starting coordinate belong to island_id=-1
                    start_row, start_col = closest_rectangle.coordinates[0], closest_rectangle.coordinates[1]
                    rect_sum = islands_rect_sum_dict[closest_rectangle.island_id]
                    # dir -> top, bottom, left, right
                    allowed_dir = [1, 1, 1, 1]
                    rectangle_coord, rect_valuation, diff = find_optimal_rectangle(square_matrix, owner_obj.old_valuation, start_row, start_col, rect_sum, closest_rectangle.island_id, islands_rect_cover_dict, allowed_dir)
                    closest_rectangle = create_rectangle(rectangle_coord, closest_rectangle.island_id, rect_valuation, owner=owner, isOccupied=True)
                assigned+=1
                # print(closest_rectangle.valuation, owner_obj.old_valuation)
                closest_rectangle.owner = owner
                closest_rectangle.isOccupied = True
                owner_obj.rectangle = closest_rectangle
                owner_obj.curr_valuation = closest_rectangle.valuation
                markRectangleOccupied(square_matrix, closest_rectangle.coordinates, closest_rectangle.owner, closest_rectangle.island_id)
    return assigned

# FIXME: Getting some unwanted overlapping small rectangles
def perform_consolidation(square_matrix, island_ids, islands_rect_cover_dict, islands_rect_sum_dict, sorted_owners, owners_dict, rect_sum, threshold=0.1):
    unassigned_owner_count = sum(value.rectangle is None for value in owners_dict.values())
    print("unassigned owner count", unassigned_owner_count)
    unassigned_rectangles = []
    while unassigned_owner_count>0:
        all_rectangles = []
        for island_id in island_ids:
            island_rectangles = get_empty_rectangles(square_matrix, island_id, islands_rect_cover_dict, islands_rect_sum_dict)
            all_rectangles.extend(island_rectangles)
        print("no of rectangles", len(all_rectangles))
        assigned_count = assign_rectangles(square_matrix, sorted_owners, all_rectangles, threshold, rect_sum, islands_rect_cover_dict, islands_rect_sum_dict)
        temp = 0.1
        while assigned_count == 0:
            print("could not assign rectangle to anyone, changing threshold to", temp+threshold)
            assigned_count = assign_rectangles(square_matrix, sorted_owners, all_rectangles, threshold+temp, rect_sum, islands_rect_cover_dict, islands_rect_sum_dict)
            temp+=0.1
            if temp>1:
                break
        print("assigned", assigned_count)
        unassigned_owner_count -= assigned_count
        print("remaining", unassigned_owner_count)
        # if no further assignment is possible or all have been assigned
        if assigned_count==0 or unassigned_owner_count==0:
            for island_id in island_ids:
                island_rectangles = get_empty_rectangles(square_matrix, island_id, islands_rect_cover_dict, islands_rect_sum_dict)
                unassigned_rectangles.extend(island_rectangles)
            break
        
    # marking remaining rectangles as parcel with no_owner
    temp = 1
    for rectangle in unassigned_rectangles:
        if rectangle.valuation!=0:
            owner_name = "no_owner_"+str(temp)
            rectangle.owner = owner_name
            owner_obj = create_owner(owner_name, 0, [], curr_valuation=0, rectangle=rectangle)
            owners_dict[owner_name] = owner_obj
            markRectangleOccupied(square_matrix, rectangle.coordinates, owner_name, rectangle.island_id)
            temp+=1
    
    return owners_dict
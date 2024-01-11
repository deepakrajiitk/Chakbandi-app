import sys
sys.path.append('..') 
from classes.rectangle_class import *

# return candidate having minimal absolute difference between org valuation and rectangle valuation
def find_min_difference_candidate(candidate_names, candidate_values, target_value):
    candidates = zip(candidate_names, candidate_values)
    min_difference_candidate = min(
        candidates,
        key=lambda candidate: abs(candidate[1] - target_value)
    )
    return min_difference_candidate

# Return the best non-overlapping rectangle candidate
def get_best_candidate(square_matrix, top_candidate, bottom_candidate, left_candidate, right_candidate, total_valuation, island_id, islands_rect_cover_dict, rect_sum):
    top_value, bottom_value, left_value, right_value = 0, 0, 0, 0
    if isValidRectangle(top_candidate, island_id, islands_rect_cover_dict) and isRectangleOccupied(top_candidate, square_matrix, "top", island_id):
        top_value = rect_sum.get_sum(top_candidate)
    if isValidRectangle(bottom_candidate, island_id, islands_rect_cover_dict) and isRectangleOccupied(bottom_candidate, square_matrix, "bottom", island_id):
        bottom_value = rect_sum.get_sum(bottom_candidate)
    if isValidRectangle(left_candidate, island_id, islands_rect_cover_dict) and isRectangleOccupied(left_candidate, square_matrix, "left", island_id):
        left_value = rect_sum.get_sum(left_candidate)
    if isValidRectangle(right_candidate, island_id, islands_rect_cover_dict) and isRectangleOccupied(right_candidate, square_matrix, "right", island_id):
        right_value = rect_sum.get_sum(right_candidate)

    min_difference_candidate = find_min_difference_candidate([top_candidate, bottom_candidate, left_candidate, right_candidate], [top_value, bottom_value, left_value, right_value], total_valuation) 

    best_candidate = min_difference_candidate[0]
    best_value = min_difference_candidate[1]   
    
    return best_candidate, best_value

# find the most optimal rectangle for the given owner in the preferred region
def find_optimal_rectangle(square_matrix, total_valuation, start_row, start_col, rect_sum, island_id, islands_rect_cover_dict, allowed_dir):
    # rectangle coordinates = x1, y1 => left top  x2, y2 => bottom right
    rectangle = [start_row, start_col, start_row, start_col]
    curr_valuaton = rect_sum.get_sum(rectangle)
    best_diff = abs(total_valuation - rect_sum.get_sum(rectangle))
    while True:
        top_candidate = (rectangle[0]-1, rectangle[1], rectangle[2], rectangle[3]) if allowed_dir[0] else None
        bottom_candidate = (rectangle[0], rectangle[1], rectangle[2]+1, rectangle[3]) if allowed_dir[1] else None
        left_candidate = (rectangle[0], rectangle[1]-1, rectangle[2], rectangle[3]) if allowed_dir[2] else None
        right_candidate = (rectangle[0], rectangle[1], rectangle[2], rectangle[3]+1) if allowed_dir[3] else None
        
        best_candidate, candidate_valuation = get_best_candidate(square_matrix, top_candidate, bottom_candidate, left_candidate, right_candidate, total_valuation, island_id, islands_rect_cover_dict, rect_sum)
        curr_diff = abs(total_valuation-candidate_valuation)
        # print(best_candidate, candidate_valuation, curr_diff, best_diff)
        if best_diff<=curr_diff:
            break
        best_diff = curr_diff
        curr_valuaton = candidate_valuation
        rectangle = best_candidate
        
    return rectangle, curr_valuaton, best_diff
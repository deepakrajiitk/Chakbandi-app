import numpy as np

class OwnerClass:
    def __init__(self, owner_name, old_valuation, preference_list, curr_valuation, rectangle):
        self.owner_name = owner_name
        self.old_valuation = old_valuation
        self.preference_list = preference_list
        self.curr_valuation = curr_valuation
        self.rectangle = rectangle

def create_owner(owner_name, old_valuation, preference_list, curr_valuation=0, rectangle=None):
    return OwnerClass(owner_name, old_valuation, preference_list, curr_valuation, rectangle)

# find total valuation of each owner according to old parcel structure
def calculate_owners_total_valuation(all_owners, owners_dict, square_matrix):
    valuation_matrix = np.array([[obj.valuation for obj in row] for row in square_matrix])
    owner_matrix = np.array([[obj.old_owner for obj in row] for row in square_matrix])
    
    for owner in all_owners:
        owner_obj = owners_dict[owner]
        indices = np.where(owner_matrix == owner)
        total_valuation = np.sum(valuation_matrix[indices])
        owner_obj.old_valuation = total_valuation
    
    return owners_dict

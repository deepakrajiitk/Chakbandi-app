from shapely.geometry import Polygon, MultiPolygon, Point, LineString

class RectangleClass:
    def __init__(self, coordinates, island_id, valuation, old_owner, owner, isOccupied):
        # coordinates are according to square matrix not actual polygon coodinates
        self.coordinates = coordinates
        self.island_id = island_id
        self.old_owner = old_owner
        self.owner = owner
        self.valuation = valuation
        self.isOccupied = isOccupied

def create_rectangle(coordinates, island_id, valuation, old_owner=None, owner=None, isOccupied=0):
    return RectangleClass(coordinates, island_id, valuation, old_owner, owner, isOccupied)

# Check if rectangle is out of bound or out of island cover
def isValidRectangle(rectangle, island_id, islands_rect_cover_dict):
    if rectangle is None:
        return False
    x1, y1, x2, y2 = rectangle[0], rectangle[1], rectangle[2], rectangle[3]
    min_row, min_col, max_row, max_col = islands_rect_cover_dict[island_id]
    if x1<min_row or x2<min_row or x1>max_row or x2>max_row or y1<min_col or y2<min_col or y1>max_col or y2>max_col:
        return False
    return True

# Check if given rectangle is overlapping with other rectangle or not
# taken part as argument to do the minimal number of cells checking
def isRectangleOccupied(rectangle, square_matrix, part, island_id):
    x1, y1, x2, y2 = rectangle[0], rectangle[1], rectangle[2], rectangle[3]
    if part == "top":
        for i in range(y1, y2+1):
            # first check if for overlapping
            # second check is to check if the cell is of other island then we can use it even if it is occupied(because its valuation will not be taken in to account, only those cell in that line which are in same island will be taken in account)
            if square_matrix[x1][i].isOccupied and square_matrix[x1][i].island_id==island_id:
                return False
    if part == "bottom":
        for i in range(y1, y2+1):
            if square_matrix[x2][i].isOccupied and square_matrix[x2][i].island_id==island_id:
                return False
    if part == "left":
        for i in range(x1, x2+1):
            if square_matrix[i][y1].isOccupied and square_matrix[i][y1].island_id==island_id:
                return False
    if part == "right":
        for i in range(x1, x2+1):
            if square_matrix[i][y2].isOccupied and square_matrix[i][y2].island_id==island_id:
                return False
    return True

# mark all the cells of assigned rectagle as occupied
# TODO: Find someway to mark rectangle optimally
def markRectangleOccupied(square_matrix, rectangle, owner, island_id):
    x1, y1, x2, y2 = rectangle[0], rectangle[1], rectangle[2], rectangle[3]
    for i in range(x1, x2+1):
        for j in range(y1, y2+1):
            # mark cell only of the same island
            if square_matrix[i][j].island_id == island_id:
                square_matrix[i][j].isOccupied = True
                # RECHECK: This might create some problems
                square_matrix[i][j].owner = owner

def save_rectangles_as_shp(owners_dict, include_island_only = None):
    rows = []
    for owner, owner_obj in owners_dict.items():
        if owner_obj.rectangle:
            rectangle = owner_obj.rectangle.coordinates
            if (include_island_only is not None) and owner_obj.rectangle.island_id!=include_island_only:
                continue
            x1, y1, x2, y2 = rectangle[0], rectangle[1], rectangle[2], rectangle[3]
            x2+=1
            y2+=1
            coordinates = [(x1*cell_size + xmin, y1*cell_size + ymin), (x2*cell_size + xmin, y1*cell_size + ymin), (x2*cell_size + xmin, y2*cell_size + ymin), (x1*cell_size + xmin, y2*cell_size + ymin)]
            is_non_finite_flag = False
            for coord in coordinates:
                if isNonFinite(coord):
                    is_non_finite_flag = True
                    break
            if not is_non_finite_flag:
                polygon = Polygon(coordinates)
                rows.append([polygon, owner_obj.curr_valuation, owner_obj.old_valuation, owner, owner_obj.rectangle.island_id])

    gdf = gpd.GeoDataFrame(rows, columns=["geometry", "valuation", "owner_old_valuation", "owner", "island_id"], crs = parcels.crs)
    return gdf
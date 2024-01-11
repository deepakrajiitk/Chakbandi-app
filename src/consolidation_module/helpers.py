import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon, Point, LineString


# Custom function to calculate intersection area
def calculate_intersection_area(row, parcels):
    if(not np.isnan(row['index_right'])):
        grid_geom = row['geometry']
        parcel_geom = parcels.iloc[int(row['index_right'])]['geometry']
        intersection = grid_geom.intersection(parcel_geom)
        return intersection.area
    return 0.0

# function return minimal rectangle coordinates which cover the given island
def find_minimal_rectangular_cover(square_matrix, island_ids):
    result = {}

    island_id_matrix = np.array([[obj.island_id for obj in row] for row in square_matrix])

    for island_id in island_ids:
        indices = np.where(island_id_matrix == island_id)
        min_row, min_col = np.min(indices, axis=1)
        max_row, max_col = np.max(indices, axis=1)
        # NOTE: I have change return format, this might create problem
        result[island_id] = (min_row, min_col, max_row, max_col)

    return result

def isNonFinite(coord):
    x, y = coord
    if x==float('inf') or x==float('-inf') or y==float('inf') or y==float('-inf'):
        return True
    return False 

def rectangles_to_gdf(owners_dict, cell_size, xmin, ymin, crs, include_island_only = None):
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

    gdf = gpd.GeoDataFrame(rows, columns=["geometry", "valuation", "owner_old_valuation", "owner", "island_id"], crs = crs)
    return gdf

# function to handle non-finite values
def replace_inf_with_zero(geom):
    if geom.has_z:
        if geom.geom_type == 'LineString':
            coordinates = list(geom.coords)
            updated_coordinate = []
            for coord in coordinates:
                x, y, z  = coord
                if z=="-inf":
                    print("--------------------")
                if z==float('-inf') or z==float('inf'):
                    z = 0
                if x==float('-inf') or x==float('inf'):
                    continue
                if y==float('-inf') or y==float('inf'):
                    continue
                updated_coordinate.append((x, y, z))
            return LineString(updated_coordinate)

        if geom.geom_type == 'Polygon':
            coordinates = list(geom.exterior.coords)
            updated_coordinate = []
            for coord in coordinates:
                x, y, z  = coord
                if z==float('-inf') or z==float('inf'):
                    z = 0
                if x==float('-inf') or x==float('inf'):
                    continue
                if y==float('-inf') or y==float('inf'):
                    continue
                updated_coordinate.append((x, y, z))
            return Polygon(updated_coordinate)

        elif geom.geom_type == 'MultiPolygon':
            polygons = [replace_inf_with_zero(polygon) for polygon in geom.geoms]
            return MultiPolygon(polygons)
    return geom

def overlay_gdf(gdf1, gdf2, island_ids):
    result_gdf = gpd.GeoDataFrame(columns=['geometry'])
    for island_id in island_ids:
        island_gdf1 = gdf1[gdf1['island_id'] == island_id]
        island_gdf2 = gdf2[gdf2['island_id'] == island_id].unary_union
        island_gdf2 = gpd.GeoDataFrame(geometry=[island_gdf2]).set_crs(island_gdf1.crs)
        overlay_gdf = gpd.overlay(island_gdf1, island_gdf2, how="intersection")
        result_gdf = result_gdf._append(overlay_gdf, ignore_index=True)
    
    result_gdf['geometry'] = result_gdf['geometry'].apply(replace_inf_with_zero)
    result_gdf = result_gdf[result_gdf['geometry'].notna()]
    
    return result_gdf

# Define a function to replace -inf with 0 in the 'z' coordinates
def replace_inf_with_zero(geom):
    if geom.has_z:
        if geom.geom_type == 'LineString':
            coordinates = list(geom.coords)
            updated_coordinate = []
            for coord in coordinates:
                x, y, z  = coord
                if z=="-inf":
                    print("--------------------")
                if z==float('-inf') or z==float('inf'):
                    z = 0
                if x==float('-inf') or x==float('inf'):
                    continue
                if y==float('-inf') or y==float('inf'):
                    continue
                updated_coordinate.append((x, y, z))
            return LineString(updated_coordinate)

        if geom.geom_type == 'Polygon':
            coordinates = list(geom.exterior.coords)
            updated_coordinate = []
            for coord in coordinates:
                x, y, z  = coord
                if z==float('-inf') or z==float('inf'):
                    z = 0
                if x==float('-inf') or x==float('inf'):
                    continue
                if y==float('-inf') or y==float('inf'):
                    continue
                updated_coordinate.append((x, y, z))
            return Polygon(updated_coordinate)

        elif geom.geom_type == 'MultiPolygon':
            polygons = [replace_inf_with_zero(polygon) for polygon in geom.geoms]
            return MultiPolygon(polygons)
    return geom
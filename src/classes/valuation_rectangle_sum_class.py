import numpy as np

# Class created to get the rectangle sum in minimum time. This work separately for each island.
# Different object will be created for each island
class ValuationRectangleSum:
    def __init__(self, square_matrix, island_id, islands_rect_cover_dict):
        self.min_row, self.min_col, self.max_row, self.max_col = islands_rect_cover_dict[island_id]
        num_rows, num_cols = self.max_row-self.min_row+1, self.max_col-self.min_col+1
        self.prefix_sum = np.zeros((num_rows+1, num_cols+1))
        for i in range(num_rows):
            for j in range(num_cols):
                valuation = 0 if square_matrix[i+self.min_row][j+self.min_col].island_id!=island_id else square_matrix[i+self.min_row][j+self.min_col].valuation
                self.prefix_sum[i+1][j+1] = (valuation + 
                    self.prefix_sum[i][j + 1] +
                    self.prefix_sum[i + 1][j] -
                    self.prefix_sum[i][j]
                )

    # Return any rectangle sum in square matrix in O(1)
    def get_sum(self, rectangle):
        x1, y1, x2, y2 = rectangle[0]-self.min_row, rectangle[1]-self.min_col, rectangle[2]-self.min_row, rectangle[3]-self.min_col
        x2 += 1
        y2 += 1
        return self.prefix_sum[x2][y2] - self.prefix_sum[x1][y2] - self.prefix_sum[x2][y1] + self.prefix_sum[x1][y1]
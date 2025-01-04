import numpy as np
import numba
from numba import njit, prange
import time

@njit
def check_aggregate(grid, i, j):
    size = grid.shape
    if grid[i, j] in (1, 2, 3):  # Check if the cell is non-aggregated A, B, or C
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:  # Check neighbors
            ni, nj = i + di, j + dj
            if 0 <= ni < size[0] and 0 <= nj < size[1]:

                if grid[ni, nj] == 0:  # Skip if the neighbor is empty
                    continue

                # If A and B are adjacent, form a C molecule
                elif (grid[i, j] == 1 and grid[ni, nj] == 2) or (grid[i, j] == 2 and grid[ni, nj] == 1):
                    # Reset both cells
                    grid[i, j] = 0
                    grid[ni, nj] = 0

                    # Create a C molecule at random between the two cells
                    if np.random.rand() > 0.5:
                        grid[i, j] = 3  # C molecule
                    else:
                        grid[ni, nj] = 3  # C molecule
                    return

                # For all aggregation scenarios: A-A, B-B, C-C, A-C, B-C, aggregate both
                elif (grid[i, j] in (1, 2, 3) and grid[ni, nj] in (1, 2, 3)):
                    # Aggregating the cells
                    grid[i, j] += 10  # Convert current to aggregated form
                    if grid[ni, nj]<=10: grid[ni, nj] += 10  # Convert neighbor to aggregated form if not already
                    return  # After aggregation, we exit since it’s handled

                







@njit(parallel=True)
def step_parallel(grid, size, ndif_ratio, stepNo):
    updated = np.zeros(size, dtype=np.bool_)

    # Define offsets for the passes
    pass_offsets = [
        (0, 0),  # Pass 1 (Top-left)
        (0, 1),  # Pass 2 (Top-center)
        (0, 2),  # Pass 3 (Top-right)
        (1, 0),  # Pass 4 (Center-left)
        (1, 1),  # Pass 5 (Center)
        (1, 2),  # Pass 6 (Center-right)
        (2, 0),  # Pass 7 (Bottom-left)
        (2, 1),  # Pass 8 (Bottom-center)
        (2, 2)   # Pass 9 (Bottom-right)
    ]

    # Precompute the list of (i, j) pairs for each pass
    indices = []
    for offset in pass_offsets:
        i_offset, j_offset = offset
        for i in range(i_offset, size[0]):
            for j in range(j_offset, size[1]):
                if (i % 3 == i_offset) and (j % 3 == j_offset):
                    indices.append((i, j))

    # Convert to NumPy array for prange
    indices = np.array(indices)

    # Parallel loop over the precomputed indices
    for index in prange(len(indices)):
        i, j = indices[index]
        cell = grid[i, j]

        # Skip if already updated in this pass or if cell is empty/aggregated
        if updated[i, j] or cell == 0 or cell >= 10:
            continue

        # Random axis selection using np.random.rand()
        if np.random.rand() < ndif_ratio:  # 0 for x-axis
            axis = 0
        else:
            axis = 1

        # Random direction
        value = 1 if np.random.rand() > 0.5 else -1

        new_i, new_j = i, j
        if axis == 0:
            new_i = i + value
        else:
            new_j = j + value

        # Ensure new_i and new_j are integers
        new_i = int(new_i)
        new_j = int(new_j)

        # Check for bounds
        if new_i < 0 or new_i >= size[0] or new_j < 0 or new_j >= size[1]:
            continue

        # Swap current cell with the new position if it's empty, else don't do anything
        if grid[new_i, new_j] == 0 and not updated[new_i, new_j]:
            grid[new_i, new_j], grid[i, j] = grid[i, j], grid[new_i, new_j]
            updated[new_i, new_j], updated[i, j] = True, True

            # Check aggregation around the new cell
            check_aggregate(grid, i, j)


@njit(parallel=True)
def num_monomers_parallel(grid):
    count = 0
    for r in prange(grid.shape[0]):
        for cell in grid[r]:
            if cell in (1, 2, 3):  # Count only non-aggregated monomers
                count += 1
    return count




@njit
def find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]  # Path compression
        x = parent[x]
    return x

@njit
def union(parent, size, rank, x, y):
    rootX = find(parent, x)
    rootY = find(parent, y)
    
    if rootX != rootY:
        if rank[rootX] > rank[rootY]:
            parent[rootY] = rootX
            size[rootX] += size[rootY]  # Update size of rootX
        elif rank[rootX] < rank[rootY]:
            parent[rootX] = rootY
            size[rootY] += size[rootX]  # Update size of rootY
        else:
            parent[rootY] = rootX
            size[rootX] += size[rootY]  # Update size of rootX
            rank[rootX] += 1

@njit(parallel=True)
def num_islands_parallel(grid):
    rows, cols = grid.shape
    parent = np.arange(rows * cols)  # Each cell is its own parent
    size = np.ones(rows * cols, dtype=np.int32)  # Initialize sizes to 1
    rank = np.zeros(rows * cols, dtype=np.int32)

    # Iterate over the grid to perform union operations
    for r in prange(rows):
        for c in range(cols):
            if grid[r, c] == 1:
                if r + 1 < rows and grid[r + 1, c] == 1:
                    union(parent, size, rank, r * cols + c, (r + 1) * cols + c)
                if c + 1 < cols and grid[r, c + 1] == 1:
                    union(parent, size, rank, r * cols + c, r * cols + (c + 1))

    # Count distinct roots (islands) and their sizes
    root_set = set()
    island_sizes = {}
    
    for r in range(rows):
        for c in range(cols):
            if grid[r, c] == 1:
                root = find(parent, r * cols + c)
                root_set.add(root)
                if root in island_sizes:
                    island_sizes[root] += 1
                else:
                    island_sizes[root] = 1

    return len(root_set), [island_sizes[root] for root in root_set]  # Return number of islands and sizes




class Simulation:
    def __init__(self, size, ndif_ratio, debug=False):
        print("Sim v 0.2 (Numba Parallelized)")

        self.size = size
        self.ndif_ratio = ndif_ratio
        self.stepNo = 1
        self.grid = np.zeros(size, dtype=np.int32)  # Integer grid

    def Step(self):
        self.stepNo += 1
        step_parallel(self.grid, self.size, self.ndif_ratio, self.stepNo)

    def Deposit(self, type, x_center=None, y_center=None, x_side=None, y_side=None):
        """ 
        Deposits a Monomer in a restricted area of the grid.
        If x_center, y_center, x_side, and y_side are not provided, the entire grid is used.
        """
        if x_center is None:
            x_center = self.grid.shape[0] // 2
        if y_center is None:
            y_center = self.grid.shape[1] // 2

        if x_side is None:
            x_side = self.grid.shape[0]
        if y_side is None:
            y_side = self.grid.shape[1]

        i = int(np.random.random() * x_side + x_center - x_side / 2)
        j = int(np.random.random() * y_side + y_center - y_side / 2)

        if self.grid[i, j] == 0:  # Only place if empty
            match type:
                case "A":
                    self.grid[i, j] = 1
                
                case "B":
                    self.grid[i, j] = 2 
                
                case "C":
                    self.grid[i, j] = 3

                case _:
                    raise ValueError 
                    
                    

    def Print(self):
        for i in range(self.size[0]):
            row = []
            for j in range(self.size[1]):
                if self.grid[i, j] == 0:
                    row.append(" ")
                elif self.grid[i, j] == 1:
                    row.append("\033[33m█\033[0m")  # Yellow for non-aggregated A
                elif self.grid[i, j] == 2:
                    row.append("\033[36m█\033[0m")  # Cyan for non-aggregated B
                elif self.grid[i, j] == 3:
                    row.append("\033[37m█\033[0m")  # IDK for non-aggregated C
                elif self.grid[i, j] == 11:
                    row.append("\033[31m█\033[0m")  # Red for aggregated A
                elif self.grid[i, j] == 12:
                    row.append("\033[32m█\033[0m")  # Green for aggregated B
                elif self.grid[i, j] == 13:
                    row.append("\033[30m█\033[0m")  # IDK for aggregated C
            print("".join(row))
    
    def NumMonomers(self):
        return num_monomers_parallel(self.grid)


    def NumIslands(self):
        # Clone the grid for the island detection
        grid = (self.grid >= 10).astype(np.int32)  # 1 for aggregated monomers, 0 otherwise
        island_count, cells_per_island = num_islands_parallel(grid)
        return island_count, cells_per_island



@njit(parallel=True)
def parallel_sum_and_count(cells):
    total_cells = 0
    count = 0
    for i in prange(len(cells)):
        if cells[i] > 0:
            total_cells += cells[i]
            count += 1
    return total_cells, count



if __name__ == "__main__":
    import os
    import time
    import matplotlib.pyplot as plt

    island_cellEvo = []
    island_numEvo = []

    size = (1000, 1000)
    sim = Simulation(size, 0.5)

    steps = 40000
    Ndif_A = 10
    Ndif_B = 15
    coverage_limit = 0.2  # stop limiter

    sim_cells = size[0] * size[1]
    i = 0
    while True:
        i += 1
        if i % Ndif_A == 0:
            sim.Deposit("A")
        if i % Ndif_B == 0:
            sim.Deposit("B")

        start_time_step = time.time()
        sim.Step()
        elapsed_time_step = (time.time() - start_time_step) * 1000  # Convert to ms

        start_time_islands = time.time()
        isl, cells = sim.NumIslands()
        # isl, cells = 0,[0]
        elapsed_time_islands = (time.time() - start_time_islands) * 1000  # Convert to ms

        # Compute the average cells per island
        if i % 1 == 0:
            start_time_sum = time.time()
            if isl > 0:
                total_cells, count = parallel_sum_and_count(cells)
                if count > 0:
                    island_cellEvo.append(total_cells / count)
                else:
                    island_cellEvo.append(0)
            else:
                island_cellEvo.append(0)
            elapsed_time_sum = (time.time() - start_time_sum) * 1000  # Convert to ms

            cells_num = sum(cells)
            fill_ratio = cells_num / sim_cells
            
            time.sleep(1/25)
            os.system("clear")
            sim.Print()
            print(f'SIM {i}\tfill {cells_num}/{sim_cells} {fill_ratio*100:.4f}%\t'
                  f'Time (Step): {elapsed_time_step:.4f} ms\t'
                  f'Time (Islands): {elapsed_time_islands:.4f} ms\t'
                  f'Time (Sum): {elapsed_time_sum:.4f} ms')
            # print(i)

            if fill_ratio >= coverage_limit:
                break

            # input()


    plt.plot(island_numEvo, label="number of islands")
    plt.plot(island_cellEvo, label="avg cells per islands")
    plt.legend()
    plt.show()

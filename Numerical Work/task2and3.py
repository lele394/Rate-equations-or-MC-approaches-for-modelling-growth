import numpy as np
import random as rd
# from numba import njit

class Monomer:
    def __init__(self, type, moves=0):

        self.type = type

        self.moves = moves
        self.aggregated = False

        self.lastMoved = 0 # Index to skip moving a monomer twice, linked to sim

    def aggregate(self, position=0):
        # print(f'Aggregating {position}')
        self.aggregated = True


class Simulation:
    def __init__(self, size, ndif_ratio, debug=False):
        """
        
        ndif_ratio : chance to go on x, inversely 1-ndif_ratio is the chance to go on y. Ratio is a wrong term, might fix it one day
        
        """
        print("Sim v 0.2")

        self.debug = debug
        
        self.size = size
        # self.rate = rate

        self.ndif_ratio = ndif_ratio

        self.stepNo = 1

        # Initialize the grid 
        self.grid = np.zeros(size, dtype=object)



    def Step(self):
        self.stepNo += 1
        # Loop through the entire grid. No caching was added so heck ye we go hard on CPU
        for (i, j), cell in np.ndenumerate(self.grid):
            if type(cell) == Monomer and not cell.aggregated and cell.lastMoved != self.stepNo: # Removed check for remaining moves on Monomers
                if self.debug: print(f'Cell {[i,j]}, aggregated {cell.aggregated} moves {cell.moves} last moved {cell.lastMoved} simStep {self.stepNo}')
                # Update Monomer position if cell is not aggregated
                axis, value = np.random.choice(["i", "j"], p=[self.ndif_ratio, 1.0-self.ndif_ratio]), np.random.choice([-1, 1])

                # Update new cell
                try:
                    match axis:
                        case "i":
                            self.grid[i, j].moves -= 1
                            self.grid[i, j].lastMoved = self.stepNo
                            self.grid[i+value, j], self.grid[i, j] = self.grid[i, j], self.grid[i+value, j]
                            new_i, new_j = i+value, j 
                            
                        case "j":
                            self.grid[i, j].moves -= 1
                            self.grid[i, j].lastMoved = self.stepNo
                            self.grid[i, j+value], self.grid[i, j] = self.grid[i, j], self.grid[i, j+value]
                            new_i, new_j = i, j+value
                            
                    # Now check for aggregation and fix both monomer if they are
                    for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                        try:# try catch for out fof bound handling
                            if type(self.grid[new_i + di, new_j + dj]) == Monomer:
                                self.grid[new_i, new_j].aggregate([new_i, new_j])
                                self.grid[new_i + di, new_j + dj].aggregate([new_i + di, new_j + dj])
                        
                        # If we're out of the box
                        except IndexError as e:
                            if self.debug:print(f'Aggregation test error {e}')
                            

                # What to do if out of bounds
                except IndexError as e:
                    if self.debug:print(f'Killing {(i,j)} since Out of bound error {e}')
                    self.grid[i, j] = 0 # just dies if goes out




    def Deposit(self, type, x_center=None, y_center=None, x_side=None, y_side=None):
        """ 
        Deposits a Monomer in a restricted area of the grid.
        If x_center, y_center, x_side, and y_side are not provided, the entire grid is used.
        """
        # If x_center or y_center are not provided, use the center of the grid
        if x_center is None:
            x_center = self.grid.shape[0] / 2
        if y_center is None:
            y_center = self.grid.shape[1] / 2

        # If x_side or y_side are not provided, use the entire grid side lengths
        if x_side is None:
            x_side = self.grid.shape[0]
        if y_side is None:
            y_side = self.grid.shape[1]

        # Calculate the bounds of the restricted area and cast to int
        # x_max = x_center + x_side / 2
        # y_max = y_center + y_side / 2

        # Randomly choose indices within the restricted area
        i = int(np.random.random() * x_side+x_center - x_side/2)
        j = int(np.random.random() * y_side+y_center - y_side/2)

        # Place a Monomer object at the randomly selected (i, j) position
        self.grid[i, j] = Monomer(type)




    def Print(self):
        """

        
        
        """
        output = []  # List to store the grid as a string
        for i in range(self.size[0]):
            row = []
            for j in range(self.size[1]):
                if self.grid[i, j] == 0:
                    row.append(" ")  # Add a space for 0
                else:
                    if self.grid[i, j].aggregated:
                        match self.grid[i, j].type:
                            case "A":
                                row.append("\033[31m█\033[0m")  # Red for aggregated
                            case "B":
                                row.append("\033[32m█\033[0m")  # Red for aggregated
                    else:
                        match self.grid[i, j].type:
                            case "A":
                                row.append("\033[33m█\033[0m")  # Red for aggregated
                            case "B":
                                row.append("\033[36m█\033[0m")  # Red for aggregated
            output.append("".join(row))  # Join the row into a string and add it to the output list
        print("\n".join(output))  # Print the entire grid as a single string


    def NumMonomers(self):
        # Clone of the grid, could use true/false instead of strings but oh well
        grid = [
                    [0 if isinstance(x, str) or x == 0 or getattr(x, 'aggregated', False) else 1 for x in row]
                    for row in self.grid
                ]
        return np.sum(grid)


    def NumIslands(self):
        # Clone of the grid, could use true/false instead of strings but oh well
        grid = [
                    ['0' if isinstance(x, str) or x == 0 or not getattr(x, 'aggregated', False) else '1' for x in row]
                    for row in self.grid
                ]
        


        
        rows, cols = len(grid), len(grid[0])
        
        def dfs(r, c):
            """Depth first search algorithm to count the number of islands"""
            # Boundary and check if it's land
            if r < 0 or c < 0 or r >= rows or c >= cols or grid[r][c] == '0':
                return 1

            
            # Mark the land as visited by setting it to '0'
            grid[r][c] = '0'
            
            s =0 #counter for cell number

            # Visit all adjacent cells (up, down, left, right)
            s += dfs(r + 1, c)  # down
            s += dfs(r - 1, c)  # up
            s += dfs(r, c + 1)  # right
            s += dfs(r, c - 1)  # left
            return s
        
        island_count = 0

        cells_per_island = []
        
        # Traverse every cell in the grid
        for r in range(rows):
            for c in range(cols):
                # Start a DFS if we find an unvisited land cell
                if grid[r][c] == '1':
                    cells_per_island.append(dfs(r, c))
                    island_count += 1  # Increase the island count after finishing the DFS
        
        return island_count, cells_per_island
    

    def NumMonomers(self):
        count = 0
        for r in self.grid:
            for cell in r:
                if cell != 0:
                    if not cell.aggregated:
                        count += 1
        return count



# ANSI escape codes for text colors and backgrounds
def print_ansi_table():
    reset = "\033[0m"
    
    print("ANSI Escape Code Table:")
    print("\nForeground Colors:")
    
    # Foreground colors
    for i in range(30, 38):  # Standard colors
        print(f"\033[{i}m\\033[{i}m  Sample Text {reset}")
    print()

    print("Bright Foreground Colors:")
    for i in range(90, 98):  # Bright colors
        print(f"\033[{i}m\\033[{i}m  Sample Text {reset}")
    print()
    
    print("Background Colors:")
    # Background colors
    for i in range(40, 48):  # Standard background colors
        print(f"\033[{i}m\\033[{i}m  Sample Text {reset}")
    print()
    
    print("Bright Background Colors:")
    for i in range(100, 108):  # Bright background colors
        print(f"\033[{i}m\\033[{i}m  Sample Text {reset}")
    print()



if __name__ == "__main__":


    # print_ansi_table()
    # quit()

    import os
    import time
    import matplotlib.pyplot as plt

    island_cellEvo = []
    island_numEvo = []

    size = (67, 112)
    sim = Simulation(size, 0.01)

    steps = 40000
    Ndif_A = 300
    Ndif_B = 300
    coverage_limit = 0.2 # stop limiter


    sim_cells = size[0]*size[1]
    i=0
    while True:
        i+=1
        # deposit a monomer every n steps
        if i%Ndif_A == 0:
            sim.Deposit("A")

        if i%Ndif_B == 0:
            sim.Deposit("B")

        sim.Step() # step

        # Compute average number of cells per island
        isl, cells = sim.NumIslands()
        island_numEvo.append(isl)


        try:
            island_cellEvo.append( sum(cells)/len(cells) )
        except ZeroDivisionError:
            island_cellEvo.append( 0 )

        # print("==========================================================")
        fill_ratio = sum(cells)/sim_cells
        if i%100 == 0:os.system("clear");print(f'SIM {i}/{steps}\t{i/steps*100}%\t{fill_ratio*100}%');sim.Print()


        #Stop condition due to coverage limit, here we take aggregated coverage
        if fill_ratio >= coverage_limit:
            break

    plt.plot(island_numEvo, label="number of islands")
    plt.plot(island_cellEvo, label="avg cells per islands")
    plt.legend()
    plt.show()



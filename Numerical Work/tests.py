# from task2and3 import Simulation, Monomer
from task2_parallel import Simulation ; _IS_PARALLEL = True  # With numba

from matplotlib.colors import ListedColormap
import os
import time
import numpy as np
import matplotlib.pyplot as plt

def yolo(Ndif_A, Ndif_B, rate, steps=40000, size=(500, 500), coverage_limit=0.2):
    island_cellEvo = []
    island_numEvo = []
    monomer_numEvo = []

    sim = Simulation(size, rate)

    sim_cells = size[0] * size[1]
    i = 0
    while True:
        i += 1
        # deposit a monomer every n steps
        if i % Ndif_A == 0:
            sim.Deposit("A")
        if i % Ndif_B == 0:
            sim.Deposit("B")

        start_step_time = time.time()  # Start time for sim.Step()
        sim.Step()  # step
        step_time = time.time() - start_step_time  # Calculate elapsed time

        start_island_time = time.time()  # Start time for sim.NumIslands()
        isl, cells = sim.NumIslands()
        island_time = time.time() - start_island_time  # Calculate elapsed time

        start_monomer_time = time.time()  # Start time for sim.NumMonomers()
        monomer_count = sim.NumMonomers()
        monomer_time = time.time() - start_monomer_time  # Calculate elapsed time

        island_numEvo.append(isl)
        monomer_numEvo.append(monomer_count)

        try:
            island_cellEvo.append(sum(cells) / len(cells))
        except ZeroDivisionError:
            island_cellEvo.append(0)

        fill_ratio = sum(cells) / sim_cells  # < Replaced below when not using parallel as the DFS implementation is bugged.
        if not _IS_PARALLEL:
            new_array = np.array([[1 if isinstance(cell, Monomer) else 0 for cell in row] for row in sim.grid])
            fill_ratio = np.sum(new_array.flatten()) / sim_cells

        if i % 1000 == 0:
            os.system("clear")
            print(f'SIM {i}\tFill Ratio: {fill_ratio*100:.2f}%\t'
                  f'Step Time: {step_time:.6f}s\t'
                  f'Num Islands Time: {island_time:.6f}s\t'
                  f'Num Monomers Time: {monomer_time:.6f}s')

        # Stop condition due to coverage limit, here we take aggregated coverage
        if fill_ratio >= coverage_limit:
            print('Reached 20% fill')
            break

    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # Plot the image on the left
    ax1.imshow([
        [0 if isinstance(x, str) or x == 0 or not getattr(x, 'aggregated', False) else 1 for x in row]
        for row in sim.grid
    ])

    # Grid for print
    new_array = np.zeros(sim.grid.shape, dtype=int)

    new_array[sim.grid == 0] = 0
    new_array[sim.grid == 1] = 1
    new_array[sim.grid == 2] = 2
    new_array[sim.grid == 11] = 3
    new_array[sim.grid == 12] = 4

    if _IS_PARALLEL:
        cbar = ax1.imshow(new_array, cmap='viridis')

    ax1.set_title(f'Evolution for Ndif = {Ndif_A}')
    ax1.axis('off')  # Hide the axis

    # Plot the graphs on the right
    ax2.plot(monomer_numEvo, label="number of monomers")
    ax2.plot(island_numEvo, label="number of islands")
    ax2.plot(island_cellEvo, label="avg cells per islands")
    ax2.set_title('Graph of Evolution')
    ax2.legend()

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Show the plot
    plt.show()

yolo(10, 20, 0.5)

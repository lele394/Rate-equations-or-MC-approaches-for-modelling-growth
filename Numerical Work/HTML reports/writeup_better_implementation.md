<style>
/* Ensure the title doesn't get numbered */
.title {
    counter-reset: none;  /* No numbering for title */
}

/* Initialize the counter after the first h1 (document title) */
h1:not(.title) {
    counter-reset: h2; /* Reset h2 counter when a new h1 starts */
    margin-left: 0;
    counter-increment: h1; /* Increment h1 counter */
}

h1:not(.title)::before {
    content: counter(h1, upper-roman) ". "; /* Roman numeral (I., II., ...) */
}

/* For H2: Decimal with h1 (e.g., I.1) */
h2 {
    counter-reset: h3; /* Reset h3 counter when a new h2 starts */
    margin-left: 0;
}

h2::before {
    counter-increment: h2;
    content: counter(h1, upper-roman) "." counter(h2) " "; /* I.1 */
}

/* For H3: Lowercase letters (e.g., I.1.a) */
h3 {
    counter-reset: h4; /* Reset h4 counter when a new h3 starts */
    margin-left: 0;
}

h3::before {
    counter-increment: h3;
    content: counter(h1, upper-roman) "." counter(h2) "." counter(h3, lower-alpha) " "; /* I.1.a */
}

/* For H4: Lowercase roman numerals (e.g., I.1.a.i) */
h4 {
    counter-reset: h5; /* Reset h5 counter when a new h4 starts */
    margin-left: 0;
}

h4::before {
    counter-increment: h4;
    content: counter(h1, upper-roman) "." counter(h2) "." counter(h3, lower-alpha) "." counter(h4, lower-roman) " "; /* I.1.a.i */
}

/* For H5: Decimal with all previous levels (e.g., I.1.a.i.1) */
h5 {
    margin-left: 0;
}

h5::before {
    counter-increment: h5;
    content: counter(h1, upper-roman) "." counter(h2) "." counter(h3, lower-alpha) "." counter(h4, lower-roman) "." counter(h5) " "; /* I.1.a.i.1 */
}

</style>







<h1 style="text-align: center;" class="title">
    Numerical Methods : <br>
    Rate equations or MC approaches for modelling growth <br>
    Write up for a better implementation
</h1>

<p style="text-align: center;">
Léo BECHET, M2 CompuPhys 2024-2025
</p>

# Introduction

In recent evaluations of our simulation, it became apparent that a more efficient implementation could be achieved. The current design relies on updating each cell of the grid, a process that requires a total of nine passes for each update cycle. However, by leveraging a list-based approach to track only the cells that need updates, we can significantly reduce computational overhead.

# Principle

While we would still maintain the 2D grid structure to represent the simulation environment, the key optimization lies in tracking only the active cells (those that require updating) rather than iterating over the entire grid. This approach enables targeted updates, leading to performance gains, particularly in sparse simulations where only a fraction of the grid is active.

## Single-Threaded Approach

In single-threaded systems, this optimization could bring substantial improvements. By using a linked list to manage the active cells dynamically, we eliminate the need for predefined arrays and reduce unnecessary cell checks. In large grids, this would greatly decrease the number of cells processed during each cycle. 

The process would involve navigating through the list of active cells, updating them as needed, and adding newly activated neighboring cells to the list while marking them to avoid duplicate updates.

### Update Challenge

One challenge arises when updating a cell that neighbors an inactive or non-aggregated cell. This would require scanning the remaining portion of the list to identify and aggregate these neighboring cells for updates. Efficiently managing this while avoiding duplicate entries or redundant updates is crucial to maintaining performance.

## Multithreaded Approach

For a multithreaded implementation, a queue-based system could replace the linked list, allowing updates to be dispatched to multiple threads in parallel. This approach would distribute the workload across cores, further increasing performance, especially for larger grids or simulations with a high number of active cells.

### Race Condition Challenge

However, the multithreaded approach introduces potential race condition issues. If two neighboring cells are updated simultaneously by different threads, their interactions could lead to unpredictable results. To mitigate this, the dispatcher would need to implement checks to prevent updates within a defined proximity (e.g., a Z-cell radius) from occurring in parallel. Ensuring thread safety in such a system would be essential to avoid data corruption or erroneous updates.

# Conclusion

While this targeted-update approach could lead to significant performance improvements, especially in scenarios with relatively few active cells, it may also introduce overhead in simulations with a large number of active cells due to the added complexity of managing the linked list or thread queue. The performance benefits will depend on the nature of the simulation, with sparse grids benefiting the most.

We encourage others in the community to explore this approach, experiment with its implementation, and share their findings, challenges, and optimizations for handling larger-scale simulations.
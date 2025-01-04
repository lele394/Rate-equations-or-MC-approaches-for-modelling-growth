import numpy as np
import matplotlib.pyplot as plt

# Parameters
steps = 1000  # Number of time steps
n = 10  # Initial value of n
n_1 = 1  # Initial value of n_1
F_d = 0  # Constant value of F_d
alpha = 0  # Constant value of alpha
dt = 0.01  # Time step
ratio_min = 0  # Minimum F_a/D_1 ratio
ratio_max = 1000  # Maximum F_a/D_1 ratio
num_ratios = 500  # Number of ratios to sweep

# Generate a range of F_a/D_1 ratios
ratios = np.linspace(ratio_min, ratio_max, num_ratios)

# Initialize arrays to store the results
n_map = np.zeros((num_ratios, steps))
n1_map = np.zeros((num_ratios, steps))

# Loop through each F_a/D_1 ratio
for r_index, ratio in enumerate(ratios):
    F_a = ratio  # F_a is now defined as the ratio * D_1
    D_1 = 1.0  # Set D_1 to 1 for simplicity
    
    # Reset initial conditions for each ratio
    n = 1.0
    n_1 = 0.5
    
    # Compute initial accelerations (forces)
    def acceleration_n1(n, n_1):
        return F_a - F_d - 2 * D_1 * n_1**2 - D_1 * n_1 * n - alpha * F_a * n

    def acceleration_n(n_1):
        return D_1 * n_1**2 + alpha * F_a * n_1

    a_n1 = acceleration_n1(n, n_1)
    a_n = acceleration_n(n_1)

    # Simulate the loop using Velocity Verlet
    for i in range(steps):
        # Update positions
        n_1 += 0.5 * a_n1 * dt**2
        n += 0.5 * a_n * dt**2
        
        # Compute new accelerations based on updated positions
        new_a_n1 = acceleration_n1(n, n_1)
        new_a_n = acceleration_n(n_1)
        
        # Update positions using the new accelerations
        n_1 += 0.5 * new_a_n1 * dt**2
        n += 0.5 * new_a_n * dt**2
        
        # Update velocities (derivatives)
        dn_1 = 0.5 * (a_n1 + new_a_n1) * dt
        dn = 0.5 * (a_n + new_a_n) * dt
        
        # Store results
        n_map[r_index, i] = n
        n1_map[r_index, i] = n_1
        
        # Update accelerations for the next step
        a_n1 = new_a_n1
        a_n = new_a_n

# Plotting the 2D maps
plt.figure(figsize=(14, 6))

# Plot for n
plt.subplot(1, 2, 1)
plt.imshow(n_map, aspect='auto', extent=[0, steps, ratio_min, ratio_max], origin='lower')
plt.colorbar(label='n')
plt.title('2D Map of n Evolution')
plt.xlabel('Step Number')
plt.ylabel('F_a/D_1 Ratio')

# Plot for n_1
plt.subplot(1, 2, 2)
plt.imshow(n1_map, aspect='auto', extent=[0, steps, ratio_min, ratio_max], origin='lower')
plt.colorbar(label='n_1')
plt.title('2D Map of n_1 Evolution')
plt.xlabel('Step Number')
plt.ylabel('F_a/D_1 Ratio')

plt.tight_layout()
plt.show()

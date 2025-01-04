import matplotlib.pyplot as plt

# Sample data for demonstration purposes
steps = 100
n = 0
n_1 = 0
F_a = 1.2
F_d = 0.8
D_1 = 0.1
alpha = 0.05

# Initialize lists
l_n = [n]
l_n1 = [n_1]
l_dn = []
l_dn1 = []

# Simulate the loop
for i in range(steps):
    dn_1 = F_a - F_d - 2 * D_1 * n_1**2 - D_1 * n_1 * n - alpha * F_a * n
    dn = D_1 * n_1**2 + alpha * F_a * n_1

    n += dn
    n_1 += dn_1

    l_n.append(n)
    l_n1.append(n_1)

    l_dn.append(dn)
    l_dn1.append(dn_1)

# Plotting the evolution of n and n_1
plt.figure(figsize=(12, 6))

# First subplot
plt.subplot(2, 1, 1)
plt.plot(l_n, label='n')
plt.plot(l_n1, label='n_1')
plt.title('Evolution of n and n_1')
plt.xlabel('Steps')
plt.ylabel('Values')
plt.legend()

# Second subplot
plt.subplot(2, 1, 2)
plt.plot(l_dn, label='dn')
plt.plot(l_dn1, label='dn_1')
plt.title('Evolution of dn and dn_1')
plt.xlabel('Steps')
plt.ylabel('Values')
plt.legend()

plt.tight_layout()
plt.show()

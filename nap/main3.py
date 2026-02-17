import numpy as np
import matplotlib.pyplot as plt

# Define nozzle parameters
r_in = 1.0      # Inlet radius
r_t = 0.5       # Throat radius
r_o = 1.2       # Outlet radius
L1 = 5.0        # Length of converging section
L2 = 5.0        # Length of diverging section
L = L1 + L2      # Total length of nozzle

# Grid parameters
n_x = 100       # Number of grid points in the x-direction (axial direction)
n_r = 50        # Number of grid points in the r-direction (radial direction)

# Generate the x-grid (uniform grid along the length of the nozzle)
x = np.linspace(0, L, n_x)

# Generate the nozzle shape (radius at each x-point)
r_shape = np.zeros_like(x)
for i, xi in enumerate(x):
    if xi <= L1:
        # Converging section
        r_shape[i] = r_in - (r_in - r_t) * (xi / L1)
    else:
        # Diverging section
        r_shape[i] = r_t + (r_o - r_t) * ((xi - L1) / L2)

# Plot the nozzle shape and grid points
plt.figure(figsize=(8, 6))

# Plot the nozzle shape
plt.plot(x, r_shape, label="Nozzle Shape", color='orange', lw=2)

# For each axial position, plot radial grid points within the nozzle radius
for i, xi in enumerate(x):
    # Generate radial points from 0 to r(x) at each axial position
    r_grid = np.linspace(0, r_shape[i], n_r)
    
    # Plot each radial grid point at xi
    plt.scatter(np.ones_like(r_grid) * xi, r_grid, color='blue', marker='o', s=5)

# Labels and plot
plt.xlabel("Position along nozzle (x)")
plt.ylabel("Radius (r)")
plt.title("Converging-Diverging Nozzle Grid")
plt.legend()
plt.grid(True)
plt.show()

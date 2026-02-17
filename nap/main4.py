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

# Generate the radial grid (non-uniform grid in the radial direction)
r = np.linspace(0, r_o, n_r)

# Create a 2D grid for the flow variables
rho = np.ones((n_x, n_r))   # Density
u = np.zeros((n_x, n_r))    # Velocity in x-direction
v = np.zeros((n_x, n_r))    # Velocity in y-direction
p = np.ones((n_x, n_r)) * 1e5  # Pressure (initially 1 bar)
E = np.ones((n_x, n_r)) * 1e5  # Energy (initially 1 bar)

# Initial conditions
rho[:,:] = 1.0  # Set initial density to 1 (can adjust for more realism)
u[:,:] = 10.0  # Set initial velocity to 10 m/s at the inlet
v[:,:] = 0.0   # Assume no velocity in the y-direction
p[:,:] = 1e5    # Initial pressure of 1 bar (100,000 Pa)
E[:,:] = 1e5    # Initial energy (could be estimated based on temperature)

# Boundary conditions at x=0 (inlet)
u[0, :] = 10.0  # Inlet velocity
p[0, :] = 1e5    # Inlet pressure
rho[0, :] = 1.0  # Inlet density

# Boundary conditions at x=L (outlet)
u[-1, :] = u[-2, :]  # Zero-gradient for velocity at the outlet
v[-1, :] = v[-2, :]  # Zero-gradient for velocity in y-direction
p[-1, :] = p[-2, :]  # Zero-gradient for pressure at the outlet
rho[-1, :] = rho[-2, :]  # Zero-gradient for density at the outlet



# Time step size and other parameters
dt = 1e-5  # Time step
dx = x[1] - x[0]  # Axial spacing
dr = r[1] - r[0]  # Radial spacing

# Define the fluxes in x and y directions (upwind scheme)
def flux_x(rho, u, p, E, i, j):
    return rho[i, j] * u[i, j]  # Mass flux in the x-direction (upwind)

def flux_y(rho, v, p, E, i, j):
    return rho[i, j] * v[i, j]  # Mass flux in the y-direction (upwind)

def flux_momentum_x(rho, u, p, E, i, j):
    return rho[i, j] * u[i, j]**2 + p[i, j]  # Momentum flux in the x-direction

def flux_momentum_y(rho, v, p, E, i, j):
    return rho[i, j] * v[i, j] * u[i, j]  # Momentum flux in the y-direction

def flux_energy(rho, u, p, E, i, j):
    return (rho[i, j] * u[i, j] * E[i, j] + p[i, j] * u[i, j])  # Energy flux in the x-direction


# Function to update the flow variables using the finite volume method
def update_fields(rho, u, v, p, E):
    # Initialize new fields for the next time step
    rho_new = np.copy(rho)
    u_new = np.copy(u)
    v_new = np.copy(v)
    p_new = np.copy(p)
    E_new = np.copy(E)
    
    # Update density using continuity equation
    for i in range(1, n_x-1):
        for j in range(1, n_r-1):
            # Calculate fluxes in the x and y directions
            flux_x_in = flux_x(rho, u, p, E, i, j)
            flux_x_out = flux_x(rho, u, p, E, i+1, j)
            flux_y_in = flux_y(rho, v, p, E, i, j)
            flux_y_out = flux_y(rho, v, p, E, i, j+1)
            
            # Update the density
            rho_new[i, j] = rho[i, j] - dt * (flux_x_out - flux_x_in) / dx - dt * (flux_y_out - flux_y_in) / dr

    # Update velocity in the x-direction using momentum equation
    for i in range(1, n_x-1):
        for j in range(1, n_r-1):
            flux_momentum_x_in = flux_momentum_x(rho, u, p, E, i, j)
            flux_momentum_x_out = flux_momentum_x(rho, u, p, E, i+1, j)
            u_new[i, j] = u[i, j] - dt * (flux_momentum_x_out - flux_momentum_x_in) / dx
            
    # Update velocity in the y-direction (similar to the x-direction)
    for i in range(1, n_x-1):
        for j in range(1, n_r-1):
            flux_momentum_y_in = flux_momentum_y(rho, v, p, E, i, j)
            flux_momentum_y_out = flux_momentum_y(rho, v, p, E, i, j+1)
            v_new[i, j] = v[i, j] - dt * (flux_momentum_y_out - flux_momentum_y_in) / dr

    # Update energy using energy equation
    for i in range(1, n_x-1):
        for j in range(1, n_r-1):
            flux_energy_in = flux_energy(rho, u, p, E, i, j)
            flux_energy_out = flux_energy(rho, u, p, E, i+1, j)
            E_new[i, j] = E[i, j] - dt * (flux_energy_out - flux_energy_in) / dx
    
    return rho_new, u_new, v_new, p_new, E_new


# Number of time steps to run the simulation
n_steps = 1000

# Time-stepping loop
for step in range(n_steps):
    # Update the flow variables
    rho, u, v, p, E = update_fields(rho, u, v, p, E)

    # Optionally, plot the results at intervals (for visualization)
    if step % 100 == 0:
        plt.clf()
        plt.contourf(x, r, rho.T)
        plt.colorbar(label='Density')
        plt.title(f'Time Step {step}')
        plt.xlabel('x')
        plt.ylabel('r')
        plt.pause(0.1)

plt.show()

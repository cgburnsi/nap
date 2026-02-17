import numpy as np
import matplotlib.pyplot as plt

class Geometry:
    def __init__(self):
        # Define nozzle parameters
        self.r_in = 1.0      # Inlet radius
        self.r_t = 0.5       # Throat radius
        self.r_o = 1.2       # Outlet radius
        self.L1 = 5.0        # Length of converging section
        self.L2 = 5.0        # Length of diverging section
        self.L = self.L1 + self.L2  # Total length of nozzle
        
        # Grid parameters
        self.n_x = 100       # Number of grid points in the x-direction (axial direction)
        self.n_r = 50        # Number of grid points in the r-direction (radial direction)
        
        # Generate the grid points
        self.x = np.linspace(0, self.L, self.n_x)
        self.r_shape = np.zeros_like(self.x)
        
        # Generate the nozzle shape (radius at each x-point)
        for i, xi in enumerate(self.x):
            if xi <= self.L1:
                # Converging section
                self.r_shape[i] = self.r_in - (self.r_in - self.r_t) * (xi / self.L1)
            else:
                # Diverging section
                self.r_shape[i] = self.r_t + (self.r_o - self.r_t) * ((xi - self.L1) / self.L2)

    def get_radius_at(self, x):
        # Find the radius at a specific axial location
        if x <= self.L1:
            return self.r_in - (self.r_in - self.r_t) * (x / self.L1)
        else:
            return self.r_t + (self.r_o - self.r_t) * ((x - self.L1) / self.L2)
    
class Grid:
    def __init__(self, geometry):
        self.geometry = geometry
        
        self.L = geometry.L
        self.r_o = geometry.r_o
        self.n_x = geometry.n_x
        self.n_r = geometry.n_r
        
        # Grid spacings in x and r directions
        self.dx = self.L / (self.n_x - 1)
        
        # Axial positions (x)
        self.x = np.linspace(0, self.L, self.n_x)
        
        # Radial positions (r) for each x (follow the nozzle geometry)
        self.r_shape = np.zeros_like(self.x)
        for i, xi in enumerate(self.x):
            if xi <= self.geometry.L1:
                self.r_shape[i] = self.geometry.r_in - (self.geometry.r_in - self.geometry.r_t) * (xi / self.geometry.L1)
            else:
                self.r_shape[i] = self.geometry.r_t + (self.geometry.r_o - self.geometry.r_t) * ((xi - self.geometry.L1) / self.geometry.L2)
        
        # Create the radial grid at each x position, ensuring the grid matches the nozzle geometry
        self.r = [np.linspace(1e-5, self.r_shape[i], self.n_r) for i in range(self.n_x)]
        
        # Create a meshgrid for r and x
        self.X, self.R = np.meshgrid(self.x, self.r_shape)
        
        # Initialize the r_shape array as a 2D array
        self.r_shape_grid = np.zeros_like(self.X)
        for i, xi in enumerate(self.x):
            self.r_shape_grid[:, i] = self.r_shape[i]

           
class Visualization:
    @staticmethod
    def plot_geometry(geometry):
        plt.figure(1, figsize=(10, 5))
        plt.clf()

        # Plot the nozzle geometry
        plt.plot(geometry.x, geometry.r_shape, label='Nozzle Wall', color='b')
        plt.fill_between(geometry.x, geometry.r_shape, 0, alpha=0.2, color='b')
        plt.title("Nozzle Geometry")
        plt.xlabel("Axial Distance (x)")
        plt.ylabel("Radial Distance (r)")
        plt.legend()
        plt.grid(True)
        plt.show()
    
    @staticmethod
    def plot_grid(grid):
        plt.figure(figsize=(10, 5))  # Set the figure size
        plt.clf()  # Clear the previous plot
        
        # Plot the grid points at the (x, r) locations
        plt.scatter(grid.X, grid.R, c='blue', s=10)  # Plot all points in the meshgrid
        
        plt.title("Computational Grid")
        plt.xlabel("Axial Distance (x)")
        plt.ylabel("Radial Distance (r)")
        plt.grid(True)
        plt.show()
    
    
    
    @staticmethod
    def plot_pressure(grid, pressure):
        # Transpose pressure to match the dimensions of X and R
        pressure_plot = plt.contourf(grid.X.T, grid.R.T, pressure.T, levels=50, cmap='viridis')  # Transpose pressure
        plt.colorbar(pressure_plot)
        plt.title("Pressure Field")
        plt.xlabel("Axial Distance (x)")
        plt.ylabel("Radial Distance (r)")
        plt.grid(True)
        plt.show()



class FlowSolver:
    def __init__(self, grid, gamma=1.4, R=287.0, p_inlet=101325, p_outlet=100000, T_inlet=300, time_step=0.01):
        self.grid = grid
        self.gamma = gamma
        self.R = R
        self.time_step = time_step  # Add time_step as an attribute

        # Initialize flow variables (ρ, u, v, p)
        self.rho = np.ones((grid.n_x, grid.n_r)) * p_inlet / (R * T_inlet)  # Ideal gas law
        self.u = np.zeros((grid.n_x, grid.n_r))  # Axial velocity
        self.v = np.zeros((grid.n_x, grid.n_r))  # Radial velocity
        self.p = np.ones((grid.n_x, grid.n_r)) * p_inlet  # Pressure
        self.e = self.p / ((gamma - 1) * self.rho)  # Internal energy per unit mass

        # Boundary conditions
        self.p_inlet = p_inlet
        self.p_outlet = p_outlet
        self.T_inlet = T_inlet
        self.rho_inlet = p_inlet / (R * T_inlet)
        self.u_inlet = 10.0  # Assumed initial subsonic inlet velocity

    def apply_boundary_conditions(self):
        """
        Apply boundary conditions at the inlet, outlet, wall, and symmetry boundaries.
        """
        # Inlet (x = 0)
        self.rho[0, :] = self.rho_inlet
        self.u[0, :] = self.u_inlet
        self.v[0, :] = 0  # No radial velocity at the inlet
        self.p[0, :] = self.p_inlet

        # Outlet (x = L)
        self.p[-1, :] = self.p_outlet  # Subsonic outlet: specify pressure
        self.rho[-1, :] = self.rho[-2, :]  # Extrapolate density
        self.u[-1, :] = self.u[-2, :]  # Extrapolate velocity
        self.v[-1, :] = self.v[-2, :]  # Extrapolate radial velocity

        # Wall boundary (no-slip condition)
        self.u[:, -1] = 0  # Axial velocity = 0 at the wall
        self.v[:, -1] = 0  # Radial velocity = 0 at the wall
        # Adiabatic wall: extrapolate energy or use wall temperature
        self.p[:, -1] = self.p[:, -2]  # Extrapolate pressure near the wall

        # Symmetry boundary (r = 0)
        self.v[:, 0] = 0  # No radial velocity at the centerline
        self.rho[:, 0] = self.rho[:, 1]  # Extrapolate density
        self.u[:, 0] = self.u[:, 1]  # Extrapolate axial velocity
        self.p[:, 0] = self.p[:, 1]  # Extrapolate pressure

    def compute_fluxes(self):
        """
        Compute the fluxes in both axial and radial directions using Lax-Friedrichs scheme.
        """
        # Initialize flux arrays
        F_rho = np.zeros_like(self.rho)
        F_u = np.zeros_like(self.u)
        F_v = np.zeros_like(self.v)
        F_p = np.zeros_like(self.p)
        F_e = np.zeros_like(self.e)

        # Loop through the grid and compute fluxes (Lax-Friedrichs scheme)
        for i in range(1, self.grid.n_x - 1):  # Axial direction
            for j in range(1, self.grid.n_r - 1):  # Radial direction
                # Axial fluxes (F_x)
                F_rho[i, j] = 0.5 * (self.rho[i+1, j] * self.u[i+1, j] + self.rho[i-1, j] * self.u[i-1, j]) \
                                - 0.5 * (self.rho[i+1, j] - self.rho[i-1, j])
                F_u[i, j] = 0.5 * (self.rho[i+1, j] * self.u[i+1, j]**2 + self.p[i+1, j] + self.rho[i-1, j] * self.u[i-1, j]**2 + self.p[i-1, j]) \
                                - 0.5 * (self.rho[i+1, j] * self.u[i+1, j] - self.rho[i-1, j] * self.u[i-1, j])
                F_v[i, j] = 0.5 * (self.rho[i+1, j] * self.u[i+1, j] * self.v[i+1, j] + self.rho[i-1, j] * self.u[i-1, j] * self.v[i-1, j]) \
                                - 0.5 * (self.rho[i+1, j] * self.v[i+1, j] - self.rho[i-1, j] * self.v[i-1, j])
                F_p[i, j] = 0.5 * (self.u[i+1, j] * (self.p[i+1, j] + self.e[i+1, j] / self.rho[i+1, j]) + self.u[i-1, j] * (self.p[i-1, j] + self.e[i-1, j] / self.rho[i-1, j])) \
                                - 0.5 * (self.u[i+1, j] - self.u[i-1, j])

                # Energy flux (F_e)
                F_e[i, j] = 0.5 * (self.u[i+1, j] * (self.e[i+1, j] + self.p[i+1, j] / self.rho[i+1, j])) \
                             + 0.5 * (self.u[i-1, j] * (self.e[i-1, j] + self.p[i-1, j] / self.rho[i-1, j]))

        return F_rho, F_u, F_v, F_p, F_e

    def runge_kutta_step(self, delta_t=None):
        """
        Perform a single Runge-Kutta time-stepping step (4th-order).
        """
        if delta_t is None:
            delta_t = self.time_step  # Use the time_step if not provided

        # Stage 1: Compute the fluxes for the initial state
        F_rho_1, F_u_1, F_v_1, F_p_1, F_e_1 = self.compute_fluxes()

        # Stage 2: Compute the fluxes for the state updated with half the previous fluxes
        # Update the flow variables by half the fluxes from stage 1
        self.rho += 0.5 * delta_t * F_rho_1
        self.u += 0.5 * delta_t * F_u_1
        self.v += 0.5 * delta_t * F_v_1
        self.p += 0.5 * delta_t * F_p_1
        self.e += 0.5 * delta_t * F_e_1

        F_rho_2, F_u_2, F_v_2, F_p_2, F_e_2 = self.compute_fluxes()

        # Stage 3: Compute the fluxes for the state updated with half the fluxes from stage 2
        self.rho += 0.5 * delta_t * F_rho_2
        self.u += 0.5 * delta_t * F_u_2
        self.v += 0.5 * delta_t * F_v_2
        self.p += 0.5 * delta_t * F_p_2
        self.e += 0.5 * delta_t * F_e_2

        F_rho_3, F_u_3, F_v_3, F_p_3, F_e_3 = self.compute_fluxes()

        # Stage 4: Compute the fluxes for the final state updated with the full fluxes from stage 3
        self.rho += delta_t * F_rho_3
        self.u += delta_t * F_u_3
        self.v += delta_t * F_v_3
        self.p += delta_t * F_p_3
        self.e += delta_t * F_e_3

        F_rho_4, F_u_4, F_v_4, F_p_4, F_e_4 = self.compute_fluxes()

        # Update the flow variables using the weighted average of the fluxes
        self.rho += (delta_t / 6.0) * (F_rho_1 + 2 * F_rho_2 + 2 * F_rho_3 + F_rho_4)
        self.u += (delta_t / 6.0) * (F_u_1 + 2 * F_u_2 + 2 * F_u_3 + F_u_4)
        self.v += (delta_t / 6.0) * (F_v_1 + 2 * F_v_2 + 2 * F_v_3 + F_v_4)
        self.p += (delta_t / 6.0) * (F_p_1 + 2 * F_p_2 + 2 * F_p_3 + F_p_4)
        self.e += (delta_t / 6.0) * (F_e_1 + 2 * F_e_2 + 2 * F_e_3 + F_e_4)

        # Apply boundary conditions after updating the flow variables
        self.apply_boundary_conditions()



if __name__ == '__main__':
    # Create instances of the Geometry, Grid, and FlowSolver classes
    geom = Geometry()
    grid = Grid(geometry=geom)
    #solver = FlowSolver(grid, time_step=0.01)

    # Apply initial boundary conditions
    #solver.apply_boundary_conditions()

    # Run multiple time steps in a loop (e.g., 100 time steps)
    #for step in range(100):
    #    solver.runge_kutta_step()  # Perform one Runge-Kutta step
        # Optionally, visualize or print progress every few steps
    #    if step % 10 == 0:
    #        print(f"Step {step}")
    #        Visualization.plot_geometry(geom)  # Visualize geometry
    #        Visualization.plot_grid(grid)      # Visualize grid
            #Visualization.plot_pressure(grid, solver.p)  # Visualize pressure field

    # Final visualization after the simulation
    Visualization.plot_geometry(geom)
    Visualization.plot_grid(grid)
    #Visualization.plot_pressure(grid, solver.p)  # Final pressure visualization




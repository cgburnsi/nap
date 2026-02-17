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
        
        # Axial positions (x)
        self.x = np.linspace(0, self.L, self.n_x)
        
        # Radial positions for each axial location (scaled by nozzle radius at each point)
        self.r = np.linspace(0, 1, self.n_r)  # Normalized radial grid
        self.r_shape = np.zeros((self.n_x, self.n_r))

        # Generate radial positions at each axial point based on nozzle shape
        for i, xi in enumerate(self.x):
            # Scale the radial grid to match the nozzle's radius at this axial position
            self.r_shape[i, :] = self.geometry.get_radius_at(xi) * self.r

        # Create meshgrid of x and r for visualization (rectangular grid)
        self.X, self.R = np.meshgrid(self.x, self.r_shape)

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

# Example Usage:
if __name__ == '__main__':
    # Create instances of the Geometry and Grid classes
    geom = Geometry()
    grid = Grid(geometry=geom)

    # Final visualization after the simulation
    Visualization.plot_geometry(geom)
    Visualization.plot_grid(grid)

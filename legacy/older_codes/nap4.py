import numpy as np
import matplotlib.pyplot as plt

   
class Sim:
    def __init__(self,l_max=10, m_max=5, max_iter=400, PT=482633, TT=299.817, THETA=0, PE=101352.9):
        self.l_max      = l_max                                         # [-] Maximum x-axis (Axial Direction) grid points
        self.m_max      = m_max                                         # [-] Maximum y-axis (Radial Direction) grid points
        self.n_max      = max_iter                                      # [-] Maximum number of time steps
        
        # Solution Surface Arrays
        self.U      = np.zeros((1, self.m_max, self.l_max))             # [m] Velocity in the x-direction at each grid point   
        self.V      = np.zeros((1, self.m_max, self.l_max))             # [m] Velocity in the y-direction at each grid point
        self.P      = np.zeros((1, self.m_max, self.l_max))             # [m] Pressure at each grid point
        self.RO     = np.zeros((1, self.m_max, self.l_max))             # [m] Density at each grid point
        
        # Boundary Conditions
        self.PT         = np.full(self.m_max, PT)                       # [Pa] Inlet boundary pressure (uniform across inlet)
        self.TT         = np.full(self.m_max, TT)                       # [K] Inlet boundary temperature (uniform across inlet)
        self.theta      = np.full(self.m_max, THETA)                    # [radians] Inlet boundary flow angle (uniform across inlet)
        self.PE         = PE                                            # [Pa] Exit boundary pressure condition

    
        
    
    
    
    
    
    
    
    def plot_surface(self):
        surface_data = self.U[0,:,:]
        rows, cols = surface_data.shape
        x, y = np.meshgrid(np.arange(cols), np.arange(rows))
    
        # Plot the points
        plt.scatter(x, y, c=surface_data.flatten(), cmap='viridis')  # Color by the value at each point
        plt.colorbar(label='Value')  # Optional colorbar to show the value scale
        plt.xlabel('Column index')
        plt.ylabel('Row index')
        plt.title(f'Iteration {1}')
        plt.gca().invert_yaxis()  # Optional: To invert the y-axis if rows increase downward
        plt.show()
            
    
    
    
    
    
        
if __name__ == '__main__':
    a = Sim()
    #a.plot_surface()
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
''' NAP flow

1. Set defaults
2. Read in Data
3. Calculate Nozzle Wall and Normals
4. Convert metric units to english
5. Convert to internal units
6. Calculate the Initial-Data Surface
7. Calculate the Initial-Data Surface Mass Flow and Thrust
8. Initialize the Time Step Integration Loop Parameters
9. Enter the Time Integration Loop
9a. Calculate Delta T
9b. Determin if the exit flow is subsonic or supersonic
9c. Calculate the Nozzle Wall and Interior Mesh Points
9d. Extrapolate the exit mesh points for supersonic flow
9e. Calculate the nozzle inlet mesh points
9f. Determine the maximum (Delta U)/U
9g. Compute the solution surface mass flowrate and thrust
9h. Calculate and print the solution surface
9i. Check for convergence of the steady state solution


1. read in data
2. set default




'''
    
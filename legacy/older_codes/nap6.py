
import numpy as np
import convert as cv
import matplotlib.pyplot as plt

# Default Values so I don't have to look them up all the time: N1D = 1, NDIM = 1, NGCB = 0, NGEOM = 2

class Solution:
    def __init__(self, l_max=15, m_max=5, n_max=400, x_i=0.007874, r_i=0.0635, r_t=0.02032, x_e=0.1143, rc_i=0.02032, rc_t=0.0127, 
                 ang_i=44.88, ang_e=15, PT=482633, TT=299.817, THETA=0, PE=101352.9, t_conv=3e5, t_stop=1.0, gamma=1.4,
                 R_gas=287.052874, FDT=1.0):
        self.l_max      = l_max                                         # [-] Maximum x-axis (Axial Direction) grid points
        self.m_max      = m_max                                         # [-] Maximum y-axis (Radial Direction) grid points
        self.n_max      = n_max                                         # [-] Maximum number of time steps
        self.x_i        = x_i                                           # [m] x-axis location for inlet point
        self.x_t        = np.nan                                        # [m] x-axis location for throat (computed later)
        self.r_i        = r_i                                           # [m] y-axis location for wall inlet point
        self.r_t        = r_t                                           # [m] y-axis location for the nozzle throat
        self.x_e        = x_e                                           # [m] x-axis location for nozzle exit
        self.r_e        = np.nan                                        # [m] y-axis location for nozzle exit (computed later)
        self.rc_i       = rc_i                                          # [m] Radius of Curvature at nozzle wall inlet contraction
        self.rc_t       = rc_t                                          # [m] Radius of Curvature at nozzle throat
        self.ang_i      = np.deg2rad(ang_i)                             # [radians] Contraction Angle
        self.ang_e      = np.deg2rad(ang_e)                             # [radians] Nozzle Exit Angle
        self.r_star     = np.nan                                        # [m] Area per unit depth where Mach = 1
        self.r_stars    = np.nan                                        # [m^2] Area per unit depth divided by np.pi
        self.a_t        = np.pi * r_t**2                                # [m^2] Throat Area
        
        self.dx         = (self.x_e - self.x_i) / (self.l_max - 1.0)    # TODO (Should this be L_max - 0) [m] distance between equally spaced x axis points   
        self.dy         = 1.0 / (self.m_max - 1.0)                      # TODO (Should this be m_max-0?) [m] distance between equally spaced y-axis points
        self.yw         = np.zeros(self.l_max)                          # [m] y-axis location of wall (equally spaced)
        self.xw         = np.zeros_like(self.yw)                        # [m] x-axis location of wall (equally spaced)
        self.yw_i       = np.zeros_like(self.yw)                        # [m] y-axis location of wall (non-equally spaced)
        self.xw_i       = np.zeros_like(self.yw)                        # [m] x-axis location of wall (non-equally spaced)
        self.nxny       = np.zeros_like(self.yw)                        # [-] negative of wall slopes corresponding to yw elements
        
        # Boundary Conditions
        self.PT         = np.full(self.m_max, PT)                       # [Pa] Inlet boundary pressure (uniform across inlet)
        self.TT         = np.full(self.m_max, TT)                       # [K] Inlet boundary temperature (uniform across inlet)
        self.theta      = np.full(self.m_max, THETA)                    # [radians] Inlet boundary flow angle (uniform across inlet)
        self.PE         = PE                                            # [Pa] Exit boundary pressure condition
        
        # Fluid Properties
        self.g          = 9.807                                         # [m/s^2] acceleration due due to gravity on Earth
        self.R_gas      = R_gas                                         # [J/kg-K] Specific gas constant
        self.gamma      = gamma                                         # [-] Ratio of specific heats for fluid
       
        # Calculation Helper Values
        self.gam_1      = self.gamma / (self.gamma -1)                  # [-] Calculation Helper
        self.gam_2      = (self.gamma - 1) / 2.0                        # [-] Calculation Helper
        self.l_t        = 0                                             # [m] Index for the throat location in the yw array
        self.l_1        = self.l_max - 1                                # [-] locator value for the edge of the grid in x-axis direction
        self.l_2        = self.l_max - 2                                # [-] locator value for the edge of the grid in x-axis direction
        self.l_3        = self.l_max - 3                                # [-] locator value for the edge of the grid in x-axis direction
        self.m_1        = self.m_max - 1                                # [-] locator value for the edge of the grid in y-axis direction
        self.m_2        = self.m_max - 2                                # [-] locator value for the edge of the grid in y-axis direction

        # Simulation Parameters
        self.t_conv     = t_conv                                        # [-] Axial velocity steady-state convergence criteria
        self.t_stop     = t_stop                                        # [sec] Time stop for simulation
        self.FDT        = FDT                                           # [-] Premultiplier on the C-F-L paramter
        
        # Solution Arrays
        self.U          = np.zeros((1, self.m_max, self.l_max))         # [m/s] Velocity x-axis (Axial Direction)
        self.V          = np.zeros((1, self.m_max, self.l_max))         # [m/s] Velocity y-axis (Radial Direction)
        self.P          = np.zeros((1, self.m_max, self.l_max))         # [Pa] Pressure Field
        self.RO         = np.zeros((1, self.m_max, self.l_max))         # [kg/m^3] Density Field
        self.M_one_dim  = np.zeros(self.l_max)                          # [-] Mach Number for one dimensional initial surface

    def build_wall_geometry(self):
        x_tan = self.x_i + self.rc_i * np.sin(self.ang_i)
        r_tan = self.r_i + self.rc_i * (np.cos(self.ang_i) - 1.0)
        r_t_1 = self.r_t - self.rc_t * (np.cos(self.ang_i) - 1.0)
        x_t_1 = x_tan + (r_tan - r_t_1) / np.tan(self.ang_i)
        
        if x_t_1 < x_tan: x_t_1 = x_tan; r_t_1 = r_tan

        self.x_t    = x_t_1 + self.rc_t * np.sin(self.ang_i)
        x_t_2       = self.x_t + self.rc_t * np.sin(self.ang_e)
        r_t_2       = self.r_t + self.rc_t * (1.0 - np.cos(self.ang_e))
        self.r_e    = r_t_2 + (self.x_e - x_t_2) * np.tan(self.ang_e)
        self.l_t    = 0

        for L in range(0, self.l_max):
            x = self.x_i + self.dx * L
            self.xw[L] = x
            if x >= self.x_i and x <= x_tan:
                self.yw[L]      = self.r_i + self.rc_i * (np.cos(np.arcsin((x - self.x_i) / self.rc_i)) - 1.0)
                self.nxny[L]    = (x - self.x_i) / (self.yw[L] - self.r_i + self.rc_i)
            elif x > x_tan and x <= x_t_1:
                self.yw[L]      = r_t_1 + (x_t_1 - x) * np.tan(self.ang_i)
                self.nxny[L]    = np.tan(self.ang_i)
            elif x > x_t_1 and x <= self.x_t:
                self.yw[L]      = self.r_t + self.rc_t * (1.0 - np.cos(np.arcsin((self.x_t - x) / self.rc_t)))
                self.nxny[L]    = (self.x_t - x) / (self.rc_t + self.r_t - self.yw[L])
            elif x > self.x_t and x <= x_t_2:
                self.yw[L]      = self.r_t + self.rc_t * (1.0-np.cos(np.arcsin((x - self.x_t) / self.rc_t)))
                self.nxny[L]    = (self.x_t - x) / (self.rc_t + self.r_t - self.yw[L])
            elif x > x_t_2 and x <= self.x_e:
                self.yw[L]      = r_t_2 + (x - x_t_2) * np.tan(self.ang_e)
                self.nxny[L]    = -np.tan(self.ang_e)
            
            # Determine where the index is for the throat in the yw (outer wall) array
            if L <= 0: 
                continue
            else:
                if self.yw[L] < self.yw[self.l_t]: 
                    self.l_t = L;

    def plot_wall_contour(self):
        plt.figure(figsize=(10, 5))
        plt.plot(self.xw, self.yw, label='Outer Wall')
        plt.xlabel("x-axis [m]")
        plt.ylabel("y-axis [m]")
        plt.grid()
        plt.legend()
    
    def plot_surface(self, iter_index, surf):
        surface_data = surf[iter_index,:,:]
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
        
    def one_dimisional(self):
        # Calculation helper values
        GRGAS = 1.0 / (self.R_gas)
        ACOEF = 2.0 / (self.gamma + 1.0)
        BCOEF = (self.gamma - 1.0) / (self.gamma + 1.0)
        CCOEF = (self.gamma + 1.0) / 2.0 / (self.gamma - 1.0)
        
        MN3 = 0.1
        
        def Newton_Raphson(MN3, A_ratio):
            for _ in range(0, 100):                         # Hardcode the number of iterations in the Newton Method.
                ABM = ACOEF + BCOEF * MN3 ** 2
                ABMC = ABM ** CCOEF
                FM = (ABMC / MN3) - A_ratio
                FPM = ABMC * (2.0 * BCOEF * CCOEF / ABM - 1.0 / MN3 ** 2)
                OMN3 = MN3
                MN3 = OMN3 - FM / FPM
            
                if OMN3 > 0.99 and OMN3 < 1.01:         MN3 = 0.5*(OMN3+MN3)
                if MN3 > 1.0 and OMN3 < 1.0:            MN3 = 0.99
                if MN3 < 1.0 and OMN3 > 1.0:            MN3 = 1.01
                if MN3 > 50.0:                          MN3 = 50.0
                if MN3 < 0.0:                           MN3 = -MN3; break
                if abs(MN3 - OMN3) / OMN3 <= 0.0005:    break
            
            return MN3
        
        # Calculate the local One Dimensional Mach Number
        for L in range(0, self.l_max):
            x = self.x_i + self.dx * (L-1)                          # [m] Current location along the x-axis direction at 'M=0'
            if x < self.x_t:                                    # Subsonic region of nozzle upstream of throat
                area = np.pi * self.yw[L]**2                    # [m**2] Area at current location
                A_ratio = area / self.a_t                       # [-] Area Ratio at current location
                MN3 = Newton_Raphson(.01, A_ratio)              # [-] Solve for the local Mach number (subsonic Mach Guess)
            elif x >= self.x_t:
                area = np.pi * self.yw[L]**2                    # [m**2] Area at current location
                A_ratio = area / self.a_t                       # [-] Area Ratio at current location
                MN3 = Newton_Raphson(5, A_ratio)              # [-] Solve for the local Mach number (supersonic Mach Guess)
            self.M_one_dim[L] = MN3
                
        
        # Calculate P, RO, U, and V for the initial surface from the local Mach number             
            for M in range(0, self.m_max):
                DEM = 1.0 + self.gam_2 * self.M_one_dim[L]**2
                DEMP = DEM**self.gam_1
                self.P[0][M][L] = self.PT[M] / DEMP
                self.RO[0][M][L] = self.P[0][M][L] * GRGAS / (self.TT[M] / DEM)
                Q = MN3 * np.sqrt(self.gamma * self.P[0][M][L] / self.RO[0][M][L])
                print(Q, MN3)
                self.U[0][M][L] = Q
                self.V[0][M][L] = 0.0

            #print(f'L = {L:<2}, x = {x*12:<3.2f}, U = {self.U[0, 0, L]:<3.2f}, V = {self.V[0, 0, L]:<3.2f}, P = {self.P[0, 0, L]/144:<3.2f}, RO = {self.RO[0, 0, L]:<3.4f} Area Ratio = {A_ratio:<3.5f}, Q = {Q:<3.2f}, Mach {MN3:<3.2f}')

    def masflo(self):
        # Helper Calculations
        LDUM = self.l_max - 1
        
        # Calculate the mass flow and thrust for the 1-D Initial Data Surface
        a_inlet     = (np.pi * self.yw[0]**2)
        a_throat    = (np.pi * self.yw[self.l_t]**2)
        a_exit      = (np.pi * self.yw[LDUM]**2)
        vm_i        = np.sqrt(self.U[0, 0, 0]**2 + self.V[0, 0, 0]**2)
        vm_t        = np.sqrt(self.U[0, 0, self.l_t]**2 + self.V[0, 0, self.l_t]**2)
        vm_e        = np.sqrt(self.U[0, 0, LDUM]**2 + self.V[0, 0, LDUM]**2)
        mass_i      = self.RO[0, 0, 0] * vm_i * a_inlet * self.g
        mass_t      = self.RO[0, 0, self.l_t] * vm_t * a_throat * self.g
        mass_e      = self.RO[0, 0, LDUM] * vm_e * a_exit * self.g
        thrust      = self.RO[0, 0, LDUM] * self.U[0, 0, LDUM]**2 * a_exit
        
        print(f'vm_e = {vm_e:<3.2f}, mass_i = {mass_i:<3.4f}, mass_t = {mass_t:<3.4f}, mass_e = {mass_e:<3.4f}, thrust = {thrust:<3.2f}')


    def print_one_dim_data_base_units(self):
        cols = {"L [-]": 5, "x [m]": 8, "U [m/s]": 10, "V [m/s]": 10, "P [Pa]": 12, "RO [kg/m3]": 12, "Mach [-]": 10}

        headers = list(cols.keys())
        # Print headers
        header_row = "|".join(f"{h:>{cols[h]}}" for h in headers)
        print(header_row)
        print("-" * len(header_row))

        # Print Data
        for idx in range(0, self.l_max):
            print(f"{idx:>{cols['L [-]']}}|{self.xw[idx]:{cols['x [m]']}.2f}|"
                  f"{self.U[0, 0, idx]:>{cols['U [m/s]']}.2f}|{self.V[0, 0, idx]:>{cols['V [m/s]']}.2f}|"
                  f"{self.P[0, 0, idx]:>{cols['P [Pa]']}.2f}|{self.RO[0, 0, idx]:>{cols['RO [kg/m3]']}.4f}|"
                  f"{self.M_one_dim[idx]:>{cols['Mach [-]']-2}.2f}|")

    def print_one_dim_data_mixed_units(self):
        cols = {"L [-]": 5, "x [m]": 8, "U [m/s]": 10, "V [m/s]": 10, "P [Pa]": 12, "P [psi]": 8, "RO [kg/m3]": 12, "Mach [-]": 10}

        headers = list(cols.keys())
        # Print headers
        header_row = "|".join(f"{h:>{cols[h]}}" for h in headers)
        print(header_row)
        print("-" * len(header_row))

        # Print Data
        for idx in range(0, self.l_max):
            print(f"{idx:>{cols['L [-]']}}|{self.xw[idx]:{cols['x [m]']}.2f}|"
                  f"{self.U[0, 0, idx]:>{cols['U [m/s]']}.2f}|{self.V[0, 0, idx]:>{cols['V [m/s]']}.2f}|"
                  f"{self.P[0, 0, idx]:>{cols['P [Pa]']}.2f}|{cv.convert(self.P[0, 0, idx], 'Pa', 'psi'):>{cols['P [psi]']}.2f}|"
                  f"{self.RO[0, 0, idx]:>{cols['RO [kg/m3]']}.4f}|"
                  f"{self.M_one_dim[idx]:>{cols['Mach [-]']-2}.2f}|")
 
    def print_one_dim_data_eng_units(self):
        cols = {"L [-]": 5, "x [in]": 8, "U [ft/s]": 10, "V [ft/s]": 10, "P [psi]": 8, "RO [lbm/ft3]": 12, "Mach [-]": 10}

        headers = list(cols.keys())
        # Print headers
        header_row = "|".join(f"{h:>{cols[h]}}" for h in headers)
        print(header_row)
        print("-" * len(header_row))

        # Print Data
        for idx in range(0, self.l_max):
            print(f"{idx:>{cols['L [-]']}}|{cv.convert(self.xw[idx], 'm', 'in'):{cols['x [in]']}.2f}|"
                  f"{cv.convert(self.U[0, 0, idx], 'm/s', 'ft/s'):>{cols['U [ft/s]']}.2f}|"
                  f"{cv.convert(self.V[0, 0, idx], 'm/s', 'ft/s'):>{cols['V [ft/s]']}.2f}|"
                  f"{cv.convert(self.P[0, 0, idx], 'Pa', 'psi'):>{cols['P [psi]']}.2f}|"
                  f"{cv.convert(self.RO[0, 0, idx], 'kg/m^3', 'lb/ft^3'):>{cols['RO [lbm/ft3]']}.4f}|"
                  f"{self.M_one_dim[idx]:>{cols['Mach [-]']-2}.2f}|")

        
        
        
             

   
if __name__ == '__main__':
    
    # NAP Document Case #1 - Converging-Diverging Nozzle (45 degree Inlet, 15 degree Exit)
    c1 = {'x_i':  cv.convert(0.31, 'in', 'm'),   'r_i':  cv.convert(2.50, 'in', 'm'),
          'r_t':  cv.convert(0.80, 'in', 'm'),   'x_e':  cv.convert(4.05, 'in', 'm'),
          'rc_i': cv.convert(0.8, 'in', 'm'),    'rc_t': cv.convert(0.5, 'in', 'm'),
          'PT':   cv.convert(70, 'psi', 'Pa'),   'TT':   cv.convert(80, 'degF', 'K'),
          'PE':   cv.convert(14.7, 'psi', 'Pa'), 'THETA': 0.0,
          'R_gas': cv.convert(53.353, 'ft*lbf/lb*degR', 'J/kg*K'), 'gamma': 1.4, 
          'FDT': 1.6, 'ang_i': 44.88, 'ang_e': 15.0, 't_conv': 3.0e5, 't_stop': 1.0}
        
    a = Solution(x_i=c1['x_i'], r_i=c1['r_i'], r_t=c1['r_t'], x_e=c1['x_e'], rc_i=c1['rc_i'], rc_t=c1['rc_t'], 
                 ang_i=c1['ang_i'], ang_e=c1['ang_e'], PT=c1['PT'], TT=c1['TT'], THETA=c1['THETA'], PE=c1['PE'], 
                 t_conv=c1['t_conv'], t_stop=c1['t_stop'], gamma=c1['gamma'], R_gas=c1['R_gas'], FDT=c1['FDT'])
       
    a.build_wall_geometry()
    #a.plot_wall_contour()
    a.one_dimisional()
    
    #a.print_one_dim_data_mixed_units()
    a.print_one_dim_data_eng_units()

    a.masflo()
    #a.plot_surface(0, a.U)

    
    gamma = 1.4
    gam_1 = gamma / (gamma -1) 
    gam_2 = (gamma - 1) / 2.0                        # [-] Calculation Helper

    M = .06
    M = .1269
    P = 70*(1-gam_2 * M**2)**gam_1
    
    print(P)
    
    
    
    
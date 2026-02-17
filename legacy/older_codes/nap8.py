import numpy as np
import convert as cv
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



class Thruster:
    def __init__(self, x_pts=80, y_pts=20, l_c=cv.convert(2, 'mm', 'm'), r_c=cv.convert(1.5, 'mm', 'm'), 
                 r_t=cv.convert(.15, 'mm', 'm'), expan=139, ang_i=44.88, ang_e=15):
        self.x_pts      = x_pts                     # [-] Number of grid points in the x-axis direction (axial along nozzle)
        self.y_pts      = y_pts                     # [-] Number of grid points in the y-axis direction (radial from centerline)
        
        # Thruster Geometry Parameters
        self.l_c        = l_c                       # [m] Combustion Chamber Length
        self.r_c        = r_c                       # [m] Combustion Chamber Radius
        self.r_t        = r_t                       # [m] Throat Radius
        self.r_e        = np.sqrt(expan*r_t**2)     # [m] Nozzle Exit Radius
        self.expan      = expan                     # [m] Nozzle Exit Expansion Ratio
        self.ang_i      = np.deg2rad(ang_i)         # [radians] Nozzle Contraction Angle
        self.ang_e      = np.deg2rad(ang_e)         # [radians] Nozzle Exit Angle
        self.points     = {}                        # [dict] x and y locations for the calculation points of the geometry
        self.arc_pts    = {}                        # [dict] x and y locations for the center of the wall geometry arcs
        
        self.yw         = np.zeros(self.x_pts)      # [m] Chamber wall points
        self.xy         = np.zeros(self.x_pts)      # [m] x-axis location for each of the wall y-coordinates 
        
        
        # Calculations
        self.calculate_wall()
    
    def calculate_wall(self):
        # Radii multiplers based on the throat used in the determination of the other points
        r1 = 10    * self.r_t
        r2 = 1.5   * self.r_t
        r3 = 0.382 * self.r_t
                
        # Helper function for circular arcs
        def calculate_arc(center_x, center_y, radius, angle_start, angle_end, num_points=50):
            angles = np.linspace(angle_start, angle_end, num_points)
            x = center_x + radius * np.cos(angles)
            y = center_y + radius * np.sin(angles)
            return x, y
        
        # Angles for the arcs
        alpha0, alphaf = np.radians(90),  np.radians(90)-self.ang_i
        beta0, betaf   = np.radians(-90)-self.ang_i, np.radians(-90)
        gamma0, gammaf = np.radians(-90), self.ang_e - np.radians(90)
        
        self.points['A'] = {'x': 0, 'y': self.r_c}          # Chamber Inlet 
        self.points['B'] = {'x': self.l_c, 'y': self.r_c}   # Chamber Contraction Start
        
        # End of the Combustion Chamber Contraction Arc
        self.arc_pts['p1'] = {'x': self.l_c, 'y': self.r_c - r1}
        Cx, Cy = self.arc_pts['p1']['x'] + r1 * np.cos(alphaf), self.arc_pts['p1']['y'] + r1 * np.sin(alphaf)
        arc_bc_x, arc_bc_y = calculate_arc(self.arc_pts['p1']['x'], self.arc_pts['p1']['y'], r1, alpha0, alphaf)
        self.points['C'] = {'x': Cx, 'y': Cy}
  
        # Start of Nozzle Inlet Arc
        m1 = np.tan(np.radians(90) - self.ang_i)
        Dy = (self.r_t + r2) - (r2 * np.sin(alphaf))
        Dx = Cx - (Dy - Cy) * m1
        self.points['D'] = {'x': Dx, 'y': Dy}
        
        # End of Nozzle Inlet Arc (Throat)
        self.arc_pts['p2'] = {'x': Dx + r2 * np.cos(alphaf), 'y': self.r_t + r2}
        arc_de_x, arc_de_y = calculate_arc(self.arc_pts['p2']['x'], self.arc_pts['p2']['y'], r2, beta0, betaf)
        self.points['E'] = {'x': self.arc_pts['p2']['x'], 'y': self.r_t}
        
        # Throat to the end of the exit arc
        self.arc_pts['p3'] = {'x': self.arc_pts['p2']['x'], 'y': self.r_t + r3}
        Fx, Fy = self.arc_pts['p3']['x'] + r3 * np.cos(gammaf), self.arc_pts['p3']['y'] + r3 * np.sin(gammaf)
        arc_ef_x, arc_ef_y = calculate_arc(self.arc_pts['p3']['x'], self.arc_pts['p3']['y'], r3, gamma0, gammaf)
        self.points['F'] = {'x': Fx, 'y': Fy}
        
        # Point G
        m2 = np.tan(self.ang_e)
        Gy = self.r_e
        Gx = Fx + (Gy - Fy) / m2
        self.points['G'] = {'x': Gx, 'y': Gy}
         
        # Combine arcs and line lengths into x and y vectors
        x_full = np.concatenate([np.array([self.points['A']['x']]),                # Point A
                                 arc_bc_x,                                    # Chamber contraction arc
                                 np.array([self.points['C']['x'], self.points['D']['x']]),  # Linear segment C-D
                                 arc_de_x,                                    # Throat entrance arc
                                 np.array([self.points['E']['x'], self.points['F']['x']]),  # Linear segment E-F
                                 arc_ef_x,                                    # Throat exit arc
                                 np.array([self.points['G']['x']])                 # Point G
                                 ])
        y_full = np.concatenate([np.array([self.points['A']['y']]),                 # Point A
                                 arc_bc_y,                                    # Chamber contraction arc
                                 np.array([self.points['C']['y'], self.points['D']['y']]),  # Linear segment C-D
                                 arc_de_y,                                    # Throat entrance arc
                                 np.array([self.points['E']['y'], self.points['F']['y']]),  # Linear segment E-F
                                 arc_ef_y,                                    # Throat exit arc
                                 np.array([self.points['G']['y']])                 # Point G
                                 ])
        
        # Compute cumulative arc length (Needed for the dx calculation later)
        distances = np.sqrt(np.diff(x_full)**2 + np.diff(y_full)**2)
        cumulative_length = np.concatenate(([0], np.cumsum(distances)))
        
        # Discretize the wall points
        num_points = 100
        interp_length = np.linspace(0, cumulative_length[-1], num_points)
        self.interp_x = interp1d(cumulative_length, x_full, kind='linear')(interp_length)
        self.interp_y = interp1d(cumulative_length, y_full, kind='linear')(interp_length)
        
        #for idx in range(0, len(self.yw)):
        #    print(idx)


    def plot_wall(self):
        
        # extract points from self.points dictonary
        x = [cv.convert(self.points[p]['x'], 'm', 'mm') for p in self.points]
        y = [cv.convert(self.points[p]['y'], 'm', 'mm') for p in self.points]
        # extract arc center points from self.arc_pts dictonary
        arcs_x = [cv.convert(self.arc_pts[arc]['x'], 'm', 'mm') for arc in self.arc_pts]
        arcs_y = [cv.convert(self.arc_pts[arc]['y'], 'm', 'mm') for arc in self.arc_pts]
        
        plt.figure(1, figsize=(10, 5))
        plt.cla()
        plt.plot(x,y, 'r*')                 # [-] Geometry Definition Points on the Wall
        plt.plot(arcs_x,arcs_y, 'g*')       # [-] Arc Definition ceter points
        
        #plt.plot(cv.convert(self.interp_x, 'm', 'mm'), cv.convert(self.interp_y, 'm', 'mm'), 'g*')
        
        plt.axis('equal')
        plt.grid(True)
        plt.title('Wall Geometry')
        plt.xlabel('x-axis [mm]')
        plt.ylabel('y-axis [mm]')








if __name__ == '__main__':
    
    a = Thruster(x_pts=80, y_pts=20)
    
    a.plot_wall()
    
    
    
    
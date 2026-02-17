import numpy as np
import convert as cv
import matplotlib.pyplot as plt
from dataclasses import dataclass
from scipy.interpolate import interp1d


class Thruster:
    def __init__(self, x_pts=80, y_pts=20, l_c=0.002, r_c=0.0015, r_t=0.00015, expan=139, ang_i=44.88, ang_e=15):
        # Grid parameters
        self.x_pts = x_pts
        self.y_pts = y_pts

        # Thruster geometry (units in meters)
        self.l_c = l_c
        self.r_c = r_c
        self.r_t = r_t
        self.expan = expan
        self.r_e = np.sqrt(expan * r_t**2)
        self.ang_i = np.deg2rad(ang_i)
        self.ang_e = np.deg2rad(ang_e)
        self.l_x = np.nan                        # [m] Total length of the thruster
        self.dx = np.nan

        # Geometry points and arc centers
        self.points = {}
        self.arc_pts = {}
        self.interp_pts = {}  # Interpolated wall points as a dictionary

        # Calculate wall geometry
        self.calculate_wall()

    @staticmethod
    def calculate_arc(center_x, center_y, radius, angle_start, angle_end, num_points=50):
        """Calculate an arc given its center, radius, and angular range."""
        angles = np.linspace(angle_start, angle_end, num_points)
        x = center_x + radius * np.cos(angles)
        y = center_y + radius * np.sin(angles)
        return x, y

    def calculate_wall(self):
        # Radii multipliers
        r1, r2, r3 = 10 * self.r_t, 1.5 * self.r_t, 0.382 * self.r_t

        # Define angles for arcs
        alpha0, alphaf = np.radians(90), np.radians(90) - self.ang_i
        beta0, betaf = np.radians(-90) - self.ang_i, np.radians(-90)
        gamma0, gammaf = np.radians(-90), self.ang_e - np.radians(90)

        # Define key points
        self.points['A'] = {'x': 0, 'y': self.r_c}
        self.points['B'] = {'x': self.l_c, 'y': self.r_c}
        self.arc_pts['p1'] = {'x': self.l_c, 'y': self.r_c - r1}

        # Arc BC
        arc_bc_x, arc_bc_y = self.calculate_arc(self.arc_pts['p1']['x'], self.arc_pts['p1']['y'], r1, alpha0, alphaf)
        self.points['C'] = {'x': arc_bc_x[-1], 'y': arc_bc_y[-1]}

        # Arc DE
        m1 = np.tan(np.radians(90) - self.ang_i)
        Dy = (self.r_t + r2) - (r2 * np.sin(alphaf))
        Dx = self.points['C']['x'] - (Dy - self.points['C']['y']) * m1
        self.points['D'] = {'x': Dx, 'y': Dy}
        self.arc_pts['p2'] = {'x': Dx + r2 * np.cos(alphaf), 'y': self.r_t + r2}
        arc_de_x, arc_de_y = self.calculate_arc(self.arc_pts['p2']['x'], self.arc_pts['p2']['y'], r2, beta0, betaf)
        self.points['E'] = {'x': arc_de_x[-1], 'y': self.r_t}

        # Arc EF
        self.arc_pts['p3'] = {'x': self.arc_pts['p2']['x'], 'y': self.r_t + r3}
        arc_ef_x, arc_ef_y = self.calculate_arc(self.arc_pts['p3']['x'], self.arc_pts['p3']['y'], r3, gamma0, gammaf)
        self.points['F'] = {'x': arc_ef_x[-1], 'y': arc_ef_y[-1]}

        # Straight line to G
        m2 = np.tan(self.ang_e)
        Gy = self.r_e
        Gx = self.points['F']['x'] + (Gy - self.points['F']['y']) / m2
        self.points['G'] = {'x': Gx, 'y': Gy}

        # Combine arcs and key points
        x_full = np.concatenate([np.array([self.points['A']['x']]),
                                 arc_bc_x,
                                 np.array([self.points['C']['x'], self.points['D']['x']]),
                                 arc_de_x,
                                 np.array([self.points['E']['x'], self.points['F']['x']]),
                                 arc_ef_x,
                                 np.array([self.points['G']['x']])])
        y_full = np.concatenate([np.array([self.points['A']['y']]),
                                 arc_bc_y,
                                 np.array([self.points['C']['y'], self.points['D']['y']]),
                                 arc_de_y,
                                 np.array([self.points['E']['y'], self.points['F']['y']]),
                                 arc_ef_y,
                                 np.array([self.points['G']['y']])])

        # Interpolate wall points
        distances = np.sqrt(np.diff(x_full)**2 + np.diff(y_full)**2)
        cumulative_length = np.concatenate(([0], np.cumsum(distances)))
        interp_length = np.linspace(0, cumulative_length[-1], self.x_pts)
        interp_x = interp1d(cumulative_length, x_full, kind='linear')(interp_length)
        interp_y = interp1d(cumulative_length, y_full, kind='linear')(interp_length)

        # Store interpolated points in a dictionary
        self.interp_pts = {'x': interp_x, 'y': interp_y}
        
        # Update the geometry data
        self.l_x = self.points['G']['x']

    def plot_wall(self):
        # Plot key points and interpolated wall
        x = [p['x'] * 1e3 for p in self.points.values()]  # Convert to mm
        y = [p['y'] * 1e3 for p in self.points.values()]
        arcs_x = [arc['x'] * 1e3 for arc in self.arc_pts.values()]
        arcs_y = [arc['y'] * 1e3 for arc in self.arc_pts.values()]

        plt.figure(1, figsize=(10, 5))
        plt.cla()
        plt.plot(x, y, 'ro', label='Key Points')
        plt.plot(arcs_x, arcs_y, 'go', label='Arc Centers')
        plt.plot(self.interp_pts['x'] * 1e3, self.interp_pts['y'] * 1e3, 'b*', label='Wall')
        plt.axis('equal')
        plt.grid(True)
        plt.title('Thruster Wall Geometry')
        plt.xlabel('x-axis [mm]')
        plt.ylabel('y-axis [mm]')
        plt.legend()
        plt.show()




if __name__ == '__main__':
    
    a = Thruster(x_pts=80, y_pts=20)
    
    a.plot_wall()
    
    
    
    
    
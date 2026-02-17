import numpy as np
from dataclasses import dataclass
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

@dataclass
class Point:
    x: float
    y: float

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
        self.l_x = np.nan

        # Geometry points
        self.points = {}
        self.interpolated_wall = None

        # Calculate wall geometry
        self.calculate_wall_geometry()

    def calculate_arc(self, center: Point, radius, angle_start, angle_end, num_points=50):
        """Calculate an arc given its center, radius, and angular range."""
        angles = np.linspace(angle_start, angle_end, num_points)
        x = center.x + radius * np.cos(angles)
        y = center.y + radius * np.sin(angles)
        return x, y

    def calculate_wall_geometry(self):
        """Calculate the geometry of the thruster wall."""
        # Radii multipliers
        r1, r2, r3 = 10 * self.r_t, 1.5 * self.r_t, 0.382 * self.r_t

        # Define arcs and key points
        self.points['A'] = Point(0, self.r_c)
        self.points['B'] = Point(self.l_c, self.r_c)
        arc_bc_center = Point(self.l_c, self.r_c - r1)

        # Arc BC
        arc_bc_x, arc_bc_y = self.calculate_arc(arc_bc_center, r1, np.pi / 2, np.pi / 2 - self.ang_i)
        self.points['C'] = Point(arc_bc_x[-1], arc_bc_y[-1])

        # Further points and arcs can be similarly calculated...
        # For brevity, the rest of the geometry computation is omitted.

        # Interpolate wall points
        self.interpolate_wall()

    def interpolate_wall(self):
        """Interpolate the points along the wall."""
        # Combine all geometry points into one array for interpolation
        x_full = np.array([p.x for p in self.points.values()])
        y_full = np.array([p.y for p in self.points.values()])

        # Calculate distances and cumulative lengths
        distances = np.sqrt(np.diff(x_full)**2 + np.diff(y_full)**2)
        cumulative_length = np.concatenate(([0], np.cumsum(distances)))
        interp_length = np.linspace(0, cumulative_length[-1], self.x_pts)

        # Interpolate x and y coordinates
        interp_x = interp1d(cumulative_length, x_full, kind='linear')(interp_length)
        interp_y = interp1d(cumulative_length, y_full, kind='linear')(interp_length)

        self.interpolated_wall = np.column_stack((interp_x, interp_y))

    def plot_wall(self):
        """Plot the thruster wall geometry."""
        plt.figure(figsize=(10, 5))
        if self.interpolated_wall is not None:
            plt.plot(self.interpolated_wall[:, 0] * 1e3, self.interpolated_wall[:, 1] * 1e3, 'b-', label='Wall')
        for label, point in self.points.items():
            plt.plot(point.x * 1e3, point.y * 1e3, 'ro')
            plt.text(point.x * 1e3, point.y * 1e3, label)
        plt.axis('equal')
        plt.grid(True)
        plt.title('Thruster Wall Geometry')
        plt.xlabel('x-axis [mm]')
        plt.ylabel('y-axis [mm]')
        plt.legend()
        plt.show()

# Example usage
if __name__ == '__main__':
    thruster = Thruster()
    thruster.plot_wall()

import numpy as np
import convert as cv
        
class NAP:
    def __init__(self):
        # Grid Parameters
        self.lmax   = 100                   # [-] (int) Number of points in axial direction
        self.mmax   = 100                    # [-] (int) Number of points in radial direction
        self.nmax   = 2                    # [-] (int) Maximum number of time steps 

        # Boundary Condition Parameters
        self.pt     = cv.convert(70, 'psi', 'Pa')    * np.ones(self.mmax)           # [Pa] (float) Inlet Pressure BC
        self.tt     = cv.convert(80, 'degF', 'degC') * np.zeros(self.mmax)          # [K] (float) Inlet Temperature BC
        self.theta  = np.zeros(self.mmax)                                           # [rad] (float) Inlet Flow Angle BC
            
        # Simulation Parameters
        self.ngeom  = 2                    # [-] Specifies Geometry Type; (2 = Circular-arc, conical nozzle)
        self.n1d    = 1                    # [-] (int) Specifies the type of initial surface (1 = subsonic-sonic-supersonic)
        self.ndim   = 1                    # [-] (int) Select flow geometry; 0-two-dimensional, planar flow; 1-axisymmetric

        # Geometry Parameters
        self.xi     = cv.convert(0.31, 'in', 'm')   # [m] (float) Inlet Axial Start Location (x-axis)
        self.ri     = cv.convert(2.5, 'in', 'm')    # [m] (float) Inlet Radial Start Location (y-axis)
        self.rt     = cv.convert(0.8, 'in', 'm')    # [m] (float) Throat Radius
        self.xe     = cv.convert(7.0, 'in', 'm')    # [m] (float) Exit Axial Location (x-axis)
        self.rci    = cv.convert(3.0, 'in', 'm')    # [m] (float) Inlet Radius of Curvature
        self.rct    = cv.convert(0.5, 'in', 'm')    # [m] (float) Throat Radius of Curvature
        self.angi   = 44.88 * np.pi/180.0           # [rad] (float) Nozzle Contraction Angle
        self.ange   = 15.00 * np.pi/180.0           # [rad] (float) Nozzle Exit Angle

        # Fluid Parameters
        self.gamma   = 1.4                              # [-] (int) Specific heat ratio for working fluid (air)
        self.rgas    = 287.0                            # [J/kg-K] (float) Gas constant for working fluid (air)
        self.g       = 9.81                             # [m/s**2] Gravitation constant 

        # Solution Arrays
        self.u       = np.zeros((self.lmax, self.mmax, self.nmax))      # [m/s] (float) Axial Fluid Velocity
        self.v       = np.zeros((self.lmax, self.mmax, self.nmax))      # [m/s] (float) Radial Fluid Velocity
        self.p       = np.zeros((self.lmax, self.mmax, self.nmax))      # [Pa] (float) Pressure field
        self.ro      = np.zeros((self.lmax, self.mmax, self.nmax))      # [kg/m**3] (float) Density field

        # Wall and Centerbody Arrays
        self.xw     = np.zeros(self.lmax)       # [m] (float) Axial Nozzle Wall Coordinates
        self.yw     = np.zeros(self.lmax)       # [m] (float) Radial Nozzle Wall Coordinates
        self.nxny   = np.zeros(self.lmax)       # [?] (float) Wall Normals at the xw,yw point? (Maybe it's interpolated between points?)
        self.dx     = 0.0                       # [m] (float) Axial grid spacing           
        self.dy     = 0.0                       # [m] (float) Radial Grid spacing
        
        # 'Mystery' Parameters 
        self.m1      = self.mmax - 1                    # [?] (int) ? 
        
        # Intermediate Calculation Values
        self.grgas   = 1.0 / (self.rgas * self.g)                       # [kg-K-s^2/J-m] Intermediate Calculation (1/Rgas*g)
        self.gam1    = self.gamma / (self.gamma - 1.0)                  # [-] Intermediate Calculation Value
        self.gam2    = (self.gamma - 1.0) / 2.0                         # [-] Intermediate Calculation Value
        self.acoef   = 2.0 / (self.gamma + 1.0)                         # [-] Intermediate Calculation Value
        self.bcoef   = (self.gamma - 1.0) / (self.gamma + 1.0)          # [-] Intermediate Calculation Value
        self.ccoef   = (self.gamma + 1.0) / 2.0 / (self.gamma - 1.0)    # [-] Intermediate Calculation Value


        # Calculation Flags
        self.nstag      = 0             # [-] (int) Stagnation Point Flag (pretty sure about this description)
        
        # Calculated Values 
        self.xt         = 0.0           # [m] (float) Axial Location of Nozzle Throat
        self.re         = 0.0           # [m] (float) Radial Location of Nozzle Exit
        self.pe         = 0.0           # [Pa] (float) Nozzle Exit Pressure
        self.massi      = 0.0           # [kg/s] (float) Nozzle Inlet Plane Flowrate
        self.masst      = 0.0           # [kg/s] (float) Nozzle Throat Plane Flowrate
        self.masse      = 0.0           # [kg/s] (float) Nozzle Exit Plane Flowrate
        self.thrust     = 0.0           # [N] (float) Nozzle Exit Plane Thrust
        self.lt         = 0             # [-] (int) Calculated nozzle throat axial index location 

    def _circular_arc_nozzle(self):
        ''' Calculate 'Typical Nozzle Type' (Circular-Arc, Conical Nozzle); self.ngeom=2 option '''
        
        if self.rci == 0.0: raise ValueError('RCI must be non-zero for circular-arc nozzle calculations.')
        if self.rct == 0.0: raise ValueError('RCT must be non-zero for circular-arc nozzle calculations.')
                
        # Calculate tangent points
        xtan = self.xi + self.rci * np.sin(self.angi)
        rtan = self.ri + self.rci * (np.cos(self.angi) - 1.0)
        print(xtan, rtan)
        rt1 = self.rt - self.rct * (np.cos(self.angi) - 1.0)
        xt1 = xtan + (rtan - rt1) / np.tan(self.angi)
        
        if xt1 < xtan:
            xt1 = xtan
            rt1 = rtan
        
        xt = xt1 + self.rct * np.sin(self.angi)
        xt2 = xt + self.rct * np.sin(self.ange)
        rt2 = self.rt + self.rct * (1.0 - np.cos(self.ange))
        re = rt2 + (self.xe - xt2) * np.tan(self.ange)
        
        # Store the calculated parameters
        self.xt = xt
        self.re = re
        
        # Find throat index (lt) - closest point to throat position
        dx = (self.xe - self.xi) / (self.lmax - 1)
        self.lt = int(round((xt - self.xi) / dx)) + 1
        
        # Show key parameters
        print(f"Parameters: xi={self.xi}, ri={self.ri}, rt={self.rt}, xe={self.xe}")
        print(f"Calculated: xt={self.xt}, re={self.re}, throat index={self.lt}")
        
        # Generate geometry points
        for l in range(self.lmax):
            x = self.xi + l * dx
            self.xw[l] = x
            
            if self.xi <= x <= xtan:
                # Circular inlet section
                y = self.ri + self.rci * (np.cos(np.arcsin((x - self.xi) / self.rci)) - 1.0)
                nx_ny = (x - self.xi) / (y - self.ri + self.rci)
                
            elif xtan < x <= xt1:
                # Conical section
                y = rt1 + (xt1 - x) * np.tan(self.angi)
                nx_ny = np.tan(self.angi)
                
            elif xt1 < x <= xt:
                # Circular throat inlet section
                y = self.rt + self.rct * (1.0 - np.cos(np.arcsin((xt - x) / self.rct)))
                nx_ny = (xt - x) / (self.rct + self.rt - y)
                
            elif xt < x <= xt2:
                # Circular throat exit section
                y = self.rt + self.rct * (1.0 - np.cos(np.arcsin((x - xt) / self.rct)))
                nx_ny = (xt - x) / (self.rct + self.rt - y)
            
            elif xt2 < x <= self.xe:
                # Conical exit section
                y = rt2 + (x - xt2) * np.tan(self.ange)
                nx_ny = -np.tan(self.ange)
            else:
                y = 0
                nx_ny = 0
                print(f"Warning: Point x={x} is outside the defined geometry range")
           
            self.yw[l] = y
            self.nxny[l] = nx_ny
            self.dx = dx
            self.dy = 1.0 / (self.mmax - 1)

    def plot_nozzle_geometry(self):
        """Plot the nozzle wall geometry with all calculated points."""
        import matplotlib.pyplot as plt
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 5))
        
        # Plot all wall points as dots
        ax.scatter(self.xw, self.yw, s=10, color='blue', zorder=2, label='Wall points')
        
        # Connect the dots with a line to show the wall profile
        ax.plot(self.xw, self.yw, '-', color='black', linewidth=1, zorder=1, alpha=0.5)
        
        # Calculate transition points for better visualization
        xtan = self.xi + self.rci * np.sin(self.angi)
        rtan = self.ri + self.rci * (np.cos(self.angi) - 1.0)
        rt1 = self.rt - self.rct * (np.cos(self.angi) - 1.0)
        xt1 = xtan + (rtan - rt1) / np.tan(self.angi)
        if xt1 < xtan:
            xt1, rt1 = xtan, rtan
        xt = xt1 + self.rct * np.sin(self.angi)
        xt2 = xt + self.rct * np.sin(self.ange)
        rt2 = self.rt + self.rct * (1.0 - np.cos(self.ange))
        
        # Highlight key transition points
        key_points = [
            (self.xi, self.ri, 'xi, ri'),
            (xtan, rtan, 'xtan, rtan'),
            (xt1, rt1, 'xt1, rt1'),
            (xt, self.rt, 'xt, rt'),
            (xt2, rt2, 'xt2, rt2'),
            (self.xe, self.re, 'xe, re')
        ]
        
        for x, y, label in key_points:
            ax.scatter(x, y, s=50, color='red', zorder=3)
            ax.annotate(label, (x, y), xytext=(5, 5), textcoords='offset points')
        
        # Color different sections of the nozzle
        section_colors = {
            "Circular inlet": "lightblue",
            "Conical section": "lightgreen",
            "Throat inlet": "lightyellow",
            "Throat exit": "mistyrose",
            "Conical exit": "lavender"
        }
        
        # Get indices for each section
        sections = [
            ("Circular inlet", [i for i, x in enumerate(self.xw) if self.xi <= x <= xtan]),
            ("Conical section", [i for i, x in enumerate(self.xw) if xtan < x <= xt1]),
            ("Throat inlet", [i for i, x in enumerate(self.xw) if xt1 < x <= xt]),
            ("Throat exit", [i for i, x in enumerate(self.xw) if xt < x <= xt2]),
            ("Conical exit", [i for i, x in enumerate(self.xw) if xt2 < x <= self.xe])
        ]
        
        # Plot each section with a different color
        for name, indices in sections:
            if indices:  # Ensure we have points in this section
                x_section = [self.xw[i] for i in indices]
                y_section = [self.yw[i] for i in indices]
                ax.scatter(x_section, y_section, s=20, color=section_colors[name], 
                          edgecolors='blue', zorder=2, label=name)
        
        # Set labels and title
        ax.set_xlabel('Axial Position')
        ax.set_ylabel('Radial Position')
        ax.set_title('Nozzle Wall Geometry')
        ax.grid(True, alpha=0.3)
        
        # Create a legend
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
        
        # Set equal aspect ratio
        ax.set_aspect('equal')
        
        # Show the plot
        plt.tight_layout()
        plt.show()
        
        return fig, ax  # Return the figure and axes for further customization if needed

    def _calculate_IVS(self):
        # Initialize constants
        nxck   = 0
        
        # Set throat area references
        self.rstar  = self.rt
        self.rstars = self.rt**2 
        
        for l in range(self.lmax):                                       # Main loop over axial stations
            x = self.xi + self.dx * l                               # Calculate the step size in x axis (axial direction)
            
            # Determine Mach number at this station
            if control.n1d < 0: 
                pass                                        # Use subsonic-to-supersonic transition logic
            elif gcb.ngcb != 0:  # With centerbody
                if l+1 < geometry.lt:
                    # Before throat - need to calculate Mach number
                    pass
                elif l+1 > geometry.lt:
                    # After throat
                    if nxck == 1:
                        # Already processed throat transition
                        pass
                    else:
                        # First point after throat - set initial Mach number
                        mn3 = 1.1 if control.n1d in (1, 3) else 0.9  # Supersonic (1.1) or subsonic (0.9) after throat
                        nxck = 1
                else:  # At throat
                    mn3 = 1.0  # Sonic at throat
            elif x < geometry.xt:
                # Before throat - need to calculate Mach number
                pass
            elif x > geometry.xt:
                # After throat
                if nxck == 1:
                    # Already processed throat transition
                    pass
                else:
                    # First point after throat - set initial Mach number
                    if control.n1d == 1 or control.n1d == 3:
                        mn3 = 1.1  # Supersonic after throat
                    else:
                        mn3 = 0.9  # Subsonic after throat
                    nxck = 1
            else:  # At throat
                mn3 = 1.0  # Sonic at throat
            
            # Calculate area ratio for this station
            if geometry.ndim == 1:  # Axisymmetric
                rads = geometry.yw[l]**2 - gcb.ycb[l]**2 if gcb.ngcb != 0 else geometry.yw[l]**2
                aratio = rads / control.rstars
            else:  # 2D planar
                rad = geometry.yw[l] - gcb.ycb[l] if gcb.ngcb != 0 else geometry.yw[l]
                aratio = rad / control.rstar
            
            
            
            # Solve for Mach number using Newton-Raphson method if not at throat
            if abs(x - geometry.xt) > 1e-10:  # Not at throat
                # Initial guess for Mach number if not already set
                if 'mn3' not in locals():
                    mn3 = 2.0 if control.n1d == -1 or control.n1d > 2 else 0.01
                
                # Newton-Raphson iteration
                for iter in range(20):
                    abm = acoef + bcoef * mn3**2
                    abmc = abm**ccoef
                    
                    # Function value: f(MN3) = ABMC/MN3 - ARATIO
                    fm = abmc / mn3 - aratio
                    
                    # Derivative: f'(MN3) = d/dMN3 [ ABM^CCOEF / MN3 ]
                    fpm = abmc * ((2.0 * bcoef * ccoef) / abm - 1.0 / (mn3**2))
                    
                    omn3 = mn3
                    mn3 = omn3 - fm / fpm
                    
                    # Prevent crossing sonic point inappropriately
                    if mn3 > 1.0 and omn3 < 1.0:
                        mn3 = 0.99
                    if mn3 < 1.0 and omn3 > 1.0:
                        mn3 = 1.01
                    
                    # Handle negative Mach numbers
                    if mn3 < 0.0:
                        mn3 = -mn3
                        continue
                    
                    # Check for convergence
                    if abs(mn3 - omn3) / omn3 <= 0.0005:
                        break
                else:
                    print(f"***** THE 1-D SOLUTION FAILED TO CONVERGE IN 20 ITERATIONS AT L={l+1} *****")
            
            # Calculate flow properties and fill in arrays
            dem = 1.0 + control.gam2 * mn3 * mn3
            demp = dem**control.gam1
            dnxny = (geometry.nxny[l] - gcb.nxnycb[l]) / control.m1 if gcb.ngcb != 0 else geometry.nxny[l] / control.m1
            
            # Loop over M index (radial points)
            for m in range(control.mmax):
                # Calculate pressure
                solution.p[l, m, 0] = bcc.pt[m] / demp
                
                # Calculate temperature and density
                temp = bcc.tt[m] / dem
                solution.ro[l, m, 0] = solution.p[l, m, 0] * grgas / temp
                
                # Calculate velocity magnitude
                q = mn3 * np.sqrt(control.gamma * solution.p[l, m, 0] / solution.ro[l, m, 0])
                
                # Calculate flow direction based on wall normal
                if gcb.ngcb != 0:
                    dn = gcb.nxnycb[l] + dnxny * m
                else:
                    dn = dnxny * m
                
                dns = dn * dn
                
                # Resolve velocity components
                if dns < 1e-10:  # Flow along axis
                    solution.u[l, m, 0] = q
                    solution.v[l, m, 0] = 0.0
                else:
                    # Calculate flow direction from wall normal
                    solution.u[l, m, 0] = q / np.sqrt(1.0 + dns)
                    solution.v[l, m, 0] = -np.sign(dn) * solution.u[l, m, 0] * np.abs(dn)
        
        return solution

    def setup(self):
        print('---- Simulation Set up ----')
        print('1  - Set up Geometry [_init_geometry()]')
        self._circular_arc_nozzle()
        self.plot_nozzle_geometry()
        print('2  - Calculate initial value surface (IVS)')
        print('3  - Calculate mass flowrate')
        print('4  - Calculate thrust')
        print('5  - Set up Loop Parameters')
        
    def run(self):
        print('')
        print('---- Main Simulation Loop ----')
        print('6  - Calculate Next Time Step')
        print('7  - Calculate Interior Points')
        print('8  - Calculate Wall Points')
        print('9  - Calculate Inlet Points (if subsonic inlet is set)')
        print('10 - Calculate Exit Points (if subsonic outlet is set)')
        print('11 - Calculate mass flowrate')
        print('12 - Calculate thrust')
        print('13 - Check for Convergence to Steady State Solution')
    

    

if __name__ == '__main__':
  
    sim = NAP()
   
    sim.setup()
    sim.run()    
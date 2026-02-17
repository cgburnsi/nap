import numpy as np
import matplotlib.pyplot as plt
import math
import convert as cv

class GeometryData:
    def __init__(self, lmax=81, nwpts=0):
        """
        Initialize geometry data with flexible array sizes
        
        Parameters:
            lmax (int): Maximum number of points in the x-direction
        """
        self.ndim   = 1                 # [-] (int) Select flow geometry; 0-two-dimensional, planar flow; 1-axisymmetric
        self.ngeom  = None              # [-] (int) Specifies nozzle/wall geometry 
        self.iint   = 1                 # [-] (int) Specifies the order of interpolation used (2 is max)
        self.idif   = 1                 # [-] (int) Specifies the order of differentiation used (5 is max)
        self.jflag  = 0                 # [-] (int) Specifies if an exhaust calculation is needed
        self.xi     = 0.0               # [m] (float) Axial location of nozzle-wall inlet
        self.ri     = 0.0               # [m] (float) Radial location of nozzle-wall inlet
        self.xt     = 0.0               # [m] (float) Axial location of nozzle-wall throat
        self.rt     = 0.0               # [m] (float) Radial location of the nozzle-wall throat
        self.xe     = 0.0               # [m] (float) Axial location of the nozzle-wall exit
        self.re     = 0.0               # [m] (float) Radial location of the nozzle-wall exit
        self.rci    = 0.0               # [m] (float) Radius of curvature of nozzle-wall inlet        
        self.rct    = 0.0               # [m] (float) Radius of curvature of the nozzle-wall throat
        self.angi   = 0.0               # [radian] (float) Angle of converging section of nozzle-wall
        self.ange   = 0.0               # [radian] (float) Angle of diverging section of nozzle-wall        

        # Set up the general nozzle wall geometry (ngeom=3)
        self.nwpts  = nwpts             # [-] (int) Specifies the number of plots in the nozzle-wall
        self.xwi    = np.zeros(nwpts)   # [m] (float) 1D array of non-equally spaced axial coordinates on nozzle-wall
        self.ywi    = np.zeros(nwpts)   # [m] (float) 1D array of non-equally spaced radial coordinates on nozzle-wall

        self.xw     = np.zeros(lmax)
        self.yw     = np.zeros(lmax)
        self.nxny   = np.zeros(lmax)
        self.lt     = 0
    
    def generate_wall(self, control):            # Generates the nozzle wall geometry
        if self.ngeom == 1:
            self.constant_area_duct(control)
        elif self.ngeom == 2:
            self.circular_arc_nozzle(control)
        elif self.ngeom == 3:
            self.general_wall_1(control)
        elif self.ngeom == 4:
            self.general_wall_2(control)
        else:
            print('Wall Geometry not specified')
            
    def constant_area_duct(self, control):           # ngeom = 1 option
        print('Constant Area_Duct')
        
    def circular_arc_nozzle(self, control):          # ngeom = 2 option
        """
        Set up a circular-arc conical nozzle geometry using pre-populated geometry object
        
        Parameters:
            geometry: GeometryData object with parameters already set
            gcb: GcbData object to populate
            params: Dictionary with any remaining parameters
            control: ControlData object
        """
        lmax = control.lmax
        
        print("Circular-Arc, Conical Nozzle Case")
        
        if self.rci == 0.0 or self.rct == 0.0:
            print("Error: RCI and RCT must be non-zero")
            return
        
        # Convert angles to radians
        ani = self.angi * np.pi / 180.0
        ane = self.ange * np.pi / 180.0
        
        # Calculate tangent points
        xtan = self.xi + self.rci * np.sin(ani)
        rtan = self.ri + self.rci * (np.cos(ani) - 1.0)
        
        rt1 = self.rt - self.rct * (np.cos(ani) - 1.0)
        xt1 = xtan + (rtan - rt1) / np.tan(ani)
        
        if xt1 < xtan:
            xt1 = xtan
            rt1 = rtan
        
        xt = xt1 + self.rct * np.sin(ani)
        xt2 = xt + self.rct * np.sin(ane)
        rt2 = self.rt + self.rct * (1.0 - np.cos(ane))
        re = rt2 + (self.xe - xt2) * np.tan(ane)
        
        # Store the calculated parameters
        self.xt = xt
        self.re = re
        
        # Find throat index (lt) - closest point to throat position
        dx = (self.xe - self.xi) / (lmax - 1)
        self.lt = int(round((xt - self.xi) / dx)) + 1
        
        # Show key parameters
        print(f"Parameters: xi={self.xi}, ri={self.ri}, rt={self.rt}, xe={self.xe}")
        print(f"Calculated: xt={self.xt}, re={self.re}, throat index={self.lt}")
        
        # Generate geometry points
        for l in range(lmax):
            x = self.xi + l * dx
            self.xw[l] = x
            
            if self.xi <= x <= xtan:
                # Circular inlet section
                y = self.ri + self.rci * (np.cos(np.arcsin((x - self.xi) / self.rci)) - 1.0)
                nx_ny = (x - self.xi) / (y - self.ri + self.rci)
                
            elif xtan < x <= xt1:
                # Conical section
                y = rt1 + (xt1 - x) * np.tan(ani)
                nx_ny = np.tan(ani)
                
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
                y = rt2 + (x - xt2) * np.tan(ane)
                nx_ny = -np.tan(ane)
            else:
                y = 0
                nx_ny = 0
                print(f"Warning: Point x={x} is outside the defined geometry range")
           
            self.yw[l] = y
            self.nxny[l] = nx_ny
            
    def general_wall_1(self):               # ngeom = 3 option
        print('General Wall Coordinate #1')
        
    def general_wall_2(self):               # geom = 4 option
        print('General Wall Coordinates #2')
        

class GcbData:
    def __init__(self, lmax=81):
        """
        Initialize center body geometry data with flexible array sizes
        
        Parameters:
            lmax (int): Maximum number of points in the x-direction
        """
        self.ngcb = 0
        self.xicb = 0.0
        self.ricb = 0.0
        self.xtcb = 0.0
        self.rtcb = 0.0
        self.xecb = 0.0
        self.recb = 0.0
        self.rcicb = 0.0
        self.rctcb = 0.0
        self.angicb = 0.0
        self.angecb = 0.0
        self.xcb = np.zeros(lmax)
        self.ycb = np.zeros(lmax)
        self.xcbi = np.zeros(lmax)
        self.ycbi = np.zeros(lmax)
        self.nxnycb = np.zeros(lmax)
        self.ncbpts = 0
        self.iintcb = 0
        self.idifcb = 0
        self.lecb = 0

class ControlData:
    def __init__(self):
        self.lmax       = 0         # [-] (int) Number of points in axial direction
        self.mmax       = 0         # [-] (int) Number of points in radial direction
        self.nmax       = 0         # [-] (int) Maximum number of time steps
        self.nprint     = 0
        self.tconv      = 0.0
        self.fdt        = 0.0
        self.gamma      = 0.0       # [-] (fluid) Specific heat ratio for working fluid
        self.rgas       = 0.0       # [J/kg-K] (float) Gas constant for working fluid
        self.gam1       = 0.0
        self.gam2       = 0.0
        self.l1         = 0
        self.l2         = 0
        self.l3         = 0
        self.m1         = 0         # [?] (int) ?
        self.m2         = 0
        self.dx         = 0.0
        self.dy         = 0.0
        self.dt         = 0.0
        self.n          = 0
        self.n1         = 0
        self.n3         = 0
        self.nasm       = 0
        self.ivel       = 0
        self.ichar      = 0
        self.n1d        = 0         # [-] (int) Specifies the type of initial surface (subsonic-sonic-supersonic)
        self.ljet       = 0
        self.jflag      = 0
        self.ierr       = 0
        self.iui        = 0
        self.iuo        = 0
        self.dxr        = 0.0
        self.dyr        = 0.0
        self.ld         = 0
        self.md         = 0
        self.lmd1       = 0
        self.lmd3       = 0
        self.ib         = 0
        self.rstar      = 0.0
        self.rstars     = 0.0
        self.nplot      = 0
        self.g          = 0.0       # [m/s**2] Gravitation constant 
        self.pc         = 0.0
        self.tc         = 0.0
        self.lc         = 0
        self.plow       = 0.0
        self.rolow      = 0.0
    
class BccData:
    def __init__(self, mmax=21):
        """
        Initialize boundary condition data with flexible array sizes
        
        Parameters:
            mmax (int): Maximum number of points in the y-direction
        """
        self.pt = np.zeros(mmax)      # Total pressure array
        self.tt = np.zeros(mmax)      # Total temperature array
        self.theta = np.zeros(mmax)   # Flow angle array
        self.pe = 0.0                 # Exit pressure
        self.masse = 0.0              # Exit mass flow
        self.massi = 0.0              # Inlet mass flow
        self.masst = 0.0              # Total mass flow
        self.thrust = 0.0             # Thrust
        self.nstag = 0                # Stagnation point flag

class SolutionData:
    def __init__(self, lmax=81, mmax=21, nmax=2):
        """
        Initialize solution data arrays with flexible dimensions
        
        Parameters:
            lmax (int): Maximum number of points in the x-direction
            mmax (int): Maximum number of points in the y-direction
            nmax (int): Maximum number of time steps
        """
        self.u = np.zeros((lmax, mmax, nmax))
        self.v = np.zeros((lmax, mmax, nmax))
        self.p = np.zeros((lmax, mmax, nmax))
        self.ro = np.zeros((lmax, mmax, nmax))

def onedim(control, geometry, gcb, bcc, solution):
    """
    Calculate the 1-D initial-data surface with subsonic flow upstream of the throat
    
    Parameters:
        control: ControlData object
        geometry: GeometryData object
        gcb: GcbData object
        bcc: BccData object
        solution: SolutionData object
    
    Returns:
        solution: Updated SolutionData object
    """
    # Initialize variables
    mn3    = 2.0 if control.n1d == -1 or control.n1d > 2 else 0.01  # Set mach number guess to subsonic or supersonic
    grgas  = 1.0 / (control.rgas * control.g)
    nxck   = 0
    acoef  = 2.0 / (control.gamma + 1.0)
    bcoef  = (control.gamma - 1.0) / (control.gamma + 1.0)
    ccoef  = (control.gamma + 1.0) / 2.0 / (control.gamma - 1.0)
    
    # Set rstar and rstars
    control.rstar  = (geometry.yw[geometry.lt - 1]    - gcb.ycb[geometry.lt - 1])    if gcb.ngcb != 0 else geometry.rt
    control.rstars = (geometry.yw[geometry.lt - 1]**2 - gcb.ycb[geometry.lt - 1]**2) if gcb.ngcb != 0 else geometry.rt**2
    
    # Store Mach numbers to improve initial guesses
    mn3_history = np.zeros(control.lmax)
    
    # Find throat location index for reference
    throat_index = 0
    if gcb.ngcb != 0:
        throat_index = geometry.lt - 1
    else:
        # Find the point closest to xt
        for i in range(control.lmax):
            x = geometry.xi + control.dx * i
            if abs(x - geometry.xt) < 1e-10:
                throat_index = i
                break
    
        def newton_raphson_solver(area_ratio, mach_guess, is_subsonic=False, near_throat=False):
            """
            Solves isentropic area-Mach relation using Newton-Raphson method
            
            Parameters:
                area_ratio: Current area ratio A/A*
                mach_guess: Initial guess for Mach number
                is_subsonic: Force solution to stay subsonic if True
                near_throat: Use tighter tolerance if near throat
                
            Returns:
                mach: Converged Mach number
                converged: Boolean indicating if solution converged
            """
            max_iterations = 50
            relaxation_factor = 0.5  # Reduced from 0.7 for better stability
            tolerance = 0.0005
            
            # Use tighter tolerance near throat
            if near_throat:
                tolerance = 0.0001
            
            # Make sure initial guess is appropriate
            if is_subsonic:
                # For subsonic flow, ensure initial guess is well within subsonic range
                mn3 = min(mach_guess, 0.8)
            else:
                # For supersonic flow
                mn3 = max(mach_guess, 1.1)
            
            for iter_count in range(max_iterations):
                abm = acoef + bcoef * mn3**2
                abmc = abm**ccoef
                fm = abmc / mn3 - area_ratio
                fpm = abmc * (2.0 * bcoef * ccoef / abm - 1.0 / mn3**2)
                
                omn3 = mn3
                
                # Apply relaxation
                delta = fm / fpm
                mn3 = omn3 - relaxation_factor * delta
                
                # Enforce constraints
                if is_subsonic:
                    # Force solution to stay subsonic for points upstream of throat
                    if mn3 >= 0.99:
                        mn3 = 0.9
                    elif mn3 < 0.01:
                        mn3 = 0.01
                else:
                    # Force solution to stay supersonic for points downstream of throat
                    if mn3 <= 1.01:
                        mn3 = 1.1
                    elif mn3 > 5.0:  # Cap maximum Mach number
                        mn3 = 5.0
                
                # Check convergence
                if abs(mn3 - omn3) / omn3 <= tolerance:
                    return mn3, True
            
            return mn3, False  # Did not converge    
        
    # Main loop over L index
    for l in range(1, control.lmax + 1):
        # Convert to 0-based indexing for arrays
        l_idx = l - 1
        
        if l == 1 and (control.n1d == -1 or control.n1d > 2):
            continue
            
        x = geometry.xi + control.dx * l_idx
        
        # Calculate area ratio
        if geometry.ndim == 1:
            # Area = π*(r_wall² - r_centerbody²)
            area = math.pi * (geometry.yw[l_idx]**2 - gcb.ycb[l_idx]**2)
            throat_area = math.pi * control.rstars
            aratio = area / throat_area
        else:
            rad = geometry.yw[l_idx] - gcb.ycb[l_idx]
            aratio = rad / control.rstar
        
        print(f"Point {l_idx}: x={x:.4f}, area_ratio={aratio:.4f}, Mach={mn3:.4f}")
        
        # Determine initial Mach number and solution approach based on location
        is_subsonic = (l_idx <= throat_index)
        near_throat = abs(aratio - 1.0) < 0.1
        
        if l_idx <= throat_index:
            # Upstream of or at throat - enforce subsonic flow
            # Initial guess for subsonic flow
            if l > 1 and mn3_history[l_idx - 1] > 0:
                # Use previous point's Mach number as initial guess
                mn3 = min(0.95, mn3_history[l_idx - 1])  # Keep subsonic
            else:
                # Make initial guess based on area ratio
                # Use isentropic flow relation for initial guess
                # For subsonic flow, larger area = lower Mach number
                mn3 = 0.3 + 0.6 * (1.0 - aratio)  # Will be between 0.3 and 0.9
            
            # At throat, set to exactly sonic condition
            if l_idx == throat_index:
                mn3 = 0.95  # Slightly subsonic at throat
        else:
            # Downstream of throat - allow supersonic flow
            if nxck == 0:
                if control.n1d == 1 or control.n1d == 3:
                    mn3 = 1.1  # Start supersonic
                elif control.n1d == 2 or control.n1d == 4:
                    mn3 = 0.9  # Start subsonic
                nxck = 1
            
            if l_idx > throat_index + 1 and mn3_history[l_idx - 1] > 0:
                mn3 = mn3_history[l_idx - 1]  # Use previous solution
        
        # Replace the existing solver call with this:
        # Call Newton-Raphson solver with clear subsonic/supersonic enforcement
        if l_idx <= throat_index:
            # Points at or upstream of throat MUST be subsonic 
            # (except exactly at throat where it's sonic)
            mn3, converged = newton_raphson_solver(aratio, mn3, is_subsonic=True, near_throat=near_throat)
            
            # At throat, set to exactly sonic condition if needed
            if near_throat and abs(aratio - 1.0) < 0.01:
                mn3 = 0.99  # Almost sonic at throat
        else:
            # Points downstream of throat can be supersonic
            mn3, converged = newton_raphson_solver(aratio, mn3, is_subsonic=False, near_throat=near_throat)        
        
        if not converged:
            print(f"***** THE 1-D SOLUTION FOR THE INITIAL-DATA SURFACE FAILED TO CONVERGE AT L={l} *****")
        
        # Store the converged Mach number for next iteration
        mn3_history[l_idx] = mn3
        
        # Fill in 2-D arrays
        dem = 1.0 + control.gam2 * mn3 * mn3
        demp = dem**control.gam1
        dnxny = (geometry.nxny[l_idx] - gcb.nxnycb[l_idx]) / control.m1
        
        for m in range(1, control.mmax + 1):
            # Convert to 0-based indexing
            m_idx = m - 1
            
            solution.p[l_idx, m_idx, 0] = bcc.pt[m_idx] / demp
            temp = bcc.tt[m_idx] / dem
            solution.ro[l_idx, m_idx, 0] = solution.p[l_idx, m_idx, 0] * grgas / temp
            
            q = mn3 * math.sqrt(control.gamma * solution.p[l_idx, m_idx, 0] / solution.ro[l_idx, m_idx, 0])
            dn = gcb.nxnycb[l_idx] + dnxny * m_idx
            dns = dn * dn
            
            if abs(dns) < 1e-10:  # Equivalent to DNS.EQ.0.0 in Fortran
                solution.u[l_idx, m_idx, 0] = q
                solution.v[l_idx, m_idx, 0] = 0.0
            else:
                sign = 1.0
                if dn > 0.0:
                    sign = -1.0
                solution.u[l_idx, m_idx, 0] = q / math.sqrt(1.0 + dns)
                solution.v[l_idx, m_idx, 0] = sign * q / math.sqrt(1.0 + 1.0 / dns)
    
    return solution

def plot_nozzle_geometry(geometry, gcb=None, ax=None):
    """
    Plot the nozzle geometry
    
    Parameters:
        geometry: GeometryData object with the nozzle geometry
        gcb: Optional GcbData object with center body data
        ax: Optional matplotlib axes to plot on
    
    Returns:
        ax: The matplotlib axes with the plot
    """
    import matplotlib.pyplot as plt
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    
    # Extract x and y coordinates from the geometry object
    x_values = [geometry.xw[l] for l in range(len(geometry.xw)) if geometry.xw[l] > 0]
    y_values = [geometry.yw[l] for l in range(len(geometry.yw)) if geometry.xw[l] > 0]
    
    # Plot the wall profile
    ax.plot(x_values, y_values, 'b-', linewidth=2, label='Wall')
    
    # Plot centerline (y=0)
    ax.plot([min(x_values), max(x_values)], [0, 0], 'k--', label='Centerline')
    
    # Mirror the profile for visualization (for axisymmetric nozzles)
    ax.plot(x_values, [-y for y in y_values], 'b-', linewidth=2)
    
    # Plot center body if provided and present
    if gcb is not None and gcb.ngcb != 0:
        cb_x = [gcb.xcb[l] for l in range(len(gcb.xcb)) if gcb.xcb[l] > 0]
        cb_y = [gcb.ycb[l] for l in range(len(gcb.ycb)) if gcb.xcb[l] > 0]
        ax.plot(cb_x, cb_y, 'r-', linewidth=2, label='Center Body')
        ax.plot(cb_x, [-y for y in cb_y], 'r-', linewidth=2)
    
    # Mark throat position
    if hasattr(geometry, 'xt') and geometry.xt > 0:
        ax.axvline(x=geometry.xt, color='r', linestyle='--', label='Throat')
    
    # Set plot properties
    ax.grid(True)
    ax.set_aspect('equal')
    ax.set_title('Nozzle Geometry')
    ax.set_xlabel('Axial Distance')
    ax.set_ylabel('Radial Distance')
    ax.legend()
    
    return ax

def create_test_data():
    ''' Creates class instances (control, geometry, gcb, bcc, and solution) '''
    
    # Create control data
    control         = ControlData()
    control.n1d     = 1                 # [-] (int) Specifies the type of initial surface (subsonic-sonic-supersonic)
    control.lmax    = 20                # [-] (int) Number of points in axial direction
    control.mmax    = 7                # [-] (int) Number of points in radial direction
    control.nmax    = 2                 # [-] (int) Maximum number of time steps 
    control.gamma   = 1.4               # [-] (int) Specific heat ratio for working fluid (air)
    control.rgas    = 287.0             # [J/kg-K] (float) Gas constant for working fluid (air)
    control.g       = 9.81              # [m/s**2] Gravitation constant 
    control.m1      = control.mmax - 1  # [?] (int) ?
    control.gam1    = control.gamma / (control.gamma - 1.0)
    control.gam2    = (control.gamma - 1.0) / 2.0
    
    # Create geometry data
    geometry        = GeometryData(control.lmax)
    gcb             = GcbData(control.lmax)        
    geometry.ndim   = 1                             # Axisymmetric flow
    geometry.ngeom  = 2                             # Circular-arc, conical nozzle
    geometry.xi     = cv.convert(0.31, 'in', 'm')   # Inlet x position
    geometry.ri     = cv.convert(2.5, 'in', 'm')    # Inlet radius
    geometry.rt     = cv.convert(0.8, 'in', 'm')    # Throat radius
    geometry.xe     = cv.convert(7.0, 'in', 'm')    # Exit x position
    geometry.rci    = cv.convert(3.0, 'in', 'm')    # Inlet radius of curvature
    geometry.rct    = cv.convert(0.5, 'in', 'm')    # Throat radius of curvature
    geometry.angi   = 44.88                         # Inlet angle (degrees)
    geometry.ange   = 15.0                          # Exit angle (degrees)

    geometry.generate_wall(control)                 # Generate the wall points
    
    # Set common control parameters based on the geometry
    control.dx = (geometry.xe - geometry.xi) / (control.lmax - 1)
    
    # Create BCC data with the correct size
    bcc = BccData(control.mmax)
    # Set uniform inlet conditions
    bcc.pt.fill(cv.convert(70, 'psi', 'Pa'))        # Set inlet gas pressure to 70 psig
    bcc.tt.fill(cv.convert(80, 'degF', 'degC'))     # Set inlet gas temperature to 80 degF
    
    # Create solution data
    solution = SolutionData(control.lmax, control.mmax, control.nmax)
    
    return control, geometry, gcb, bcc, solution

def main():
    # Create test data with a circular-arc conical nozzle
    control, geometry, gcb, bcc, solution = create_test_data()
    
    # Run the onedim function
    solution = onedim(control, geometry, gcb, bcc, solution)
    
    # Create figure for plots
    fig = plt.figure(figsize=(16, 12))
    
    # Plot 1: Nozzle geometry
    ax1 = fig.add_subplot(221)
    plot_nozzle_geometry(geometry, gcb, ax=ax1)
    
    # Extract results for visualization
    lmax = control.lmax
    mmax = control.mmax
    
    # Create proper structured grid for plotting
    x_1d = np.array([geometry.xw[l] for l in range(lmax)])
    y_max = np.array([geometry.yw[l] for l in range(lmax)])
    
    # Create evenly spaced y-coordinates for each x position
    y_1d = np.linspace(0, 1, mmax)  # Normalized y coordinates
    
    # Create meshgrid for plotting
    X, Y_norm = np.meshgrid(x_1d, y_1d, indexing='ij')
    
    # Scale Y by the wall profile to create the actual y coordinates
    Y = Y_norm * y_max[:, np.newaxis]
    
    # Extract flow variables
    u = solution.u[:, :, 0]
    v = solution.v[:, :, 0]
    p = solution.p[:, :, 0]
    rho = solution.ro[:, :, 0]
    
    # Calculate Mach number
    speed = np.sqrt(u**2 + v**2)
    sound_speed = np.sqrt(control.gamma * p / rho)
    mach = np.zeros_like(speed)
    # Avoid division by zero
    mask = sound_speed > 1e-10
    mach[mask] = speed[mask] / sound_speed[mask]
    
    # Plot 2: Pressure distribution
    ax2 = fig.add_subplot(222)
    c2 = ax2.contourf(X, Y, p, 50, cmap=plt.cm.viridis)
    ax2.set_title('Pressure Distribution')
    ax2.set_xlabel('X Position')
    ax2.set_ylabel('Y Position')
    plt.colorbar(c2, ax=ax2, label='Pressure (Pa)')
    
    # Plot 3: Mach number
    ax3 = fig.add_subplot(223)
    c3 = ax3.contourf(X, Y, mach, 50, cmap=plt.cm.viridis)
    ax3.set_title('Mach Number Distribution')
    ax3.set_xlabel('X Position')
    ax3.set_ylabel('Y Position')
    plt.colorbar(c3, ax=ax3, label='Mach Number')
    
    # Plot 4: Flow vectors
    ax4 = fig.add_subplot(224)
    c4 = ax4.contourf(X, Y, mach, 20, cmap=plt.cm.viridis, alpha=0.6)
    
    # Use quiver to show flow direction
    skip = (slice(None, None, 3), slice(None, None, 3))
    ax4.quiver(X[skip], Y[skip], u[skip], v[skip], 
               color='black', scale=50, width=0.002, headwidth=4)
    
    plt.colorbar(c4, ax=ax4, label='Mach Number')
    ax4.set_title('Flow Vectors')
    ax4.set_xlabel('X Position')
    ax4.set_ylabel('Y Position')
    
    plt.tight_layout()
    plt.savefig('nozzle_flow_results.png', dpi=300)
    plt.show()
    
    # Print some flow statistics
    throat_index = 0
    for i in range(lmax):
        x = geometry.xi + control.dx * i
        if abs(x - geometry.xt) < control.dx:
            throat_index = i
            break
    
    inlet_mach = np.mean(mach[0, :])
    throat_mach = np.mean(mach[throat_index, :])
    exit_mach = np.mean(mach[-1, :])
    
    print("\nFlow Statistics:")
    print(f"Inlet Mach number: {inlet_mach:.3f}")
    print(f"Throat Mach number: {throat_mach:.3f}")
    print(f"Exit Mach number: {exit_mach:.3f}")
    
    print(f"\nInlet pressure: {np.mean(p[0, :]):.1f} Pa")
    print(f"Throat pressure: {np.mean(p[throat_index, :]):.1f} Pa")
    print(f"Exit pressure: {np.mean(p[-1, :]):.1f} Pa")
    
    print("\nCalculation complete. Results visualized.")


if __name__ == "__main__":
    main()
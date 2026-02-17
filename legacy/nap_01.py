import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
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
        self.rt     = 0.0               # [m] (float) Radial location of the nozzle-wall throat
        self.xe     = 0.0               # [m] (float) Axial location of the nozzle-wall exit
        self.rci    = 0.0               # [m] (float) Radius of curvature of nozzle-wall inlet        
        self.rct    = 0.0               # [m] (float) Radius of curvature of the nozzle-wall throat
        self.angi   = 0.0               # [radian] (float) Angle of converging section of nozzle-wall
        self.ange   = 0.0               # [radian] (float) Angle of diverging section of nozzle-wall        
        self.re     = 0.0

        # Set up the general nozzle wall geometry (ngeom=3)
        self.nwpts  = nwpts             # [-] (int) Specifies the number of plots in the nozzle-wall
        self.xwi    = np.zeros(nwpts)   # [m] (float) 1D array of non-equally spaced axial coordinates on nozzle-wall
        self.ywi    = np.zeros(nwpts)   # [m] (float) 1D array of non-equally spaced radial coordinates on nozzle-wall

        self.xw     = np.zeros(lmax)
        self.yw     = np.zeros(lmax)
        self.nxny   = np.zeros(lmax)
        self.lt     = 0
        self.xt     = 0.0

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
        self.lmax = 0
        self.mmax = 0
        self.nmax = 0
        self.nprint = 0
        self.tconv = 0.0
        self.fdt = 0.0
        self.gamma = 0.0
        self.rgas = 0.0
        self.gam1 = 0.0
        self.gam2 = 0.0
        self.l1 = 0
        self.l2 = 0
        self.l3 = 0
        self.m1 = 0
        self.m2 = 0
        self.dx = 0.0
        self.dy = 0.0
        self.dt = 0.0
        self.n = 0
        self.n1 = 0
        self.n3 = 0
        self.nasm = 0
        self.ivel = 0
        self.ichar = 0
        self.n1d = 0
        self.ljet = 0
        self.jflag = 0
        self.ierr = 0
        self.iui = 0
        self.iuo = 0
        self.dxr = 0.0
        self.dyr = 0.0
        self.ld = 0
        self.md = 0
        self.lmd1 = 0
        self.lmd3 = 0
        self.ib = 0
        self.rstar = 0.0
        self.rstars = 0.0
        self.nplot = 0
        self.g = 0.0
        self.pc = 0.0
        self.tc = 0.0
        self.lc = 0
        self.plow = 0.0
        self.rolow = 0.0

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
    mn3 = 0.01
    if control.n1d == -1 or control.n1d > 2:
        mn3 = 2.0
    
    grgas = 1.0 / (control.rgas * control.g)
    nxck = 0
    acoef = 2.0 / (control.gamma + 1.0)
    bcoef = (control.gamma - 1.0) / (control.gamma + 1.0)
    ccoef = (control.gamma + 1.0) / 2.0 / (control.gamma - 1.0)
    
    # Set rstar and rstars
    if gcb.ngcb != 0:
        control.rstar = geometry.yw[geometry.lt - 1] - gcb.ycb[geometry.lt - 1]
        control.rstars = geometry.yw[geometry.lt - 1]**2 - gcb.ycb[geometry.lt - 1]**2
    else:
        control.rstar = geometry.rt
        control.rstars = geometry.rt**2
    
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
    
    # Main loop over L index
    for l in range(1, control.lmax + 1):
        # Convert to 0-based indexing for arrays
        l_idx = l - 1
        
        if l == 1 and (control.n1d == -1 or control.n1d > 2):
            continue
            
        x = geometry.xi + control.dx * l_idx
        
        # Calculate area ratio
        if geometry.ndim == 1:
            rads = geometry.yw[l_idx]**2 - gcb.ycb[l_idx]**2
            aratio = rads / control.rstars
        else:
            rad = geometry.yw[l_idx] - gcb.ycb[l_idx]
            aratio = rad / control.rstar
        
        # Determine initial Mach number and solution approach based on location
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
        
        # Newton-Raphson iteration
        max_iterations = 50
        relaxation_factor = 0.7
        
        for iter_count in range(max_iterations):
            abm = acoef + bcoef * mn3**2
            abmc = abm**ccoef
            fm = abmc / mn3 - aratio
            fpm = abmc * (2.0 * bcoef * ccoef / abm - 1.0 / mn3**2)
            
            omn3 = mn3
            
            # Apply relaxation
            delta = fm / fpm
            mn3 = omn3 - relaxation_factor * delta
            
            # For points upstream of throat, force subsonic flow
            if l_idx < throat_index and mn3 >= 1.0:
                mn3 = 0.95  # Keep subsonic
            
            # Ensure Mach number remains positive and reasonable
            if mn3 < 0.01:
                mn3 = 0.01
            if mn3 > 5.0:
                mn3 = 5.0  # Cap maximum Mach number
            
            # Check convergence
            tolerance = 0.0005
            if abs(aratio - 1.0) < 0.1:
                tolerance = 0.0001  # Tighter tolerance near throat
                
            if abs(mn3 - omn3) / omn3 <= tolerance:
                break
        else:
            print(f"***** THE 1-D SOLUTION FOR THE INITIAL-DATA SURFACE FAILED TO CONVERGE IN {max_iterations} ITERATIONS AT L={l} *****")
            # Use last value anyway and continue
        
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

# Now let's create an example that demonstrates the use of this function

def create_test_data():
    """Create test data for the onedim function"""
    
    # Create geometry data
    geometry = GeometryData()
    geometry.ndim       = 1             # [-] (int) Set flow geometry to axisymmetric
    geometry.ngeom      = 2                                 # [-] (int) Set geometry for circular-arc, conical nozzle wall
    geometry.xi         = cv.convert(0.31, 'in', 'm')       # [m] (float) Axial location of nozzle-wall inlet
    geometry.ri         = cv.convert(2.5, 'in', 'm')        # [m] (float) Radial location of nozzle-wall inlet
    geometry.rt         = cv.convert(0.8, 'in', 'm')        # [m] (float) Radial location of nozzle-wall throat        
    geometry.xt = 1.0      # Throat x position
    geometry.lt = 20       # Index of throat position
    
    # Create control data
    control = ControlData()
    control.lmax = 40      # Number of points in x-direction
    control.mmax = 15      # Number of points in y-direction
    control.nmax = 2       # Number of time steps (only using first one here)
    control.gamma = 1.4    # Specific heat ratio for air
    control.rgas = 287.0   # Gas constant for air (J/kg·K)
    control.g = 9.81       # Gravity (m/s²)
    control.n1d = 1        # Flag for 1D solution
    control.gam1 = control.gamma / (control.gamma - 1.0)
    control.gam2 = (control.gamma - 1.0) / 2.0
    control.dx = 0.05      # Grid spacing
    control.m1 = control.mmax - 1
    

    
    # Create wall geometry
    x_range = np.linspace(0.0, 2.0, control.lmax)
    
    # Simple converging-diverging nozzle profile
    def nozzle_profile(x):
        if x < geometry.xt:
            # Converging section - cosine profile
            return 1.0 - 0.5 * (1.0 - np.cos((x / geometry.xt) * np.pi/2))
        else:
            # Diverging section - parabolic
            return 0.5 + 0.3 * ((x - geometry.xt) / geometry.xt)**2
    
    for i, x in enumerate(x_range):
        geometry.xw[i] = x
        geometry.yw[i] = nozzle_profile(x)
        # Simple wall slope
        if i > 0:
            geometry.nxny[i] = (geometry.yw[i] - geometry.yw[i-1]) / (geometry.xw[i] - geometry.xw[i-1])
        else:
            geometry.nxny[i] = 0.0
    
    # Create GCB data (center body)
    gcb = GcbData()
    gcb.ngcb = 0
    
    # Create BCC data with the correct size
    bcc = BccData(control.mmax)
    # Set uniform inlet conditions
    bcc.pt.fill(101325.0)  # Standard atmospheric pressure (Pa)
    bcc.tt.fill(300.0)     # Standard temperature (K)
    
    # Create solution data
    solution = SolutionData(control.lmax, control.mmax, control.nmax)
    
    return control, geometry, gcb, bcc, solution


def main():
    # Create test data
    control, geometry, gcb, bcc, solution = create_test_data()
    # Run the onedim function
    solution = onedim(control, geometry, gcb, bcc, solution)
    
    # Extract results for visualization
    lmax = control.lmax
    mmax = control.mmax
    
    # Create proper structured grid for plotting
    x_1d = np.array([geometry.xw[l] for l in range(lmax)])
    y_max = np.array([geometry.yw[l] for l in range(lmax)])
    
    # Create evenly spaced y-coordinates for each x position
    y_1d = np.linspace(0, 1, mmax)  # Normalized y coordinates
    
    # Create meshgrid for plotting - note the indexing for correct orientation
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
    
    # Create plots
    plt.figure(figsize=(14, 10))
    
    # Plot 1: Pressure
    plt.subplot(2, 2, 1)
    plt.contourf(X, Y, p, 50, cmap=cm.viridis)
    plt.colorbar(label='Pressure (Pa)')
    plt.title('Pressure Distribution')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    
    # Plot 2: Velocity magnitude
    plt.subplot(2, 2, 2)
    plt.contourf(X, Y, speed, 50, cmap=cm.viridis)
    plt.colorbar(label='Velocity (m/s)')
    plt.title('Velocity Distribution')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    
    # Plot 3: Mach number
    plt.subplot(2, 2, 3)
    plt.contourf(X, Y, mach, 50, cmap=cm.viridis)
    plt.colorbar(label='Mach Number')
    plt.title('Mach Number Distribution')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    
    # Plot 4: For streamplot, we'll use quiver instead which is more flexible
    plt.subplot(2, 2, 4)
    plt.contourf(X, Y, mach, 20, cmap=cm.viridis, alpha=0.6)
    
    # Use quiver to show flow direction instead of streamplot
    # Sample the grid to avoid too many arrows
    skip = (slice(None, None, 3), slice(None, None, 3))
    plt.quiver(X[skip], Y[skip], u[skip], v[skip], 
               color='black', scale=50, width=0.002, headwidth=4)
    
    plt.colorbar(label='Mach Number')
    plt.title('Flow Vectors')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    
    plt.tight_layout()
    plt.savefig('onedim_results.png', dpi=300)
    plt.show()
    
    print("Calculation complete. Results visualized.")
    
if __name__ == "__main__":
    main()
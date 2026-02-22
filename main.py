import numpy as np
import matplotlib.pyplot as plt

def plot_wall_grid(geometry):
    """
    Plots only the top wall geometry and its specific grid points.
    """   
    # Plot the individual grid points on the wall
    plt.plot(geometry['xw'], geometry['yw'], 'ro', markersize=2, label='Wall Grid Points')

    # Plot the wall line for continuity
    plt.plot(geometry['xw'], geometry['yw'], color='blue', linewidth=1.5, label='Wall Line')
   
    plt.xlabel('X (Axial Distance)')
    plt.ylabel('Y (Radius)')
    plt.title('Nozzle Wall Geometry - Top Half Grid Points')
    plt.legend()
    plt.grid(True)
    plt.axis('equal') # Ensures physical proportions are maintained
    plt.show()
   

def plot_isentropic_data(u, v, p, d, geometry, params):
    """
    Generates a 2x2 grid plot of the initial 1D isentropic flow properties.
    """
    gamma = params['gamma']
    g = params['g']  # 32.174
    x = geometry['xw']
    
    # Extract values for the first radial grid point (m=0) and first time step (n=0)
    pres = p[:, 0, 0]
    dens = d[:, 0, 0]
    vel_u = u[:, 0, 0]
    
    # Corrected speed of sound calculation (psia -> psf conversion)
    # a = sqrt(gamma * P_psia * 144 * g / rho)
    a_sound = np.sqrt(gamma * pres * 144.0 * g / dens)
    
    # Corrected Mach calculation
    mach = vel_u / a_sound

    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    
    # Subplot 1: Mach Number
    axs[0, 0].plot(x, mach, 'b-', label=r'$M$')
    axs[0, 0].set_title('Mach Number')
    axs[0, 0].set_ylabel(r'$M$')
    axs[0, 0].grid(True)
    
    # Subplot 2: Pressure
    axs[0, 1].plot(x, pres, 'r-', label=r'$P$')
    axs[0, 1].set_title('Pressure (psia)')
    axs[0, 1].set_ylabel(r'$P$ (psia)')
    axs[0, 1].grid(True)
    
    # Subplot 3: Density
    axs[1, 0].plot(x, dens, 'g-', label=r'$d$')
    axs[1, 0].set_title('Density ($slug/ft^3$)')
    axs[1, 0].set_ylabel(r'$d$ ($slug/ft^3$)')
    axs[1, 0].grid(True)
    
    # Subplot 4: Axial Velocity
    axs[1, 1].plot(x, vel_u, 'k-', label=r'$U$')
    axs[1, 1].set_title('Axial Velocity (ft/s)')
    axs[1, 1].set_ylabel(r'$U$ (ft/s)')
    axs[1, 1].grid(True)
    
    for ax in axs.flat:
        ax.set_xlabel('X (Axial Distance)')
        ax.legend()
    
    plt.tight_layout()
    plt.show()




def plot_computational_domain(geometry, params):
    """
    Plots the actual calculated grid nodes in both physical and 
    computational space using the geometry dictionary.
    """
    lmax = params['lmax']
    mmax = params['mmax']
    dy = 1.0 / (mmax - 1)
    
    # 1. Create Physical Space Plot (X, Y)
    plt.figure(figsize=(12, 10))
    plt.subplot(2, 1, 1)
    
    for l in range(lmax):
        x = geometry['xw'][l]
        be = geometry['be'][l]
        y_cb = geometry['ycb'][l]
        
        # Calculate physical Y for every radial point at this axial station
        for m in range(mmax):
            eta = m * dy
            y_phys = (eta / be) + y_cb
            
            # Plot the node
            plt.plot(x, y_phys, 'ro', markersize=1)
            
        # Draw vertical grid lines (constant xi)
        y_wall = geometry['yw'][l]
        plt.plot([x, x], [y_cb, y_wall], 'k-', linewidth=0.2)

    # Draw horizontal grid lines (constant eta)
    for m in range(mmax):
        eta = m * dy
        y_line = []
        for l in range(lmax):
            y_line.append((eta / geometry['be'][l]) + geometry['ycb'][l])
        plt.plot(geometry['xw'], y_line, 'b-', linewidth=0.5)

    plt.title('Physical Domain Mesh (Actual Calculated Nodes)')
    plt.xlabel('X (Axial)')
    plt.ylabel('Y (Radial)')
    plt.axis('equal')
    plt.grid(False)

    # 2. Create Computational Space Plot (L, M)
    plt.subplot(2, 1, 2)
    for l in range(lmax):
        for m in range(mmax):
            plt.plot(l, m, 'go', markersize=1)
            
    plt.title('Computational Domain Mesh (Indices L, M)')
    plt.xlabel('Axial Index (L)')
    plt.ylabel('Radial Index (M)')
    plt.grid(True, linestyle=':', alpha=0.5)
    
    plt.tight_layout()
    plt.show()
    
    
    
    
    
    
    


def setup_geometry(params):
    # Step 1c - Geometry Setup GEOM (Define nozzle wall)
    # This routine calculates the physical coordinates and slopes based on geom.f.
    lmax = params['lmax']
    xi = params['xi']
    xe = params['xe']
    ngeom = params['ngeom']
    
    # Initialize arrays for wall coordinates (xw, yw), wall slopes (nxny), 
    # and x-derivatives (xwi).
    geometry = {
        'xw': np.zeros(lmax),
        'yw': np.zeros(lmax),
        'nxny': np.zeros(lmax),
        'xwi': np.ones(lmax), 
        'ycb': np.zeros(lmax),
        'nxnycb': np.zeros(lmax)
    }
    
    # Calculate grid spacing 
    dx = (xe - xi) / (lmax - 1)
    params['dx'] = dx

    # Case Selection based on NGEOM 
    if ngeom == 1:
        # Case 10: Constant Area Duct 
        ri = params['ri']
        for l in range(lmax):
            geometry['xw'][l] = xi + l * dx
            geometry['yw'][l] = ri
            geometry['nxny'][l] = 0.0
            
    elif ngeom == 2:
        # Case 30: Circular-Arc, Conical Nozzle 
        ri = params['ri']
        rt = params['rt']
        rci = params['rci']
        rct = params['rct']
        # Convert angles to radians 
        ani = params['angi'] * np.pi / 180.0
        ane = params['ange'] * np.pi / 180.0
        
        # Calculate transition (tangent) points for the geometry sections 
        xtan = xi + rci * np.sin(ani)
        rtan = ri + rci * (np.cos(ani) - 1.0)
        rt1 = rt - rct * (np.cos(ani) - 1.0)
        xt1 = xtan + (rtan - rt1) / np.tan(ani)
        
        if xt1 < xtan: 
            xt1 = xtan
            rt1 = rtan
            
        xt_throat = xt1 + rct * np.sin(ani)
        xt2 = xt_throat + rct * np.sin(ane)
        rt2 = rt + rct * (1.0 - np.cos(ane))
        
        params['xt'] = xt_throat
        
        # Loop through axial grid points 
        for l in range(lmax):
            x = xi + l * dx
            geometry['xw'][l] = x
            
            # Identify nozzle section 
            if x <= xtan:
                # Section 50: Inlet Circular Arc 
                # YW calculation using circular arc formula
                geometry['yw'][l] = ri + rci * (np.cos(np.arcsin((x - xi) / rci)) - 1.0)
                
                # CORRECTED Slope: NXNY(L) = (X - XI) / (YW(L) - RI + RCI)
                # This ensures the normal vector is calculated correctly for the arc.
                geometry['nxny'][l] = (x - xi) / (geometry['yw'][l] - ri + rci)
                
            elif x <= xt1:
                # Section 60: Conical Converging 
                geometry['yw'][l] = rt1 + (xt1 - x) * np.tan(ani)
                geometry['nxny'][l] = np.tan(ani)
                
            elif x <= xt_throat:
                # Section 70: Converging Throat Arc
                geometry['yw'][l] = rt - rct * (np.cos(np.arcsin((x - xt_throat) / rct)) - 1.0)
                geometry['nxny'][l] = (x - xt_throat) / (geometry['yw'][l] - rt + rct)
                
            elif x <= xt2:
                # Section 80: Diverging Throat Arc
                geometry['yw'][l] = rt + rct * (1.0 - np.cos(np.arcsin((x - xt_throat) / rct)))
                geometry['nxny'][l] = (xt_throat - x) / (geometry['yw'][l] - rt - rct)
                
            else:
                # Section 90: Conical Expansion
                geometry['yw'][l] = rt2 + (x - xt2) * np.tan(ane)
                geometry['nxny'][l] = -np.tan(ane)  
                
    return geometry


def solve_mach_newton(area_ratio, gamma, is_supersonic, mn_guess):
    """
    Newton-Raphson solver for the Area-Mach relation.
    Optimized constants from the NAP source code.
    """
    # Optimized isentropic constants
    acoef = 2.0 / (gamma + 1.0)
    bcoef = (gamma - 1.0) / (gamma + 1.0)
    ccoef = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    
    aratio = max(area_ratio, 1.0)
    if aratio == 1.0:
        return 1.0
        
    mn = mn_guess
    for _ in range(30):
        # f(m) and f'(m) exactly as defined in the 1977 report
        abm = acoef + bcoef * mn**2
        if abm <= 0: break
        
        abmc = abm**ccoef
        f_m = abmc / mn - aratio
        df_m = abmc * (2.0 * bcoef * ccoef * mn**2 / abm - 1.0) / mn**2
        
        mn_new = mn - f_m / df_m
        
        # Branch enforcement: prevents jumping across the sonic point
        if is_supersonic and mn_new < 1.0: mn_new = 1.01
        if not is_supersonic and mn_new > 1.0: mn_new = 0.99
            
        if abs(mn_new - mn) < 0.00001:
            return mn_new
        mn = mn_new
        
    return mn

def set_initial_conditions(u, v, p, d, params, geometry):
    """
    Step 1d: Fills the grid with 1D isentropic initial values using 
    English Engineering units.
    """
    gamma = params['gamma']
    rgas = params['rgas']
    g = params['g']
    pc = params['pc']
    tc = params['tc']
    ndim = params.get('ndim', 1) # 1 for Axisymmetric, 0 for Planar
    
    # Unit Scaling: Convert pc from psia to psf to get density in lbm/ft^3
    pc_psf = pc * 144.0
    dc = pc_psf / (rgas * tc)
    
    r_star = np.min(geometry['yw'])
    throat_idx = np.argmin(geometry['yw'])
    
    # Initialize Mach guess (subsonic at inlet)
    mn_guess = 0.1
    
    for l in range(params['lmax']):
        # Area Ratio logic: Axisymmetric uses Y^2, Planar uses Y
        if ndim == 1:
            aratio = (geometry['yw'][l]**2) / (r_star**2)
        else:
            aratio = geometry['yw'][l] / r_star
            
        is_supersonic = l > throat_idx
        
        if l == throat_idx:
            m_l = 1.0
            mn_guess = 1.1 # Force expansion branch immediately after throat
        else:
            # Shift guess if we just crossed the throat
            if is_supersonic and mn_guess < 1.0:
                mn_guess = 1.1
                
            m_l = solve_mach_newton(aratio, gamma, is_supersonic, mn_guess)
            mn_guess = m_l # March the guess to the next station for stability
        
        # Isentropic Relations for P and d
        tr = 1.0 + (gamma - 1.0) / 2.0 * m_l**2
        
        # Store Pressure in psia, Density in lbm/ft^3
        p[l, :, 0] = pc / (tr**(gamma / (gamma - 1.0)))
        d[l, :, 0] = dc / (tr**(1.0 / (gamma - 1.0)))
        
        # Velocity calculation: a = sqrt(gamma * R * g * T)
        t_l = tc / tr
        u_l = m_l * np.sqrt(gamma * rgas * g * t_l)
        
        u[l, :, 0] = u_l
        v[l, :, 0] = 0.0
        
    return u, v, p, d


def calculate_grid_mapping(geometry, params):
    # Step 1e - Mapping (Calculate transformation functions for grid)
    lmax = params['lmax']
    mmax = params['mmax']
    dy = 1.0 / (mmax - 1)
    
    # Initialize the mapping coefficient arrays
    geometry['be'] = np.zeros(lmax)
    geometry['al'] = np.zeros((lmax, mmax))
    geometry['de'] = np.zeros((lmax, mmax))
    
    for l in range(lmax):
        # Beta calculation: 1 / local height
        height = geometry['yw'][l] - geometry['ycb'][l]
        be = 1.0 / height
        geometry['be'][l] = be
        
        for m in range(mmax):
            y_comp = m * dy # Normalized radial coordinate (0 to 1)
            
            # Alpha: accounts for wall and centerbody slopes
            geometry['al'][l, m] = be * (geometry['nxnycb'][l] + y_comp * (geometry['nxny'][l] - geometry['nxnycb'][l]))
            
            # Delta: accounts for axial grid stretching
            geometry['de'][l, m] = -be * y_comp * geometry['xwi'][l]
            
    return geometry







def initalize_simulation():
    ''' Initalization of simulation. '''    
    # Step 1 - Initialization (MAIN, GEOM, ONEDIM)
    # Step 1a - Input Data (grid size, gas constants, boundary conditions)
    params = {'gamma': 1.4,             # Ratio of Specific Heats
              'rgas': 53.35,           # Gas Constant
              'g': 32.174,              # gravitation acceleration constant
              'lmax': 81,               # Axial Grid Points
              'mmax': 21,               # Radial Grid Points
              'tconv': 0.001,           # Convergence Tolerance
              'nmax': 1000,             # Maximum Time Steps
              'ngeom': 2,               # Case 2 corresponds to Circular-Arc/Conical Nozzle
              'xi': 0.31,               # Inlet X
              'xe': 4.05,               # Exit X
              'ri': 2.5,                # Inlet Radius
              'rt': 0.8,                # Throat Radius
              'rci': 0.8,               # Inlet Arc Radius
              'rct': 0.5,               # Throat Arc Radius
              'angi': 44.88,            # Inlet Angle (degrees)
              'ange': 15.0,             # Exit Angle (degrees)
              'ndim': 1,                # 1 for Axisymmetric, 0 for Planar
              'pc': 100.0,              # Stagnation Pressure (psia)
              'tc': 530.0,              # Stagnation Temperature (degR)
              }
                
    
    # Step 1b - Initialize solution arrays: [axial, radial, time_step_index]
    u = np.zeros((params['lmax'], params['mmax'], 2))
    v = np.zeros((params['lmax'], params['mmax'], 2))
    p = np.zeros((params['lmax'], params['mmax'], 2))
    d = np.zeros((params['lmax'], params['mmax'], 2))
    
    # Step 1c - Geometry Setup GEOM (Define nozzle wall) 
    geometry = setup_geometry(params)
    
    # Step 1d - Initial Values (Initial flow conditions using 1D isentropic equations)
    u, v, p, d = set_initial_conditions(u, v, p, d, params, geometry)
    
    # Step 1e - Mapping (Calculate transformation functions for grid)
    geometry = calculate_grid_mapping(geometry, params)
    
    return u, v, p, d, params, geometry
    


if __name__ == '__main__':
    
    u, v, p, d, params, geometry = initalize_simulation()
    
    # Step 2 - Time Integration Loop (MAIN)
    # Step 2a - Calculate Time Step (CFL condition)
    # Step 2b - Predictor Step (INTER - first MacCormack Step)
    # Step 2c - Boundary Conditions (INLET, WALL, EXITT)
    # Step 2d - Corrector Step (INTER - second MacCormack Step)
    # Step 2e - ConvergencyTCONV (check delu/u)
    # Step 2f - Calculate mass flowrate and thrust (MASFLO)
    
    # Step 3 - Plotting and Reporting
    plot_wall_grid(geometry)
    plot_isentropic_data(u, v, p, d, geometry, params)
    plot_computational_domain(geometry, params)
    
    
    
'''# Step 2 - Time Integration Loop (MAIN)
def run_simulation(u, v, p, ro, params, geometry):
    # The main time loop from Fortran: DO 650 N=1,NMAX.
    for n in range(params['nmax']):
        # Indices for current (n1) and next (n3) time steps
        n1 = n % 2
        n3 = (n + 1) % 2
        
        # Step 2a - Calculate Time Step (CFL condition)
        # Calculates DT based on the Courant-Friedrichs-Lewy stability criterion.
        dt = calculate_time_step(u, v, p, ro, n1, params)

        # Step 2b - Predictor Step (INTER - first MacCormack Step)
        # Uses backward finite differences to predict flow at N+1.
        u, v, p, ro = call_inter(u, v, p, ro, n1, n3, params, step_type='predictor')

        # Step 2c - Boundary Conditions (INLET, WALL, EXITT)
        # Updates the predicted values at the boundaries using characteristic methods[cite: 1, 2].
        u, v, p, ro = apply_boundary_conditions(u, v, p, ro, n3, params, geometry)

        # Step 2d - Corrector Step (INTER - second MacCormack Step)
        # Uses forward differences and averages with predictor values for 2nd order accuracy.
        u, v, p, ro = call_inter(u, v, p, ro, n1, n3, params, step_type='corrector')
        
        # Re-apply boundaries to the corrected solution.
        u, v, p, ro = apply_boundary_conditions(u, v, p, ro, n3, params, geometry)

        # Step 2e - ConvergencyTCONV (check delu/u)
        # Checks if the maximum change in axial velocity is below tolerance.
        if check_convergence(u, n1, n3, params):
            print(f"Converged at step {n}")
            break

# Step 3 - Plotting and Reporting
# Equivalent to CALL MASFLO and CALL PLOT in the Fortran code.
def report_results(u, v, p, ro, n3, params):
    # Calculate final mass flow and thrust.
    # Generate output tables and visualization plots.
    pass
'''
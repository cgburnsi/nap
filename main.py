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
    
    
def plot_simulation_results(u, v, p, d, geometry, params):
    """
    Generates a 2x2 grid plot of the 2D results comparing Centerline vs Wall.
    Similar to the isentropic plot but showing 2D data variation.
    """
    gamma, g, x = params['gamma'], params['g'], geometry['xw']
    # Final results are usually in the most recently updated n-index (n3)
    # If converged, n1 and n3 are very similar.
    n_idx = (params['iteration_count'] + 1) % 2
    
    # 1. Properties at Centerline (m=0)
    p_cl, d_cl, u_cl, v_cl = p[:, 0, n_idx], d[:, 0, n_idx], u[:, 0, n_idx], v[:, 0, n_idx]
    m_cl = np.sqrt(u_cl**2 + v_cl**2) / np.sqrt(gamma * p_cl * 144.0 * g / d_cl)
    
    # 2. Properties at Wall (m=MMAX-1)
    p_w, d_w, u_w, v_w = p[:, -1, n_idx], d[:, -1, n_idx], u[:, -1, n_idx], v[:, -1, n_idx]
    m_w = np.sqrt(u_w**2 + v_w**2) / np.sqrt(gamma * p_w * 144.0 * g / d_w)
    
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    
    # Mach Number
    axs[0, 0].plot(x, m_cl, 'b-', label='Centerline')
    axs[0, 0].plot(x, m_w, 'b--', label='Wall')
    axs[0, 0].set_title('Mach Number')
    axs[0, 0].set_ylabel('M')
    
    # Pressure (psia)
    axs[0, 1].plot(x, p_cl, 'r-', label='Centerline')
    axs[0, 1].plot(x, p_w, 'r--', label='Wall')
    axs[0, 1].set_title('Pressure (psia)')
    axs[0, 1].set_ylabel('P (psia)')
    
    # Density (slug/ft^3)
    axs[1, 0].plot(x, d_cl, 'g-', label='Centerline')
    axs[1, 0].plot(x, d_w, 'g--', label='Wall')
    axs[1, 0].set_title('Density (slug/ft^3)')
    axs[1, 0].set_ylabel('d')
    
    # Axial Velocity (ft/s)
    axs[1, 1].plot(x, u_cl, 'k-', label='Centerline')
    axs[1, 1].plot(x, u_w, 'k--', label='Wall')
    axs[1, 1].set_title('Axial Velocity (ft/s)')
    axs[1, 1].set_ylabel('U (ft/s)')
    
    for ax in axs.flat:
        ax.set_xlabel('X (Axial Distance)')
        ax.legend()
        ax.grid(True)
    
    plt.tight_layout()
    plt.savefig('simulation_output_2d.png')
    plt.show()
    
def plot_pressure_contours(u, v, p, d, geometry, params):
    """
    Generates a 2D contour plot of the pressure field for the top half of the nozzle.
    """
    lmax, mmax = params['lmax'], params['mmax']
    dy = 1.0 / (mmax - 1)
    n_idx = (params.get('iteration_count', 0) + 1) % 2
    
    X = np.zeros((lmax, mmax))
    Y = np.zeros((lmax, mmax))
    P_data = p[:, :, n_idx]

    for l in range(lmax):
        be = geometry['be'][l]
        y_cb = geometry['ycb'][l]
        for m in range(mmax):
            X[l, m] = geometry['xw'][l]
            Y[l, m] = (m * dy / be) + y_cb

    plt.figure(figsize=(14, 5))
    
    # Plot only the top half
    cp = plt.contourf(X, Y, P_data, levels=50, cmap='jet')
    plt.colorbar(cp, label='Pressure (psia)')
    
    # Plot the boundaries
    plt.plot(geometry['xw'], geometry['yw'], 'k-', linewidth=2, label='Nozzle Wall')
    plt.plot(geometry['xw'], geometry['ycb'], 'k--', alpha=0.8, label='Centerline')
    
    plt.title('2D Pressure Distribution (Centerline to Wall)')
    plt.xlabel('X (Axial Distance)')
    plt.ylabel('Y (Radius)')
    plt.axis('equal')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.3)
    
    plt.savefig('pressure_contours_top.png')
    plt.show()
    
    
def plot_mach_contours(u, v, p, d, geometry, params):
    """
    Generates a 2D contour plot of the Mach Number field for the top half.
    """
    lmax, mmax = params['lmax'], params['mmax']
    gamma, g = params['gamma'], params['g']
    dy = 1.0 / (mmax - 1)
    n_idx = (params.get('iteration_count', 0) + 1) % 2
    
    X = np.zeros((lmax, mmax))
    Y = np.zeros((lmax, mmax))
    Mach = np.zeros((lmax, mmax))

    for l in range(lmax):
        be = geometry['be'][l]
        y_cb = geometry['ycb'][l]
        for m in range(mmax):
            X[l, m] = geometry['xw'][l]
            Y[l, m] = (m * dy / be) + y_cb
            
            ul, vl, pl, dl = u[l, m, n_idx], v[l, m, n_idx], p[l, m, n_idx], d[l, m, n_idx]
            a_local = np.sqrt(gamma * pl * 144.0 * g / dl)
            Mach[l, m] = np.sqrt(ul**2 + vl**2) / a_local

    plt.figure(figsize=(14, 5))
    
    # Plot only the top half
    cp = plt.contourf(X, Y, Mach, levels=50, cmap='viridis')
    plt.colorbar(cp, label='Mach Number (M)')
    
    # Explicitly highlight the Sonic Line (M=1.0)
    plt.contour(X, Y, Mach, levels=[1.0], colors='white', linestyles='dashed', linewidths=1.5)

    plt.plot(geometry['xw'], geometry['yw'], 'k-', linewidth=2, label='Nozzle Wall')
    plt.plot(geometry['xw'], geometry['ycb'], 'k--', alpha=0.8, label='Centerline')
    
    plt.title('2D Mach Number Field (Centerline to Wall)')
    plt.xlabel('X (Axial Distance)')
    plt.ylabel('Y (Radius)')
    plt.axis('equal')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.3)
    
    plt.savefig('mach_contours_top.png')
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
        'xwi': np.zeros(lmax),
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
              'pe': 5.5,                # Static exit pressure for pressure BC mode (psia)
              'exit_bc_mode': 'pressure', # 'pressure' or 'extrapolate'
              'fdt': 0.15,              # Fractional Delta Time (time-step factor)
              # --- New Shock/AV Parameters from shock.f ---
              'iav': 0,                 # Artificial Viscosity Switch (1=On)
              'cav': 1.0,               # Viscosity Coefficient 
              'xla': 1.0,               # Lambda coefficient 
              'xmu': 1.0,               # Mu coefficient 
              'rkmu': 1.0,              # Energy/RK coefficient 
              'cta': 0.9,               # Relaxation Factor 
              'lss': 1,                 # Starting L-index for AV calculation [cite: 44]
              }
    
    
    # Step 1b - Initialize solution arrays: [axial, radial, time_step_index]
    u = np.zeros((params['lmax'], params['mmax'], 2))
    v = np.zeros((params['lmax'], params['mmax'], 2))
    p = np.zeros((params['lmax'], params['mmax'], 2))
    d = np.zeros((params['lmax'], params['mmax'], 2))
    
    # Initialize Q (Viscosity) arrays in params for persistence 
    params['q_old'] = {
        'u': np.zeros((params['lmax'], params['mmax'])),
        'v': np.zeros((params['lmax'], params['mmax'])),
        'p': np.zeros((params['lmax'], params['mmax']))
    }
    
    
    # Step 1c - Geometry Setup GEOM (Define nozzle wall) 
    geometry = setup_geometry(params)
    
    # Step 1d - Initial Values (Initial flow conditions using 1D isentropic equations)
    u, v, p, d = set_initial_conditions(u, v, p, d, params, geometry)
    
    # Step 1e - Mapping (Calculate transformation functions for grid)
    geometry = calculate_grid_mapping(geometry, params)
    
    return u, v, p, d, params, geometry
  





def calculate_shock_viscosity(u, v, p, d, n1, params, geometry, q_old):
    """
    Step 2g - Local Artificial Viscosity (SHOCK IPASS=1).
    Fixed NameError by rigidly defining all intermediate gradient terms from shock.f.
    """
    lmax, mmax = params['lmax'], params['mmax']
    gamma, rgas, g = params['gamma'], params['rgas'], params['g']
    dxr, dyr = 1.0 / params['dx'], 1.0 / (1.0 / (mmax - 1))
    rg = rgas * g
    
    cav, xla, xmu, rkmu = params['cav'], params['xla'], params['xmu'], params['rkmu']
    cta, lss = params['cta'], params.get('lss', 1)
    gam_m1 = gamma - 1.0
    rk_const = params.get('gam1', gamma / gam_m1) 
    
    dt, dx, dy = params['dt'], params['dx'], 1.0/(mmax-1)
    rdum = cav * dt * dx * dy * 2.0 

    q_u, q_v, q_p = np.zeros((lmax, mmax)), np.zeros((lmax, mmax)), np.zeros((lmax, mmax))

    for l in range(lss, lmax - 1):
        be = geometry['be'][l]
        be_m = geometry['be'][l-1]
        be_p = geometry['be'][l+1]
        
        for m in range(mmax):
            al = geometry['al'][l, m]
            
            # 1. Divergence Check (DIV)
            ux = 0.5 * (u[l+1, m, n1] - u[l-1, m, n1]) * dxr
            if m == 0:
                uy, vy = (u[l, 1, n1] - u[l, 0, n1]) * dyr, (v[l, 1, n1] - v[l, 0, n1]) * dyr
            elif m == mmax - 1:
                uy, vy = (u[l, m, n1] - u[l, m-1, n1]) * dyr, (v[l, m, n1] - v[l, m-1, n1]) * dyr
            else:
                uy, vy = 0.5 * (u[l, m+1, n1] - u[l, m-1, n1]) * dyr, 0.5 * (v[l, m+1, n1] - v[l, m-1, n1]) * dyr
            
            div_check = ux + al*uy + be*vy
            if div_check >= 0.0:
                continue

            # 2. Axial Gradients and Mapping
            t = p[l, m, n1] / (d[l, m, n1] * rg)
            ux1, ux2 = (u[l, m, n1] - u[l-1, m, n1]) * dxr, (u[l+1, m, n1] - u[l, m, n1]) * dxr
            vx1, vx2 = (v[l, m, n1] - v[l-1, m, n1]) * dxr, (v[l+1, m, n1] - v[l, m, n1]) * dxr
            tx1, tx2 = (t - p[l-1, m, n1]/(d[l-1, m, n1]*rg)) * dxr, (p[l+1, m, n1]/(d[l+1, m, n1]*rg) - t) * dxr
            
            be1, be2 = 0.5*(be_m + be), 0.5*(be_p + be)
            al1, al2 = 0.5*(geometry['al'][l-1, m] + al), 0.5*(geometry['al'][l+1, m] + al)

            # 3. Radial Gradients (M-Branching)
            ty1 = ty2 = tx3 = tx4 = 0.0
            if m == 0: # Centerline Point
                uy1 = 0.5*(u[l,1,n1] + u[l-1,1,n1] - u[l,0,n1] - u[l-1,0,n1]) * dyr
                vy1 = 0.5*(v[l,1,n1] + v[l-1,1,n1] - v[l,0,n1] - v[l-1,0,n1]) * dyr
                uy2 = 0.5*(u[l+1,1,n1] + u[l,1,n1] - u[l+1,0,n1] - u[l,0,n1]) * dyr
                vy2 = 0.5*(v[l+1,1,n1] + v[l,1,n1] - v[l+1,0,n1] - v[l,0,n1]) * dyr
                ux3 = ux4 = vx3 = vx4 = 0.0
                theta_w = np.arctan(-geometry['nxnycb'][l])
                theta_int = np.arctan(v[l, 1, n1] / u[l, 1, n1])
                vm_int = np.sqrt(u[l, 1, n1]**2 + v[l, 1, n1]**2)
                r_theta = 2.0 * theta_w - theta_int
                ur, vr = vm_int * np.cos(r_theta), vm_int * np.sin(r_theta)
                uy3, vy3 = (u[l, 0, n1] - ur) * dyr, (v[l, 0, n1] - vr) * dyr
                uy4, vy4 = (u[l, 1, n1] - u[l, 0, n1]) * dyr, (v[l, 1, n1] - v[l, 0, n1]) * dyr
                ty4 = (p[l, 1, n1]/(d[l, 1, n1]*rg) - t) * dyr
                ty3 = -ty4
                be3, be4, al3, al4 = be, be, al, geometry['al'][l, 1]
            elif m == mmax - 1: # Wall Point
                uy1 = 0.5*(u[l,m,n1] + u[l-1,m,n1] - u[l,m-1,n1] - u[l-1,m-1,n1]) * dyr
                vy1 = 0.5*(v[l,m,n1] + v[l-1,m,n1] - v[l,m-1,n1] - v[l-1,m-1,n1]) * dyr
                uy2 = 0.5*(u[l+1,m,n1] + u[l,m,n1] - u[l+1,m-1,n1] - u[l,m-1,n1]) * dyr
                vy2 = 0.5*(v[l+1,m,n1] + v[l,m,n1] - v[l+1,m-1,n1] - v[l,m-1,n1]) * dyr
                ux3 = ux4 = vx3 = vx4 = 0.0
                uy3, vy3 = (u[l,m,n1]-u[l,m-1,n1])*dyr, (v[l,m,n1]-v[l,m-1,n1])*dyr
                theta_w = np.arctan(-geometry['nxny'][l])
                theta_int = np.arctan(v[l, m-1, n1] / u[l, m-1, n1])
                vm_int = np.sqrt(u[l, m-1, n1]**2 + v[l, m-1, n1]**2)
                r_theta = 2.0 * theta_w - theta_int
                ur, vr = vm_int * np.cos(r_theta), vm_int * np.sin(r_theta)
                uy4, vy4 = (ur - u[l, m, n1]) * dyr, (vr - v[l, m, n1]) * dyr
                ty3 = (t - p[l, m-1, n1]/(d[l, m-1, n1]*rg)) * dyr
                ty4 = -ty3
                be3, be4, al3, al4 = be, be, geometry['al'][l, m-1], al
            else: # Interior Point
                uy1 = 0.25*(u[l,m+1,n1] + u[l-1,m+1,n1] - u[l,m-1,n1] - u[l-1,m-1,n1]) * dyr
                uy2 = 0.25*(u[l,m+1,n1] + u[l+1,m+1,n1] - u[l,m-1,n1] - u[l+1,m-1,n1]) * dyr
                vy1 = 0.25*(v[l,m+1,n1] + v[l-1,m+1,n1] - v[l,m-1,n1] - v[l-1,m-1,n1]) * dyr
                vy2 = 0.25*(v[l,m+1,n1] + v[l+1,m+1,n1] - v[l,m-1,n1] - v[l+1,m-1,n1]) * dyr
                ux3 = 0.25*(u[l+1,m,n1] + u[l+1,m-1,n1] - u[l-1,m,n1] - u[l-1,m-1,n1]) * dxr
                ux4 = 0.25*(u[l+1,m,n1] + u[l+1,m+1,n1] - u[l-1,m,n1] - u[l-1,m+1,n1]) * dxr
                vx3 = 0.25*(v[l+1,m,n1] + v[l+1,m-1,n1] - v[l-1,m,n1] - v[l-1,m-1,n1]) * dxr
                vx4 = 0.25*(v[l+1,m,n1] + v[l+1,m+1,n1] - v[l-1,m,n1] - v[l-1,m+1,n1]) * dxr
                uy3, uy4 = (u[l,m,n1]-u[l,m-1,n1])*dyr, (u[l,m+1,n1]-u[l,m,n1])*dyr
                vy3, vy4 = (v[l,m,n1]-v[l,m-1,n1])*dyr, (v[l,m+1,n1]-v[l,m,n1])*dyr
                ty3, ty4 = (t - p[l, m-1, n1]/(d[l, m-1, n1]*rg)) * dyr, (p[l, m+1, n1]/(d[l, m+1, n1]*rg) - t) * dyr
                be3, be4 = be, be
                al3, al4 = geometry['al'][l, m-1], geometry['al'][l, m+1]

            # 4. Composite Gradients
            uxy1, uxy2 = ux1 + al1*uy1, ux2 + al2*uy2
            uxy3, uxy4 = ux3 + al3*uy3, ux4 + al4*uy4
            vxy1, vxy2 = vx1 + al1*vy1, vx2 + al2*vy2
            vxy3, vxy4 = vx3 + al3*vy3, vx4 + al4*vy4
            
            # FIXED: Explicitly define all BUY and BVY terms
            buy1, buy2, buy3, buy4 = be1*uy1, be2*uy2, be3*uy3, be4*uy4
            bvy1, bvy2, bvy3, bvy4 = be1*vy1, be2*vy2, be3*vy3, be4*vy4
            
            uxy12 = 0.5 * (ux1 + ux2 + al3*uy3 + al4*uy4)
            vxy12 = 0.5 * (vx1 + vx2 + al3*vy3 + al4*vy4)
            buy34, bvy34 = 0.5*(buy3 + buy4), 0.5*(bvy3 + bvy4)
            txy1, txy2 = tx1 + al1*ty1, tx2 + al2*ty2
            txy3, txy4 = tx3 + al3*ty3, tx4 + al4*ty4
            
            # 5. Viscosity Coefficients
            div = uxy12 + bvy34
            vid = vxy12 - buy34
            rla_term, rmu_term = xla * rdum * abs(div)/be, xmu * rdum * abs(vid)/be
            rk = rmu_term * rk_const * rg / rkmu
            rlp2m = rla_term + 2.0 * rmu_term

            # 6. Axisymmetric Terms
            uvta = vvta = pvta = pcta = 0.0
            if params['ndim'] == 1:
                y_phys = m * dy / be + geometry['ycb'][l]
                if not (m == 0 and geometry['ycb'][l] == 0.0):
                    rlpm, rla2 = rla_term + rmu_term, 2.0 * rla_term
                    uvta = (rlpm * vxy12 + rmu_term * buy34) / y_phys
                    vvta = rlp2m * (bvy34 + v[l, m, n1] / y_phys) / y_phys
                    pvta = (rlp2m * v[l, m, n1]**2 / y_phys + rla2 * v[l, m, n1] * (bvy34 + uxy12)) / y_phys
                    pcta = rk * 0.5 * (be4*ty4 + be3*ty3) / y_phys
                    duvta = rlpm * be * (vxy4 - vxy3) * dyr + rmu_term * be * (buy4 - buy3) * dyr
                    dvvta = rlp2m * 0.5 * be * (bvy4 - bvy3) * dyr
                    dpvta = (rlp2m + rla2) * bvy34**2 + rla2 * bvy34 * uxy12
                    dpcta = rk * be * (be4*ty4 - be3*ty3) * dyr
                    if abs(uvta) > abs(duvta): uvta = duvta
                    if abs(vvta) > abs(dvvta): vvta = dvvta
                    if abs(pvta) > abs(dpvta): pvta = dpvta
                    if abs(pcta) > abs(dpcta): pcta = dpcta
                else: # Axis limit label 120
                    uvta = (rla_term+rmu_term)*be*(vxy4-vxy3)*dyr + rmu_term*be*(buy4-buy3)*dyr
                    vvta = rlp2m * 0.5 * be * (bvy4 - bvy3) * dyr
                    pvta = (rlp2m + 2.0*rla_term) * bvy34**2 + 2.0*rla_term * bvy34 * uxy12
                    pcta = rk * be * (be4*ty4 - be3*ty3) * dyr

            # 7. Final Accumulation
            q_u[l, m] = rlp2m*(uxy2-uxy1)*dxr + al*(rlp2m*(uxy4-uxy3) + rla_term*(bvy4-bvy3))*dyr + rmu_term*be*(vxy4-buy4-vxy3+buy3)*dyr + uvta
            q_v[l, m] = rmu_term*(vxy2-buy2-vxy1+buy1)*dxr + rmu_term*al*(vxy4-buy4-vxy3+buy3)*dyr + be*(rla_term*(uxy4-uxy3) + rlp2m*(bvy4-bvy3))*dyr + vvta
            
            work_diss = rlp2m*(uxy12**2 + bvy34**2) + rmu_term*(vxy12**2 + buy34**2) + 2.0*rla_term*uxy12*bvy34 + 2.0*rmu_term*buy34*vxy12
            heat_diss = rk * ((txy2 - txy1) * dxr + al * (txy4 - txy3) * dyr + be * (be4*ty4 - be3*ty3) * dyr)
            q_p[l, m] = d[l, m, n1] * gam_m1 * (work_diss + heat_diss + pvta + pcta)

            # 8. Relaxation
            q_u[l, m] = cta * q_u[l, m] + (1.0 - cta) * q_old['u'][l, m]
            q_v[l, m] = cta * q_v[l, m] + (1.0 - cta) * q_old['v'][l, m]
            q_p[l, m] = cta * q_p[l, m] + (1.0 - cta) * q_old['p'][l, m]

    return q_u, q_v, q_p




def smooth_flow_field(u, v, p, d, n3, params, geometry):
    """
    Step 2h - Global Flow Smoothing (SHOCK IPASS=2).
    Applies a weighted average to the solution to prevent numerical ringing.
    """
    smp = params.get('smp', 0.1) # Smoothing factor
    
    # If SMP is out of range, the routine returns without action
    if smp < 0.0 or smp >= 1.0:
        return u, v, p, d

    l1, m1 = params['lmax'] - 2, params['mmax'] - 2
    lmax, mmax = params['lmax'], params['mmax']
    smp4 = 0.25 * (1.0 - smp) # Smoothing coefficient

    # 1. Smooth the Wall (M=MMAX)
    m_wall = mmax - 1
    for l in range(1, lmax - 1):
        u[l, m_wall, n3] = smp4 * (u[l-1, m_wall, n3] + u[l+1, m_wall, n3] + 2.0 * u[l, m_wall-1, n3]) + smp * u[l, m_wall, n3]
        v[l, m_wall, n3] = -u[l, m_wall, n3] * geometry['nxny'][l] # Maintain tangency
        p[l, m_wall, n3] = smp4 * (p[l-1, m_wall, n3] + p[l+1, m_wall, n3] + 2.0 * p[l, m_wall-1, n3]) + smp * p[l, m_wall, n3]
        d[l, m_wall, n3] = smp4 * (d[l-1, m_wall, n3] + d[l+1, m_wall, n3] + 2.0 * d[l, m_wall-1, n3]) + smp * d[l, m_wall, n3]

    # 2. Smooth the Axis (M=0)
    for l in range(1, lmax - 1):
        u[l, 0, n3] = smp4 * (u[l-1, 0, n3] + u[l+1, 0, n3] + 2.0 * u[l, 1, n3]) + smp * u[l, 0, n3]
        v[l, 0, n3] = -u[l, 0, n3] * geometry['nxnycb'][l] # Centerbody tangency
        p[l, 0, n3] = smp4 * (p[l-1, 0, n3] + p[l+1, 0, n3] + 2.0 * p[l, 1, n3]) + smp * p[l, 0, n3]
        d[l, 0, n3] = smp4 * (d[l-1, 0, n3] + d[l+1, 0, n3] + 2.0 * d[l, 1, n3]) + smp * d[l, 0, n3]

    # 3. Smooth the Interior Points
    for l in range(1, lmax - 1):
        for m in range(1, mmax - 1):
            # Weighted average of the 4 cardinal neighbors
            u[l, m, n3] = smp4 * (u[l-1, m, n3] + u[l+1, m, n3] + u[l, m-1, n3] + u[l, m+1, n3]) + smp * u[l, m, n3]
            v[l, m, n3] = smp4 * (v[l-1, m, n3] + v[l+1, m, n3] + v[l, m-1, n3] + v[l, m+1, n3]) + smp * v[l, m, n3]
            p[l, m, n3] = smp4 * (p[l-1, m, n3] + p[l+1, m, n3] + p[l, m-1, n3] + p[l, m+1, n3]) + smp * p[l, m, n3]
            d[l, m, n3] = smp4 * (d[l-1, m, n3] + d[l+1, m, n3] + d[l, m-1, n3] + d[l, m+1, n3]) + smp * d[l, m, n3]

    return u, v, p, d






def calculate_time_step(u, v, p, d, n, params):
    """
    Step 2a - Calculate Time Step using NAP standard FDT.
    """
    gamma = params['gamma']
    rgas = params['rgas']
    g = params['g']
    dx = params['dx']
    fdt = params.get('fdt', 1.6) # Updated per NAP p.36

    # Calculate local sound speed
    # a = sqrt(gamma * P * 144 * g / d)
    a_sound = np.sqrt(gamma * p[:, :, n] * 144.0 * g / d[:, :, n])

    # Wave speeds in the axial direction
    wave_speed = np.abs(u[:, :, n]) + a_sound
    
    # Global minimum delta t
    dt_min = np.min(dx / wave_speed)
    
    return fdt * dt_min
  

def call_inter(u, v, p, d, n1, n3, dt, params, geometry, step_type='predictor'):
    """
    Step 2b/2d - Interior Solver (MacCormack Step).
    Strictly follows axis/interior branching from inter.f to prevent 1/0 errors.
    """
    gamma, ndim = params['gamma'], params['ndim']
    dxr = 1.0 / params['dx']
    dyr = 1.0 / (1.0 / (params['mmax'] - 1))
    lmax, mmax = params['lmax'], params['mmax']
    
    ichar = 1 if step_type == 'predictor' else 2
    l_start = 1
    # Predictor loops to LMAX, Corrector to LMAX-1
    l_end = (lmax - 1) if ichar == 1 else (lmax - 2)
    m_end = mmax - 2 # Interior points exclude the wall boundary

    for l in range(l_start, l_end + 1):
        for m in range(0, m_end + 1):
            n_idx = n1 if ichar == 1 else n3
            l_ref = l-1 if ichar == 1 else l+1
            m_ref = m-1 if ichar == 1 else m+1

            ub, vb, pb, rob = u[l,m,n_idx], v[l,m,n_idx], p[l,m,n_idx], d[l,m,n_idx]
            asb, be = gamma * pb / rob, geometry['be'][l]

            if m == 0: # Symmetry Axis (M=1 branch in inter.f)
                dudx = (ub - u[l_ref, m, n_idx]) * dxr if ichar == 1 else (u[l_ref, m, n_idx] - ub) * dxr
                dpdx = (pb - p[l_ref, m, n_idx]) * dxr if ichar == 1 else (p[l_ref, m, n_idx] - pb) * dxr
                drodx = (rob - d[l_ref, m, n_idx]) * dxr if ichar == 1 else (d[l_ref, m, n_idx] - rob) * dxr
                # Special radial derivative for axis
                dvdy = (4.0 * v[l, 1, n_idx] - v[l, 2, n_idx]) * 0.5 * dyr
                
                urhs = -ub * dudx - dpdx / rob
                # Axisymmetric RORHS includes (1+NDIM) factor
                rorhs = -ub * drodx - rob * dudx - (1 + ndim) * rob * be * dvdy
                prhs = -ub * dpdx + asb * (rorhs + ub * drodx)
                v[l, m, n3] = 0.0 # Force zero radial velocity on axis
            else: # Interior Points (M > 1 branch in inter.f)
                al, de = geometry['al'][l, m], geometry['de'][l, m]
                uvb = ub * al + vb * be + de
                
                dudx = (ub - u[l_ref, m, n_idx]) * dxr if ichar == 1 else (u[l_ref, m, n_idx] - ub) * dxr
                dvdx = (vb - v[l_ref, m, n_idx]) * dxr if ichar == 1 else (v[l_ref, m, n_idx] - vb) * dxr
                dpdx = (pb - p[l_ref, m, n_idx]) * dxr if ichar == 1 else (p[l_ref, m, n_idx] - pb) * dxr
                drodx = (rob - d[l_ref, m, n_idx]) * dxr if ichar == 1 else (d[l_ref, m, n_idx] - rob) * dxr
                
                dudy = (ub - u[l, m_ref, n_idx]) * dyr if ichar == 1 else (u[l, m_ref, n_idx] - ub) * dyr
                dvdy = (vb - v[l, m_ref, n_idx]) * dyr if ichar == 1 else (v[l, m_ref, n_idx] - vb) * dyr
                dpdy = (pb - p[l, m_ref, n_idx]) * dyr if ichar == 1 else (p[l, m_ref, n_idx] - pb) * dyr
                drody = (rob - d[l, m_ref, n_idx]) * dyr if ichar == 1 else (d[l, m_ref, n_idx] - rob) * dyr
                
                # ATERM bypasses at m=0 to prevent division by zero
                radius = m * (1.0 / (mmax - 1)) / be + geometry['ycb'][l]
                aterm = (rob * vb / radius) if ndim == 1 else 0.0
                
                urhs = -ub * dudx - uvb * dudy - (dpdx + al * dpdy) / rob
                vrhs = -ub * dvdx - uvb * dvdy - be * dpdy / rob
                rorhs = -ub * drodx - uvb * drody - rob * (dudx + al * dudy + be * dvdy) - aterm
                prhs = -ub * dpdx - uvb * dpdy + asb * (rorhs + ub * drodx + uvb * drody)
                
                if ichar == 1: v[l, m, n3] = v[l, m, n1] + vrhs * dt
                else: v[l, m, n3] = (v[l, m, n1] + v[l, m, n3] + vrhs * dt) * 0.5
            
            # Step forward for U, P, and Density
            if ichar == 1:
                u[l, m, n3] = u[l, m, n1] + urhs * dt
                p[l, m, n3] = p[l, m, n1] + prhs * dt
                d[l, m, n3] = d[l, m, n1] + rorhs * dt
            else:
                u[l, m, n3] = (u[l, m, n1] + u[l, m, n3] + urhs * dt) * 0.5
                p[l, m, n3] = (p[l, m, n1] + p[l, m, n3] + prhs * dt) * 0.5
                d[l, m, n3] = (d[l, m, n1] + d[l, m, n3] + rorhs * dt) * 0.5
            
            # Clamp to prevent non-physical negatives during transients
            if p[l, m, n3] <= 0: p[l, m, n3] = 0.001 * params['pc']
            if d[l, m, n3] <= 0: d[l, m, n3] = 0.001 / params['g']

    return u, v, p, d


def apply_inlet_bc(u, v, p, d, n1, n3, dt, params, geometry, step_type):
    """
    Step 2c - Inlet Boundary Condition.
    Solves for the inlet state using characteristic variables and total properties.
    Includes stabilization for the iterative solver.
    """
    gamma = params['gamma']
    rgas = params['rgas']
    g = params['g']
    mmax = params['mmax']
    dxr = 1.0 / params['dx']
    
    gam1 = gamma / (gamma - 1.0)
    gam2 = (gamma - 1.0) / 2.0
    grgb = gamma * rgas * g

    n_state = n1 if step_type == 'predictor' else n3

    for m in range(mmax):
        # 1. Characteristic Origin (X2) Interpolation
        # Current state at the boundary
        u0 = u[0, m, n_state]
        p0 = p[0, m, n_state]
        d0 = d[0, m, n_state]
        a0 = np.sqrt(gamma * p0 * 144.0 * g / d0)
        
        # Characteristic distance back into the interior
        dist = (u0 - a0) * dt
        weight = -dist * dxr # Weighting factor for interpolation

        # Interpolate properties at the characteristic origin (X2)
        u2 = u0 + (u[1, m, n_state] - u0) * weight
        p2 = p0 + (p[1, m, n_state] - p0) * weight
        d2 = d0 + (d[1, m, n_state] - d0) * weight
        a2 = np.sqrt(gamma * p2 * 144.0 * g / d2)

        # 2. Iterative Solve for Boundary State (MN3)
        # Initial guess for Mach number from previous time step
        mn3 = u[0, m, n3] / a0 if u[0, m, n3] > 1e-9 else 0.05
        
        for _ in range(20): # Iterate to converge on the new Mach number
            # Limit Mach number to prevent overflow during transients
            mn3 = max(1e-5, min(mn3, 0.95))
            
            dem = 1.0 + gam2 * mn3**2
            
            # Calculate new total pressure and temperature based on guess
            p_new = params['pc'] / (dem**gam1)
            t_new = params['tc'] / dem
            
            # Prevent division by zero if T_new gets too small
            if t_new < 1e-9: t_new = 1e-9

            # Compatibility Equation: Connects boundary to interior state
            # This updates the velocity based on the pressure change
            u[0, m, n3] = u2 + (p_new - p2) / (d2 * a2)
            
            # Calculate the new Mach number from the updated velocity and temp
            a_new = np.sqrt(grgb * t_new)
            mn3_new = u[0, m, n3] / a_new

            # Check for convergence
            if abs(mn3_new - mn3) < 1e-4:
                break
            
            # Update guess for the next iteration (with simple under-relaxation)
            mn3 = 0.5 * (mn3 + mn3_new)

        # 3. Final Update of Boundary Properties
        v[0, m, n3] = 0.0 # Inlet is purely axial
        p[0, m, n3] = p_new
        d[0, m, n3] = p_new * 144.0 / (rgas * t_new * g)

    return u, v, p, d




def apply_wall_bc(u, v, p, d, n1, n3, dt, params, geometry, step_type='predictor'):
    """
    Step 2c - Wall Boundary Conditions (MoC).
    Strict translation of wall.f characteristic coupling logic.
    """
    gamma, g = params['gamma'], params['g']
    mmax, lmax = params['mmax'], params['lmax']
    dxr = 1.0 / params['dx']
    dyr = 1.0 / (1.0 / (mmax - 1))
    ichar = 1 if step_type == 'predictor' else 2
    ndim = params.get('ndim', 1)

    m = mmax - 1
    m_int = m - 1
    sign = 1.0
    dys = sign * dyr
    y_b = 1.0
    l_end = (lmax - 1) if ichar == 1 else (lmax - 2)

    for l in range(1, l_end + 1):
        u1, v1, p1, ro1 = u[l, m, n1], v[l, m, n1], p[l, m, n1], d[l, m, n1]
        a1 = np.sqrt(max(gamma * p1 * 144.0 * g / ro1, 1e-12))

        if ichar == 2:
            u3, v3, p3, ro3 = u[l, m, n3], v[l, m, n3], p[l, m, n3], d[l, m, n3]
            a3 = np.sqrt(max(gamma * p3 * 144.0 * g / ro3, 1e-12))
        else:
            u3, v3, p3, ro3, a3 = u1, v1, p1, ro1, a1

        bu = (u1 - u[l, m_int, n1]) * dys
        bv = (v1 - v[l, m_int, n1]) * dys
        bp = (p1 - p[l, m_int, n1]) * dys
        bro = (ro1 - d[l, m_int, n1]) * dys
        cu, cv, cp, cro = u1 - bu * y_b, v1 - bv * y_b, p1 - bp * y_b, ro1 - bro * y_b

        du = (u1 - u[l-1, m, n1]) * dxr
        dv = (v1 - v[l-1, m, n1]) * dxr
        dp = (p1 - p[l-1, m, n1]) * dxr
        dro = (ro1 - d[l-1, m, n1]) * dxr
        du_m = (u[l, m_int, n1] - u[l-1, m_int, n1]) * dxr
        dv_m = (v[l, m_int, n1] - v[l-1, m_int, n1]) * dxr
        dp_m = (p[l, m_int, n1] - p[l-1, m_int, n1]) * dxr
        dro_m = (d[l, m_int, n1] - d[l-1, m_int, n1]) * dxr
        bdu = (du - du_m) * dys
        bdv = (dv - dv_m) * dys
        bdp = (dp - dp_m) * dys
        bdro = (dro - dro_m) * dys
        cdu = du - bdu * y_b
        cdv = dv - bdv * y_b
        cdp = dp - bdp * y_b
        cdro = dro - bdro * y_b

        al, be, de = geometry['al'][l, m], geometry['be'][l], geometry['de'][l, m]
        als = np.sqrt(al * al + be * be)
        uv3 = u3 * al + v3 * be + de
        al2 = al

        y2 = y_b
        u2, v2, p2, ro2, a2 = u1, v1, p1, ro1, a1
        for _ in range(3):
            uv2 = u2 * al2 + v2 * be + de
            y2 = y_b - (uv2 + sign * al * als * a2 + uv3 + sign * als * a3) * dt * 0.5
            u2 = bu * y2 + cu
            v2 = bv * y2 + cv
            p2 = bp * y2 + cp
            ro2 = bro * y2 + cro
            al2 = y2 * al
            ad = gamma * p2 * 144.0 * g / max(ro2, 1e-12)
            a2 = np.sqrt(max(ad, 1e-12))

        du1 = du
        dv1 = dv
        dp1 = dp
        dro1 = dro
        du2 = bdu * y2 + cdu
        dv2 = bdv * y2 + cdv
        dp2 = bdp * y2 + cdp
        dro2 = bdro * y2 + cdro

        aterm2 = 0.0
        if ndim == 1:
            denom2 = geometry['ycb'][l] + y2 / be
            if abs(denom2) > 1e-12:
                aterm2 = ro2 * v2 / denom2

        psi21 = -u1 * du1 - dp1 / ro1
        psi31 = -u1 * dv1
        psi41 = -u1 * dp1 + a1 * a1 * u1 * dro1
        psi12 = -u2 * dro2 - ro2 * du2 - aterm2
        psi22 = -u2 * du2 - dp2 / ro2
        psi32 = -u2 * dv2
        psi42 = -u2 * dp2 + a2 * a2 * u2 * dro2

        if ichar == 2:
            du3 = (u[l+1, m, n3] - u3) * dxr
            dv3 = (v[l+1, m, n3] - v3) * dxr
            dp3 = (p[l+1, m, n3] - p3) * dxr
            dro3 = (d[l+1, m, n3] - ro3) * dxr
            aterm3 = 0.0
            if ndim == 1:
                denom3 = geometry['ycb'][l] + 1.0 / be
                if abs(denom3) > 1e-12:
                    aterm3 = ro3 * v2 / denom3
            psi13 = -u3 * dro3 - ro3 * du3 - aterm3
            psi23 = -u3 * du3 - dp3 / ro3
            psi33 = -u3 * dv3
            psi43 = -u3 * dp3 + a3 * a3 * u3 * dro3

            psi21b = 0.5 * (psi21 + psi23)
            psi31b = 0.5 * (psi31 + psi33)
            psi41b = 0.5 * (psi41 + psi43)
            psi12b = 0.5 * (psi12 + psi13)
            psi22b = 0.5 * (psi22 + psi23)
            psi32b = 0.5 * (psi32 + psi33)
            psi42b = 0.5 * (psi42 + psi43)
        else:
            psi21b, psi31b, psi41b = psi21, psi31, psi41
            psi12b, psi22b, psi32b, psi42b = psi12, psi22, psi32, psi42

        abr = geometry['nxny'][l]
        alb = 0.5 * (al2 + al) / als
        beb = be / als
        a1b = 0.5 * (a1 + a3)
        a2b = 0.5 * (a2 + a3)
        ro2b = 0.5 * (ro2 + ro3)

        u_new = (u1 - abr * v1 + (psi21b - abr * psi31b) * dt) / (1.0 + abr * abr)
        v_new = -u_new * abr
        p_rhs = -sign * ro2b * a2b * (alb * (u_new - u2) + beb * (v_new - v2))
        p_rhs += (psi42b + a2b * a2b * psi12b + sign * ro2b * a2b * (alb * psi22b + beb * psi32b)) * dt
        p_new = p2 + p_rhs / (144.0 * g)
        if p_new <= 0.0:
            p_new = params.get('plow', 0.001) * params['pc']

        ro_new = ro1 + (p_new - p1 - psi41b * dt) / max(a1b * a1b, 1e-12)
        if ro_new <= 0.0:
            ro_new = params.get('rolow', 0.001) / g

        u[l, m, n3], v[l, m, n3], p[l, m, n3], d[l, m, n3] = u_new, v_new, p_new, ro_new

    return u, v, p, d




def apply_exit_bc(u, v, p, d, n1, n3, dt, params, geometry):
    l_exit = params['lmax'] - 1
    mmax = params['mmax']
    mode = params.get('exit_bc_mode', 'pressure')

    for m in range(mmax):
        if mode == 'pressure':
            # Pressure-driven exit boundary (closer to EXITT behavior).
            p[l_exit, m, n3] = params['pe']
            u[l_exit, m, n3] = 2.0 * u[l_exit-1, m, n3] - u[l_exit-2, m, n3]
            v[l_exit, m, n3] = 2.0 * v[l_exit-1, m, n3] - v[l_exit-2, m, n3]
            d[l_exit, m, n3] = 2.0 * d[l_exit-1, m, n3] - d[l_exit-2, m, n3]
            if d[l_exit, m, n3] <= 0:
                d[l_exit, m, n3] = max(d[l_exit-1, m, n3], 1e-8)
        else:
            # Supersonic outflow extrapolation.
            u[l_exit, m, n3] = 2.0 * u[l_exit-1, m, n3] - u[l_exit-2, m, n3]
            v[l_exit, m, n3] = 2.0 * v[l_exit-1, m, n3] - v[l_exit-2, m, n3]
            p[l_exit, m, n3] = 2.0 * p[l_exit-1, m, n3] - p[l_exit-2, m, n3]
            d[l_exit, m, n3] = 2.0 * d[l_exit-1, m, n3] - d[l_exit-2, m, n3]

    return u, v, p, d


def apply_boundary_conditions(u, v, p, d, n1, n3, dt, params, geometry, step_type='predictor'):
    """
    Step 2c - Boundary Conditions wrapper.
    Passes the step_type to ensure MoC logic uses correct indices.
    """
    u, v, p, d = apply_inlet_bc(u, v, p, d, n1, n3, dt, params, geometry, step_type)
    u, v, p, d = apply_wall_bc(u, v, p, d, n1, n3, dt, params, geometry, step_type)
    u, v, p, d = apply_exit_bc(u, v, p, d, n1, n3, dt, params, geometry) # Exit usually 0-order
    
    return u, v, p, d


def check_convergence(u, n1, n3, params):
    """
    Step 2e - Check for steady-state convergence.
    Monitors the maximum relative change in axial velocity (u).
    """
    lmax = params['lmax']
    mmax = params['mmax']
    tconv = params['tconv']
    
    # Initialize maximum change tracker
    max_delu = 0.0
    
    # Iterate through all grid points (interior and boundaries)
    for l in range(lmax):
        for m in range(mmax):
            u_old = u[l, m, n1]
            u_new = u[l, m, n3]
            
            # Avoid division by zero at the start of the simulation
            if abs(u_old) > 1e-6:
                # Relative change: |(u_new - u_old) / u_old|
                delu = abs((u_new - u_old) / u_old)
                
                if delu > max_delu:
                    max_delu = delu
                    
    # The simulation has converged if the largest change is below TCONV
    converged = max_delu < tconv
    
    # Optional: Print progress every 100 steps to monitor stability
    if params.get('iteration_count', 0) % 100 == 0:
        print(f"Iteration {params.get('iteration_count', 0)}: Max Delta-U = {max_delu:.6e}")
        
    return converged

def calculate_performance(u, v, p, d, n, params, geometry):
    """
    Step 2f - Performance Metrics (MASFLO).
    Integrates flow properties and stores results in the params dictionary.
    """
    g = params['g']
    lmax = params['lmax']
    mmax = params['mmax']
    dy = 1.0 / (mmax - 1)  # Computational dy
    
    # Identify indices for Inlet (0), Throat (LT), and Exit (LDUM)
    # LDUM logic from masflo.f line 280
    lt = np.argmin(geometry['yw'])
    ldum = lmax - 1
    if params.get('jflag', 0) != 0:
        ldum = params.get('ljet', lmax) - 1
    
    # Initialize accumulators
    mass_i, mass_t, mass_e, thrust = 0.0, 0.0, 0.0, 0.0
    
    # Station heights (Wall - Centerbody)
    hi = geometry['yw'][0] - geometry['ycb'][0]
    ht = geometry['yw'][lt] - geometry['ycb'][lt]
    he = geometry['yw'][ldum] - geometry['ycb'][ldum]
    
    ndim = params.get('ndim', 1) 
    
    # Loop across the radial mesh (M=1 to M1 in Fortran)
    for m in range(mmax - 1):
        # Radial coordinates for the integration segments
        radi = m * dy * hi + geometry['ycb'][0]
        radt = m * dy * ht + geometry['ycb'][lt]
        rade = m * dy * he + geometry['ycb'][ldum]
        
        # Calculate area increments (AREAI, AREAT, AREAE)
        if ndim == 1: # Axisymmetric integration (PI * (R_outer^2 - R_inner^2))
            areai = np.pi * ((radi + dy * hi)**2 - radi**2)
            areat = np.pi * ((radt + dy * ht)**2 - radt**2)
            areae = np.pi * ((rade + dy * he)**2 - rade**2)
        else: # Planar integration (DY * Height)
            areai, areat, areae = dy * hi, dy * ht, dy * he
            
        # Average Momentum (rho * u) for the segment
        roui = 0.5 * (d[0, m, n] * u[0, m, n] + d[0, m+1, n] * u[0, m+1, n])
        rout = 0.5 * (d[lt, m, n] * u[lt, m, n] + d[lt, m+1, n] * u[lt, m+1, n])
        roue = 0.5 * (d[ldum, m, n] * u[ldum, m, n] + d[ldum, m+1, n] * u[ldum, m+1, n])
        
        # Momentum flux squared (rho * u^2) for thrust
        # Logic matches line 860: THRUST = THRUST + ROUE2 * AREAE
        roue2 = 0.5 * (d[ldum, m, n] * u[ldum, m, n]**2 + d[ldum, m+1, n] * u[ldum, m+1, n]**2)
        
        # Summation for mass flow (multiplied by G)
        mass_i += roui * areai * g
        mass_t += rout * areat * g
        mass_e += roue * areae * g
        thrust += roue2 * areae
        
    # Add results directly to the params dictionary
    params['mass_flow_inlet'] = mass_i
    params['mass_flow_throat'] = mass_t
    params['mass_flow_exit'] = mass_e
    params['thrust'] = thrust
    
    return params


def run_simulation(u, v, p, d, params, geometry):
    # The main time loop from Fortran: DO 650 N=1,NMAX.
    for n in range(params['nmax']):
        # Indices for current (n1) and next (n3) time steps
        n1, n3 = n % 2, (n + 1) % 2
        params['iteration_count'] = n
        
        # Step 2a - Calculate Time Step (CFL condition)
        # Calculates DT based on the Courant-Friedrichs-Lewy stability criterion.
        dt = calculate_time_step(u, v, p, d, n1, params)
        params['dt'] = dt # Required by calculate_shock_viscosity

        q_u = np.zeros((params['lmax'], params['mmax']))
        q_v = np.zeros((params['lmax'], params['mmax']))
        q_p = np.zeros((params['lmax'], params['mmax']))
        if params.get('iav', 0) == 1:
            # Step 2b - Calculate Shock Viscosity (SHOCK IPASS=1)
            q_u, q_v, q_p = calculate_shock_viscosity(u, v, p, d, n1, params, geometry, params['q_old'])
            # Save current Q for the next step's relaxation (CTA)
            params['q_old'] = {'u': q_u, 'v': q_v, 'p': q_p}

        # Step 2c - Predictor Step (INTER - first MacCormack Step)
        # Uses backward finite differences to predict flow at N+1.
        u, v, p, d = call_inter(u, v, p, d, n1, n3, dt, params, geometry, step_type='predictor')
        u, v, p, d = apply_boundary_conditions(u, v, p, d, n1, n3, dt, params, geometry, step_type='predictor')
        
        # Step 2d - Corrector Step (INTER - second MacCormack Step)
        # Uses forward differences and averages with predictor values for 2nd order accuracy.
        u, v, p, d = call_inter(u, v, p, d, n1, n3, dt, params, geometry, step_type='corrector')
        u, v, p, d = apply_boundary_conditions(u, v, p, d, n1, n3, dt, params, geometry, step_type='corrector')
        
        # Step 2e - Apply Shock Terms to the Corrected Solution
        # U(LMN3) = U(LMN3) + QUT(L,M)
        u[:, :, n3] += q_u
        v[:, :, n3] += q_v
        p[:, :, n3] += q_p
        
        # Step 2f - Global Flow Smoothing (SHOCK IPASS=2)
        u, v, p, d = smooth_flow_field(u, v, p, d, n3, params, geometry)
        
        # Step 2g - ConvergencyTCONV (check delu/u)
        # Checks if the maximum change in axial velocity is below tolerance.
        if check_convergence(u, n1, n3, params):
            print(f"Converged at step {n}")
            break    
        
        # Step 2h - Calculate mass flowrate and thrust (MASFLO)        
        calculate_performance(u, v, p, d, n3, params, geometry)            
        
    return u, v, p, d, params, geometry






if __name__ == '__main__':
    
    # Step 1 - Initialization (MAIN, GEOM, ONEDIM)
    u, v, p, d, params, geometry = initalize_simulation()
    
    # Step 2 - Time Integration Loop (MAIN)
    u, v, p, d, params, geometry = run_simulation(u, v, p, d, params, geometry)
    
    
    # Step 3 - Plotting and Reporting
    #plot_wall_grid(geometry)
    #plot_isentropic_data(u, v, p, d, geometry, params)
    #plot_computational_domain(geometry, params)
    plot_simulation_results(u, v, p, d, geometry, params)
    plot_pressure_contours(u, v, p, d, geometry, params)    
    plot_mach_contours(u, v, p, d, geometry, params)

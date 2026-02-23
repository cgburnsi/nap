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
              'ftd': 1.6,               # Fractional Delta Time (maybe?)
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
    Translates the logic from inter.f[cite: 4, 11].
    """
    gamma = params['gamma']
    dxr = 1.0 / params['dx']
    dyr = 1.0 / (1.0 / (params['mmax'] - 1)) # DYR = 1/DY
    ndim = params['ndim']
    lmax = params['lmax']
    mmax = params['mmax']
    
    # ICHAR logic from inter.f: Predictor = 1, Corrector = 2 [cite: 10]
    ichar = 1 if step_type == 'predictor' else 2
    
    # Start/End indices for loops 
    # Fortran DO 30 L=2,LMAX (Predictor) or DO 70 L=2,L1 (Corrector)
    l_start = 1
    l_end = (lmax - 1) if ichar == 1 else (lmax - 2)
    m_start = 0 # Symmetry axis (M=1 in Fortran)
    m_end = mmax - 1 # Wall (M1 in Fortran)

    for l in range(l_start, l_end + 1):
        for m in range(m_start, m_end): # Interior and axis points
            # Mapping values from Step 1e 
            al = geometry['al'][l, m]
            be = geometry['be'][l]
            de = geometry['de'][l, m]
            
            # Current values 
            ub = u[l, m, n1] if ichar == 1 else u[l, m, n3]
            vb = v[l, m, n1] if ichar == 1 else v[l, m, n3]
            pb = p[l, m, n1] if ichar == 1 else p[l, m, n3]
            rob = d[l, m, n1] if ichar == 1 else d[l, m, n3]
            asb = gamma * pb / rob
            
            # 1. Spatial Derivatives (Backward for Predictor, Forward for Corrector) 
            if ichar == 1: # Predictor (Backward differences)
                dudx = (ub - u[l-1, m, n1]) * dxr
                dpdx = (pb - p[l-1, m, n1]) * dxr
                drodx = (rob - d[l-1, m, n1]) * dxr
                dvdx = (vb - v[l-1, m, n1]) * dxr
                if m == 0: # Symmetry Axis
                    dvdy = (4.0 * v[l, 1, n1] - v[l, 2, n1]) * 0.5 * dyr
                    dudy = dvdy = dpdy = drody = 0.0
                else: # Interior
                    dudy = (ub - u[l, m-1, n1]) * dyr
                    dvdy = (vb - v[l, m-1, n1]) * dyr
                    dpdy = (pb - p[l, m-1, n1]) * dyr
                    drody = (rob - d[l, m-1, n1]) * dyr
            else: # Corrector (Forward differences)
                dudx = (u[l+1, m, n3] - ub) * dxr
                dpdx = (p[l+1, m, n3] - pb) * dxr
                drodx = (d[l+1, m, n3] - rob) * dxr
                dvdx = (v[l+1, m, n3] - vb) * dxr
                if m == 0: # Symmetry Axis
                    dvdy = (4.0 * v[l, 1, n3] - v[l, 2, n3]) * 0.5 * dyr
                    dudy = dvdy = dpdy = drody = 0.0
                else: # Interior
                    dudy = (u[l, m+1, n3] - ub) * dyr
                    dvdy = (v[l, m+1, n3] - vb) * dyr
                    dpdy = (p[l, m+1, n3] - pb) * dyr
                    drody = (d[l, m+1, n3] - rob) * dyr

            # 2. RHS Calculations 
            uvb = ub * al + vb * be + de
            aterm = (rob * vb / (m * (1.0/(mmax-1)) / be + geometry['ycb'][l])) if ndim == 1 else 0.0
            
            urhs = -ub * dudx - uvb * dudy - (dpdx + al * dpdy) / rob
            vrhs = -ub * dvdx - uvb * dvdy - be * dpdy / rob
            rorhs = -ub * drodx - uvb * drody - rob * (dudx + al * dudy + be * dvdy) - aterm
            prhs = -ub * dpdx - uvb * dpdy + asb * (rorhs + ub * drodx + uvb * drody)

            # 3. Step Forward 
            if ichar == 1: # Predictor
                u[l, m, n3] = u[l, m, n1] + urhs * dt
                v[l, m, n3] = v[l, m, n1] + vrhs * dt if m != 0 else 0.0
                p[l, m, n3] = p[l, m, n1] + prhs * dt
                d[l, m, n3] = d[l, m, n1] + rorhs * dt
            else: # Corrector (Average)
                u[l, m, n3] = (u[l, m, n1] + u[l, m, n3] + urhs * dt) * 0.5
                v[l, m, n3] = (v[l, m, n1] + v[l, m, n3] + vrhs * dt) * 0.5 if m != 0 else 0.0
                p[l, m, n3] = (p[l, m, n1] + p[l, m, n3] + prhs * dt) * 0.5
                d[l, m, n3] = (d[l, m, n1] + d[l, m, n3] + rorhs * dt) * 0.5

    return u, v, p, d


def apply_inlet_bc(u, v, p, d, n1, n3, dt, params, geometry):
    # Logic from inlet.f 
    gamma, rgas, g = params['gamma'], params['rgas'], params['g']
    mmax, dxr = params['mmax'], 1.0 / params['dx']
    dy = 1.0 / (mmax - 1)
    dyr = 1.0 / dy
    gam1 = gamma / (gamma - 1.0)
    gam2 = (gamma - 1.0) / 2.0

    for m in range(mmax):
        # Local mapping and initial values
        al, be, de = geometry['al'][0, m], geometry['be'][0], geometry['de'][0, m]
        u2, p2, d2 = u[0, m, n1], p[0, m, n1], d[0, m, n1]
        a2 = np.sqrt(gamma * p2 * 144.0 * g / d2)
        
        # Characteristic backstep (X2) calculation
        # x2 is the location upstream from where the characteristic originated
        x2 = (u2 - a2) * 0.5 * dt
        
        # Linear interpolation for properties at X2
        # Here we simplify the BU/BP logic from the Fortran loop
        u_interp = u2 + (u[1, m, n1] - u2) * dxr * x2
        p_interp = p2 + (p[1, m, n1] - p2) * dxr * x2
        d_interp = d2 + (d[1, m, n1] - d2) * dxr * x2
        
        # Solve for boundary Mach (MN3) via iteration
        mn3 = np.sqrt(u[0, m, n3]**2 + v[0, m, n3]**2) / a2 if u[0, m, n3] != 0 else 0.01
        for _ in range(10):
            dem = 1.0 + gam2 * mn3**2
            p_new = params['pc'] / (dem**gam1)
            t_new = params['tc'] / dem
            
            # Compatibility equation for U
            # U = U2 + (P_new - P2) / (rho * a) ...
            u[0, m, n3] = u_interp + (p_new - p_interp) / (d_interp * a2)
            v[0, m, n3] = 0.0 # Standard axial inlet assumption
            
            mn3_new = u[0, m, n3] / np.sqrt(gamma * rgas * g * t_new)
            if abs(mn3_new - mn3) < 0.001: break
            mn3 = mn3_new
            
        d[0, m, n3] = p_new * 144.0 / (rgas * t_new * g)
        p[0, m, n3] = p_new
        
    return u, v, p, d


def apply_wall_bc(u, v, p, d, n1, n3, dt, params, geometry):
    # Logic from wall.f [cite: 27]
    # IB=1 for Wall (m=MMAX), IB=2 for Symmetry Axis (m=0)
    for ib in [1, 2]:
        m_idx = (params['mmax'] - 1) if ib == 1 else 0
        sign = 1.0 if ib == 1 else -1.0
        
        for l in range(1, params['lmax']):
            # Tangency slope [cite: 32]
            slope = geometry['nxny'][l] if ib == 1 else geometry['nxnycb'][l]
            
            # Simplified Wall Characteristic solve
            # U = (U_old - slope * V_old) / (1 + slope^2)
            u_wall = (u[l, m_idx, n1] - slope * v[l, m_idx, n1]) / (1.0 + slope**2)
            v_wall = -u_wall * slope
            
            u[l, m_idx, n3] = u_wall
            v[l, m_idx, n3] = v_wall
            # P and RO are typically updated via compatibility equations
            # For now, we use a simple zero-gradient extrapolation for pressure
            p[l, m_idx, n3] = p[l, (m_idx-1 if ib==1 else m_idx+1), n3]
            d[l, m_idx, n3] = d[l, (m_idx-1 if ib==1 else m_idx+1), n3]
            
    return u, v, p, d


def apply_exit_bc(u, v, p, d, n1, n3, dt, params, geometry):
    # Logic from exitt.f [cite: 12]
    l_exit = params['lmax'] - 1
    
    for m in range(params['mmax']):
        # Simple Supersonic Outflow (Extrapolation) [cite: 12]
        # In a fully supersonic exit, we extrapolate from the interior
        u[l_exit, m, n3] = 2.0 * u[l_exit-1, m, n3] - u[l_exit-2, m, n3]
        v[l_exit, m, n3] = 2.0 * v[l_exit-1, m, n3] - v[l_exit-2, m, n3]
        p[l_exit, m, n3] = 2.0 * p[l_exit-1, m, n3] - p[l_exit-2, m, n3]
        d[l_exit, m, n3] = 2.0 * d[l_exit-1, m, n3] - d[l_exit-2, m, n3]
        
    return u, v, p, d

def apply_boundary_conditions(u, v, p, d, n1, n3, dt, params, geometry):
    # Step 2c - Boundary Conditions wrapper
    u, v, p, d = apply_inlet_bc(u, v, p, d, n1, n3, dt, params, geometry)
    u, v, p, d = apply_wall_bc(u, v, p, d, n1, n3, dt, params, geometry)
    u, v, p, d = apply_exit_bc(u, v, p, d, n1, n3, dt, params, geometry)
    
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
        n1 = n % 2
        n3 = (n + 1) % 2
        
        # Step 2a - Calculate Time Step (CFL condition)
        # Calculates DT based on the Courant-Friedrichs-Lewy stability criterion.
        dt = calculate_time_step(u, v, p, d, n1, params)

        # Step 2b - Predictor Step (INTER - first MacCormack Step)
        # Uses backward finite differences to predict flow at N+1.
        u, v, p, d = call_inter(u, v, p, d, n1, n3, dt, params, geometry, step_type='predictor')

        # Step 2c - Boundary Conditions (INLET, WALL, EXITT)
        # Updates the predicted values at the boundaries using characteristic methods[cite: 1, 2].
        u, v, p, d = apply_boundary_conditions(u, v, p, d, n1, n3, dt, params, geometry)

        # Step 2d - Corrector Step (INTER - second MacCormack Step)
        # Uses forward differences and averages with predictor values for 2nd order accuracy.
        u, v, p, d = call_inter(u, v, p, d, n1, n3, dt, params, geometry, step_type='corrector')
        
        # Re-apply boundaries to the corrected solution.
        u, v, p, d = apply_boundary_conditions(u, v, p, d, n1, n3, dt, params, geometry)

        # Step 2e - ConvergencyTCONV (check delu/u)
        # Checks if the maximum change in axial velocity is below tolerance.
        if check_convergence(u, n1, n3, params):
            print(f"Converged at step {n}")
            break    
        
        # Step 2f - Calculate mass flowrate and thrust (MASFLO)        
        calculate_performance(u, v, p, d, n, params, geometry)
            
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
    
    




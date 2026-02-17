import numpy as np
import convert as cv
import matplotlib.pyplot as plt

''' The following are default values from the Fortran Code
N1D = 1, NDIM = 1, NGCB = 0, NGEOM = 2

IUI = 1, IU0 = 1

'''



# Geometry Variables
l_max      = 20                                     # [-] Maximum x-axis (Axial Direction) grid points
m_max      = 8                                      # [-] Maximum y-axis (Radial Direction) grid points
n_max      = 400                                    # [-] Maximum number of time steps
x_i        = cv.convert(0.31, 'in', 'm')            # [m] x-axis location for inlet point
x_t        = np.nan                                 # [m] x-axis location for throat (computed later)
r_i        = cv.convert(2.5, 'in', 'm')             # [m] y-axis location for wall inlet point
r_t        = cv.convert(0.8, 'in', 'm')             # [m] y-axis location for the nozzle throat
x_e        = cv.convert(4.05, 'in', 'm')             # [m] x-axis location for nozzle exit
r_e        = np.nan                                 # [m] y-axis location for nozzle exit (computed later)
rc_i       = cv.convert(0.8, 'in', 'm')             # [m] Radius of Curvature at nozzle wall inlet contraction
rc_t       = cv.convert(0.5, 'in', 'm')             # [m] Radius of Curvature at nozzle throat
ang_i      = np.deg2rad(44.88)                      # [radians] Contraction Angle
ang_e      = np.deg2rad(15.0)                       # [radians] Nozzle Exit Angle
r_star     = np.nan                                 # [m] Area per unit depth where Mach = 1
r_stars    = np.nan                                 # [m^2] Area per unit depth divided by np.pi
a_i         = np.pi * r_i**2                        # [m^2] Inlet Cross Sectional Area
a_t        = np.pi * r_t**2                         # [m^2] Throat Area
a_e        = np.nan                                 # [m^2] Exit Area

# Grid Variables
dx         = (x_e - x_i) / (l_max - 1.0)            # TODO (Should this be L_max - 0) [m] distance between equally spaced x axis points   
dy         = 1.0 / (m_max - 1.0)                    # TODO (Should this be m_max-0?) [m] distance between equally spaced y-axis points
yw         = np.zeros(l_max)                        # [m] y-axis location of wall (equally spaced)
xw         = np.zeros(l_max)                               # [m] x-axis location of wall (equally spaced)
yw_i       = np.zeros(l_max)                                 # [m] y-axis location of wall (non-equally spaced)
xw_i       = np.zeros(l_max)                                # [m] x-axis location of wall (non-equally spaced)
nxny       = np.zeros(l_max)                                # [-] negative of wall slopes corresponding to yw elements

# Boundary Conditions
PT         = np.full(m_max, cv.convert(70, 'psi', 'Pa'))                     # [Pa] Inlet boundary pressure (uniform across inlet)
TT         = np.full(m_max, cv.convert(80, 'degF', 'K'))                     # [K] Inlet boundary temperature (uniform across inlet)
theta      = np.full(m_max, np.deg2rad(0))                  # [radians] Inlet boundary flow angle (uniform across inlet)
PE         = 14.7                                     # [Pa] Exit boundary pressure condition

# Fluid Properties
g          = 9.807                                         # [m/s^2] acceleration due due to gravity on Earth
R_gas      = cv.convert(53.3523, 'ft*lbf/lb*degR', 'J/kg*K', )               # [J/kg-K] Specific gas constant
gamma      = 1.4                                         # [-] Ratio of specific heats for fluid

# Calculation Helper Values
gam_1      = gamma / (gamma -1)                  # [-] Calculation Helper
gam_2      = (gamma - 1) / 2.0                        # [-] Calculation Helper
l_t        = 0                                             # [m] Index for the throat location in the yw array
l_1        = l_max - 1                                # [-] locator value for the edge of the grid in x-axis direction
l_2        = l_max - 2                                # [-] locator value for the edge of the grid in x-axis direction
l_3        = l_max - 3                                # [-] locator value for the edge of the grid in x-axis direction
m_1        = m_max - 1                                # [-] locator value for the edge of the grid in y-axis direction
m_2        = m_max - 2                                # [-] locator value for the edge of the grid in y-axis direction

# Simulation Parameters
t_conv     = 3e-5                                          # [-] Axial velocity steady-state convergence criteria
t_stop     = 1.0                                        # [sec] Time stop for simulation
FDT        = 1.6                                           # [-] Premultiplier on the C-F-L paramter

# Solution Arrays
U          = np.zeros((1, m_max, l_max))         # [m/s] Velocity x-axis (Axial Direction)
V          = np.zeros((1, m_max, l_max))         # [m/s] Velocity y-axis (Radial Direction)
P          = np.zeros((1, m_max, l_max))         # [Pa] Pressure Field
RO         = np.zeros((1, m_max, l_max))         # [kg/m^3] Density Field
T          = np.zeros((1, m_max, l_max))         # [K] Temperature Field
Mach       = np.zeros((1, m_max, l_max))         # [-] Mach Number f


def plot_wall_contour(xw, yw):
    plt.figure(1, figsize=(10, 5))
    plt.cla()
    plt.plot(xw, yw, label='Outer Wall')
    plt.xlabel("x-axis [m]")
    plt.ylabel("y-axis [m]")
    plt.grid()
    plt.legend()

def calculate_wall_geometry(l_max, dx, xw, yx, nxny, x_i, r_i, r_t, rc_i, rc_t, ang_i, ang_e):
    x_tan = x_i + rc_i * np.sin(ang_i)
    r_tan = r_i + rc_i * (np.cos(ang_i) - 1.0)
    r_t_1 = r_t - rc_t * (np.cos(ang_i) - 1.0)
    x_t_1 = x_tan + (r_tan - r_t_1) / np.tan(ang_i)

    if x_t_1 < x_tan: x_t_1 = x_tan; r_t_1 = r_tan
    
    x_t    = x_t_1 + rc_t * np.sin(ang_i)
    x_t_2  = x_t + rc_t * np.sin(ang_e)
    r_t_2  = r_t + rc_t * (1.0 - np.cos(ang_e))
    r_e    = r_t_2 + (x_e - x_t_2) * np.tan(ang_e)
    l_t    = 0
    a_e    = np.pi * r_e**2

    for L in range(0, l_max):
        xw[L] = x_i + dx * L
        if xw[L] >= x_i and xw[L] <= x_tan:
            yw[L]      = r_i + rc_i * (np.cos(np.arcsin((xw[L] - x_i) / rc_i)) - 1.0)
            nxny[L]    = (xw[L] - x_i) / (yw[L] - r_i + rc_i)
        elif xw[L] > x_tan and xw[L] <= x_t_1:
            yw[L]      = r_t_1 + (x_t_1 - xw[L]) * np.tan(ang_i)
            nxny[L]    = np.tan(ang_i)
        elif xw[L] > x_t_1 and xw[L] <= x_t:
            yw[L]      = r_t + rc_t * (1.0 - np.cos(np.arcsin((x_t - xw[L]) / rc_t)))
            nxny[L]    = (x_t - xw[L]) / (rc_t + r_t - yw[L])
        elif xw[L] > x_t and xw[L] <= x_t_2:
            yw[L]      = r_t + rc_t * (1.0-np.cos(np.arcsin((xw[L] - x_t) / rc_t)))
            nxny[L]    = (x_t - xw[L]) / (rc_t + r_t - yw[L])
        elif xw[L] > x_t_2 and xw[L] <= x_e:
            yw[L]      = r_t_2 + (xw[L] - x_t_2) * np.tan(ang_e)
            nxny[L]    = -np.tan(ang_e)
        
        # Determine where the index is for the throat in the yw (outer wall) array
        if L <= 0: 
            continue
        else:
            if yw[L] < yw[l_t]: 
                l_t = L;

    return x_t, l_t, r_e, a_e, xw, yw, nxny

def Mach_f(guess, k, expan):
    a1 = k+1; a2 = k-1; a3 = a1/a2; a4 = 2/a1; a5 = a2/2; a6 = a3/2
    return (expan - (1/guess)*(a4*(1+a5*guess**2))**a6)

def Mach_df(guess, k):
    a1 = k+1; a2 = k-1; a3 = a1/a2; a4 = 2/a1; a5 = a2/2; a6 = a3/2 
    return ((a4*(a5*guess**2 + 1))**a6/guess**2 - a3*a4*a5*(a4*(a5*guess**2 + 1))**(a6 - 1))
    
def Newton_Raphson(M_guess, A_ratio, gamma):
    #Useful Parameters
    tol = 1e-4          # [-] Calculation tolerance for convergence
    MaxIter = 200       # [#] Maximum number of iterations
    j = 1               # [#] Start at 1 because first step is manual calc
    
    # Calculate the initial point
    f = Mach_f(M_guess, gamma, A_ratio)
    df = Mach_df(M_guess, gamma)
    M = np.zeros(MaxIter+1)    
    M[0] = M_guess - (f/df)
    
    # Iterate using Newton-Raphson to determine the actual mach number    
    while (M[j-1] >=tol) and (j<=MaxIter):
        f = Mach_f(M[j-1], gamma, A_ratio)
        df = Mach_df(M[j-1], gamma)
        M[j] = M[j-1] - (f/df)
        j = j+1
        
    M_last = M[-1]
    return M_last

def one_dimisional(l_max, m_max, R_gas, gamma, gam_1, gam_2, x_t, a_t, yw, Mach, PT, TT, T, U, V, RO, P):
    # Calculation helper values
    GRGAS = 1.0 / (R_gas)
      
    # Calculate the local One-Dimensional Mach Number
    for L in range(0, l_max):
        x = x_i + dx * (L)   # (Maybe this should be L and not L-1)??                         # [m] Current location along the x-axis direction at 'M=0'
        if x < x_t:                                    # Subsonic region of nozzle upstream of throat
            area = np.pi * yw[L]**2                    # [m**2] Area at current location
            A_ratio = area / a_t                       # [-] Area Ratio at current location
            MN3 = Newton_Raphson(0.1, A_ratio, gamma)          # [-] Solve for the local Mach number (subsonic Mach Guess)
        elif x >= x_t:
            area = np.pi * yw[L]**2                    # [m**2] Area at current location
            A_ratio = area / a_t                       # [-] Area Ratio at current location
            MN3 = Newton_Raphson(5.0, A_ratio, gamma)              # [-] Solve for the local Mach number (supersonic Mach Guess)
    
    # Calculate P, RO, U, and V for the initial surface from the local Mach number             
        for M in range(0, m_max):
            DEM = 1.0 + gam_2 * MN3**2
            DEMP = DEM**gam_1
            Mach[0, M, L] = MN3
            P[0][M][L] = PT[M] / DEMP
            RO[0][M][L] = P[0][M][L] * GRGAS / (TT[M] / DEM)
            T[0][M][L] = P[0][M][L] / (RO[0][M][L] * R_gas);
            Q = MN3 * np.sqrt(gamma * P[0][M][L] / RO[0][M][L])
            U[0][M][L] = Q
            V[0][M][L] = 0.0

        #print(f'L = {L:<2}, x = {x:<3.2f}, U = {U[0, 0, L]:<3.2f}, V = {V[0, 0, L]:<3.2f}, P = {P[0, 0, L]:<3.2f}, RO = {RO[0, 0, L]:<3.4f} Area Ratio = {A_ratio:<3.5f}, Q = {Q:<3.2f}, Mach {MN3:<3.2f}')

    return Mach, U, V, RO, P, T

def masflo(l_max, yw, U, V, RO, a_t, a_e):
    # Helper Calculations
    LDUM = l_max-1
    
    # Calculate the mass flow and thrust for the 1-D Initial Data Surface
    a_inlet     = (np.pi * yw[0]**2)
    #a_throat    = (np.pi * self.yw[self.l_t]**2)
    a_exit      = (np.pi * yw[LDUM]**2)
    vm_i        = np.sqrt(U[0, 0, 0]**2 + V[0, 0, 0]**2)
    #vm_t        = np.sqrt(U[0, 0, l_t]**2 + V[0, 0, l_t]**2)
    vm_e        = np.sqrt(U[0, 0, LDUM]**2 + V[0, 0, LDUM]**2)
    mass_i      = RO[0, 0, 0] * vm_i * a_inlet; print(mass_i, a_inlet)
    #mass_t      = self.RO[0, 0, self.l_t] * vm_t * a_throat 
    #mass_e      = self.RO[0, 0, LDUM] * vm_e * a_exit 
    thrust      = RO[0, 0, LDUM] * U[0, 0, LDUM]**2 * a_exit
    
    #print(f'vm_e = {vm_e:<3.2f}, mass_i = {mass_i:<3.4f}, mass_t = {mass_t:<3.4f}, mass_e = {mass_e:<3.4f}, thrust = {thrust:<3.2f}')
    print(f'Expansion Ratio = {a_e / a_t:<3.2f} [-], vel_exit = {vm_e:<3.2f} [m/s], mdot = {cv.convert(mass_i, "kg/s", "g/s"):<3.5f} [g/s], thrust = {cv.convert(thrust, "N", "mN"):<3.2f} [mN]')
    print(f'Expansion Ratio = {a_e / a_t:<3.2f} [-], vel_exit = {cv.convert(vm_e, "m/s", "ft/s"):<3.2f} [m/s], mdot = {cv.convert(mass_i, "kg/s", "lb/s"):<3.5f} [lbm/s], thrust = {cv.convert(thrust, "N", "lbf"):<3.2f} [lbf]')

def area_ratio_equation(M, gamma):
    term1 = 1 / M
    term2 = (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
    term3 = (1 + (gamma - 1) / 2 * M**2) ** ((gamma + 1) / (2 * (gamma - 1)))
    
    return term1 * term2 * term3

def print_one_dim_data_eng_units(l_max, xw, U, V, P, RO, Mach, T):
    cols = {"L [-]": 5, "x [in]": 8, "U [ft/s]": 10, "V [ft/s]": 10, "P [psi]": 8, "RO [lbm/ft3]": 12, "Mach [-]": 10, "T [degF]": 10}

    headers = list(cols.keys())
    # Print headers
    header_row = "|".join(f"{h:>{cols[h]}}" for h in headers)
    print(header_row)
    print("-" * len(header_row))

    # Print Data
    for idx in range(0, l_max):
        print(f"{idx:>{cols['L [-]']}}|{cv.convert(xw[idx], 'm', 'in'):{cols['x [in]']}.2f}|"
              f"{cv.convert(U[0, 0, idx], 'm/s', 'ft/s'):>{cols['U [ft/s]']}.2f}|"
              f"{cv.convert(V[0, 0, idx], 'm/s', 'ft/s'):>{cols['V [ft/s]']}.2f}|"
              f"{cv.convert(P[0, 0, idx], 'Pa', 'psi'):>{cols['P [psi]']}.2f}|"
              f"{cv.convert(RO[0, 0, idx], 'kg/m^3', 'lb/ft^3'):>{cols['RO [lbm/ft3]']}.4f}|"
              f"{Mach[0,0, idx]:>{cols['Mach [-]']-2}.2f}|{cv.convert(T[0, 0, idx], 'K', 'degF'):>{cols['T [degF]']}.2f}")



def inter(i_char, MDUM): # Calculate the Interior Mesh Points
    ATERM = 0.0
    
    if i_char == 2:
        ...
        # Go to 40
    else:
        # Compute the tentative Solution at T + DT
        MDUM = 1




def map(IP, L, M, yw, dy, nxny, xw_i, LD1):
    # Mapping the physical plane to the rectangular computational plane
    BE = 1 / yw[L]
    if IP == 0: return
    
    y = (M-1)*dy
    AL=BE * y*nxny[L]
    DE = -BE*y*xw_i[L]                                                
    if IP == 1: return
    
    BE1 = 1.0 / yw[LD1]                                        
    AL1 = BE1 * y*nxny[LD1]
    DE1 = -BE1 * y*xw_i[LD1]                                               
    
    return BE, BE1
    
    
    
    

if __name__ == '__main__':
    x_t, l_t, r_e, a_e, xw, yw, nxny = calculate_wall_geometry(l_max, dx, xw, yw, nxny, x_i, r_i, r_t, rc_i, rc_t, ang_i, ang_e)
    
    Mach, U, V, RO, P, T = one_dimisional(l_max, m_max, R_gas, gamma, gam_1, gam_2, x_t, a_t, yw, Mach, PT, TT, T, U, V, RO, P)
    
    print_one_dim_data_eng_units(l_max, xw, U, V, P, RO, Mach, T)
    
    masflo(l_max, yw, U, V, RO, a_t, a_e)
        
    #plot_wall_contour(xw, yw)
 
    LD1 = 0         # I don't know where this is defined in the fortran code.

    
    # Calculate the Initial Value Surface
    for IU in range(0,2):
        x = x_i - dx
        for L in range(0, l_max):
            x = x + dx
            AL, AL1, BE, BE1 = map(0, L, 0, yw,)




    # Things needed for inter()
    i_char = 1
    i_b = 1














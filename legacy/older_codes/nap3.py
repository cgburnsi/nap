
import numpy as np
import matplotlib.pyplot as plt


class Solution:
    def __init__(self, l_max=81, m_max=21, n_max=400, x_i=0.007874, r_i=0.0635, r_t=0.02032, x_e=0.1143, rc_i=0.02032, rc_t=0.0127, 
                 ang_i=44.88, ang_e=15, PT=482633, TT=299.817, THETA=0, PE=101352.9, t_conv=3e5, t_stop=1.0, gamma=1.4,
                 R_gas=287.052874, g=9.807):
        self.l_max      = l_max                                         # [-] Maximum x-axis (Axial Direction) grid points
        self.m_max      = m_max                                         # [-] Maximum y-axis (Radial Direction) grid points
        self.n_max      = n_max                                         # [-] Maximum number of time steps
        self.x_i        = x_i                                           # [m] x-axis location for inlet point
        self.x_t        = np.nan                                        # [m] x-axis location for throat (computed later)
        self.r_i        = r_i                                           # [m] y-axis location for wall inlet point
        self.r_t        = r_t                                           # [m] y-axis location for the nozzle throat
        self.x_e        = x_e                                           # [m] x-axis location for nozzle exit
        self.r_e        = np.nan                                        # [m] y-axis location for nozzle exit (computed later)
        self.rc_i       = rc_i                                          # [m] Radius of Curvature at nozzle wall inlet contraction
        self.rc_t       = rc_t                                          # [m] Radius of Curvature at nozzle throat
        self.ang_i      = np.deg2rad(ang_i)                             # [radians] Contraction Angle
        self.ang_e      = np.deg2rad(ang_e)                             # [radians] Nozzle Exit Angle
        self.r_star     = np.nan                                        # [m] Area per unit depth where Mach = 1
        self.r_stars    = np.nan                                        # [m^2] Area per unit depth divided by np.pi
        
        self.dx         = (self.x_e - self.x_i) / (self.l_max - 1.0)    # TODO (Should this be L_max - 0) [m] distance between equally spaced x axis points   
        self.dy         = 1.0 / (self.m_max - 1.0)                      # TODO (Should this be m_max-0?) [m] distance between equally spaced y-axis points
        self.yw         = np.zeros(self.l_max)                          # [m] y-axis location of wall (equally spaced)
        self.xw         = np.zeros_like(self.yw)                        # [m] x-axis location of wall (equally spaced)
        self.yw_i       = np.zeros_like(self.yw)                        # [m] y-axis location of wall (non-equally spaced)
        self.xw_i       = np.zeros_like(self.yw)                        # [m] x-axis location of wall (non-equally spaced)
        self.nxny       = np.zeros_like(self.yw)                        # [-] negative of wall slopes corresponding to yw elements
        
        # Boundary Conditions
        self.PT         = np.full(self.m_max, PT)                       # [Pa] Inlet boundary pressure (uniform across inlet)
        self.TT         = np.full(self.m_max, TT)                       # [K] Inlet boundary temperature (uniform across inlet)
        self.theta      = np.full(self.m_max, THETA)                    # [radians] Inlet boundary flow angle (uniform across inlet)
        self.PE         = PE                                            # [Pa] Exit boundary pressure condition
        
        # Fluid Properties
        self.g          = g                                             # [m/s^2] acceleration due due to gravity
        self.R_gas      = R_gas                                         # [J/kg-K] Specific gas constant
        self.gamma      = gamma                                         # [-] Ratio of specific heats for fluid
       
        # Calculation Helper Values
        self.gam_1      = self.gamma / (self.gamma -1)                  # [-] Calculation Helper
        self.gam_2      = (self.gamma - 1) / 2.0                        # [-] Calculation Helper
        self.l_t        = 0                                             # [?] Not really sure what this does yet (integer?)
        self.l_1        = self.l_max - 1                                # [-] locator value for the edge of the grid in x-axis direction
        self.l_2        = self.l_max - 2                                # [-] locator value for the edge of the grid in x-axis direction
        self.l_3        = self.l_max - 3                                # [-] locator value for the edge of the grid in x-axis direction
        self.m_1        = self.m_max - 1                                # [-] locator value for the edge of the grid in y-axis direction
        self.m_2        = self.m_max - 2                                # [-] locator value for the edge of the grid in y-axis direction

        # Simulation Parameters
        self.t_conv     = t_conv                                        # [-] Axial velocity steady-state convergence criteria
        self.t_stop     = t_stop                                        # [sec] Time stop for simulation
        
        # Solution Arrays
        self.U          = np.zeros((self.l_max, self.m_max, 2))         # [m/s] Velocity x-axis (Axial Direction)
        self.V          = np.zeros_like(self.U)                         # [m/s] Velocity y-axis (Radial Direction)
        self.P          = np.zeros_like(self.U)                         # [Pa] Pressure Field
        self.RO         = np.zeros_like(self.U)                         # [kg/m^3] Density Field

    def build_wall_geometry(self):
        x_tan = self.x_i + self.rc_i * np.sin(self.ang_i)
        r_tan = self.r_i + self.rc_i * (np.cos(self.ang_i) - 1.0)
        r_t_1 = self.r_t - self.rc_t * (np.cos(self.ang_i) - 1.0)
        x_t_1 = x_tan + (r_tan - r_t_1) / np.tan(self.ang_i)
        
        if x_t_1 < x_tan: x_t_1 = x_tan; r_t_1 = r_tan

        self.x_t    = x_t_1 + self.rc_t * np.sin(self.ang_i)
        x_t_2       = self.x_t + self.rc_t * np.sin(self.ang_e)
        r_t_2       = self.r_t + self.rc_t * (1.0 - np.cos(self.ang_e))
        self.r_e    = r_t_2 + (self.x_e - x_t_2) * np.tan(self.ang_e)

        for L in range(0, self.l_max):
            x = self.x_i + self.dx * L
            self.xw[L] = x
            if x >= self.x_i and x <= x_tan:
                self.yw[L]      = self.r_i + self.rc_i * (np.cos(np.arcsin((x - self.x_i) / self.rc_i)) - 1.0)
                self.nxny[L]    = (x - self.x_i) / (self.yw[L] - self.r_i + self.rc_i)
            elif x > x_tan and x <= x_t_1:
                self.yw[L]      = r_t_1 + (x_t_1 - x) * np.tan(self.ang_i)
                self.nxny[L]    = np.tan(self.ang_i)
            elif x > x_t_1 and x <= self.x_t:
                self.yw[L]      = self.r_t + self.rc_t * (1.0 - np.cos(np.arcsin((self.x_t - x) / self.rc_t)))
                self.nxny[L]    = (self.x_t - x) / (self.rc_t + self.r_t - self.yw[L])
            elif x > self.x_t and x <= x_t_2:
                self.yw[L]      = self.r_t + self.rc_t * (1.0-np.cos(np.arcsin((x - self.x_t) / self.rc_t)))
                self.nxny[L]    = (self.x_t - x) / (self.rc_t + self.r_t - self.yw[L])
            elif x > x_t_2 and x <= self.x_e:
                self.yw[L]      = r_t_2 + (x - x_t_2) * np.tan(self.ang_e)
                self.nxny[L]    = -np.tan(self.ang_e)
            
            # This is somewhat of a mystery still.  I'm not sure what it is doing.
            if L > 0:
                if self.yw[L] > self.yw[self.l_t]: self.l_t = L;

    def plot_wall_contour(self):
        plt.figure(figsize=(10, 5))
        plt.plot(self.xw, self.yw, label='Outer Wall')
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (Torr)")
        plt.grid()
        plt.legend()

    def one_dimisional(self):
        # Calculation helper values
        GRGAS = 1.0 / (self.R_gas * self.g)
        ACOEF = 2.0 / (self.gamma + 1.0)
        BCOEF = (self.gamma - 1.0) / (self.gamma + 1.0)
        CCOEF = (self.gamma + 1.0) / 2.0 / (self.gamma - 1.0)
        
        MN3 = 0.01
        NXDK = 0.0

        self.r_star = self.r_t; self.r_stars = self.r_t**2
        
        for L in range(0, self.l_max):
            x = self.x_i + self.dx * L
            if x < self.x_t:
                RADS = self.yw[L]**2 
                ARATIO = RADS / self.r_stars
                
                # Newton Raphson Iteration Loop
                for _ in range(0, 19):
                    ABM = ACOEF + BCOEF * MN3**2
                    ABMC = ABM**CCOEF
                    FM = ABMC / MN3 - ARATIO
                    FPM = ABMC * (2.0 * BCOEF * CCOEF/ABM-1.0/MN3**2)
                    OMN3 = MN3
                    MN3 = OMN3 - FM/FPM
                    if MN3 > 1.0 and OMN3 < 1.0:  MN3 = 0.99
                    if MN3 < 1.0 and OMN3 > 1.0:  MN3 = 1.01
                    if MN3 < 0.0: MN3 = -MN3
                    if (np.abs(MN3 - OMN3) / OMN3) <= 0.0005:
                        print(L, 'hello from Newton Raphson Solver')
                        # Fill in 2-D arrays
                        DEM = 1.0 + self.gam_2 * MN3**2
                        DEMP = DEM**self.gam_1
                        DNXNY = self.nxny[L] / self.m_1
                        for M in range(0, self.m_max):
                            self.P[L][M][0] = self.PT[M] / DEMP
                            TEMP = self.TT[M] / DEM
                            self.RO[L][M][0] = self.P[L][M][0] * GRGAS / TEMP
                            Q = MN3 * np.sqrt(self.gamma * self.P[L][M][0] / self.RO[L][M][0])
                            DN = DNXNY * (M-1)   
                            DNS = DN**2
                            if DNS == 0:
                                self.U[L][M][0] = Q
                                self.V[L][M][0] = 0.0
                            else:
                                SIGN = 1.0
                                if DN > 0.0: 
                                    SIGN = -1.0
                                self.U[L][M][0] = Q / np.sqrt(1.0 + DNS)
                                self.V[L][M][0] = SIGN * Q / np.sqrt(1.0 + (1.0 / DNS))
                                









def get_default_dicts():
    DESC = {'TITLE': ''}
    
    # Control Parameters
    CNTRLC = {'LMAX': 80, 'MMAX': 20, 'NMAX': 0, 'NPRINT': 0, 'TCONV': 0, 'FDT': 1.0, 'TSTOP': 1.0, 
              'GAMMA': 1.4, 'RGAS': 53.35, 'GAM1': 0, 'GAM2': 0, 'L1': 0, 'L2': 0, 'L3': 0, 
              'M1': 0, 'M2': 0, 'DX': 0, 'DY':0, 'DT':0, 'N': 0, 'N1':0, 'N3':0, 'NASM': 1, 
              'IVEL':0, 'ICHAR':0, 'N1D': 1, 'LJET': 0, 'JFLAG': 0, 'IERR': 0, 'IUI': 1, 'IUO': 1, 
              'DXR':0, 'DYR':0, 'LD':0, 'MD':0, 'LMD1':0, 'LMD3':0, 'IB':0, 'RSTAR': 0, 
              'RSTARS': 0, 'NPLOT': -1, 'G': 32.174, 'PC': 144.0, 'TC': 460.0, 'LC': 12.00, 'PLOW':0.01, 
              'ROLOW': 0.0001, 'NAME': 0, 'IPUNCH': 0, 'IEX': 1, 'NCONVI': 1, 'IUNIT': 0}

    AV = {'IAV': 0, 'CAV': 4, 'NST': 0, 'SMP': 0.95, 'LSS': 2, 'CTA': 0.5, 'XMU': 0.2, 
          'XLA': 1.0, 'RKMU': 0.7, 'QUT': np.zeros((80, 20)), 'QVT': np.zeros((80, 20)), 
          'QPT': np.zeros((80, 20))}

    ONESID = {'UD':4, 'VD':4, 'PD': 4, 'ROD':4}

    SOLUTN = {'NSTART': 0, 'TSTART': 0, 
              'U': np.zeros((80, 20, 2)),
              'V': np.zeros((80, 20, 2)),
              'P': np.zeros((80, 20, 2)),
              'RO':np.zeros((80, 20, 2))}

    GEMTRYC = {'NGEOM': 0, 'XI': 0.31, 'RI': 2.5, 'XT': 0, 'RT': 0.8, 'XE': 4.05, 'RE': 0, 'RCI': 0.8,
               'RCT': 0.5, 'ANGI': 44.88, 'ANGE': 15, 'XW': np.zeros(80), 
               'YW': np.zeros(80), 'XWI':np.zeros(80), 
               'YWI':np.zeros(80), 'NXNY': np.zeros(80), 'NWPTS': 0, 
               'IINT': 1, 'IDIF': 1, 'LT': 0, 'NDIM': 1, 'DX': 0}
    
    GCB = {'NGCB':0, 'XICB':0, 'RICB':0, 'XTCB':0, 'RTCB':0, 'XECB':0, 'RECB':0, 'RCICB':0, 'RCTCB':0,
           'ANGICB':0, 'ANGECB':0, 'XCB':np.zeros(80), 'YCB':np.zeros(80),
           'XCBI':np.zeros(80), 'YCBI':np.zeros(80), 
           'NXNYCB':np.zeros(80), 'NCBPTS':0, 'IINTCB': 1, 'IDIFCB': 1, 'LECB':0}

    BCC = {'PT': np.zeros(20), 'TT': np.zeros(20),
           'THETA':np.zeros(20), 'PE': 14.7, 'MASSE':0, 'MASSI':0, 'MASST':0,
           'THRUST':0, 'NSTAG': 0, 'ISUPER': 0}
    
    return DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC

def case1_dicts(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC):
    DESC['TITLE'] = 'CASE NO. 1 - CONVERGING-DIVERGING NOZZLE (45 DEG INLET, 15 DEG EXIT)'
    CNTRLC['NMAX'] = 400
    CNTRLC['TCONV'] = 0.003
    CNTRLC['FDT'] = 1.6
    CNTRLC['JFLAG'] = 0
    
    GEMTRYC['NGEOM'] = 2
    GEMTRYC['XI']=0.31
    GEMTRYC['RI']=2.5
    GEMTRYC['RT']=0.8
    GEMTRYC['XE']=4.05
    GEMTRYC['RCI']=0.8
    GEMTRYC['RCT']=0.5
    GEMTRYC['ANGI']=44.88
    GEMTRYC['ANGE']=15.0
    
    BCC['PT'] = np.full(20, 70.0)
    BCC['TT'] = np.full(20, 80.0)
           
    return DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC

def print_title(items):
    print(items['TITLE'])

def print_controls(variables, columns=5, column_width=20):
    print('Control Parameters\n')
    items = [f"{name}={value}" for name, value in variables.items()]
    for i, item in enumerate(items):
        print(f"{item:<{column_width}}", end="")  # Left-align each item with fixed width
        if (i + 1) % columns == 0:  # Check if the row is complete
            print()  # Move to the next line
    if len(items) % columns != 0:  # If items don't perfectly fill the grid, add a newline
        print()    

def print_boundary_conditions(dictionary, array_keys, columns=4, column_width=20):
    print('Boundary Conditions\n')
    # Print array keys as single lines
    for key in array_keys:
        if key in dictionary:
            print(f"{key}=", end="")
            print(" ".join(str(value) for value in dictionary[key]))  # Print the array on a single line

    # Filter out keys with arrays
    non_array_items = {k: v for k, v in dictionary.items() if k not in array_keys}

    # Print the remaining items in a grid layout
    items = [f"{name}={value}" for name, value in non_array_items.items()]
    for i, item in enumerate(items):
        print(f"{item:<{column_width}}", end="")
        if (i + 1) % columns == 0:
            print()  # Move to the next line
    if len(items) % columns != 0:
        print()  # Final newline if needed


def geom(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC):
    if GEMTRYC['NGEOM'] == 1:                   # Constant Area Duct
        print('\nNOZZLE GEOMETRY\n')
        if CNTRLC['IUI'] == 1:
            print(f'CONSTANT AREA DUCT SPECIFIED BY XI={GEMTRYC["XI"]} (IN), RI={GEMTRYC["RI"]} (IN), and XE={GEMTRYC["XE"]} (IN)')
        if CNTRLC['IUI'] == 2:
            print(f'CONSTANT AREA DUCT SPECIFIED BY XI={GEMTRYC["XI"]} (CM), RI={GEMTRYC["RI"]} (CM), and XE={GEMTRYC["XE"]} (CM)')
    
        LT = CNTRLC['LMAX']
        GEMTRYC['DX'] = (GEMTRYC['XE'] - GEMTRYC['XI']) / (CNTRLC['LMAX']-1)
        GEMTRYC['XT'] = GEMTRYC['XE']
        GEMTRYC['RT'] = GEMTRYC['RI']
        GEMTRYC['RE'] = GEMTRYC['RI']

        GEMTRYC['YW'] = np.full(CNTRLC['LMAX'], GEMTRYC['RI'])
        GEMTRYC['NXNY'] = np.full(CNTRLC['LMAX'], 0)
        
        if CNTRLC['JFLAG'] == 1:
            XWL = GEMTRYC['XI'] + (CNTRLC['LJET']-2) * GEMTRYC['DX']
            if CNTRLC['IUI'] == 1:
                print(f'EXHAUST JET CALCULATION HAS BEEN REQUESTED, THE NOZZLE ENDS AT X={XWL} (IN).  THE MESH POINTS L={CNTRLC["LJET"]} to L={CNTRLC["LMAX"]} ARE AN INITIAL APPROXIMATION TO THE EXHAUST JET BOUNDARY.')
            if CNTRLC['IUI'] == 2:
                print(f'EXHAUST JET CALCULATION HAS BEEN REQUESTED, THE NOZZLE ENDS AT X={XWL} (CM)).  THE MESH POINTS L={CNTRLC["LJET"]} to L={CNTRLC["LMAX"]} ARE AN INITIAL APPROXIMATION TO THE EXHAUST JET BOUNDARY.')
 
        if CNTRLC['IUI'] == 2:
            GEMTRYC['YW'] = GEMTRYC['YW']/2.54           
    
    elif GEMTRYC['NGEOM'] == 2:                 # Circular-Arc, Conical Nozzle Case
        print('\nNOZZLE GEOMETRY\n')
        if GEMTRYC['RCI'] == 0 or GEMTRYC['RCT'] == 0:
            print('RCI OR RCT WAS SPECIFIED AS ZERO')
            CNTRLC['IERR'] = 1
            return DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC
        ANI = np.deg2rad(GEMTRYC['ANGI'])
        ANE = np.deg2rad(GEMTRYC['ANGE'])
        XTAN = GEMTRYC['XI'] + GEMTRYC['RCI'] * np.sin(ANI)
        RTAN = GEMTRYC['RI'] + GEMTRYC['RCI'] * (np.cos(ANI)-1.0)
        RT1 = GEMTRYC['RT'] - GEMTRYC['RCT'] * (np.cos(ANI)-1.0)
        XT1 = XTAN + (RTAN - RT1) / np.tan(ANI)
        if XT1 < XTAN:
            XT1 = XTAN
            RT1 = RTAN
        GEMTRYC['XT'] = XT1 + GEMTRYC['RCT'] * np.sin(ANI)
        XT2 = GEMTRYC['XT'] + GEMTRYC['RCT'] * np.sin(ANE)
        RT2 = GEMTRYC['RT'] + GEMTRYC['RCT'] * (1.0 - np.cos(ANE))
        GEMTRYC['RE'] = RT2 + (GEMTRYC['XE'] - XT2) * np.tan(ANE)
        LT = 0
        GEMTRYC['DX'] = (GEMTRYC['XE'] - GEMTRYC['XI']) / (CNTRLC['LMAX']-1)
        if CNTRLC['IUI'] == 1:
            print(f'A CIRCULAR-ARC, CONICAL NOZZLE HAS BEEN SPECIFIED BY XI={GEMTRYC["XI"]} (IN), RI={GEMTRYC["RI"]} (IN), RT={GEMTRYC["RT"]} (IN), XE={GEMTRYC["XE"]} (IN), RCI={GEMTRYC["RCI"]} (IN), RCT={GEMTRYC["RCT"]} (IN), ANGI={GEMTRYC["ANGI"]} (DEG), AND ANGE={GEMTRYC["ANGE"]} (DEG).  THE COMPUTED VALUES ARE XT={GEMTRYC["XT"]} (IN) AND RE={GEMTRYC["RE"]} (IN)')
        elif CNTRLC['IUI'] == 2:
            print(f'A CIRCULAR-ARC, CONICAL NOZZLE HAS BEEN SPECIFIED BY XI={GEMTRYC["XI"]} (CM), RI={GEMTRYC["RI"]} (CM), RT={GEMTRYC["RT"]} (CM), XE={GEMTRYC["XE"]} (CM), RCI={GEMTRYC["RCI"]} (CM), RCT={GEMTRYC["RCT"]} (CM), ANGI={GEMTRYC["ANGI"]} (DEG), AND ANGE={GEMTRYC["ANGE"]} (DEG).  THE COMPUTED VALUES ARE XT={GEMTRYC["XT"]} (CM) AND RE={GEMTRYC["RE"]} (CM)')

        for L in range(0, CNTRLC['LMAX']):
            X = GEMTRYC['XI'] + (L-0)*GEMTRYC['DX']
            GEMTRYC['XW'][L] = GEMTRYC['XI'] + GEMTRYC['DX'] * L
            if X >= GEMTRYC['XI'] and X <= XTAN:
                GEMTRYC['YW'][L] = GEMTRYC['RI']+GEMTRYC['RCI']*(np.cos(np.arcsin((X-GEMTRYC['XI'])/GEMTRYC['RCI']))-1.0)
                GEMTRYC['NXNY'][L] = (X-GEMTRYC['XI'])/(GEMTRYC['YW'][L]-GEMTRYC['RI']+GEMTRYC['RCI'])
            elif X > XTAN and X <= XT1:
                GEMTRYC['YW'][L] = RT1 + (XT1-X)*np.tan(ANI)
                GEMTRYC['NXNY'][L] = np.tan(ANI)
            elif X > XT1 and X <= GEMTRYC['XT']:
                GEMTRYC['YW'][L] = GEMTRYC['RT']+GEMTRYC['RCT']*(1.0-np.cos(np.arcsin((GEMTRYC['XT'] - X)/GEMTRYC['RCT'])))
                GEMTRYC['NXNY'][L] = (GEMTRYC['XT']-X)/(GEMTRYC['RCT']+GEMTRYC['RT'] - GEMTRYC['YW'][L])
            elif X > GEMTRYC['XT'] and X <= XT2:
                GEMTRYC['YW'][L] = GEMTRYC['RT']+GEMTRYC['RCT']*(1.0-np.cos(np.arcsin((X - GEMTRYC['XT'])/GEMTRYC['RCT'])))
                GEMTRYC['NXNY'][L] = (GEMTRYC['XT']-X)/(GEMTRYC['RCT']+GEMTRYC['RT'] - GEMTRYC['YW'][L])
            elif X > XT2 and X <= GEMTRYC['XE']:
                GEMTRYC['YW'][L] = RT2 + (X - XT2) * np.tan(ANE)
                GEMTRYC['NXNY'][L] = -np.tan(ANE)
            
            if L > 0:
                if GEMTRYC['YW'][L] > GEMTRYC['YW'][LT]: LT = L
                
        if CNTRLC['JFLAG'] == 1:
            XWL = GEMTRYC['XI'] + (CNTRLC['LJET']-2) * GEMTRYC['DX']
            if CNTRLC['IUI'] == 1:
                print(f'EXHAUST JET CALCULATION HAS BEEN REQUESTED, THE NOZZLE ENDS AT X={XWL} (IN).  THE MESH POINTS L={CNTRLC["LJET"]} to L={CNTRLC["LMAX"]} ARE AN INITIAL APPROXIMATION TO THE EXHAUST JET BOUNDARY.')
            if CNTRLC['IUI'] == 2:
                print(f'EXHAUST JET CALCULATION HAS BEEN REQUESTED, THE NOZZLE ENDS AT X={XWL} (CM)).  THE MESH POINTS L={CNTRLC["LJET"]} to L={CNTRLC["LMAX"]} ARE AN INITIAL APPROXIMATION TO THE EXHAUST JET BOUNDARY.')

        if CNTRLC['IUI'] == 2:
            GEMTRYC['YW'] = GEMTRYC['YW']/2.54     
            GEMTRYC['XT'] = GEMTRYC['XT']/2.54
            GEMTRYC['RT'] = GEMTRYC['RT']/2.54
            
            if GCB['NGCB'] == 0:
                GEMTRYC['XI'] = GEMTRYC['XI']/2.54
                GEMTRYC['XE'] = GEMTRYC['XE']/2.54
                GEMTRYC['DX'] = GEMTRYC['DX']/2.54
    elif GEMTRYC['NGEOM'] == 3:
        print('Not Implemented')
    elif GEMTRYC['NGEOM'] == 4:
        print('Not Implemented')
    else:
        print('NGEOM out of range (1-4 only)')
   
             
    return DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC
    


def plot_wall_contour(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC):
    plt.figure(figsize=(10, 5))
    plt.plot(GEMTRYC['XW'], GEMTRYC['YW'], label="Outer Wall")
    plt.xlabel('x-axis [m]')
    plt.ylabel('y-axis [m]')
    plt.grid()
    plt.legend()

    return DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC

    

def ONEDIM(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC):
    # Set up before the loop
    MN3 = 0.01
    if CNTRLC['N1D'] == -1 or CNTRLC['N1D'] > 2:    MN3 = 2.0
    GRGAS = 1.0 / (CNTRLC['RGAS'] * CNTRLC['G'])
    NXDK = 0.0
    ACOEF = 2.0/(CNTRLC['GAMMA'] + 1.0)
    BCOEF = (CNTRLC['GAMMA'] - 1.0) / (CNTRLC['GAMMA'] + 1.0)
    CCOEF = (CNTRLC['GAMMA'] + 1.0) / (2.0 * (CNTRLC['GAMMA'] - 1.0))
    
    RSTAR = GEMTRYC['RT']
    RSTARS = GEMTRYC['RT']**2
    
    for L in range(0, CNTRLC['LMAX']):
        X = GEMTRYC['XI'] + (L-0)*GEMTRYC['DX']
        if X < GEMTRYC['XT']:
            RADS = GEMTRYC['YW'][L]**2 - GCB['YCB'][L]**2
            ARATIO = RADS / RSTARS
            
            # Newton Raphson Iteration Loop
            for _ in range(0, 19):
                ABM = ACOEF + BCOEF * MN3**2
                ABMC = ABM**CCOEF
                FM = ABMC / MN3 - ARATIO
                FPM = ABMC * (2.0 * BCOEF * CCOEF/ABM-1.0/MN3**2)
                OMN3 = MN3
                MN3 = OMN3 - FM/FPM
                if MN3 > 1.0 and OMN3 < 1.0:  MN3 = 0.99
                if MN3 < 1.0 and OMN3 > 1.0:  MN3 = 1.01
                if MN3 < 0.0: MN3 = -MN3
                if (np.abs(MN3 - OMN3) / OMN3) <= 0.0005:
                    # Fill in 2-D arrays
                    DEM = 1.0 + CNTRLC['GAM2'] * MN3**2
                    DEMP = DEM**CNTRLC['GAM1']
                    DNXNY = (GEMTRYC['NXNY'][L] - GCB['NXNYCB'][L] / CNTRLC['M1'])
                    for M in range(0, CNTRLC['MMAX']):
                        SOLUTN['P'][L][M][1] = BCC['PT'][M] / DEMP
                        TEMP = BCC['TT'][M] / DEM
                        SOLUTN['RO'][L][M][1] = SOLUTN['P'][L][M][1] * GRGAS / TEMP
                        Q = MN3 * np.sqrt(CNTRLC['GAMMA'] * SOLUTN['P'][L][M][1] / SOLUTN['RO'][L][M][1])
                        DN = GCB['NXNYCB'][L] + DNXNY * (M-1)   
                        DNS = DN**2
                        if DNS == 0:
                            SOLUTN['U'][L][M][1] = Q
                            SOLUTN['V'][L][M][1] = 0.0
                        else:
                            SIGN = 1.0
                            if DN > 0.0: SIGN = -1.0
                            SOLUTN['U'][L][M][1] = Q/np.sqrt(1.0 * DNS)
                            SOLUTN['V'][L][M][1] = SIGN * Q / np.sqrt(1.0 + (1.0 / DNS))

                        
                        
            
            
    return DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC

    





    
if __name__ == '__main__':

    DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC = get_default_dicts()
    (DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC) = case1_dicts(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC)

    print_title(DESC)
    print_controls(CNTRLC)
    print_boundary_conditions(BCC, ['PT', 'TT', 'THETA'], columns=4, column_width=20)

    # Calculates the Nozzle Radius and Outer Normal
    (DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC) = geom(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC)
    
    CNTRLC['DY'] = 1.0 / (CNTRLC['MMAX'] - 1)

    for M in range(1, CNTRLC['MMAX']):
        BCC['PT'][M] = BCC['PT'][0] 
        BCC['TT'][M] = BCC['TT'][0] 
        BCC['THETA'][M] = BCC['THETA'][0]
    

    if CNTRLC['IUI'] == 1:
        print(f'Boundary Conditions - PT={BCC["PT"][0]} (PSI), TT={BCC["TT"][0]} (F), THETA={BCC["THETA"][0]} (DEG), AND PE={BCC["PE"]} (PSI)')
    if CNTRLC['IUI'] == 2:
        print(f'Boundary Conditions - PT={BCC["PT"][0]} (KPA), TT={BCC["TT"][0]} (C), THETA={BCC["THETA"][0]} (DEG), AND PE={BCC["PE"]} (PSI)')

    # Convert input data units to internal units
    CNTRLC['TCONV'] = CNTRLC['TCONV']/100.0
    T = SOLUTN['TSTART'] * CNTRLC['LC']
    CNTRLC['TSTOP'] = CNTRLC['TSTOP'] * CNTRLC['LC'] 
    for L in range(0, CNTRLC['LMAX']):
        GEMTRYC['XWI'] = 0.0
    for M in range(0, CNTRLC['MMAX']):
        BCC['PT'][M] = BCC['PT'][M] * CNTRLC['PC']
        BCC['TT'][M] = BCC['TT'][M] * CNTRLC['TC']
        BCC['THETA'][M] = BCC['THETA'][M] * 0.0174533
    
    BCC['PE'] = BCC['PE'] * CNTRLC['PC']
    CNTRLC['GAM1'] = CNTRLC['GAMMA'] / (CNTRLC['GAMMA'] - 1.0)
    CNTRLC['GAM2'] = (CNTRLC['GAMMA'] - 1.0) / 2.0

    CNTRLC['L1'] = CNTRLC['LMAX'] - 1
    CNTRLC['L2'] = CNTRLC['LMAX'] - 2
    CNTRLC['L3'] = CNTRLC['LMAX'] - 3
    CNTRLC['M1'] = CNTRLC['MMAX'] - 1
    CNTRLC['M1'] = CNTRLC['MMAX'] - 2
    
    
    #Get the 1D inital surface
    (DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC) = ONEDIM(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC)


    #(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC) = plot_wall_contour(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC)
    
     #a = Solution(x_i=0.31, r_i=2.5, r_t=0.8, x_e=4.05, rc_i=0.8, rc_t=0.5, ang_i=44.88, ang_e=15.0, PT=70, TT=80, THETA=0, PE=14.7, t_conv=3e5, t_stop=1, gamma=1.4, R_gas=53.3523, g=32.2)
   
    
    a = Solution()  #l_max=5, x_i=0, x_e=1)    
    a.build_wall_geometry()
    a.one_dimisional()
    
    
    
    
    
    
''' Pick up with the a.one_dimensional() stuff.  I'm frustrated by that function call and I'm going to tak a break from this code for a while.'''    
    
    
    
    
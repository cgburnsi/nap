import numpy as np

class Data:
    def __init__(self):
        # Program Control Block Data Variables        
        self.LMAX       = 21
        self.MMAX       = 8
        self.NMAX       = 400
        self.NPRINT     = 0
        self.NASM       = 1
        self.IVEL       = 0
        self.ICHAR      = 0
        self.N1D        = 1
        self.LJET       = 0
        self.JFLAG      = 0
        self.IERR       = 0
        self.IUI        = 1
        self.IUO        = 1
        self.LD         = 0
        self.MD         = 0
        self.LMD1       = 0
        self.LMD3       = 0
        self.IB         = 0
        self.NPLOT      = -1
        self.TCONV      = 0.003
        self.FDT        = 1.6
        self.GAMMA      = 1.4
        self.RGAS       = 53.35
        self.GAM1       = 0.0
        self.GAM2       = 0.0
        self.RSTAR      = 0.0
        self.RSTARS     = 0.0
        self.DX         = 0.1
        self.DY         = 0.1
        self.DT         = 0.01
        self.G          = 31.174
        self.PC         = 144.0
        self.TC         = 0.0
        self.PLOW       = 0.01
        self.ROLOW      = 0.0001
    
        # Artificial Viscosity Variables
        self.IAV = 0
        self.NST = 0
        self.SMP = 0.95
        self.LSS = 0
        self.CAV = 4.0
        self.CTA = 0.0
        self.XMU = 0.0
        self.XLA = 1.0
        self.RKMU = 0.7
        self.QUT = np.zeros((81, 21))
        self.QVT = np.zeros((81, 21))
        self.QPT = np.zeros((81, 21))

        # Variables in ONESID common block
        self.UD = np.zeros(4)
        self.VD = np.zeros(4)
        self.PD = np.zeros(4)
        self.ROD = np.zeros(4)

        # Variables in SOLUTN common block
        self.U = np.zeros((81, 21, 2))
        self.V = np.zeros((81, 21, 2))
        self.P = np.zeros((81, 21, 2))
        self.RO = np.zeros((81, 21, 2))

        # Variables in GEMTRYC common block
        self.NGEOM = 2
        self.NXNY = np.zeros(81, dtype=int)
        self.NWPTS = 0
        self.IINT = 1
        self.IDIF = 1
        self.LT = 0
        self.NDIM = 1
        self.XI = 0.31
        self.RI = 2.5
        self.XT = 0.0
        self.RT = 0.8
        self.XE = 4.05
        self.RE = 0.0
        self.RCI = 0.8
        self.RCT = 0.5
        self.ANGI = 44.88
        self.ANGE = 15.0
        self.XW = np.zeros(81)
        self.YW = np.zeros(81)
        self.XWI = np.zeros(81)
        self.YWI = np.zeros(81)

        # Variables in GCB common block
        self.NGCB = 0
        self.NXNYCB = np.zeros(81, dtype=int)
        self.NCBPTS = 0
        self.IINTCB = 1
        self.IDIFCB = 1
        self.LECB = 0
        self.XICB = 0.0
        self.RICB = 0.0
        self.XTCB = 0.0
        self.RTCB = 0.0
        self.XECB = 0.0
        self.RECB = 0.0
        self.RCICB = 0.0
        self.RCTCB = 0.0
        self.ANGICB = 0.0
        self.ANGECB = 0.0
        self.XCB = np.zeros(81)
        self.YCB = np.zeros(81)
        self.XCBI = np.zeros(81)
        self.YCBI = np.zeros(81)

        # Variables in BCC common block
        self.PT = 70
        self.TT = 80
        self.THETA = np.zeros(21)                       # The Fortran code calls out specifically 'THETA(1)=0.0.
        self.PE = 14.7
        self.MASSE = 0.0
        self.MASSI = 0.0
        self.MASST = 0.0
        self.THRUST = 0.0
        self.NSTAG = 0
        
        # Additional Variables Not Part of the Fortran IV Common Blocks
        self.TSTART = 0.0
        self.TSTOP = 1.0
        self.NAME = 0
        self.IPUNCH = 0
        self.NSTART = 0
        self.N = 0 
        self.IEX = 1
        self.NCONVI = 1
        self.ISUPER = 0
        self.IUNIT = 0
        self.TITLE = 'Converging-Diverging Nozzle (45 Deg Inlet, 15 Deg Exit)'

        self.printed_param = {'LMAX': self.LMAX, 'MMAX': self.MMAX, 'NMAX': self.NMAX, 'NPRINT': self.NPRINT,
                              'TCONV': self.TCONV, 'FDT': self.FDT, 'NSTAG': self.NSTAG, 'NASM': self.NASM,
                              'IUNIT': self.IUNIT, 'IUI': self.IUI, 'IUO': self.IUO, 'IEX': self.IEX,
                              'NCONVI': self.NCONVI, 'TSTOP': self.TSTOP, 'N1D': self.N1D, 
                              'NPLOT': self.NPLOT, 'IPUNCH': self.IPUNCH, 'ISUPER': self.ISUPER, 
                              'IAV': self.IAV, 'CAV': self.CAV, 'XMU': self.XMU, 'XLA': self.XLA, 
                              'RKMU': self.RKMU, 'CTA': self.CTA, 'LSS': self.LSS, 'SMP': self.SMP,
                              'NST:': self.NST}
        


def print_230_250(XI, RI, XE):
    print(f'Nozzle Geometry -- Constant Area Duct: XI={XI:.2} [in], RI={RI:.2} [in], and XE={XE:.2f} [in]')
    
def print_230_260(XI, RI, XE):
    print(f'Nozzle Geometry -- Constant Area Duct: XI={XI:.2} [cm], RI={RI:.2} [cm], and XE={XE:.2f} [cm]')
 
def print_370(d):
    print(f'Exhaust Jet Calculation Requested, Nozzle End = {d.XWL:.2f} [in].')
    print(f'Exhaust Boundary Estimate: Mesh Points L={d.LJET:} to L={d.LMAX}.')

def print_390():
    print('ERROR: RCI or RCT was specified as zero.')
    
def print_690():                        
    print(f'{"-" * 90}')
    print(f'{"NAP - TWO-DIMENSIONAL, TIME-DEPENENT, INVISID NOZZLE FLOW":^90}')
    print(f'{"-" * 90}')
    print('')

def print_710(title):        
    print(f'Job Title -- {title}')
   
def print_720_730(data):
    columns = 6
    # Adjust the column width to fit the longest key and value
    key_width = max(len(key) for key in data.keys())
    value_width = max(len(f"{value:.3f}") if isinstance(value, float) else len(str(value)) for value in data.values())
    
    print('\nControl Parameters --')
    for i in range(0, len(data), columns):
        row = list(data.items())[i:i + columns]
        print(" ".join([f"{key: <{key_width}}: {value: <{value_width}}" for key, value in row]))

def print_740(gamma, rgas):
    print(f'\nFluid Model -- GAMMA = {gamma:.2f} and, R={rgas:.3f} [ft-lbf/lmb-degR]')

def print_750(gamma, rgas):
    print(f'\nFluid Model -- GAMMA = {gamma:.2f} and, R={rgas:.3f} [J/kg-K]')
 
def print_790():
    print('Flow Geometry -- Two-Dimensional, Planar Flow Specified')

def print_800():
    print('Flow Geometry -- Axisymmetric Flow Specified')
    
def print_270(d):
    print('Nozzle Geometry -- Circular Arc, Conical Nozzle')
    print(f'      XI= {d.XI:.2f} [in], RI= {d.RI:.2f} [in], RT=   {d.RT:0.2f} [in],  XE=  {d.XE:.2f} [in]')
    print(f'      RCI={d.RCI:.2f} [in], RCT={d.RCT:.2f} [in], ANGE={d.ANGI:0.2f} [deg], ANGE={d.ANGE:.2f} [deg]')
    print(f'      XT= {d.XT:.2f} [in], RE= {d.RE:.2f} [in]')
 
def print_280(d):
    print('Nozzle Geometry -- Circular Arc, Conical Nozzle')
    print(f'      XI= {d.XI:.2f} [cm], RI= {d.RI:.2f} [cm], RT=   {d.RT:0.2f} [cm],  XE=  {d.XE:.2f} [cm]')
    print(f'      RCI={d.RCI:.2f} [cm], RCT={d.RCT:.2f} [cm], ANGE={d.ANGI:0.2f} [deg], ANGE={d.ANGE:.2f} [deg]')
    print(f'      XT= {d.XT:.2f} [cm], RE= {d.RE:.2f} [cm]')
    
def print_xx():
    print('Boundary Conditions --')
   
def GEOM(data):
    if data.NGEOM   == 1: GEOM_Constant_Area(data)
    elif data.NGEOM == 2: GEOM_Conical_Nozzle(data)
    elif data.NGEOM == 3: GEOM_General_Nozzle_Coords(data)
    elif data.NGEOM == 4: GEOM_General_Nozzle_Coords_Slopes(data)

def GEOM_Constant_Area(data):
    if   data.IUI == 1: print_230_250(data.XI, data.RI, data.XE)
    elif data.IUI == 2: print_230_260(data.XI, data.RI, data.XE)
    data.LT = data.LMAX
    data.DX = (data.XE-data.XI/(data.LMAX-1))
    data.XT = data.XE
    data.RT = data.RI
    data.RE = data.RI
    data.YW = np.full(data.LMAX, data.RI)
    data.NXNY = np.full(data.LMAX, 0.0)
    
    if data.JFLAG == 0 and data.IUI == 1:
        data.YW = data.YW / 2.54    
        data.XWL = data.XI + (data.LJET-2) * data.DX
    
        if   data.IUI == 1: print_370(data)
        elif data.IUI == 2: print_370(data)
    
    return data
    
def GEOM_Conical_Nozzle(data):
    if d.RCI == 0.0 or d.RCT == 0.0: 
        print_390()
        d.IERR = 1
        return
    d.ANI = d.ANGI * 3.141593/180.0
    d.ANE = d.ANGE * 3.141593/180.0
    d.XTAN  = d.XI + d.RCI * np.sin(d.ANI)
    d.RTAN  = d.RI + d.RCI * (np.cos(d.ANE) - 1.0)
    d.RT1  = d.RT - d.RCT * (np.cos(d.ANI) - 1.0)
    d.XT1  = d.XTAN + (d.RTAN-d.RT1) / np.tan(d.ANI)
        
    if d.XT1 < d.XTAN: d.XT1, d.RT1 = d.XTAN, d.RTAN
    d.XT = d.XT1 + d.RCT * np.sin(d.ANI)
    d.XT2  = d.XT + d.RCT * np.sin(d.ANE)
    d.RT2  = d.RT + d.RCT * (1.0 - np.cos(d.ANE))
    d.RE   = d.RT2 + (d.XE-d.XT2) * np.tan(d.ANE)
    d.LT = 0
    d.DX   = (d.XE -d.XI) / (d.LMAX-1)
    if   d.IUI == 1: print_270(d)
    elif d.IUI == 2: print_280(d)
    for L in range(0, d.LMAX):
        X = d.XI + d.DX * L
        #self.xw[L] = x
        if X >= d.XI and X <= d.XTAN:
            d.YW[L]      = d.RI + d.RCI * (np.cos(np.arcsin((X - d.XI) / d.RCI)) - 1.0)
            d.NXNY[L]    = (X - d.XI) / (d.YW[L] - d.RI + d.RCI)
        elif X > d.XTAN and X <= d.XT1:
            d.YW[L]      = d.RT1 + (d.XT1 - X) * np.tan(d.ANI)
            d.NXNY[L]    = np.tan(d.ANI)
        elif X > d.XT1 and X <= d.XT:
            print('hi')
            d.YW[L]      = d.RT + d.RCT * (1.0 - np.cos(np.arcsin((d.XT - X) / d.RCT)))
            d.NXNY[L]    = (d.XT - X) / (d.RCT + d.RT - d.YW[L])
        elif X > d.XT and X <= d.XT2:
            print('hi hi')

            d.YW[L]      = d.RT + d.RCT * (1.0-np.cos(np.arcsin((X - d.XT) / d.RCT)))
            d.NXNY[L]    = (d.XT - X) / (d.RCT + d.RT - d.YW[L])
        elif X > d.XT2 and X <= d.XE:
            print('hi hi hi ')
            d.YW[L]      = d.RT2 + (X - d.XT2) * np.tan(d.ANE)
            d.NXNY[L]    = -np.tan(d.ANE)
        
        # Determine where the index is for the throat in the yw (outer wall) array
        if L <= 0: 
            continue
        else:
            if d.YW[L] < d.YW[d.LT]: 
                d.LT = L;


def GEOM_General_Nozzle_Coords(d):
    ...
    
def GEOM_General_Nozzle_Coords_Slopes(data):
    ...


if __name__ == "__main__":
    d = Data()          # NAP program Data (Class constructor is set to example case #1 as defaults)

    # Output data about the problem
    print_690() 
    print_710(d.TITLE) 
    NPRIND=np.abs(float(d.NPRINT))
    print_720_730(d.printed_param)
    if   d.IUI == 1: print_740(d.GAMMA, d.RGAS)
    elif d.IUI == 2: print_750(d.GAMMA, d.RGAS)
    if   d.NDIM == 0: print_790()
    elif d.NDIM == 1: print_800()

    # Calculate the Nozzle Radius and Normal
    # Stuff goes to GEOM: NGEOM, XI, RI, XE
    a = GEOM(d)






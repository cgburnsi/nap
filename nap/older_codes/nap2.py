import numpy as np
import matplotlib.pyplot as plt

import onedim as od

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
    CNTRLC['N1D'] = 1
    
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
    plt.xlabel("Time (s)")
    plt.ylabel("Pressure (Torr)")
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
    
    
    (DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC) = od.onedim(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC)
    print(SOLUTN['U'].shape)
    #Get the 1D inital surface
   # (DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC) = ONEDIM(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC)


    #(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC) = plot_wall_contour(DESC, CNTRLC, AV, ONESID, SOLUTN, GEMTRYC, GCB, BCC)
    
    
    
    
    
    
    
    
    
    
    
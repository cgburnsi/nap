import numpy as np
import convert as cv

from geom1 import GEOM
from map2 import MAP
from onedim1 import ONEDIM
from masflo1 import MASFLO
from shock1 import SHOCK
from inter1 import INTER
#from calcivs1 import CALCIVS



def print_title():                        
    print(f'{"-" * 90}')
    print(f'{"NAP - TWO-DIMENSIONAL, TIME-DEPENENT, INVISID NOZZLE FLOW":^90}')
    print(f'{"-" * 90}')
    print('')
    print('Job Title -- C-D Nozzle (45 deg Inlet, 15 deg Exit')






if __name__ == '__main__':
       
    # Program Controls
    LMAX, MMAX, NMAX = 5, 2, 20
    #Simulation Properties
    FDT, TSTOP, IAV = 1.6, 1.0, 0
    # Boundary Conditions
    PT              = np.full(MMAX, cv.convert(70, 'psi', 'Pa'))
    TT              = np.full(MMAX, cv.convert(80, 'degF', 'K'))
    PE              = cv.convert(14.7, 'psi', 'Pa')   
    # Fluid Properties
    GAMMA           = 1.4
    RGAS            = cv.convert(53.353, 'ft*lbf/lb*degR', 'J/kg*K')    
    # Geometry Variables
    XI, RI, RT      = cv.convert(0.31, 'in', 'm'), cv.convert(2.50, 'in', 'm'), cv.convert(0.80, 'in', 'm') 
    XE, RCI, RCT    = cv.convert(4.05, 'in', 'm'), cv.convert(0.80, 'in', 'm'), cv.convert(0.50, 'in', 'm') 
    ANI, ANE        = np.deg2rad(44.88), np.deg2rad(15.0)
    
    
    # Calculate the wall geometry
    XT, RE, LT, YW, NXNY = GEOM(LMAX, XI, RI, RT, XE, RCI, RCT, ANI, ANE)
    # Calculate One Dimensional Surface
    U, V, P, RO, T, MN = ONEDIM(LMAX, MMAX, RGAS, GAMMA, RT, XT, XI, XE, YW, PT, TT)
    # Calculate mass flowrate and thrust
    MASSI, THRUST = MASFLO(LMAX, RT, YW, LT, U, V, RO)
    







    # Here are things that need to be set up before I start moving the following to a separate function 
    XWI = np.zeros(LMAX)
    
    # Calculate the Initial-Value Surface
    DX  = (XE - XI) / (LMAX - 1)
    DY  = 1.0 / (MMAX - 1)
    X   = XI - DX
    LD1 = 1 # I don't know what this is.  I can't find any 'LD1' in the Fortran IV code.


    # Calculate the Initial Value Surface
    #CALCIVS(LMAX, MMAX, X, DX, DY, XWI, YW, NXNY, U, V, P, RO, GAMMA, RGAS, LD1)
     
    # These are 
    L1=LMAX-1
    L2=LMAX-2
    L3=LMAX-3
    M1=MMAX-1
    M2=MMAX-1       
    
    # Initialize the Time Step Integration Loop Parameters
    N1, N3, DQM, NS, NCONV, NC, LDUM, NPC, NPD = 0, 1, 0, 0, 0, 0, 0, 0, 0
    DXR, DYR = 1.0 / DX, 1.0 / DY
    DXRS, DYRS = DXR**2, DYR**2
    LD, MD = LMAX, MMAX
    LMD = LD * MD
    NSTART = 0 
    T = 0
    IVEL = 0
    NDIM = 1

    # Adjust LDUM based on Throat Location)    
    if LT != 1: LDUM = LT - 1


    # Enter the time step integration loop
    for N in range(0, NMAX):
        NPD = NPD + 1
        # Print every 10 steps
        if NPD == 10:
            NP = N + NSTART
            print(f"Time step: {NP}")
            NPD = 0
    
        # Update memory offsets
        LMD1 = LMD * N1
        LMD3 = LMD * N3
        
        # Calculate delta T
        UPAM = 0                     # initialize the max combined speed of sound    
        for L in range(LMAX):
            BE, *_ = MAP(0, L, MD, XWI, YW, NXNY, DY, LD1)

            # Calculate DXDY
            DXDY = DXRS + BE**2 * DYRS
        
            for M in range(MMAX):
                # Calculate 1D index for the current time step
                LMN1 = L + LD * (M) + LMD1  # Zero-based indexing

                # Velocity magnitude squared
                QS=U(LMN1)**2 + V(LMN1)**2
                
                QS = U[0, M, L]**2 + V[0, M, L]**2
                print(L, M, LMN1, U[0, M, L], V[0, M, L], QS)
        
                # Speed of sound squared
                AS = GAMMA * P[0, M, L] / RO[0, M, L]
        
                # Combined speed
                UPA = np.sqrt(QS * DXDY) + np.sqrt(AS * DXDY)
        
                # Track maximum combined speed
                if UPAM == 0 or UPA > UPAM:
                    UPAM = UPA
        
    # Calculate the time step
    DT = FDT / UPAM
    
    T += DT
    
    # Adjust time if it exceeds TSTOP
    if T > TSTOP:
        T -= DT
        DT = TSTOP - T
        T = TSTOP  
        
    # Determine if the Exit Flow Is Subsonic or Supersonic
    if QS >= AS: IVEL=1

    # Calculate the Nozzle Wall and Interior Mesh Points
    if IAV != 0: SHOCK()
    ICHAR, IB = 1, 1
    # Calculate Interior Points
    #INTER(LMAX, ICHAR, M1, XWI, YW, NXNY, DY, LD1, LD, LMD1, LMD3, U, V, P, RO, GAMMA, DXR, DYR, N1, NDIM)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
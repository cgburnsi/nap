import numpy as np
import convert as cv



def GEOM2(LMAX, XI, RI, RT, XE, RCI, RCT, ANI, ANE):
    
    YW, NXNY = np.zeros(LMAX), np.zeros(LMAX)
    XT, RE, LT = 0, 0, 0
    
    XTAN = XI + RCI * np.sin(ANI)
    RTAN = RI + RCI * (np.cos(ANI) - 1.0)
    RT1  = RT - RCT * (np.cos(ANI) - 1.0)
    XT1  = XTAN + (RTAN - RT1) / np.tan(ANI)
    if XT1 < XTAN:
        XT1, RT1 = XTAN, RTAN
    XT   = XT1 + RCT * np.sin(ANI)
    XT2  = XT + RCT * np.sin(ANE)
    RT2  = RT + RCT * (1.0 - np.cos(ANE))
    RE   = RT2 + (XE - XT2) * np.tan(ANE)
    DX   = (XE - XI) / (LMAX - 1)

    for L in range(0, LMAX):
        X = XI + DX * L
        if X >= XI and X <= XTAN:
            YW[L]   = RI + RCI * (np.cos(np.arcsin((X - XI) / RCI)) - 1.0)
            NXNY[L] = (X - XI) / (YW[L] - RI + RCI)
        elif X > XTAN and X <= XT1:
            YW[L]   = RT1 + (XT1 - X) * np.tan(ANI)
            NXNY[L] = np.tan(ANI)
        elif X > XT1 and X <= XT:
            YW[L]   = RT + RCT * (1.0 - np.cos(np.arcsin((XT - X) / RCT)))
            NXNY[L] = (XT - X) / (RCT + RT - YW[L])
        elif X > XT and X <= XT2:
            YW[L]   = RT + RCT * (1.0-np.cos(np.arcsin((X - XT) / RCT)))
            NXNY[L] = (XT - X) / (RCT + RT - YW[L])
        elif X > XT2 and X <= XE:
            YW[L]   = RT2 + (X - XT2) * np.tan(ANE)
            NXNY[L] = -np.tan(ANE)
        
        # Determine where the index is for the throat in the yw (outer wall) array
        if L <= 0: 
            continue
        else:
            if YW[L] < YW[LT]: 
                LT = L
                
    return XT, RE, LT, YW, NXNY


if __name__ == '__main__':
   
    # Program Controls
    LMAX = 21

    # Geometry Variables
    XI, RI, RT   = cv.convert(0.31, 'in', 'm'), cv.convert(2.50, 'in', 'm'), cv.convert(0.80, 'in', 'm') 
    XE, RCI, RCT = cv.convert(4.05, 'in', 'm'), cv.convert(0.80, 'in', 'm'), cv.convert(0.50, 'in', 'm') 
    ANI, ANE = np.deg2rad(44.88), np.deg2rad(15.0)
    
    # Calculate the wall geometry
    XT, RE, LT, YW, NXNY = GEOM2(LMAX, XI, RI, RT, XE, RCI, RCT, ANI, ANE)
    
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
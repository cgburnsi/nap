import numpy as np
from map1 import MAP



def print_header():
    print(f"{'L':>3} {'M':>3} {'XP':>7} {'YP':>7} {'UP':>7} {'VP':>7} {'PRES':>10} {'RHO':>7} {'VELMAG':>7} {'XMACH':>7} {'TEMP':>7}")

def print_row(L, M, XP, YP, UP, VP, PRES, RHO, VELMAG, XMACH, TEMP):
    print(f"{L:3} {M:3} {XP:7.3f} {YP:7.3f} {UP:7.3f} {VP:7.3f} {PRES:10.3f} {RHO:7.3f} {VELMAG:7.3f} {XMACH:7.3f} {TEMP:7.3f}")



def CALCIVS(LMAX, MMAX, X, DX, DY, XWI, YW, NXNY, U, V, P, RO, GAMMA, RGAS, LD1):
    print_header()
    
    for L in range(0, LMAX):
        X       = X + DX
        BE      =  MAP(0, L, 0, XWI, YW, NXNY, DY, LD1)
        DYIO    = DY / BE
        Y       = -DYIO
        for M in range(0, MMAX):
            Y = Y+DYIO
            VELMAG = np.sqrt(U[0,M,L]**2 + V[0,M,L]**2)
            XMACH = VELMAG / np.sqrt(GAMMA * P[0,M,L] / RO[0,M,L])
            PRES = P[0,M,L]
            RHO = RO[0,M,L]
            TEMP = P[0,M,L]/(RGAS * RHO)
            XP, YP = X, Y
            UP, VP = U[0,M,L], V[0,M,L]
            
        print_row(L, M, XP, YP, UP, VP, PRES, RHO, VELMAG, XMACH, TEMP)
    
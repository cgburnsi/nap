import numpy as np

def ONEDIM(LMAX, MMAX, RGAS, GAMMA, RT, XT, XI, XE, YW, PT, TT):
    # Calculation helper values
    GRGAS = 1.0 / (RGAS)
    GAM1  = GAMMA / (GAMMA -1)
    GAM2  = (GAMMA - 1) / 2.0 
    DX    = (XE - XI) / (LMAX - 1)
    AT    = np.pi * RT**2 
    # Solution Surfaces
    U     = np.zeros((1, MMAX, LMAX))
    V     = np.zeros((1, MMAX, LMAX))
    P     = np.zeros((1, MMAX, LMAX))
    RO    = np.zeros((1, MMAX, LMAX))
    T     = np.zeros((1, MMAX, LMAX))
    MN    = np.zeros((1, MMAX, LMAX))
    # Calculate the local One-Dimensional Mach Number
    for L in range(0, LMAX):
        X = XI + DX * L
        if X < XT:                                      # Subsonic region of nozzle upstream of throat
            CI  = np.pi * YW[L]**2                      # [m**2] Cros Sectional Area at Current Location
            CIT = CI / AT                               # [-] Area Ratio at current location
            MN3 = Newton_Raphson(0.1, CIT, GAMMA)       # [-] Solve for the local Mach number (subsonic Mach Guess)
        elif X >= XT:
            CE  = np.pi * YW[L]**2                      # [m**2] Area at current location
            CET = CE / AT                               # [-] Area Ratio at current location
            MN3 = Newton_Raphson(5.0, CET, GAMMA)       # [-] Solve for the local Mach number (supersonic Mach Guess)    
    # Calculate P, RO, U, and V for the initial surface from the local Mach number             
        for M in range(0, MMAX):
            DEM = 1.0 + GAM2 * MN3**2
            DEMP = DEM**GAM1
            MN[0, M, L] = MN3
            P[0, M, L]    = PT[M] / DEMP
            RO[0, M, L]   = P[0, M, L] * GRGAS / (TT[M] / DEM)
            T[0, M, L]    = P[0, M, L] / (RO[0, M, L] * RGAS);
            U[0, M, L]    = MN3 * np.sqrt(GAMMA * P[0, M, L] / RO[0, M, L])
            V[0, M, L]    = 0.0

        #print(f'L = {L:<2}, x = {x:<3.2f}, U = {U[0, 0, L]:<3.2f}, V = {V[0, 0, L]:<3.2f}, P = {P[0, 0, L]:<3.2f}, RO = {RO[0, 0, L]:<3.4f} Area Ratio = {A_ratio:<3.5f}, Q = {Q:<3.2f}, Mach {MN3:<3.2f}')

    return U, V, P, RO, T, MN

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

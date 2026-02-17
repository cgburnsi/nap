import numpy as np

def MASFLO(LMAX, RT, YW, LT, U, V, RO):
    # Helper Calculations
    LDUM = LMAX-1
    
    # Calculate the mass flow and thrust for the 1-D Initial Data Surface
    AREAI    = np.pi * YW[0]**2
    AREAT    = np.pi * RT**2 
    #a_throat    = (np.pi * self.yw[self.l_t]**2)
    AREAE    = np.pi * YW[LDUM]**2
    VMI      = np.sqrt(U[0, 0, 0]**2 + V[0, 0, 0]**2)
    VMT      = np.sqrt(U[0, 0, LT]**2 + V[0, 0, LT]**2)
    VME      = np.sqrt(U[0, 0, LDUM]**2 + V[0, 0, LDUM]**2)
    MASSI    = RO[0, 0, 0] * VMI * AREAI
    #mass_t      = self.RO[0, 0, self.l_t] * vm_t * a_throat 
    #mass_e      = self.RO[0, 0, LDUM] * vm_e * a_exit 
    THRUST   = RO[0, 0, LDUM] * U[0, 0, LDUM]**2 * AREAE
    
    #print(f'vm_e = {vm_e:<3.2f}, mass_i = {mass_i:<3.4f}, mass_t = {mass_t:<3.4f}, mass_e = {mass_e:<3.4f}, thrust = {thrust:<3.2f}')
    #print(f'Expansion Ratio = {AREAE/AREAT:<3.2f} [-], vel_exit = {VME:<3.2f} [m/s], mdot = {cv.convert(MASSI, "kg/s", "g/s"):<3.5f} [g/s], thrust = {cv.convert(thrust, "N", "mN"):<3.2f} [mN]')
    return MASSI, THRUST

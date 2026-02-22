import numpy as np
from mapping import map_coordinates

def inlet_boundary(state, params, geometry):
    ''' Calculates the boundary mesh points at the nozzle inlet for subsonic flow.
    
    state:    Dictionary or object containing arrays for U, V, P, D (density)
    params:   Control parameters (gamma, rgas, dt, dx, dy, etc.)
    geometry: Geometry data (ycb, mapping coefficients)
    '''
    
    # Step 1 - Function organization/set up
    mmax = params['mmax']
    
    # Step 2 - Characterization Passes (ICHAR loop) Run Two Times
    for ichar in range(1,3):
        # Step 2a - Radial Mesh Point Loop (Iterate through grip point (M) at the inlet plane.)
        for m in range(mmax):
            # Step 3 - Geometry Mapping and Indicies
            al, be, de, al1, be1, de1 = map_coordinates(2, 0, m, geometry, params)
            
            
            # Step 4 - Property Interpolation at Characteristic Origin
            # Step 5 - Calculate Cross-Derivative Source Terms (psi terms)
            # Step 6 - Solve Compatibility Equations
            # Step 7 - Stagnation Condition Iteration
            for iter_count in range(20):
                # Step 7a - Update Pressure and Temperature based on Local Mach
                # Step 7b - Calculate new velocity (u) using compatibility relations
                # Step 7c - Check Convergency
                pass
            
            # Step 8 - Final Density Calculation          






    return state













if __name__ == '__main__':
    
    params = {'mmax': 8}
    state = {}
    geometry = {}
    
    
    inlet = inlet_boundary(state, params, geometry)



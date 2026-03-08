''' src/dtypes/solver_dtypes.py 
    ------------------
    This module defines the NumPy data types (np.dtype) for organizing the NAP code:

    - Fluid States
        1: U - Conserved Variables (State Vector)
        2: W - Physical Variables (Primitive Vector)
         
    Key Concepts:
    - Numba-compatible: Uses fixed-width floats (f8) to enable nopython mode acceleration.
    - The data structure layout for this code uses the AoS (Array of Structures)
      approach to organizing the data.  It's not the most efficient for numba, but
      is useful as a code organizational strategy.
'''
import numpy as np

# CONSERVATIVE STATE (U): The variables actually solved by the Euler equations.
# Used directly in the HLLC Riemann solver loop.
CONSERVATIVE_DTYPE = np.dtype([
    ('rho',   'f8'),    # Mass density
    ('rhou',  'f8'),    # Axial momentum
    ('rhov',  'f8'),    # Radial momentum
    ('rhoE',  'f8')     # Total energy density
])


# PRIMITIVE STATE (W): The variables we usually think about physically.
# Used for Boundary Conditions, source terms (e.g., gravity/reactions), and post-processing.
# PRIMITIVE STATE: Includes T for ease of use in BCs and Viz
PRIMITIVE_DTYPE = np.dtype([
    ('rho', 'f8'),      # Density
    ('u',   'f8'),      # Axial Velocity
    ('v',   'f8'),      # Radial Velocity
    ('p',   'f8'),      # Pressure
    ('a',   'f8'),      # Speed of Sound
    ('T',   'f8')       # Temperature
])




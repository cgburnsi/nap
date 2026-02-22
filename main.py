


if __name__ == '__main__':
    
    pass

    # Step 1 - Initialization (MAIN, GEOM, ONEDIM)
    # Step 1a - Input Data (grid size, gas constants, boundary conditions)
    # Step 1b - Geometry Setup (Define nozzle wall)
    # Step 1c - Initial Values (Initial flow conditions using 1D isentropic equations)
    # Step 1d - Mapping (Calculate transformation functions for grid)
    
    # Step 2 - Time Integration Loop (MAIN)
    # Step 2a - Calculate Time Step (CFL condition)
    # Step 2b - Predictor Step (INTER - first MacCormack Step)
    # Step 2c - Boundary Conditions (INLET, WALL, EXITT)
    # Step 2d - Corrector Step (INTER - second MacCormack Step)
    # Step 2e - ConvergencyTCONV (check delu/u)
    
    # Step 3 - Plotting and Reporting
    
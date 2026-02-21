import numpy as np
from dataclasses import dataclass

@dataclass
class Solution:
    d: np.ndarray       # [lmax, mmax] Density Array
    u: np.ndarray       # [lmax, mmax] Velocity (x-axis) Array
    v: np.ndarray       # [lmax, mmax] Velocity (y-axis) Array
    p: np.ndarray       # [lmax, mmax] Pressure Array

@dataclass
class History:
    d: np.ndarray       # [n_steps+1, lmax, mmax] Density Array 
    u: np.ndarray       # [n_steps+1, lmax, mmax] Velocity (x-axis) Array
    v: np.ndarray       # [n_steps+1, lmax, mmax] Velocity (y-axis) Array
    p: np.ndarray       # [n_steps+1, lmax, mmax] Pressure Array
    
    res: np.ndarray     # [n_steps,] Residual for a given solution array
    dt:  np.ndarray     # [n_steps,] Time Step for a given solution array



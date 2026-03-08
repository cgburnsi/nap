"""
main.py  —  1D Finite Volume Method Solver
==========================================
Everything in one file. Read top to bottom and it tells the whole story.

Physics:   1D Euler equations (inviscid compressible flow)
Scheme:    Godunov / MUSCL-Hancock with HLLC Riemann solver
Time:      Explicit 3-stage Runge-Kutta (SSP-RK3)
Test case: Sod shock tube (analytical solution included for verification)

Run in Spyder with F5. All results land in the variable explorer.

Sections:
  1.  Constants and parameters          — edit this to change the case
  2.  Equation of state                 — calorically perfect gas
  3.  Conserved <-> primitive variables — packing and unpacking U
  4.  Physical flux                     — F(U) for 1D Euler
  5.  HLLC Riemann solver               — numerical flux at interfaces
  6.  Slope limiters                    — for MUSCL reconstruction
  7.  MUSCL reconstruction              — second-order face states
  8.  Boundary conditions               — ghost cells at domain edges
  9.  Residual assembly                 — sweeps all faces, builds dU/dt
  10. SSP-RK3 time integrator           — advances U by one time step
  11. CFL condition                     — computes stable dt
  12. Initial conditions                — Sod shock tube setup
  13. Analytical solution               — exact Sod solution for comparison
  14. Main solver loop                  — runs the simulation
  15. Post-processing                   — plots results
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# 1. CONSTANTS AND PARAMETERS
#    Edit this section to change the case. Nothing else needs to change.
# =============================================================================

# -- Gas properties --
GAMMA = 1.4          # ratio of specific heats (1.4 for diatomic air)
R_GAS = 287.058      # specific gas constant [J/(kg·K)]

# -- Domain --
X_START = 0.0        # left boundary [m]
X_END   = 1.0        # right boundary [m]
NCELLS  = 400        # number of cells (more = sharper features, slower)

# -- Sod shock tube initial conditions --
# Left state  [rho, u, p]
W_LEFT  = np.array([1.000, 0.0, 1.0])
# Right state [rho, u, p]
W_RIGHT = np.array([0.125, 0.0, 0.1])
# Diaphragm location (fraction of domain)
DIAPHRAGM = 0.5

# -- Numerics --
CFL         = 0.5          # Courant number (keep < 1 for stability)
T_END       = 0.20         # simulation end time [s]
SCHEME      = 'muscl'      # 'first_order' or 'muscl'
LIMITER     = 'van_leer'   # 'minmod', 'van_leer', 'superbee' (only used if muscl)

# -- Boundary conditions --
# 'transmissive' — zero-gradient (good for shock tube, waves leave domain)
# 'reflective'   — mirror velocity (solid wall)
BC_LEFT  = 'transmissive'
BC_RIGHT = 'transmissive'

# -- Output --
PLOT_ANALYTICAL = True      # overlay exact Sod solution
PRINT_INTERVAL  = 100       # print progress every N steps


# =============================================================================
# 2. EQUATION OF STATE
#    Calorically perfect gas: constant gamma, ideal gas law.
#    All thermodynamic relations live here.
# =============================================================================

def pressure_from_internal_energy(rho, e):
    """p = (gamma - 1) * rho * e"""
    return (GAMMA - 1.0) * rho * e

def sound_speed(rho, p):
    """a = sqrt(gamma * p / rho)"""
    return np.sqrt(GAMMA * p / rho)

def temperature(rho, p):
    """T = p / (rho * R)"""
    return p / (rho * R_GAS)

def total_energy(rho, u, p):
    """
    Total energy per unit volume.
    E = rho * (cv * T  +  0.5 * u^2)
      = p / (gamma - 1)  +  0.5 * rho * u^2
    """
    return p / (GAMMA - 1.0) + 0.5 * rho * u**2


# =============================================================================
# 3. CONSERVED <-> PRIMITIVE VARIABLE CONVERSIONS
#
#    Conserved variables U = [rho, rho*u, E]
#    Primitive variables W = [rho, u, p]
#
#    The solver marches U forward in time (conservation law form).
#    Primitive variables are used for the Riemann solver, reconstruction,
#    and output — they're easier to think about physically.
# =============================================================================

def prim_to_cons(W):
    """
    W = [rho, u, p]  ->  U = [rho, rho*u, E]

    Accepts shape (3,) for a single cell or (ncells, 3) for all cells.
    Returns same leading shape.
    """
    batch = (W.ndim == 2)
    if not batch:
        W = W[np.newaxis, :]

    rho = W[:, 0]
    u   = W[:, 1]
    p   = W[:, 2]
    E   = total_energy(rho, u, p)

    U = np.column_stack([rho, rho * u, E])
    return U if batch else U[0]


def cons_to_prim(U):
    """
    U = [rho, rho*u, E]  ->  W = [rho, u, p]

    Accepts shape (3,) or (ncells, 3).
    """
    batch = (U.ndim == 2)
    if not batch:
        U = U[np.newaxis, :]

    rho  = U[:, 0]
    rhou = U[:, 1]
    E    = U[:, 2]

    u = rhou / rho
    # Pressure from total energy: p = (gamma-1) * (E - 0.5 * rho * u^2)
    p = (GAMMA - 1.0) * (E - 0.5 * rhou**2 / rho)

    W = np.column_stack([rho, u, p])
    return W if batch else W[0]


def check_physical(U, location=""):
    """
    Sanity check: density and pressure must be positive.
    Negative values indicate a blown-up simulation.
    """
    W = cons_to_prim(U)
    if np.any(W[:, 0] <= 0) or np.any(W[:, 2] <= 0):
        raise ValueError(
            f"Non-physical state detected {location}: "
            f"min(rho)={W[:,0].min():.3e}, min(p)={W[:,2].min():.3e}"
        )


# =============================================================================
# 4. PHYSICAL FLUX
#    F(U) = [rho*u,  rho*u^2 + p,  u*(E+p)]
#
#    This is what appears in the conservation law:
#        dU/dt + dF/dx = 0
#
#    The Riemann solver produces a *numerical* flux F* at each interface
#    that approximates F(U) in a way that respects wave propagation.
# =============================================================================

def euler_flux(U):
    """
    1D Euler physical flux.
    U shape: (3,) or (n, 3)
    F shape: same as U
    """
    batch = (U.ndim == 2)
    if not batch:
        U = U[np.newaxis, :]

    rho  = U[:, 0]
    rhou = U[:, 1]
    E    = U[:, 2]
    u    = rhou / rho
    p    = (GAMMA - 1.0) * (E - 0.5 * rhou**2 / rho)

    F = np.column_stack([
        rhou,            # mass flux
        rhou * u + p,    # momentum flux (dynamic pressure + static pressure)
        u * (E + p)      # energy flux
    ])
    return F if batch else F[0]


def max_wave_speed(U):
    """
    Maximum signal speed across all cells: max(|u| + a).
    Used for the CFL condition.
    U shape: (ncells, 3)
    """
    W    = cons_to_prim(U)
    rho  = W[:, 0]
    u    = W[:, 1]
    p    = W[:, 2]
    a    = sound_speed(rho, p)
    return float(np.max(np.abs(u) + a))


# =============================================================================
# 5. HLLC RIEMANN SOLVER
#
#    Given left and right states at an interface, HLLC computes the
#    numerical flux F* that correctly propagates all three wave families:
#      - Left-running acoustic wave  (speed S_L)
#      - Contact/entropy wave         (speed S*)
#      - Right-running acoustic wave (speed S_R)
#
#    In 1D there's no rotation needed — everything is already aligned.
# =============================================================================

def hllc_flux(U_L, U_R):
    """
    HLLC numerical flux for a single interface.

    Parameters
    ----------
    U_L, U_R : shape (3,) — conserved states left and right of interface

    Returns
    -------
    F : shape (3,) — numerical flux through the interface
    """
    # Unpack primitive variables
    rho_L, u_L, p_L = cons_to_prim(U_L)
    rho_R, u_R, p_R = cons_to_prim(U_R)
    a_L = sound_speed(rho_L, p_L)
    a_R = sound_speed(rho_R, p_R)
    E_L = U_L[2]
    E_R = U_R[2]

    # -- Wave speed estimates (pressure-based, Batten et al.) ---------------
    # These are more robust near sonic points than simple min/max estimates.
    # The p_star estimate comes from the linearized Riemann problem.
    rho_bar = 0.5 * (rho_L + rho_R)
    a_bar   = 0.5 * (a_L   + a_R)
    p_star  = max(0.5*(p_L + p_R) - 0.5*(u_R - u_L)*rho_bar*a_bar, 0.0)

    # q factors: = 1 in rarefactions, > 1 across shocks
    q_L = 1.0 if p_star <= p_L else np.sqrt(1 + (GAMMA+1)/(2*GAMMA) * (p_star/p_L - 1))
    q_R = 1.0 if p_star <= p_R else np.sqrt(1 + (GAMMA+1)/(2*GAMMA) * (p_star/p_R - 1))

    S_L = u_L - a_L * q_L   # left wave speed (moves left, so usually negative)
    S_R = u_R + a_R * q_R   # right wave speed (moves right, so usually positive)

    # -- Pure upwind: entire Riemann fan on one side -------------------------
    if S_L >= 0.0:
        return euler_flux(U_L)    # supersonic flow to the right
    if S_R <= 0.0:
        return euler_flux(U_R)    # supersonic flow to the left

    # -- Contact wave speed S* ----------------------------------------------
    # S* is the speed of the middle wave (contact discontinuity / entropy wave).
    # Derived by requiring pressure and normal velocity to be equal
    # across the contact wave.
    num   = p_R - p_L + rho_L*u_L*(S_L - u_L) - rho_R*u_R*(S_R - u_R)
    denom = rho_L*(S_L - u_L) - rho_R*(S_R - u_R)
    S_star = num / (denom + 1e-300)   # 1e-300 guards against degenerate cases

    # -- HLLC star states ---------------------------------------------------
    # The star state is the constant state between S_L (or S_R) and S*.
    # It satisfies the Rankine-Hugoniot conditions across the bounding wave.
    def star_state(U, rho, u, E, p, S):
        """
        Compute U* (the state between the outer wave S and the contact S*).
        From Toro (1994) eq. 10.73.
        """
        factor = rho * (S - u) / (S - S_star + 1e-300)
        e_star = E/rho + (S_star - u) * (S_star + p / (rho*(S - u + 1e-300)))
        return np.array([factor,
                         factor * S_star,
                         factor * e_star])

    if S_star >= 0.0:
        # Interface is between S_L and S*: use left star state
        U_star_L = star_state(U_L, rho_L, u_L, E_L, p_L, S_L)
        F_L      = euler_flux(U_L)
        # F* = F_L + S_L * (U*_L - U_L)   [from Rankine-Hugoniot]
        return F_L + S_L * (U_star_L - U_L)
    else:
        # Interface is between S* and S_R: use right star state
        U_star_R = star_state(U_R, rho_R, u_R, E_R, p_R, S_R)
        F_R      = euler_flux(U_R)
        return F_R + S_R * (U_star_R - U_R)


def hllc_flux_batch(U_L, U_R):
    """
    Vectorized HLLC over all interfaces simultaneously.
    Much faster than looping over hllc_flux() for large meshes.

    U_L, U_R : shape (nfaces, 3)
    Returns  : shape (nfaces, 3)
    """
    rho_L = U_L[:,0]; rhou_L = U_L[:,1]; E_L = U_L[:,2]
    rho_R = U_R[:,0]; rhou_R = U_R[:,1]; E_R = U_R[:,2]

    u_L = rhou_L / rho_L
    u_R = rhou_R / rho_R
    p_L = (GAMMA-1) * (E_L - 0.5*rhou_L**2/rho_L)
    p_R = (GAMMA-1) * (E_R - 0.5*rhou_R**2/rho_R)
    a_L = np.sqrt(GAMMA * p_L / rho_L)
    a_R = np.sqrt(GAMMA * p_R / rho_R)

    # Wave speeds
    rho_bar = 0.5*(rho_L + rho_R)
    a_bar   = 0.5*(a_L   + a_R)
    p_star  = np.maximum(0.5*(p_L+p_R) - 0.5*(u_R-u_L)*rho_bar*a_bar, 0.0)

    q_L = np.where(p_star <= p_L, 1.0,
                   np.sqrt(1 + (GAMMA+1)/(2*GAMMA)*(p_star/p_L - 1)))
    q_R = np.where(p_star <= p_R, 1.0,
                   np.sqrt(1 + (GAMMA+1)/(2*GAMMA)*(p_star/p_R - 1)))

    S_L = u_L - a_L*q_L
    S_R = u_R + a_R*q_R

    num    = p_R - p_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R)
    denom  = rho_L*(S_L-u_L) - rho_R*(S_R-u_R)
    S_star = num / (denom + 1e-300)

    # Star states (vectorized)
    def vstar(rho, u, E, p, S):
        f      = rho*(S-u) / (S-S_star+1e-300)
        e_star = E/rho + (S_star-u)*(S_star + p/(rho*(S-u+1e-300)))
        return np.column_stack([f, f*S_star, f*e_star])

    U_starL = vstar(rho_L, u_L, E_L, p_L, S_L)
    U_starR = vstar(rho_R, u_R, E_R, p_R, S_R)

    # Physical fluxes
    def flux_batch(U, rho, rhou, E, u, p):
        return np.column_stack([rhou, rhou*u+p, u*(E+p)])

    F_L = flux_batch(U_L, rho_L, rhou_L, E_L, u_L, p_L)
    F_R = flux_batch(U_R, rho_R, rhou_R, E_R, u_R, p_R)

    F_starL = F_L + S_L[:,None] * (U_starL - U_L)
    F_starR = F_R + S_R[:,None] * (U_starR - U_R)

    # Select flux based on wave speed signs
    F = np.where(S_L[:,None] >= 0,    F_L,
        np.where(S_star[:,None] >= 0,  F_starL,
        np.where(S_R[:,None] >= 0,     F_starR,
                                        F_R)))
    return F


# =============================================================================
# 6. SLOPE LIMITERS
#
#    Slope limiters prevent MUSCL reconstruction from creating new extrema
#    (which would cause oscillations near shocks). Each limiter is a
#    different tradeoff between accuracy in smooth regions and sharpness
#    at discontinuities.
#
#    All limiters take the ratio r = dW_left / dW_right and return phi(r).
#    The limited slope is: dW_limited = phi(r) * dW_right
# =============================================================================

def minmod(r):
    """
    Most diffusive limiter. Very stable, first-order near shocks.
    phi = max(0, min(1, r))
    """
    return np.maximum(0.0, np.minimum(1.0, r))


def van_leer(r):
    """
    Smooth, second-order in smooth regions, first-order at extrema.
    Good general-purpose choice.
    phi = (r + |r|) / (1 + |r|)
    """
    return (r + np.abs(r)) / (1.0 + np.abs(r) + 1e-300)


def superbee(r):
    """
    Least diffusive standard limiter. Sharpest shocks, most oscillation-prone.
    phi = max(0, min(2r, 1), min(r, 2))
    """
    return np.maximum(0.0,
           np.maximum(np.minimum(2.0*r, 1.0),
                      np.minimum(r, 2.0)))


# Map limiter name string to function
LIMITERS = {
    'minmod':   minmod,
    'van_leer': van_leer,
    'superbee': superbee,
}


# =============================================================================
# 7. MUSCL RECONSTRUCTION
#
#    MUSCL (Monotone Upstream-Centered Schemes for Conservation Laws)
#    upgrades from first-order (piecewise constant) to second-order
#    (piecewise linear) accuracy by extrapolating cell-centered values
#    to face values using limited slopes.
#
#    For cell i with neighbors i-1 and i+1:
#      slope     = limiter(dW_left, dW_right)
#      W_face_R  = W_i + 0.5 * slope    (right face of cell i)
#      W_face_L  = W_i - 0.5 * slope    (left face of cell i)
#
#    We reconstruct in PRIMITIVE variables (rho, u, p) for better
#    behavior near shocks, then convert back to conserved for the
#    Riemann solver.
# =============================================================================

def muscl_reconstruct(U, limiter_fn):
    """
    Reconstruct left and right states at all interior faces.

    Returns
    -------
    U_L : shape (ncells+1, 3) — left state at each face
    U_R : shape (ncells+1, 3) — right state at each face

    Face indexing: face i separates cell i-1 (left) and cell i (right).
    Faces 0 and ncells are the boundary faces.
    """
    ncells = U.shape[0]
    nfaces = ncells + 1

    # Work in primitive variables for reconstruction
    W = cons_to_prim(U)   # shape (ncells, 3)

    # Initialize face states — will be filled by reconstruction + BCs
    W_L = np.zeros((nfaces, 3))
    W_R = np.zeros((nfaces, 3))

    # -- Interior faces (1 through ncells-1) --------------------------------
    # For face i (between cell i-1 and cell i):
    #   Left state comes from cell i-1 extrapolated rightward
    #   Right state comes from cell i   extrapolated leftward
    for i in range(1, nfaces - 1):
        # Slopes for cell i-1 (contributing left state at face i)
        if i >= 2:
            dW_left_im1  = W[i-1] - W[i-2]
        else:
            dW_left_im1  = np.zeros(3)       # ghost: zero slope at left boundary
        dW_right_im1 = W[i]   - W[i-1]

        r_im1 = dW_left_im1 / (dW_right_im1 + 1e-300 * np.sign(dW_right_im1 + 1e-300))
        phi_im1 = limiter_fn(r_im1)
        slope_im1 = phi_im1 * dW_right_im1
        W_L[i] = W[i-1] + 0.5 * slope_im1   # extrapolate rightward to face

        # Slopes for cell i (contributing right state at face i)
        dW_left_i  = W[i]   - W[i-1]
        if i < nfaces - 2:
            dW_right_i = W[i+1] - W[i]
        else:
            dW_right_i = np.zeros(3)          # ghost: zero slope at right boundary

        r_i = dW_left_i / (dW_right_i + 1e-300 * np.sign(dW_right_i + 1e-300))
        phi_i = limiter_fn(r_i)
        slope_i = phi_i * dW_right_i
        W_R[i] = W[i] - 0.5 * slope_i       # extrapolate leftward to face

    # Boundary faces handled separately in apply_boundary_conditions()

    # Convert back to conserved variables
    U_L = prim_to_cons(W_L)
    U_R = prim_to_cons(W_R)

    return U_L, U_R


def first_order_states(U):
    """
    First-order (piecewise constant) face states.
    Simplest possible: just copy cell value to neighboring faces.
    """
    ncells = U.shape[0]
    nfaces = ncells + 1

    U_L = np.zeros((nfaces, 3))
    U_R = np.zeros((nfaces, 3))

    # Interior faces
    U_L[1:-1] = U[:-1]   # left state at face i = cell i-1
    U_R[1:-1] = U[1:]    # right state at face i = cell i

    return U_L, U_R


# =============================================================================
# 8. BOUNDARY CONDITIONS
#
#    Boundary conditions are applied by setting ghost cell states at the
#    two boundary faces (face 0 at left, face ncells at right).
#
#    The ghost cell approach means boundary faces are treated identically
#    to interior faces by the Riemann solver — it just sees two states.
#
#    'transmissive': copies interior state to ghost (waves pass out cleanly)
#    'reflective':   mirrors normal velocity (solid wall, no penetration)
# =============================================================================

def apply_boundary_conditions(U_L, U_R, U, bc_left, bc_right):
    """
    Fill boundary face states (face 0 and face ncells).

    U_L, U_R : shape (nfaces, 3) — modified in-place
    U        : shape (ncells, 3) — interior solution
    """
    # -- Left boundary (face 0): ghost is to the LEFT of the domain ----------
    if bc_left == 'transmissive':
        # Ghost cell = copy of first interior cell
        # Wave leaves domain without reflection
        U_L[0] = U[0]
        U_R[0] = U[0]
    elif bc_left == 'reflective':
        # Ghost cell mirrors velocity: u_ghost = -u_interior
        U_ghost    = U[0].copy()
        U_ghost[1] = -U[0, 1]   # negate rho*u
        U_L[0] = U_ghost
        U_R[0] = U[0]

    # -- Right boundary (face ncells): ghost is to the RIGHT of the domain ---
    if bc_right == 'transmissive':
        U_L[-1] = U[-1]
        U_R[-1] = U[-1]
    elif bc_right == 'reflective':
        U_ghost    = U[-1].copy()
        U_ghost[1] = -U[-1, 1]
        U_L[-1] = U[-1]
        U_R[-1] = U_ghost

    return U_L, U_R


# =============================================================================
# 9. RESIDUAL ASSEMBLY
#
#    The residual R is the spatial operator: dU/dt = -R(U)
#
#    In FVM, R is assembled by sweeping all faces and accumulating
#    flux contributions into each cell:
#
#        R[i] = (F*[i+1] - F*[i]) / dx
#
#    where F*[i] is the numerical flux at the left face of cell i
#    and   F*[i+1] is the numerical flux at the right face of cell i.
#
#    This is the discrete divergence theorem.
# =============================================================================

def compute_residual(U, dx, scheme, limiter_name):
    """
    Compute the spatial residual R such that dU/dt = -R(U).

    Parameters
    ----------
    U            : shape (ncells, 3) — current conserved state
    dx           : float — uniform cell width
    scheme       : 'first_order' or 'muscl'
    limiter_name : str — name of slope limiter (only used if muscl)

    Returns
    -------
    R : shape (ncells, 3) — residual
    """
    limiter_fn = LIMITERS.get(limiter_name, van_leer)

    # Step 1: Reconstruct face states
    if scheme == 'muscl':
        U_L, U_R = muscl_reconstruct(U, limiter_fn)
    else:
        U_L, U_R = first_order_states(U)

    # Step 2: Apply boundary conditions to boundary faces
    U_L, U_R = apply_boundary_conditions(U_L, U_R, U, BC_LEFT, BC_RIGHT)

    # Step 3: Compute numerical flux at all faces using HLLC
    # Use vectorized batch version for speed
    F_face = hllc_flux_batch(U_L, U_R)   # shape (nfaces, 3)

    # Step 4: Accumulate flux divergence into each cell
    # Cell i receives flux from its right face minus its left face
    R = (F_face[1:] - F_face[:-1]) / dx   # shape (ncells, 3)

    return R


# =============================================================================
# 10. SSP-RK3 TIME INTEGRATOR
#
#     Strong Stability Preserving 3-stage Runge-Kutta (Shu & Osher 1988).
#     Third-order accurate in time, maintains TVD property of spatial scheme.
#
#     U^(1) = U^n + dt * L(U^n)
#     U^(2) = 3/4 * U^n + 1/4 * (U^(1) + dt * L(U^(1)))
#     U^(n+1) = 1/3 * U^n + 2/3 * (U^(2) + dt * L(U^(2)))
#
#     where L(U) = -R(U) is the spatial operator (negative residual).
# =============================================================================

def ssp_rk3_step(U, dt, dx, scheme, limiter_name):
    """
    Advance U by one time step using SSP-RK3.

    Parameters
    ----------
    U            : shape (ncells, 3) — solution at time n
    dt           : float — time step
    dx           : float — cell width
    scheme       : 'first_order' or 'muscl'
    limiter_name : str

    Returns
    -------
    U_new : shape (ncells, 3) — solution at time n+1
    """
    def L(U_):
        return -compute_residual(U_, dx, scheme, limiter_name)

    # Stage 1
    U1 = U + dt * L(U)

    # Stage 2
    U2 = 0.75 * U + 0.25 * (U1 + dt * L(U1))

    # Stage 3
    U_new = (1.0/3.0) * U + (2.0/3.0) * (U2 + dt * L(U2))

    return U_new


# =============================================================================
# 11. CFL CONDITION
#
#     The time step must be small enough that information doesn't travel
#     more than one cell width per step. The CFL condition is:
#
#         dt = CFL * dx / max(|u| + a)
#
#     where (|u| + a) is the maximum signal speed in any cell.
#     CFL < 1 guarantees stability for explicit schemes.
# =============================================================================

def compute_dt(U, dx, cfl, t, t_end):
    """
    Compute stable time step.
    Also clips dt so we land exactly on t_end.
    """
    s_max = max_wave_speed(U)
    dt    = cfl * dx / s_max
    dt    = min(dt, t_end - t)    # don't overshoot end time
    return dt


# =============================================================================
# 12. INITIAL CONDITIONS — SOD SHOCK TUBE
#
#     The Sod shock tube is the classic test for 1D Euler solvers.
#     A diaphragm at x=0.5 separates two gas states. When removed at t=0,
#     three waves emerge:
#       - Left-running rarefaction wave
#       - Stationary contact discontinuity
#       - Right-running shock wave
#
#     All three wave families are excited, making it a comprehensive test.
# =============================================================================

def setup_mesh():
    """
    Build a uniform 1D mesh.
    Returns cell centers and cell width.
    """
    dx      = (X_END - X_START) / NCELLS
    x_faces = np.linspace(X_START, X_END, NCELLS + 1)
    x_cells = 0.5 * (x_faces[:-1] + x_faces[1:])   # cell centers
    return x_cells, dx


def sod_initial_condition(x_cells):
    """
    Sod shock tube: left/right states separated at diaphragm.
    Returns U shape (ncells, 3).
    """
    ncells = len(x_cells)
    U      = np.zeros((ncells, 3))

    x_diaphragm = X_START + DIAPHRAGM * (X_END - X_START)

    for i, x in enumerate(x_cells):
        W      = W_LEFT if x < x_diaphragm else W_RIGHT
        U[i]   = prim_to_cons(W)

    return U


# =============================================================================
# 13. ANALYTICAL SOLUTION FOR SOD SHOCK TUBE
#
#     The exact solution consists of four constant states separated by
#     three waves. We solve for the intermediate states iteratively
#     using the Rankine-Hugoniot and isentropic relations.
#
#     This is used to verify the numerical solution.
# =============================================================================

def sod_analytical(x_cells, t):
    """
    Exact solution for the Sod shock tube problem at time t.

    Uses an iterative Newton solve for the intermediate pressure p*.
    Valid for any left/right initial conditions with the same gamma.

    Returns W_exact shape (ncells, 3) — [rho, u, p] at each cell.
    """
    from scipy.optimize import brentq

    rho_L, u_L, p_L = W_LEFT
    rho_R, u_R, p_R = W_RIGHT
    a_L = sound_speed(rho_L, p_L)
    a_R = sound_speed(rho_R, p_R)

    g  = GAMMA
    g1 = (g - 1) / (2 * g)
    g2 = (g + 1) / (2 * g)
    g3 = 2 * g / (g - 1)
    g4 = 2 / (g - 1)
    g5 = 2 / (g + 1)
    g6 = (g - 1) / (g + 1)

    def f_L(p):
        """Relation across left wave (rarefaction if p < p_L, shock if p > p_L)."""
        if p <= p_L:
            return g4 * a_L * ((p/p_L)**g1 - 1)
        else:
            A = g5 / rho_L
            B = g6 * p_L
            return (p - p_L) * np.sqrt(A / (p + B))

    def f_R(p):
        """Relation across right wave."""
        if p <= p_R:
            return g4 * a_R * ((p/p_R)**g1 - 1)
        else:
            A = g5 / rho_R
            B = g6 * p_R
            return (p - p_R) * np.sqrt(A / (p + B))

    def f(p):
        return f_L(p) + f_R(p) + (u_R - u_L)

    # Solve for intermediate pressure p*
    p_min = min(p_L, p_R) * 1e-6
    p_max = max(p_L, p_R) * 100
    p_star = brentq(f, p_min, p_max, xtol=1e-12)
    u_star = 0.5*(u_L + u_R) + 0.5*(f_R(p_star) - f_L(p_star))

    # Intermediate densities from isentropic / Rankine-Hugoniot
    if p_star <= p_L:
        rho_starL = rho_L * (p_star/p_L)**(1/g)   # isentropic
    else:
        rho_starL = rho_L * (p_star/p_L + g6) / (g6 * p_star/p_L + 1)

    if p_star <= p_R:
        rho_starR = rho_R * (p_star/p_R)**(1/g)
    else:
        rho_starR = rho_R * (p_star/p_R + g6) / (g6 * p_star/p_R + 1)

    # Wave speeds
    if p_star <= p_L:
        # Left rarefaction: head and tail speeds
        S_HL = u_L - a_L
        a_starL = a_L * (p_star/p_L)**g1
        S_TL = u_star - a_starL
    else:
        # Left shock speed
        S_L = u_L - a_L * np.sqrt(g2 * p_star/p_L + g1)

    if p_star <= p_R:
        # Right rarefaction
        S_HR = u_R + a_R
        a_starR = a_R * (p_star/p_R)**g1
        S_TR = u_star + a_starR
    else:
        # Right shock speed
        S_R = u_R + a_R * np.sqrt(g2 * p_star/p_R + g1)

    # Contact discontinuity speed
    S_contact = u_star

    # Evaluate exact solution at each cell center
    x_diaphragm = X_START + DIAPHRAGM * (X_END - X_START)
    xi = (x_cells - x_diaphragm) / (t + 1e-300)   # similarity variable

    W_exact = np.zeros((len(x_cells), 3))

    for i, s in enumerate(xi):
        if p_star <= p_L:
            # Left rarefaction fan
            if s <= S_HL:
                W_exact[i] = [rho_L, u_L, p_L]
            elif s <= S_TL:
                # Inside rarefaction
                u_fan   = g5 * (a_L + (s - u_L)*(g-1)/2 / 1 + u_L)
                # More precisely:
                u_fan   = (2/(g+1)) * (a_L + (g-1)/2 * u_L + s)
                a_fan   = u_fan - s
                rho_fan = rho_L * (a_fan/a_L)**g4
                p_fan   = p_L   * (a_fan/a_L)**g3
                W_exact[i] = [rho_fan, u_fan, p_fan]
            elif s <= S_contact:
                W_exact[i] = [rho_starL, u_star, p_star]
            elif s <= (S_R if p_star > p_R else S_TR):
                W_exact[i] = [rho_starR, u_star, p_star]
            else:
                W_exact[i] = [rho_R, u_R, p_R]
        else:
            # Left shock
            if s <= S_L:
                W_exact[i] = [rho_L, u_L, p_L]
            elif s <= S_contact:
                W_exact[i] = [rho_starL, u_star, p_star]
            elif s <= (S_R if p_star > p_R else S_TR):
                W_exact[i] = [rho_starR, u_star, p_star]
            else:
                W_exact[i] = [rho_R, u_R, p_R]

    return W_exact


# =============================================================================
# 14. MAIN SOLVER LOOP
#
#     This is the time-marching loop. Each iteration:
#       1. Compute stable dt from CFL condition
#       2. Advance U one time step with SSP-RK3
#       3. Check physical validity
#       4. Print progress
#
#     All solution data is stored in plain numpy arrays — inspect freely
#     in the Spyder variable explorer after the run.
# =============================================================================

def run_solver():
    """
    Run the simulation. Returns all solution data for post-processing.
    """
    print("=" * 60)
    print("  PyFVM — 1D Euler Solver")
    print("=" * 60)
    print(f"  Cells:   {NCELLS}")
    print(f"  Scheme:  {SCHEME}" + (f" + {LIMITER}" if SCHEME == 'muscl' else ""))
    print(f"  CFL:     {CFL}")
    print(f"  t_end:   {T_END}")
    print(f"  BC:      left={BC_LEFT}, right={BC_RIGHT}")
    print("=" * 60)

    # -- Setup ---------------------------------------------------------------
    x_cells, dx = setup_mesh()
    U            = sod_initial_condition(x_cells)

    # Storage for convergence monitoring
    t          = 0.0
    iteration  = 0
    time_hist  = [t]
    mass_hist  = [np.sum(U[:, 0]) * dx]   # total mass (should be conserved)

    # -- Time loop -----------------------------------------------------------
    while t < T_END:
        dt = compute_dt(U, dx, CFL, t, T_END)

        U  = ssp_rk3_step(U, dt, dx, SCHEME, LIMITER)

        t         += dt
        iteration += 1

        # Conservation check (mass should be exactly conserved)
        mass = np.sum(U[:, 0]) * dx
        mass_hist.append(mass)
        time_hist.append(t)

        if iteration % PRINT_INTERVAL == 0:
            W     = cons_to_prim(U)
            s_max = max_wave_speed(U)
            print(f"  iter={iteration:5d}  t={t:.4f}  dt={dt:.2e}  "
                  f"max_speed={s_max:.3f}  "
                  f"mass={mass:.6f}")

        # Validity check (catches blow-ups early)
        check_physical(U, f"at iteration {iteration}")

    print(f"\n  Done. {iteration} iterations, t={t:.4f}")
    print(f"  Mass conservation error: {abs(mass_hist[-1] - mass_hist[0]):.2e}")

    return x_cells, dx, U, np.array(time_hist), np.array(mass_hist)


# =============================================================================
# 15. POST-PROCESSING
#
#     Plot the solution and compare to the analytical result.
#     All arrays are available in the Spyder variable explorer.
# =============================================================================

def plot_results(x_cells, U, t_final):
    """
    Plot rho, u, p with optional analytical solution overlay.
    """
    W = cons_to_prim(U)   # shape (ncells, 3) — rho, u, p

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle(
        f"Sod Shock Tube  |  t={t_final:.3f}  |  "
        f"N={NCELLS}  |  {SCHEME}"
        + (f"+{LIMITER}" if SCHEME == 'muscl' else ""),
        fontsize=12
    )

    labels   = ['Density  ρ', 'Velocity  u', 'Pressure  p']
    colors   = ['steelblue', 'darkorange', 'forestgreen']

    for ax, col, label, var_idx in zip(axes, colors, labels, range(3)):
        ax.plot(x_cells, W[:, var_idx],
                color=col, linewidth=1.8, label='Numerical')

        if PLOT_ANALYTICAL:
            W_exact = sod_analytical(x_cells, t_final)
            ax.plot(x_cells, W_exact[:, var_idx],
                    'k--', linewidth=1.2, label='Exact', alpha=0.7)

        ax.set_xlabel('x')
        ax.set_ylabel(label)
        ax.set_xlim(X_START, X_END)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
    return fig


def plot_conservation(time_hist, mass_hist):
    """
    Plot mass conservation over time. Should be flat to machine precision.
    """
    fig, ax = plt.subplots(figsize=(7, 3))
    error = np.abs(mass_hist - mass_hist[0])
    ax.semilogy(time_hist, error + 1e-16)
    ax.set_xlabel('Time')
    ax.set_ylabel('|Mass error|')
    ax.set_title('Mass Conservation')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    return fig


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == '__main__' or True:
    # The `or True` means this runs whether you execute the file directly
    # OR run it with F5 in Spyder. Remove `or True` if you want to import
    # this file from another script without auto-running.

    # Run the solver
    x, dx, U, time_hist, mass_hist = run_solver()

    # Primitive variables — available in variable explorer after run
    W = cons_to_prim(U)
    rho = W[:, 0]
    u   = W[:, 1]
    p   = W[:, 2]
    T   = temperature(rho, p)
    a   = sound_speed(rho, p)
    M   = u / a   # Mach number

    # Plots
    fig_solution     = plot_results(x, U, T_END)
    fig_conservation = plot_conservation(time_hist, mass_hist)